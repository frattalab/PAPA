#!/usr/bin/env python3

#     Script to combine a series of PAPA predicted last exons into a single GTF file with expression information appended
#     Copyright (C) 2024  Sam Bryce-Smith samuel.bryce-smith.19@ucl.ac.uk

#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.

#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.

#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <https://www.gnu.org/licenses/>.

import pyranges as pr
import pandas as pd
import numpy as np
import os
import sys
import argparse

# Let's loop over all 3 datasets and perform same pre-processing
# goal - add PPAU, gene_name, le_id, experiment_id to each GTF of novel last exons
# Combine all LEs into a single df and annotate overlaps etc.
# Output a combined GTF for downstream analysis


gtf_name = "all_conditions.merged_last_exons.3p_end_filtered.gtf"
tx2le_name = "novel_ref_combined.tx2le.tsv"
tx2gene_name = "novel_ref_combined.le2gene.tsv"
tx2genename_name = "novel_ref_combined.le2genename.tsv"
ppau_name = "summarised_pas_quantification.ppau.tsv"

gtf_cols_to_keep = ["gene_id", "transcript_id", "event_type", "3p_extension_length", "atlas_filter", "nearest_atlas_distance","Name", "motif_filter", "pas_motifs", "condition_id"]
ppau_cols_to_keep = [0,-3,-2,-1]
ppau_cols_general = ["mean_PPAU_base", "mean_PPAU_treatment", "delta_PPAU_treatment_control"]



def main(in_dir,
         out_gtf):
    processed = []
    for exp in [os.path.join(in_dir, d) for d in os.listdir(in_dir)]:
        exp_name = os.path.basename(exp)
        print(f"Inferred experiment name/ID - {exp_name}")
        
        gtf = pr.read_gtf(os.path.join(exp, gtf_name))
        # keep minimal cols
        gtf = gtf[gtf_cols_to_keep]
        
        tx2le = pd.read_csv(os.path.join(exp, tx2le_name), sep="\t")
        le2name = pd.read_csv(os.path.join(exp, tx2genename_name), sep="\t")
        le2gene = pd.read_csv(os.path.join(exp, tx2gene_name), sep="\t")
        
        # join le_id & other ref annotation IDs
        gtf = gtf.apply(lambda df: df.merge(tx2le, on="transcript_id", how="left"))
        gtf = gtf.apply(lambda df: df.merge(le2name, on="le_id", how="left"))
        gtf = gtf.apply(lambda df: df.merge(le2gene, on="le_id", how="left", suffixes=[None, "_ref"]))
        
        # rename gene_name column to gene_name_ref for consistency with joined gene_id_ref
        gtf = gtf.apply(lambda df: df.rename(columns={"gene_name": "gene_name_ref"}))
        
        # subset to minimal PPAU info (le_id, means in each condition)
        ppau = pd.read_csv(os.path.join(exp, ppau_name), sep="\t")
        ppau_cols = ppau.columns[ppau_cols_to_keep]
        ppau = ppau.loc[:, ppau_cols]
        
        # rename condition-key specific PPAU cols to 'base' & 'treat' so diff dfs have general
        # le_id (first in list) stays the same
        ppau_cols_gen = {col: ppau_cols_general[i] for i, col in enumerate(ppau_cols[1:])}
        print(ppau_cols_gen)
        ppau = ppau.rename(columns=ppau_cols_gen)
        
        
        gtf = gtf.apply(lambda df: df.merge(ppau, on="le_id", how="left"))
        
        # add an 'experiment_id' to track where isoform was discovered
        gtf = gtf.assign("experiment_id",
                        lambda df: pd.Series([exp_name]*len(df), index=df.index))
        
        processed.append(gtf)
        

    processed[0]


    combined = pr.concat(processed)

    # some identified isoforms seem to have NA PPAUs after salmon (various possible reasons)
    quant_na_mask = combined.delta_PPAU_treatment_control.isna()
    print(f"Number of isoforms with NA quantifcation/PPAU values - {quant_na_mask.sum()}")
    # combined[quant_na_mask]


    # remove NA quants (very few, don't know how arise...)
    combined = combined[~quant_na_mask]
    # combined

    # Find last exons that overlap on the same strand (any overlap = reported, not considering 3'end consistency)
    combined = combined.cluster(strand=True)
    # combined


    # count number of experiments which have a given last exon (any overlap)
    le_counts = combined.as_df()[["Cluster", "experiment_id"]].drop_duplicates().Cluster.value_counts()
    le_counts = le_counts.reset_index().rename(columns={"Cluster": "experiment_count", "index": "Cluster"})
    print(le_counts.experiment_count.describe(percentiles=[i * 0.1 for i in range(0,11,1)]))

    combined = combined.apply(lambda df: df.merge(le_counts, on="Cluster", how="left"))
    # combined

    # Of the events only found in a single dataset, which ones do they come from?
    print("Events discovered in a single dataset")
    print(combined.subset(lambda df: df.experiment_count == 1).experiment_id.value_counts())


    combined = combined.assign("multiple_datasets",
                    lambda df: pd.Series(np.where(df["experiment_count"].gt(1), 1, 0)))


    # combined

    # events discovered in a single dataset only - are they biased towards lower/higher relative expression? (use max PPAU of either condition)

    combined = combined.assign("max_mean_PPAU",
                            lambda df: df[["mean_PPAU_base", "mean_PPAU_treatment"]].max(axis="columns"))

    print("Distribution of condition-wise max mean PPAU for each event")
    print((combined
    .as_df()
    .drop_duplicates(subset=["Cluster"])
    .groupby("multiple_datasets")
    ["max_mean_PPAU"]
    .describe(percentiles=[i * 0.1 for i in range(0,11,2)])
    ))

    print("Discovery in multiple datasets as a function of event type & motif/atlas discovery filter")
    print((combined
    .as_df()
    .assign(simple_event_type=lambda df: np.where(df["event_type"].str.contains("extension"), "extension", "spliced"))
    .drop_duplicates(subset=["Cluster"])
    .groupby("multiple_datasets")
    [["simple_event_type", "atlas_filter", "motif_filter"]]
    .value_counts(normalize=False).to_frame()
    ))

    print((combined
    .as_df()
    .assign(simple_event_type=lambda df: np.where(df["event_type"].str.contains("extension"), "extension", "spliced"))
    .drop_duplicates(subset=["Cluster"])
    .groupby("multiple_datasets")
    [["simple_event_type", "atlas_filter", "motif_filter"]]
    .value_counts(normalize=True).to_frame()
    ))

    # For shared events, what are the most common combinations of datasets?

    print("Number of events found in different combinations of datasets")
    
    print((combined.subset(lambda df: df.multiple_datasets.eq(1))
    .as_df()
    [["Cluster","experiment_id"]].
    groupby("Cluster")
    .agg(lambda col: ",".join(sorted(set(col.astype(str))))) # collapses to unique and sorts so same order.loc[lambda x: x.str.contains(",")]
    ["experiment_id"]
    # .loc[lambda x: x.str.contains(",")]
    .value_counts(normalize=False).to_frame()
    ))

    # Datasets an event identified in (fraction)
    print((combined.subset(lambda df: df.multiple_datasets.eq(1))
    .as_df()
    [["Cluster","experiment_id"]].
    groupby("Cluster")
    .agg(lambda col: ",".join(sorted(set(col.astype(str))))) # collapses to unique and sorts so same order
    ["experiment_id"]
    # .loc[lambda x: x.str.contains(",")]
    .value_counts(normalize=True).to_frame()
    ))

    # MAKE SURE HAVE 'exon' IN FEATURE FIELD SO GFFREAD WILL EXTRACT THE SEQUENCES!!!!!!!!!!!!!!!
    combined.Feature = "exon"

    # 1 event has 0 length interval (which is invalid)
    print(f"Writing combined GTF of last exons to - {out_gtf}")
    combined.subset(lambda df: (df.End - df.Start) > 0).to_gtf(out_gtf)


if __name__ == '__main__':
    
    descrpn = "Combine PAPA predicted novel last exons from multiple datasets into a single GTF"
    
    parser = argparse.ArgumentParser(description=descrpn,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter, # Add defaults to end of help strings
                                     )
    
    parser.add_argument("-i",
                        "--input-dir",
                        required=True,
                        default=argparse.SUPPRESS,
                        type=str,
                        help="Path to directory containing experiment-wise subdirectories of GTFs of predicted novel last exons to merge, along with tx2le, le2genename, le2gene & PPAU tables to additional metadata information. If multiple directories, pass paths consecutively and space-separated")
    
    parser.add_argument("-o",
                        "--output-gtf",
                        required=True,
                        type=str,
                        help="Name of output GTF file of combined last exons"
                        )

    if len(sys.argv) == 1:
        parser.print_help()
        parser.exit()

    args = parser.parse_args()
    
    main(args.input_dir,
         args.output_gtf)
