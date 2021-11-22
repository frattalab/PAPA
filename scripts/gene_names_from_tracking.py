#!/usr/bin/env python3

import pandas as pd
from papa_helpers import eprint
import sys
import os
from timeit import default_timer as timer


def main(tracking_path, tx2le_path, out_path):
    '''
    '''

    tracking_default_cols = ["post_merge_tx_id",
                             "post_merge_gene_id",
                             "ref_gene_id",
                             "overlap_class_code",
                             "novel_gtf",
                             "ref_gtf"]

    eprint(f"Reading in GFFcompare '<prefix>.tracking' file - {tracking_path}")

    tracking = pd.read_csv(tracking_path,
                           sep="\t",
                           names=tracking_default_cols,
                           usecols=[0,1,5])


    #extract gene_name from each 'ref_gtf' transcript string
    # e.g. q1:XLOC_000056|PAPA_00000154|3|0.000000|0.000000|0.000000|2251
    # want XLOC_000056
    eprint("Extracting reference gene names...")
    tracking["ref_gene_name"] = (tracking["ref_gtf"].str.split(":", expand=True)[1]
                                                    .str.split("\\|", expand=True)[0])

    # eprint(tracking)

    # Generate 2 col df of post_merge_gene_id | ref_gene_names
    # Where ref_gene_names = ref_gene_name sep by ','
    eprint("Collapsing reference gene names by gene ID...")

    gene2name = (tracking.dropna(subset=["ref_gene_name"])
                         .groupby("post_merge_gene_id")["ref_gene_name"]
                         .agg(ref_gene_names=lambda names: ",".join(names.unique()))
                         .reset_index())

    # eprint(gene2name)

    tracking = tracking.merge(gene2name, how="left", on="post_merge_gene_id")
    # eprint(tracking)

    # Merge tx2le with gene2name to get df of
    # tx_id | le_id | gene_id (post_merge) | ref_gene_names

    eprint(f"Reading in '<prefix>.tx2le.tsv' file - {tx2le_path}")
    tx2le = pd.read_csv(tx2le_path, sep="\t")

    eprint("Generating combined df of tx_id | le_id | gene_id | ref_gene_names ...")
    combined = (tx2le.merge(tracking[["post_merge_tx_id", "post_merge_gene_id", "ref_gene_names"]],
                            how="left",
                            left_on="transcript_id",
                            right_on="post_merge_tx_id")
                     .drop("post_merge_tx_id", axis=1)
                     .rename({"post_merge_gene_id": "gene_id"}, axis=1))

    # eprint(combined)

    eprint(f"Number of genes with multiple assigned reference gene names - {combined.loc[combined['ref_gene_names'].str.contains(','), 'gene_id'].nunique()}")
    eprint(f"Number of genes with single assigned reference gene name - {combined.loc[~combined['ref_gene_names'].str.contains(','), 'gene_id'].nunique()}")

    eprint(f"Writing combined TSV file to {out_path}")
    combined.to_csv(out_path,
                    sep="\t",
                    header=True,
                    index=False)

if __name__ == '__main__':

    start = timer()

    if len(sys.argv) == 1:
        eprint("Usage information")
        eprint("python gene_names_from_tracking.py <prefix>.tracking <prefix>.tx2le.tsv <output_path>")
        sys.exit(0)

    tracking_path = sys.argv[1]
    tx2le_path = sys.argv[2]
    out_path = sys.argv[3]

    main(tracking_path, tx2le_path, out_path)

    end = timer()

    eprint(f"Complete: took {round(end - start, 3)} s / {round((end - start) / 60, 3)} min (3 dp)")
