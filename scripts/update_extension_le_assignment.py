import pandas as pd
import numpy as np
from papa_helpers import eprint
import os
import argparse
import sys
from timeit import default_timer as timer


'''
If merge expression of transcripts according to shared last exons, 3'UTR extension events will be masked/aggregated along with the proximal event
Since the purpose of the 'expression_merge_by' == "last_exon" pipeline is to simplify interpretation of differential usage of novel events
Need to consider extension events as a 'different last exon' for them to be quantified separately to isoforms using the reference exon
This script will modify the 'le' output of assign_tx_to_pas.py to consider extensions as separate events.

This script takes as input:
- '<prefix>.pas_assignment.tsv' output by assign_tx_to_pas.py
- '<prefix>.tx2le.tsv' output by assign_tx_to_pas.py
- '<prefix>.le2gene.tsv' output by assign_tx_to_pas.py
- '<prefix>.match_stats.tsv' output by filter_tx_by_three_end.py
- '<prefix>'.loci file output by GFFcompare (rule merge_filtered_with_ref)

Script outputs:
- <output_prefix>.tx2le.tsv,
- <output_prefix>.pas_assignment.tsv',
- <output_prefix>.le2gene.tsv with modified/extra le_ids for extension events

1. Read in .loci file and get novel IDs into a column of lists
.loci file has the general structure
gene_id | locus/region | unknown | file1_tr | file2_tr |file..n_tr

where file numbers are in the order they were provided to GFFcompare call
each 'transcript' column (file1_tr, file2_tr) contains comma separated string of transcript IDs assigned to that locus
e.g. PAPA_00000377,PAPA_00000378,PAPA_00004797
e.g. ENST00000420770.7,ENST00000375375.7,ENST00000400661.3

since rule merge_filtered_with_ref always passes the filtered novel events first
can safely assume 4th column == col containing novel transcripts

2. Extract valid extension event IDs from '<prefix>.match_stats.tsv'

3. Identify extension event IDs from .loci df
    - Usually expect a single event, but if gene has multiple extension events, group as part of same?
    - Or load in GTF and reclassify?
    - Or check where in current gene order they appear, update all if necessary
        - E.g. extract out old le number, add one to all ds events, slot in proximal there, slot in distal at end

4. Update le_id to be one greater than max for gene in all le files

5. Output modified dfs
'''

def update_extension_le_id(df):
    '''
    '''

    # assert
    out = {}
    for name, group in df.groupby("gene_id"):
        # Track how many extension IDs for gene starting from 5'end
        # Once encounter an ext, all IDs downstream need to be shifted down by n of preceeding extensions
        ext_count = 0

        # Track le_id of most recently found ext event
        # If gene has alt extension isoforms of the same last exon
        # Will group together for simplification's sake
        prev_ext_le_id = ""

        for row_n, row in group.iterrows():
            if row["is_extension"] == 1:
                if row["le_id"] == prev_ext_le_id:
                    # Tx is a diff extension of same le - group together
                    out[row["transcript_id"]] = row["pre_le_number"] + ext_count

                else:
                    # New distinct extension event
                    ext_count += 1
                    out[row["transcript_id"]] = row["pre_le_number"] + ext_count

            elif row["is_extension"] == 0:
                out[row["transcript_id"]] = row["pre_le_number"] + ext_count


    df = (pd.Series(out, name="post_le_number")
            .to_frame()
            .reset_index()
            .rename({"index": "transcript_id"},
                    axis=1)
          )

    return df


def main(pas_assignment_path, tx2le_path, le2gene_path, match_stats_path, loci_path, tracking_path, output_prefix):
    '''
    '''

    #1. Read in GFFcompare .loci file and convert the two Tx columns into cols of lists
    loci_cols = ["gene_id", "locus", "unknown", "novel_ids", "ref_ids"]
    loci = pd.read_csv(loci_path, sep="\t", names=loci_cols)


    # cols with '-' in novel_ids have no novel tx contributing to that gene
    # Let's remove to avoid unnecessary processing
    loci = loci.loc[loci["novel_ids"] != "-", :]

    # Split comma separated strings into lists of transcript IDs
    loci[["novel_ids", "ref_ids"]] = loci[["novel_ids", "ref_ids"]].apply(lambda col: col.str.split(","), axis=0)

    #2. get a set of IDs of UTR extension events
    match_stats = pd.read_csv(match_stats_path, sep="\t")

    pre_merge_ext_ids = set(match_stats.loc[(match_stats["match_class"] == "valid") &
                                            (match_stats["isoform_class"] == "ds_3utr_extension"),
                                            "transcript_id_novel"]
                            )

    # eprint(len(pre_merge_ext_ids))

    # Note that IDs from match stats come from pre-reference merging
    # i.e. IDs in pas_assignment, tx2le etc. will be different
    # Transcript ID tracking is stored in '.tracking' file
    # Process in .tracking file to df of 'pre_merge_tx_id' | post_merge_tx_id

    tracking_default_cols = ["post_merge_tx_id",
                             "post_merge_gene_id",
                             "ref_gene_id",
                             "overlap_class_code",
                             "novel_gtf",
                             "ref_gtf"]

    tracking = pd.read_csv(tracking_path,
                           sep="\t",
                           names=tracking_default_cols,
                           usecols=[0,1,4])

    #Extract pre_merge_tx_id from the transcript info string
    # e.g. q1:XLOC_000056|PAPA_00000154|3|0.000000|0.000000|0.000000|2251

    tracking["pre_merge_tx_id"] = tracking["novel_gtf"].str.split("\\|", expand=True)[1]

    tracking = tracking[["pre_merge_tx_id", "post_merge_tx_id", "post_merge_gene_id"]]

    # None are reference transcripts, don't need to carry them around
    tracking = tracking[~tracking["pre_merge_tx_id"].isna()]
    # eprint(tracking)


    #3. Add a column with the matching extension ID for a given gene
    loci["ext_ids"] = loci["novel_ids"].apply(lambda x: [novel_tx for novel_tx in x if novel_tx in pre_merge_ext_ids])

    loci["n_extension_ids"] = loci["ext_ids"].apply(len)

    loci = loci.loc[loci["n_extension_ids"] != 0, :]

    # eprint(loci["ext_ids"].apply(len).describe())


    #4. Read in tx2le & le2gene - create df of tx_id | le_id | gene_id
    tx2le = pd.read_csv(tx2le_path, sep="\t")
    le2gene = pd.read_csv(le2gene_path, sep="\t")

    tx2le2gene = tx2le.merge(le2gene, on="le_id")

    # eprint(tx2le2gene)

    #5. Extract genes with extension IDs

    ext_gene_ids = set(loci.loc[loci["n_extension_ids"] != 0, "gene_id"])
    ext_tx2le2gene = tx2le2gene.loc[tx2le2gene["gene_id"].isin(ext_gene_ids), :]


    #For each gene, count the total number of last exons
    ext_tx2le2gene = ext_tx2le2gene.merge((ext_tx2le2gene.groupby("gene_id")
                                           ["le_id"].nunique()),
                                           on="gene_id",
                                           suffixes=[None, "_count"])

    # eprint(ext_tx2le2gene["le_id_count"].describe())

    # Pull out 'assigned LE number' from le_id string
    # XLOC_000005_1 - '1' = first last exon in gene
    # Subset to Txs where extension is most 3' in the gene
    ext_tx2le2gene["pre_le_number"] = ext_tx2le2gene["le_id"].str.split("_", expand=True).iloc[:, -1].astype(int)

    # Add column marking whether a given transcript ID is an extension event
    # tx2le transcript_ids are based on post-merge IDs
    # Need to grab 'updated' extension IDs from tracking
    post_merge_ext_ids = set(tracking.loc[tracking["pre_merge_tx_id"].isin(pre_merge_ext_ids),
                                          "post_merge_tx_id"])

    # eprint(len(post_merge_ext_ids))

    ext_tx2le2gene["is_extension"] = np.where(ext_tx2le2gene["transcript_id"].isin(post_merge_ext_ids),
                                              1,
                                              0)

    # eprint(ext_tx2le2gene["is_extension"].value_counts())

    # Update assigned last exon number for genes with extension events
    # Extra sort so extension events come after non-extensions of the same le_id
    ext_tx2le2gene = ext_tx2le2gene.sort_values(by=["gene_id", "pre_le_number", "is_extension"])

    # eprint(ext_tx2le2gene.loc[ext_tx2le2gene.gene_id == "XLOC_056316", :])


    # eprint(ext_tx2le2gene)

    upd_le_numbers = update_extension_le_id(ext_tx2le2gene)

    ext_tx2le2gene = ext_tx2le2gene.merge(upd_le_numbers, on="transcript_id")

    # eprint(ext_tx2le2gene)



    # Create new le_id
    ext_tx2le2gene["new_le_id"] = ext_tx2le2gene["gene_id"] + "_" + ext_tx2le2gene["post_le_number"].astype(str)
    ext_tx2le2gene = (ext_tx2le2gene.drop(["le_id",
                                           "le_id_count",
                                           "is_extension",
                                           "pre_le_number"],
                                          axis=1)
                                    .rename({"new_le_id": "le_id",
                                             "post_le_number": "le_number"},
                                            axis=1)
                      )

    # eprint(ext_tx2le2gene)

    # Create updated '.pas_assignment.tsv'
    # First - return pas_number & pas_id to ext_tx2le2gene
    pas_assignment = pd.read_csv(pas_assignment_path, sep="\t")
    ext_tx2le2gene = ext_tx2le2gene.merge(pas_assignment[["transcript_id",
                                                          "pas_number",
                                                          "pas_id"]],
                                                          on="transcript_id")

    # Remove genes with extensions from pas_assignment, concat with ext_tx2le2gene
    upd_pas_assignment = pd.concat([pas_assignment[~pas_assignment["gene_id"].isin(ext_gene_ids)],
                                    ext_tx2le2gene],
                                    ignore_index=True)

    upd_pas_assignment.to_csv(output_prefix + ".pas_assignment.tsv",
                              header=True,
                              index=False,
                              sep="\t")

    ## Create updated 'tx2le' file
    not_ext_tx2le2gene = tx2le2gene[~tx2le2gene["gene_id"].isin(ext_gene_ids)]

    upd_tx2le = pd.concat([not_ext_tx2le2gene[["le_id","transcript_id"]],
                           ext_tx2le2gene[["le_id", "transcript_id"]]
                           ], ignore_index=True)

    upd_tx2le.drop_duplicates().to_csv(output_prefix + ".tx2le.tsv",
                                       header=True,
                                       index=False,
                                       sep="\t")

    ## Create updated 'le2gene' file
    upd_le2gene = pd.concat([not_ext_tx2le2gene[["le_id","gene_id"]],
                             ext_tx2le2gene[["le_id","gene_id"]]],
                            ignore_index=True)

    upd_le2gene.drop_duplicates().to_csv(output_prefix + ".le2gene.tsv",
                                         header=True,
                                         index=False,
                                         sep="\t")


if __name__ == '__main__':

    start = timer()

    descrpn = """Update transcript to shared 'last exon' assignment such that novel extension events are given a distinct ID for differential usage"""

    parser = argparse.ArgumentParser(description=descrpn,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter, # Add defaults to end of help string
                                     )


    parser.add_argument("-a",
                        "--pas-assignment",
                        type=str,
                        required=True,
                        default=argparse.SUPPRESS,
                        help="Path to '<prefix>.pas_assignment.tsv' file output by assign_tx_to_pas.py")

    parser.add_argument("-t",
                        "--tx2le",
                        type=str,
                        required=True,
                        default=argparse.SUPPRESS,
                        help="Path to '<prefix>.tx2le.tsv' output by assign_tx_to_pas.py"
                        )

    parser.add_argument("-g",
                        "--le2gene",
                        type=str,
                        required=True,
                        default=argparse.SUPPRESS,
                        help="Path to '<prefix>.le2gene.tsv' output by assign_tx_to_pas.py")

    parser.add_argument("-m",
                        "--match-stats",
                        type=str,
                        required=True,
                        default=argparse.SUPPRESS,
                        help="Path to '<prefix>.match_stats.tsv' output by filter_tx_by_three_end.py"
                        )

    parser.add_argument("-l",
                        "--loci",
                        type=str,
                        required=True,
                        default=argparse.SUPPRESS,
                        help="Path to '<prefix>'.loci file output by GFFcompare (rule merge_filtered_with_ref)")

    parser.add_argument("--tracking",
                        type=str,
                        required=True,
                        default=argparse.SUPPRESS,
                        help="Path to '<prefix>'.tracking file output by GFFcompare (rule merge_filtered_with_ref)")


    parser.add_argument("-o",
                        "--output_prefix",
                        type=str,
                        default="upd_extension_ids",
                        help="path to/prefix for output files. '.pas_assignment.tsv' is suffixed for updated --pas-assignment TSV file, '.tx2le.tsv' for updated --tx2le input & '.le2gene.tsv' for updated --le2gene input.")


    if len(sys.argv) == 1:
        parser.print_help()
        parser.exit()


    args = parser.parse_args()

    main(args.pas_assignment,
         args.tx2le,
         args.le2gene,
         args.match_stats,
         args.loci,
         args.tracking,
         args.output_prefix)

    end = timer()

    eprint(f"Complete: took {round(end - start, 3)} s / {round((end - start) / 60, 3)} min (3 dp)")
