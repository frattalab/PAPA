#!/usr/bin/env python3

import pyranges as pr
import pandas as pd
import numpy as np
from papa_helpers import eprint, get_terminal_regions, _df_add_region_number
from filter_tx_by_intron_chain import _pd_merge_gr
import argparse
import os
import sys
from timeit import default_timer as timer

def cluster_to_region_number(gr, group_id_col, out_col="pas_number", cluster_col="Cluster"):
    '''
    Returns gr with 'out_col' column added
    where out_col is leftmost to rightmost cluster_col converted to a
    strand-aware 1..n order by group_id_col
    1 = most 5' site in group_id_col
    n = most 3' site in group_id_col
    '''

    c2p = gr.assign(out_col,
                    lambda df: _df_add_region_number(df, group_id_col, cluster_col))

    return c2p


def main(gtf_path, window_size, group_id_col, tx_id_col, output_prefix):
    '''
    '''

    eprint(f"Reading in input GTF file... - {gtf_path}")
    gtf = pr.read_gtf(gtf_path)

    assert group_id_col in gtf.columns
    assert gtf.stranded

    eprint(f"Extracting last exons for each transcript...")
    le = get_terminal_regions(gtf.subset(lambda df: df.Feature == "exon"),
                              source="stringtie")

    eprint("Extracting 3'ends for each transcript...")
    pas = le.three_end()

    # Extend either side of each PAS by window_size/2
    # Allows to merge closely spaced PAS as a single site
    eprint(f"Extending 3'ends by {(window_size -1)/2} nt either side to group together closely spaced poly(A) sites...")

    pas = pas.slack({"5": (window_size - 1)/2,
                     "3": (window_size - 1)/2}
                    )

    #Assign overlapping polyAsites a common ID with pr.cluster()
    #Returned 1..n 'leftmost to rightmost' in each df
    #i.e. regardless of strand, n+1 is the next rightmost 'grouped interval' from n
    # gr = pr.from_dict({"Chromosome": ["chr1"]*8, "Start": [1,1,1]*2 + [10,10], "End": [2,2,2]*2 + [11,11], "Strand": ["+"]*3+["-"]*5})
    # gr.cluster()

    eprint("Assigning 5'-3' pas_number for transcripts/poly(A) sites of each gene...")
    pas = pas.cluster(by=group_id_col)

    # Classify polyAsite_numbers 1..n for each gene
    # Where 1 = most 5', n = most 3' regardless of strand
    # adds 'pas_number' column
    pas = cluster_to_region_number(pas, group_id_col)

    pas = pas[["gene_id", "transcript_id", "pas_number"]]

    #Also perform a last exon classification
    # Principle the same (i.e 1..n in 5'-3' order in gene)
    # If last exons have any overlap then consider them the same le_number
    eprint("Assigning 5'-3' last exon number (le_number) for transcripts of each gene...")
    le = le.cluster(by=group_id_col)

    le = cluster_to_region_number(le, group_id_col, out_col="le_number")

    # Merge le with pas so that le_number and pas_number are in the same gr
    eprint("Merging pas_number & le_number for each transcript into common dataframe...")

    pas_cols = pas.columns.tolist()

    le = le.apply_pair(pas, lambda df, df_to_merge: _pd_merge_gr(df,
                                                                 df_to_merge,how="left",
                                                                 on=tx_id_col,
                                                                 suffixes=[None, "_match"],
                                                                 to_merge_cols=pas_cols)
                       )


    eprint("assigning 'pas_id' for each gene...")
    le = le.assign("pas_id",
                   lambda df: df[group_id_col] + "_" + df["pas_number"].astype(int).astype(str))

    tx_to_pas = le.as_df()[[group_id_col, tx_id_col, "pas_number", "le_number", "pas_id"]]



    eprint(f"Writing PAS and LE assignment dataframe to TSV... - {output_prefix + '.pas_assignment.tsv'}")

    tx_to_pas.sort_values(by=[group_id_col,
                              "pas_number"]).to_csv(output_prefix + ".pas_assignment.tsv",
                                                    sep="\t",
                                                    index=False,
                                                    header=True)

    eprint(f"Writing 'tx2pas' (transcript_id | pas_id) to TSV... - {output_prefix + '.tx2pas.tsv'}")
    tx_to_pas[[tx_id_col,
               "pas_id"]].drop_duplicates().to_csv(output_prefix + ".tx2pas.tsv",
                                                   sep="\t",
                                                   index=False,
                                                   header=True)


if __name__ == '__main__':

    start = timer()

    descrpn = """Script to assign transcripts of a gene according to the polyA site which they use/share"""

    parser = argparse.ArgumentParser(description=descrpn,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter, # Add defaults to end of help strings
                                     )

    parser.add_argument("-i",
                        "--input-gtf",
                        dest="input_gtf", type=str,
                        default=argparse.SUPPRESS,
                        required=True,
                        help="Path to input GTF file containing transcripts to group by used poly(A) site")

    parser.add_argument("-w",
                        "--merge-window-size",
                        dest="merge_window_size",
                        type=int,
                        default=25,
                        help="Size of window (nt) centred on polyA site to merge/consider overlapping poly(A) sites as a a single cleavage event. MUST BE ODD (PAS is 1nt interval, extend equally upstream & downstream)")

    parser.add_argument("-g","--group-attribute-key",
                        dest="group_id_key",
                        type=str,
                        default="gene_id",
                        help="Name of GTF attribute key that groups together transcripts. PolyA sites will be assigned according to values in this key")

    parser.add_argument("-t","--transcript-attribute-key",
                        dest="tx_id_key",
                        type=str,
                        default="transcript_id",
                        help="Name of GTF attribute key that defines transcripts. Used to return polyA site assignments for these features")

    parser.add_argument("-o","--output-prefix",
                        dest="output_prefix",
                        type=str,
                        default=argparse.SUPPRESS,
                        required=True,
                        help="path to/prefix for output files. '.pas_assignment.tsv' is suffixed for TSV file storing polyA site assignment per transcript/gene (tx_id | gene_id| pas_number | pas_id)")



    if len(sys.argv) == 1:
        parser.print_help()
        parser.exit()


    args = parser.parse_args()

    if (args.merge_window_size % 2) == 0:
        raise ValueError(f"-w/--merge-window-size value must be an odd number - you provided {args.merge_window_size}. This is so can extend equally in either direction of the poly(A) site (1nt length interval)")

    main(args.input_gtf, args.merge_window_size, args.group_id_key, args.tx_id_key, args.output_prefix)

    end = timer()

    eprint(f"Complete: took {round(end - start, 3)} s / {round((end - start) / 60, 3)} min (3 dp)")
