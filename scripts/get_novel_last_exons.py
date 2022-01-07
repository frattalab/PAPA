#!/usr/bin/env python3

from __future__ import print_function
import pyranges as pr
import numpy as np
import pandas as pd
from papa_helpers import eprint, add_region_number, get_terminal_regions
from timeit import default_timer as timer
import argparse
import sys

'''
Script to extract candidate novel last exons from a (filtered) assembled transcriptome

Novel last exons (spliced, bleedthrough) are extracted with the following criteria:
- 3'end does not overlap with any reference exon
    - 3'ends very unlikely to be called correctly from overlapping internal exon coverage profile)
    - 3'ends overlapping/inside reference last exons are not extensions
- Extension events are above a minimum length threshold (default: 100nt)
    - Internal and last exon extension events exactly match overlapping exon at 5'end (i.e. identical 3'ss)
    - First exon extension events match the reference 5'end within a given window tolerance (precise TSS difficult to call from coverage)
    - Annotate as first_exon_extension, internal_exon_extension, last_exon_extension
- 5'ss of putative spliced events exactly match a reference 5'ss
    - Annotate as first_intron_spliced, internal_intron_spliced or last_intron_spliced
    - Sub-classify last_intron_spliced (downstream, proximal ALE, novel_3ss etc.)

Takes as input:
- GTF(s) of predicted reference transcripts
- GTF of reference transcripts

Outputs:
- GTF of combined candidate novel last exons from all input GTFs
    - Additional attributes:
        - sample_name - input GTF name
        - last_exon_ID - unique identifier shared between events with the same 5'ss
        - papa_class - Event classification code (e.g. first_exon_extension,internal_intron_spliced etc.)
        - ref_gene_id - gene_id of matching reference gene(s) (separated by ';' if multiple)
        - ref_transcript_id - transcript_ids of matching transcripts (separated by ';' if multiple)
- 'Match status' TSV storing sample, match status, isoform class information for each transcript
- 'general stats' TSV showing number of events of each class/ after each filtering step for each sample
'''

# Che

# Pyranges default GTF columns are named as below

#https://github.com/biocore-ntnu/pyranges/blob/1ee215c645f7dbec3282555fcd0ccec610236614/pyranges/out.py#L44
pyranges_gtf_cols = "Chromosome   Source   Feature    Start     End       Score    Strand   Frame".split()

# Minimal cols needed from GTF for processing (no need to carry loads of cols)
processing_gtf_cols_en = pyranges_gtf_cols + ["gene_id",
                                              "gene_name",
                                              "transcript_id",
                                              "exon_number"]

processing_gtf_cols_n_en = pyranges_gtf_cols + ["gene_id",
                                                "gene_name",
                                                "transcript_id"]


def _df_add_region_rank(df,
                        id_col,
                        region_number_col,
                        first_key,
                        internal_key,
                        last_key):
    '''
    '''

    conditions = [df[region_number_col] == 1,
    # keep="last" sets last in ID to 'False' and all others true (negate to keep last only)
    # Safe as gr is sorted by tx_id and region_number_col prior
                  ~df.duplicated(subset=[id_col], keep="last")]

    choices = [first_key, last_key]

    decisions = np.select(conditions, choices, default=internal_key)

    return pd.Series(decisions)


def add_region_rank(gr,
                    id_col="transcript_id",
                    region_number_col="exon_number",
                    out_col="region_rank",
                    first_key="first",
                    internal_key="internal",
                    last_key="last"):
    '''
    Add a column specifying whether region corresponds to the most 3' region in the group (e.g. transcript) or not (1/0)
    '''

    # Sort regions by id_col & 5'-3' position (region_number_col must be strand_aware)

    gr = gr.apply(lambda df: df.sort_values(by=[id_col,
                                                region_number_col],
                                            ascending=True),
                  nb_cpu=1)

    # keep="last" sets last row by ID (last exon) to False and all others True
    gr = gr.assign(out_col,
                   lambda df: _df_annotate_region_rank(df,
                                                       id_col,
                                                       region_number_col,
                                                       first_key,
                                                       internal_key,
                                                       last_key)
                   )

    return gr


def extract_format_exons_introns(gr,
                                 exons_choices={"extract_regions": True,
                                                "add_region_number": True,
                                                "add_region_rank": True,
                                                "extract_last_regions": False},
                                 introns_choices={"extract_regions": True,
                                                  "add_region_number": True,
                                                  "add_region_rank": True,
                                                  "extract_last_regions": False},

                                                 ):
    '''
    Wrapper function to extract exons, introns from gr
    and add region_number/region_rank columns
    Returns 4 length tuple of (exons, last_exons, introns, last_introns)
    If choices were 'False' (i.e. do not create this object), None is returned in its place

    gr: target PyRanges object, ideally obtained for pr.read_gtf(). Must contain 'exon' & 'transcript' in 'Feature' column
    '''

    if exons["extract_regions"] or exons["extract_last_regions"]:

        exons = gr.subset(lambda df: df["Feature"] == "exon")

        if exons["add_region_number"] or exons["add_region_rank"] or exons["extract_last_regions"]:

            exons = add_region_number(exons,
                                      feature_key="exon",
                                      out_col="exon_number")

            if exons["add_region_rank"] or exons["extract_last_regions"]:
                exons = add_region_rank(exons)

                if exons["extract_last_regions"]:
                    last_exons = get_terminal_regions(exons)

                else:
                    last_exons = None

            else:
                last_exons = None
        else:
            last_exons = None

    elif not exons["extract_regions"]:
        exons = None

    else:
        exons = None
        last_exons = None



    if introns["extract_regions"] or introns["extract_last_regions"]:

        introns = gr.features.introns(by="transcript")

        if introns["add_region_number"] or introns["add_region_rank"] or introns["extract_last_regions"]:

            introns = add_region_number(introns,
                                        feature_key="intron",
                                        out_col="intron_number")

            if introns["add_region_rank"] or introns["extract_last_regions"]:
                introns = add_region_rank(introns,
                                          region_number_col="intron_number")

                if introns["extract_last_regions"]:
                    last_introns = get_terminal_regions(introns,
                                                        feature_key="intron",
                                                        region_number_col="intron_number")

                else:
                    last_introns = None

            else:
                last_introns = None
        else:
            last_introns = None

    elif not introns["extract_regions"]:
        # Return last introns only
        introns = None

    else:
        #Return neither
        introns = None
        last_introns = None


    return exons, last_exons, introns, last_introns


def no_overlap_3p_ids(gr, overlaps_gr, id_col="transcript_id"):
    '''
    Return set of IDs (str) from gr for 3'ends of regions that do not overlap at all with regions in overlaps_gr
    '''

    gr_3p = gr.three_end()

    gr_3p = gr_3p.overlap(overlaps_gr, invert=True)

    return set(gr_3p.as_df()[id_col])


def main(input_gtf_list,
         ref_gtf_path,
         min_extension_length,
         first_exon_5p_tolerance,
         trust_exon_number_ref,
         trust_exon_number_input,
         output_prefix):
    '''
    '''

    # 1. Prepare reference GTF for analysis - need:
    # - Reference introns
    # - Reference exons
    # Add columns denoting region_number by tx in 5'-3' order (<exon/intron>_number)
    # Add columns denoting whether exon/intron is first, internal or last in tx (region_rank)

    eprint("Reading in reference GTF, this can take a while...")

    start = timer()
    ref_gtf = pr.read_gtf(ref_gtf_path)
    end = timer()

    eprint(f"Complete - took {end - start} s")

    eprint(ref_gtf.as_df().info(memory_usage="deep"))

    # Ref GTFs often have a load of columns (attributes) don't need
    # Subset to essential to cut down on memory
    if trust_exon_number_ref:
        ref_gtf = ref_gtf[processing_gtf_cols_en]

    else:
        ref_gtf = ref_gtf[processing_gtf_cols_n_en]

    eprint(ref_gtf.as_df().info(memory_usage="deep"))

    eprint(" ".join(["Extracting exons & introns from reference GTF,"
                     "numbering by 5'-3' order along the transcript &",
                     "classifying as 'first', 'internal or 'last' region..."]))

    start = timer()
    # Don't need first and last ref exons/introns as separate objects
    ref_exons, _, ref_introns, _ = extract_format_exons_introns(ref_gtf,
                                                                exons_choices={"extract_regions": True,
                                                                               "add_region_number": trust_exon_number_ref,
                                                                               "add_region_rank": True,
                                                                               "extract_last_regions": False}
                                                                )
    end = timer()

    eprint(f"Complete - took {end - start} s")

    eprint(ref_exons.as_df().info(memory_usage="deep"))
    eprint(ref_introns.as_df().info(memory_usage="deep"))

    for sample_gtf in input_gtf_list:
        eprint(f"Reading in GTF at {sample_gtf}")

        start = timer()
        novel_gtf = pr.read_gtf(sample_gtf)

        end = timer()

        eprint(f"Complete - took {end - start} s")

        eprint(novel_gtf.as_df().info(memory_usage="deep"))

        eprint(" ".join(["Extracting last exons & last introns",
                         "from input GTF & numbering by",
                         "5'-3' order along the transcript"]))

        start = timer()
        _, novel_le, _, novel_li = extract_format_exons_introns(novel_gtf,
                                                                exons_choices={"extract_regions": False,
                                                                               "add_region_number": trust_exon_number_input,
                                                                               "add_region_rank": True,
                                                                               "extract_last_regions": True},
                                                                introns_choices={"extract_regions": False,
                                                                                 "add_region_number": True,
                                                                                 "add_region_rank": True,
                                                                                 "extract_last_regions": True}
                                                                )

        end = timer()

        eprint(f"Complete - took {end - start} s")

        eprint(novel_le.as_df().info(memory_usage="deep"))
        eprint(novel_li.as_df().info(memory_usage="deep"))

        # Remove tx if 3'ends of last exons overlap with any reference exon
        eprint("Getting last exons with 3'ends not overlapping any known exons")

        start = timer()
        no_exon_overlap_3p_ids = no_overlap_3p_ids(novel_le, ref_exons)
        end = timer()

        eprint(f"Complete - took {end - start} s")

        novel_le = novel_le.subset(lambda df: df["transcript_id"].isin(no_exon_overlap_3p_ids))
        novel_li = novel_li.subset(lambda df: df["transcript_id"].isin(no_exon_overlap_3p_ids))


        # Identify extension events



if __name__ == '__main__':



    main()
