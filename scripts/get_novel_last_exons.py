#!/usr/bin/env python3

from __future__ import print_function
import pyranges as pr
import numpy as np
import pandas as pd
from papa_helpers import eprint, add_region_number, get_terminal_regions
from pyranges.readers import read_gtf_restricted
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
- GTF of predicted novel transcripts
- GTF of reference transcripts

Outputs:
- GTF of candidate novel last exons from input GTF
    - Additional attributes:
        - sample_name - input GTF name
        - isoform_class - Event classification code (e.g. first_exon_extension,internal_intron_spliced etc.)
        - ref_gene_id - gene_id of matching reference gene(s) (separated by ';' if multiple)
        - ref_transcript_id - transcript_ids of matching transcripts (separated by ';' if multiple)
- 'Match status' TSV storing sample, match status, isoform class information for each transcript
- 'general stats' TSV showing summary counts of number of events of each class
'''


# Che

# Pyranges default GTF columns are named as below

#https://github.com/biocore-ntnu/pyranges/blob/1ee215c645f7dbec3282555fcd0ccec610236614/pyranges/out.py#L44
pyranges_gtf_cols = "Chromosome   Source   Feature    Start     End       Score    Strand   Frame".split()

# Minimal cols needed from GTF for processing (no need to carry loads of cols)
processing_gtf_cols_en = pyranges_gtf_cols + ["gene_id",
                                              "transcript_id",
                                              "exon_number"]

processing_gtf_cols_n_en = pyranges_gtf_cols + ["gene_id",
                                                "transcript_id"]


def _n_ids(gr, id_col):

    assert id_col in gr.columns

    return len(set(gr.as_df()[id_col]))


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

    return pd.Series(decisions, index=df.index)


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

    gr = gr.assign(out_col,
                   lambda df: _df_add_region_rank(df,
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

    if any([exons_choices["extract_regions"],
            exons_choices["extract_last_regions"]
            ]):

        exons = gr.subset(lambda df: df["Feature"] == "exon")

        if any([exons_choices["add_region_number"],
                exons_choices["add_region_rank"],
                exons_choices["extract_last_regions"]
                ]):

            # Other two depend on region_number so needs to be added
            exons = add_region_number(exons,
                                      feature_key="exon",
                                      out_col="exon_number")

            if exons_choices["add_region_rank"]:
                exons = add_region_rank(exons)

            if exons_choices["extract_last_regions"]:
                last_exons = get_terminal_regions(exons)

            else:
                last_exons = None

        if not exons_choices["extract_regions"]:
            # Return last_exons only
            exons = None

    else:
        exons = None
        last_exons = None


    if any([introns_choices["extract_regions"],
            introns_choices["extract_last_regions"]
            ]):

        introns = gr.features.introns(by="transcript")

        if any([introns_choices["add_region_number"],
                introns_choices["add_region_rank"],
                introns_choices["extract_last_regions"]
                ]):

            # Other two depend on region_number so this must always be added
            introns = add_region_number(introns,
                                        feature_key="intron",
                                        out_col="intron_number")

            if introns_choices["add_region_rank"]:
                introns = add_region_rank(introns,
                                          region_number_col="intron_number")

            if introns_choices["extract_last_regions"]:
                last_introns = get_terminal_regions(introns,
                                                    feature_key="intron",
                                                    region_number_col="intron_number")

            else:
                last_introns = None

        if not introns_choices["extract_regions"]:
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


def add_3p_extension_length(gr,
                            ref_gr,
                            id_col="transcript_id",
                            out_col="3p_extension_length",
                            suffix="_b",
                            ref_cols_to_keep=["gene_id",
                                              "transcript_id",
                                              "exon_id",
                                              "Start",
                                              "End",
                                              "region_rank"],
                            nb_cpu=1):
    '''
    Add column '3p_extension_length' reporting 3'end extension of overlapping regions in gr & ref_gr (distance relative to gr)
    Note that for each unique ID (id_col) in gr, the smallest extension will be reported
    Avoids cases where overlaps with short & long isoforms, but tx is just reassembly of long isoform
    (so it's just reassembly of longer isoform, but is an extension relative to shorter ref isoform)
    '''

    # Find columns unique to ref_gr, so can drop at the end (no suffix added)
    not_cols = gr.columns.tolist()



    # Set of reference cols want to drop at end (either suffixed or not)
    # Note that Chromosome is never copied + suffixed in a join
    ref_cols = list(set(ref_gr.columns) - set(["Chromosome"]) - set(ref_cols_to_keep))

    ref_to_drop = [col if col not in not_cols else col + suffix for col in ref_cols]

    # Only overlapping intervals kept
    joined = gr.join(ref_gr,
                     strandedness="same",
                     how=None,
                     suffix=suffix,
                     nb_cpu=nb_cpu)

    joined = joined.assign(out_col,
                           lambda df: df["End"] - df["End_b"] if (df["Strand"] == "+").all() else
                           df["Start_b"] - df["Start"],
                           nb_cpu=nb_cpu)


    # To avoid capturing extensions of shorter isoforms (that really just are the longer known isoform)
    # Pick the smallest extension for each transcripts
    joined = joined.apply(lambda df: df.sort_values([id_col, out_col],
                                                    ascending=True).drop_duplicates(subset=[id_col],
                                                                                    keep="first"),
                          nb_cpu=nb_cpu)


    # eprint(joined.columns)

    return joined.drop(ref_to_drop)


def _df_5p_end_tolerance(df, rank_col, first_key, first_5p_tolerance, other_5p_tolerance):
    '''
    '''

    if (df["Strand"] == "+").all():
        # 5'end = Start
        # conditions = [(df[rank_col] == first_key) & (abs(df["Start"] - df["Start_b"]) <= first_5p_tolerance),
        #               (df[rank_col] != first_key) & (abs(df["Start"] - df["Start_b"]) <= other_5p_tolerance)
        #               ]
        # choices = [True, True]
        decisions = np.where((df[rank_col] == first_key) & (abs(df["Start"] - df["Start_b"]) <= first_5p_tolerance) |
                             (df[rank_col] != first_key) & (abs(df["Start"] - df["Start_b"]) <= other_5p_tolerance),
                             True,
                             False
                             )

    elif (df["Strand"] == "-").all():
        # 5'end = End
        decisions = np.where((df[rank_col] == first_key) & (abs(df["End"] - df["End_b"]) <= first_5p_tolerance) |
                             (df[rank_col] != first_key) & (abs(df["End"] - df["End_b"]) <= other_5p_tolerance),
                             True,
                             False
                             )


    return pd.Series(decisions, index=df.index)


def _df_add_event_type(df, id_col, rank_col, rkey2key, collapse_by_id=True):
    '''
    '''

    assert isinstance(rkey2key, dict)

    # eprint(df[rank_col].drop_duplicates())
    # eprint(rkey2key)

    assert any([True if key in set(df[rank_col]) else False for key in rkey2key.keys()])

    conditions = [df[rank_col] == key for key in rkey2key.keys()]

    choices = list(rkey2key.values())

    decisions = np.select(conditions, choices)

    if not collapse_by_id:
        return pd.Series(decisions, index=df.index)

    else:
        # Some IDs could have multiple assignments based on overlapping exon
        # Collapse to comma-separated string if multiple

        if not df.empty:
            to_join = df.loc[:, [id_col]]
            to_join["event_type_decision"] = pd.Series(decisions, index=df.index)

            to_join = (to_join.groupby(id_col)
                       ["event_type_decision"]
                       .agg(lambda x: ",".join(list(set(x))))
                       .reset_index() #return tx_id to column
                       )

            # Now join by tx_id to get collapsed ID in correct (original) index
            df2 = df.merge(to_join, on=id_col)

            return df2["event_type_decision"]

        else:
            return pd.Series(decisions, index=df.index)



def find_extension_events(novel_le,
                          ref_exons,
                          min_extension_length,
                          first_5p_tolerance,
                          other_5p_tolerance,
                          id_col="transcript_id",
                          rank_col="region_rank",
                          event_type_outcol="event_type",
                          first_key="first_exon_extension",
                          internal_key="internal_exon_extension",
                          last_key="last_exon_extension"):
    '''
    Return gr with last exons that extend a reference exon
    Criteria:
    1 - Overlap with reference exon, novel last exon 3'end is downstream of ref_exon 3'end
    2 - Smallest extension length is >= min_extension_length
        - Smallest to avoid reassembly of longer isofrm being called as extension of shorter isoform
    3 - 5'ends of novel_le & ref_exons match within given tolerance (nt)
        - first exons should be more lenient (imprecision of TSS annotation, reassembly not nt precise)
    '''

    assert isinstance(min_extension_length, int)
    assert isinstance(first_5p_tolerance, int)
    assert isinstance(other_5p_tolerance, int)

    # Find events extending 3'ends of ref exon, add 3p_extension_length (nt) column
    # 1 length per le returned (smallest)
    novel_le_ext = add_3p_extension_length(novel_le, ref_exons)
    eprint(f"Number of events with any 3' ref exon extension - {_n_ids(novel_le_ext, id_col)}")

    ext_len_dist = (novel_le_ext.as_df()
                    ["3p_extension_length"]
                    .describe(percentiles=[i * 0.1 for i in range(1,11,1)])
                    )

    eprint(f"3' extension length distribution\n{ext_len_dist}")

    # Subset for extensions of min length
    novel_le_ext = novel_le_ext.subset(lambda df: df["3p_extension_length"] >= min_extension_length)

    eprint(f"After minimum length filter - {min_extension_length} - number of extension events - {_n_ids(novel_le_ext, id_col)}")

    # Check 5'ends overlap within given tolerance
    novel_le_ext = novel_le_ext.subset(lambda df: _df_5p_end_tolerance(df,
                                                                       rank_col,
                                                                       first_key,
                                                                       first_5p_tolerance,
                                                                       other_5p_tolerance
                                                                       )
                                       )

    eprint(f"After 5'end match tolerance filter, number of events - {_n_ids(novel_le_ext, id_col)}")


    # assign a 'event_type' column based on overlapping exon 'rank'
    novel_le_ext = novel_le_ext.assign(event_type_outcol,
                                       lambda df: _df_add_event_type(df,
                                                                     id_col,
                                                                     rank_col,
                                                                     rkey2key={"first": first_key,
                                                                               "internal": internal_key,
                                                                               "last": last_key}
                                                                     )
                                       )

    ev_types = novel_le_ext.as_df()[[id_col,
                                      event_type_outcol]
                                    ].drop_duplicates().value_counts(subset=[event_type_outcol])

    eprint(f"Number of events of each type- {ev_types}")

    return novel_le_ext


def find_spliced_events(novel_li,
                        ref_introns,
                        suffix="_b",
                        ref_cols_to_keep=["gene_id",
                                          "transcript_id",
                                          "exon_id",
                                          "Start",
                                          "End",
                                          "region_rank"],
                        id_col="transcript_id",
                        rank_col="region_rank",
                        event_type_outcol="event_type",
                        first_key="first_exon_spliced",
                        internal_key="internal_exon_spliced",
                        last_key="last_exon_spliced"):
    '''
    '''

    # Find columns unique to ref_introns, so can drop at the end (no suffix added)
    not_cols = ref_introns.columns.tolist()

    # Set of reference cols want to drop at end (either suffixed or not)
    # Note that Chromosome is never copied + suffixed in a join
    ref_cols = list(set(ref_introns.columns.tolist()) - set(ref_cols_to_keep) - set(["Chromosome"]))
    ref_to_drop = [col if col not in not_cols else col + suffix for col in ref_cols]

    # Only novel last SJs overlapping ref SJs kept
    novel_spliced = novel_li.join(ref_introns, suffix=suffix)

    eprint(f"Number of putative novel spliced events - {_n_ids(novel_spliced, id_col)}")
    # eprint(f"ref exon ranks\n {novel_spliced.as_df()[rank_col].drop_duplicates()}")

    # Subset for exact matches at 5'end
    novel_spliced = novel_spliced.subset(lambda df: (df["Strand"] == "+") & (df["Start"] == df["Start_b"]) |
                                                    (df["Strand"] == "-") & (df["End"] == df["End_b"])
                                         )

    eprint(f"After filtering for exact reference 5'ss matches, number of novel spliced events - {_n_ids(novel_spliced, id_col)}")
    # eprint(f"exact matching exon_ranks\n {novel_spliced.as_df()[rank_col].drop_duplicates()}")

    # Add event_type col based on overlapping exon 'rank'
    novel_spliced = novel_spliced.assign(event_type_outcol,
                                         lambda df: _df_add_event_type(df,
                                                                       id_col,
                                                                       rank_col,
                                                                       rkey2key={"first": first_key,
                                                                                 "internal": internal_key,
                                                                                 "last": last_key}
                                                                       )
                                         )


    ev_types = novel_spliced.as_df()[[id_col,
                                      event_type_outcol]
                                     ].drop_duplicates().value_counts(subset=[event_type_outcol])

    eprint(f"Number of events of each type - {ev_types}")

    return novel_spliced.drop(ref_to_drop)


def _filter_gr_for_tx(df,df2):

    if df.empty:
        return df

    elif df2.empty:
        return df

    return df[~df["transcript_id"].isin(set(df2["transcript_id"]))]


def main(input_gtf_path,
         ref_gtf_path,
         min_extension_length,
         first_exon_5p_tolerance,
         other_exon_5p_tolerance,
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
    # ref_gtf = pr.read_gtf(ref_gtf_path)
    ref_gtf = read_gtf_restricted(ref_gtf_path, skiprows=None)
    end = timer()

    eprint(f"Complete - took {end - start} s")

    # eprint(ref_gtf.as_df().info(memory_usage="deep"))

    # Ref GTFs often have a load of columns (attributes) don't need
    # Subset to essential to cut down on memory
    if trust_exon_number_ref:
        try:
            ref_gtf = ref_gtf[processing_gtf_cols_en]
        except KeyError:
            # using gtf_restricted (No source or Score)
            ref_gtf = ref_gtf[list(set(processing_gtf_cols_en) - set(["Source",
                                                                      "Score",
                                                                      "Frame"]
                                                                     )
                                   )
                              ]
    else:
        try:
            ref_gtf = ref_gtf[processing_gtf_cols_n_en]
        except KeyError:
            # using gtf_restricted (no Source or Score)
            ref_gtf = ref_gtf[list(set(processing_gtf_cols_n_en) - set(["Source",
                                                                        "Score",
                                                                        "Frame"]
                                                                       )
                                   )
                              ]

    # eprint(ref_gtf.as_df().info(memory_usage="deep"))

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

    # eprint(ref_exons.as_df().info(memory_usage="deep"))
    # eprint(ref_introns.as_df().info(memory_usage="deep"))

    # eprint(ref_exons.region_rank.drop_duplicates())
    # eprint(ref_introns.region_rank.drop_duplicates())

    eprint("Reading in input GTF...")

    start = timer()
    novel_gtf = pr.read_gtf(input_gtf_path)

    end = timer()

    eprint(f"Complete - took {end - start} s")

    # eprint(novel_gtf.as_df().info(memory_usage="deep"))

    eprint(" ".join(["Extracting last exons & last introns",
                     "from input GTF & numbering by",
                     "5'-3' order along the transcript"]))

    start = timer()
    _, novel_le, _, novel_li = extract_format_exons_introns(novel_gtf,
                                                            exons_choices={"extract_regions": False,
                                                                           "add_region_number": trust_exon_number_input,
                                                                           "add_region_rank": False,
                                                                           "extract_last_regions": True},
                                                            introns_choices={"extract_regions": False,
                                                                             "add_region_number": True,
                                                                             "add_region_rank": False,
                                                                             "extract_last_regions": True}
                                                            )

    end = timer()

    eprint(f"Complete - took {end - start} s")

    # eprint(novel_le.as_df().info(memory_usage="deep"))
    # eprint(novel_li.as_df().info(memory_usage="deep"))

    # Remove tx if 3'ends of last exons overlap with any reference exon
    eprint("Getting last exons with 3'ends not overlapping any known exons")

    start = timer()
    no_exon_overlap_3p_ids = no_overlap_3p_ids(novel_le, ref_exons)
    end = timer()

    eprint(f"Complete - took {end - start} s")

    eprint(f"Number of putative novel 3'ends - {len(no_exon_overlap_3p_ids)}")

    novel_le = novel_le.subset(lambda df: df["transcript_id"].isin(no_exon_overlap_3p_ids))
    novel_li = novel_li.subset(lambda df: df["transcript_id"].isin(no_exon_overlap_3p_ids))


    # Identify extension events
    # Overlap with known exons, 3'end is downstream of overlapping exon end
    # Retain extensions above a minimum length threshold
    # Note: smallest extension for each novel event selected as representative
    ## Avoids e.g. reassembly of long 3'UTR isoform being called as extension of proximal isoform
    eprint("Finding novel extension events...")

    start = timer()
    extensions = find_extension_events(novel_le,
                                       ref_exons,
                                       min_extension_length,
                                       first_exon_5p_tolerance,
                                       other_exon_5p_tolerance)
    end = timer()

    eprint(f"Complete - took {end - start} s")

    # eprint(f"extensions cols \n{extensions.columns}")

    # Identify spliced events
    # First filter novel events for non-extensions
    eprint("Finding novel spliced events - filtering for non-extensions...")

    start = timer()

    # TODO: is apply_pair any quicker than just making a big set & pr.subset?
    novel_li = novel_li.apply_pair(extensions,
                                   lambda df, df2:
                                   _filter_gr_for_tx(df, df2),
                                   as_pyranges=True)
    end = timer()

    eprint(f"Complete - took {end - start} s")

    # Valid events overlap reference intron/junction and exactly share a 5'ss
    eprint("Finding novel spliced events...")

    start = timer()
    spliced = find_spliced_events(novel_li, ref_introns)
    end = timer()

    eprint(f"Complete - took {end - start} s")

    # eprint(f"spliced columns {spliced.columns}")

    combined = pr.concat([extensions, spliced])

    # Output
    combined.to_gtf(output_prefix + "_last_exons.gtf")



if __name__ == '__main__':

    start = timer()

    descrpn = """Extract putative novel last exons from a GTF file of predicted transcripts"""

    parser = argparse.ArgumentParser(description=descrpn,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter, # Add defaults to end of help strings
                                     )

    parser.add_argument("-i",
                        "--input-gtf",
                        required=True,
                        type=str,
                        default=argparse.SUPPRESS,
                        help="Path to input GTF file of predicted transcripts (from which putative novel last exons will be identified)")

    parser.add_argument("-r",
                        "--reference-gtf",
                        required=True,
                        type=str,
                        default=argparse.SUPPRESS,
                        help="Path to GTF file containing reference transcripts")

    parser.add_argument("-m",
                        "--min-extension-length",
                        type=int,
                        default=100,
                        help="Minimum length (nt) of a putative extension last exon event for it to be retained")

    parser.add_argument("-f",
                        "--first-exon-5p-tolerance",
                        type=int,
                        default=100,
                        help="Tolerance window (nt) to match 5'ends of reference first exons & first exon extension events")

    parser.add_argument("-t",
                        "--other-exon-5p-tolerance",
                        type=int,
                        default=0,
                        help="Tolerance window (nt) to match 5'ends of reference exons & overlapping putative last exon")

    parser.add_argument("-o",
                        "--output-prefix",
                        type=str,
                        default="putative_novel",
                        help="Name of prefix for output files - GTF  (suffixed with '_last_exons.gtf'), 'match stats' ('_match_stats.tsv')")

    parser.add_argument("--trust-input-exon-number",
                        action="store_true",
                        default=False,
                        help="Whether to 'trust' the exon number attribute in input GTF as being strand-aware (not the case for StringTie)")

    parser.add_argument("--trust-ref-exon-number",
                        action="store_true",
                        default=False,
                        help="Whether to 'trust' the exon number attribute in reference GTF as being strand-aware")

    if len(sys.argv) == 1:
        parser.print_help()
        parser.exit()

    args = parser.parse_args()

    main(args.input_gtf,
         args.reference_gtf,
         args.min_extension_length,
         args.first_exon_5p_tolerance,
         args.other_exon_5p_tolerance,
         args.trust_ref_exon_number,
         args.trust_input_exon_number,
         args.output_prefix)

    end = timer()

    eprint(f"Script complete: took {round(end - start, 3)} s / {round((end - start) / 60, 3)} min (3 dp)")
