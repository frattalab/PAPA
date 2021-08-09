#! /usr/bin/env python3
from __future__ import print_function
import pyranges as pr
import numpy as np
import pandas as pd
import argparse
import os
import sys
import logging
from timeit import default_timer as timer

'''
Filter assembled transcripts for those that match reference transcripts in their intron chain up until their penultimate intron
This allows identification of novel last exons (conservatively), but limits addition of false positive transcripts
'''


def eprint(*args, **kwargs):
    '''
    Nice lightweight function to print to STDERR (saves typing, I'm lazy)
    Credit: https://stackoverflow.com/questions/5574702/how-to-print-to-stderr-in-python (MarcH)
    '''
    print(*args, file=sys.stderr, **kwargs)




def introns_from_df(df):
    '''
    '''

    n_exons = len(df)

    if n_exons < 2:
        return None
        #print(df)
        #raise Exception("at least two exons are required for transcript to have an intron")
    # n exons = n-1 introns

    strand = df["Strand"].drop_duplicates().tolist()[0]
#     print(strand)
    chrom = df["Chromosome"].drop_duplicates().tolist()[0]
    gene_id = df["gene_id"].drop_duplicates().tolist()[0]
    tx_id = df["transcript_id"].drop_duplicates().tolist()[0]
    feature = "intron"
    introns = {}
    for i in range(0, n_exons - 1):
        if strand == "+":
            intron_start = df.iloc[i, lambda x: x.columns.get_loc("End")]
            intron_end = df.iloc[i+1, lambda x: x.columns.get_loc("Start")]
            introns[str(i)] = {"Chromosome": chrom,
                               "Start": intron_start,
                               "End": intron_end,
                               "Strand": strand,
                               "Feature": feature,
                               "gene_id": gene_id,
                               "transcript_id": tx_id}
        elif strand == "-":
            intron_start = df.iloc[i, lambda x: x.columns.get_loc("End")]
            intron_end = df.iloc[i+1, lambda x: x.columns.get_loc("Start")]
            introns[str(i)] = {"Chromosome": chrom,
                               "Start": intron_start,
                               "End": intron_end,
                               "Strand": strand,
                               "Feature": feature,
                               "gene_id": gene_id,
                               "transcript_id": tx_id}
    return pd.DataFrame.from_dict(introns, orient = "index")



def introns_by_tx(gr, by="transcript_id", nb_cpu=1):
    '''
    '''
    # Sort by position (for safety)
    gr = gr.sort()

    return gr.apply(lambda df: df.groupby(by).apply(introns_from_df), nb_cpu=nb_cpu)


def rle(inarray):
        """
        run length encoding. Partial credit to R rle function.
        Multi datatype arrays catered for including non Numpy
        returns: tuple (runlengths, startpositions, values)
        https://stackoverflow.com/questions/1066758/find-length-of-sequences-of-identical-values-in-a-numpy-array-run-length-encoding
        Thomas Browne
        """
        ia = np.asarray(inarray)                # force numpy
        n = len(ia)
        if n == 0:
            return (None, None, None)
        else:
            y = ia[1:] != ia[:-1]               # pairwise unequal (string safe)
            i = np.append(np.where(y), n - 1)   # must include last element posi
            z = np.diff(np.append(-1, i))       # run lengths
            p = np.cumsum(np.append(0, z))[:-1] # positions
            return(z, p, ia[i])


def intron_id(gr):

    return gr.assign("intron_id",
                                         lambda df: pd.Series([":".join([tx_id, str(start), str(end)])
                                                     for tx_id, start, end in zip(df["transcript_id"],
                                                                                  df["Start"],
                                                                                  df["End"])
                                                    ])
                                        )

def sort_introns_by_strand(df):
    '''
    '''
    # first reset_index call removes the original index of the group (e.g. row 4005 in df)
    # second reset_index call adds the sorted index as a column to the dataframe (the order along exon in each transcript)
    if (df.Strand == '+').all():
        return df.sort_values(by=['End']).reset_index(drop=True).reset_index()
    elif (df.Strand == '-').all():
        return df.sort_values(by=['Start'], ascending=False).reset_index(drop=True).reset_index()



def validate_matching_chain(df, max_terminal_non_match=2):
    '''
    apply to grouped df
    '''

    runs, starts, vals = rle(df["match"])

    # Valid matches are:
    # All introns match (e.g. bleedthrough event)
    # All but last x introns match (usually one) (i.e. runs = 1,0)

    if np.all(vals):
        # All introns match (e.g. bleedthrough event)
        return True

    elif np.array_equal(vals, [1,0]) and runs[-1] <= max_terminal_non_match:
        # all but last x introns match (x = max_terminal_non_match) (i.e. runs = 1,0)
        return True

    else:
        return False

def stie_groupby_last_exon(df, exon_n_col, which="last"):
    '''
    '''

    if (df["Strand"] == "+").all():
        if which == "last":
            return df[exon_n_col].idxmax()
        elif which == "first":
            return df[exon_n_col].idxmin()

    elif (df["Strand"] == "-").all():
        if which == "last":
            return df[exon_n_col].idxmin()
        elif which == "first":
            return df[exon_n_col].idxmax()


def filter_multi_exon(df, exon_n_col):
    '''
    Want transcripts with > 1 exon
    '''
    if df[exon_n_col].nunique() > 1:
        return True
    else:
        return False


def get_terminal_regions(gr,
                   feature_col = "Feature",
                   feature_key = "exon",
                   id_col = "transcript_id",
                   region_number_col = "exon_number",
                   source = None,
                   which_region="last",
                   filter_single = False,
                   nb_cpu = 1):
    '''
    Return gr of last exons for each transcript_id
    In process, region_number_col will be converted to type 'int'
    StringTie merged GTFs (or whatever tool single_steps/stringtie_longreads.smk is using)
    reports exon_number that DOES NOT RESPECT STRAND (from browsing in IGV)
    i.e. for minus-strand - largest exon_number for transcript corresponds to FIRST EXON, not last
    Annotated (i.e. Ensembl) reported exon_numbers DO RESPECT STRAND (i.e. max always = last exon)

    if Do respect strand, put source = None (default)
    if Don't respect strand, put source = "stringtie" (i.e. plus strand = max, minus strand = min)
    '''

    assert source in [None, "stringtie"]
    assert which_region in ["first", "last"]
    assert region_number_col in gr.columns.tolist()
    assert feature_col in gr.columns.tolist()
    assert id_col in gr.columns.tolist()

    # Make sure only 'exon' features are in the gr
    assert gr.as_df()[feature_col].drop_duplicates().tolist() == [feature_key], "only {} entries should be present in gr".format(feature_key)

    # Make sure region_number_col is int
    mod_gr = (gr.assign(region_number_col,
                      lambda df: df[region_number_col].astype(int),
                      nb_cpu = nb_cpu)
             )


    # Filter out single-exon transcripts
    if filter_single:
        eprint("Filtering for multi-exon transcripts...")
        eprint("Before: {}".format(len(set(mod_gr.as_df()[id_col].tolist()))))

        mod_gr = (mod_gr.apply(lambda df: (df.groupby(id_col)
                                       .filter(lambda x: filter_multi_exon(df, region_number_col))
                                      )
                           ,
                           nb_cpu=nb_cpu
                          )
                 )
        eprint("After: {}".format(len(set(mod_gr.as_df()[id_col].tolist()))))




    if source is None:
        # source = None means that 1 = first region of group regardless of strand
        # Pick last region entry by max region number for each transcript (id_col)
        # Pick first region entry by min region number for each transcript (id_col)

        if which_region == "last":
            out_gr = mod_gr.apply(lambda df: df.iloc[df.groupby(id_col)[region_number_col].idxmax(),], nb_cpu = nb_cpu)

        elif which_region == "first":
            out_gr = mod_gr.apply(lambda df: df.iloc[df.groupby(id_col)[region_number_col].idxmin(),], nb_cpu = nb_cpu)

    elif source == "stringtie":
        # Numbering Doesn't respect strand - pick min if Minus strand, max if plus strand
        out_gr = (mod_gr.apply(lambda df: df.iloc[(df.groupby(id_col)
                                                  .apply(lambda df: stie_groupby_last_exon(df, region_number_col, which_region)
                                                        )
                                                 ),],
                               nb_cpu = nb_cpu
                              )
                 )


    return out_gr


def filter_transcripts_by_chain(novel_exons, novel_introns, ref_exons, ref_introns, match_type = "transcript", max_terminal_non_match=2, nb_cpu = 1):
    '''
    '''

    # Only keep essential columns (or those with potentially useful info) to save on memory
    # TO DO: this should really happen in main

    novel_cols_to_keep = [col for col in ["Feature","transcript_id"] if col in novel_exons.columns.tolist()]
    ref_cols_to_keep = [col for col in ["Feature", "transcript_id", "gene_id", "gene_name"] if col in ref_exons.columns.tolist()]

    assert match_type in ["transcript", "any"], "match_type must be one of 'transcript' or 'any'. value passed - {}".format(str(match_type))

    #1. Find introns by transcript & give each intron a unique ID
    # print("finding introns...")
    # t1 = timer()

    #novel_introns = introns_by_tx(novel_exons, nb_cpu=nb_cpu).sort()
    #ref_introns = introns_by_tx(ref_exons, nb_cpu=nb_cpu).sort()

    # t2 = timer()

    # print("took {} (s)".format(t2 - t1))
    eprint("filtering transcripts by intron chain matching...")

    eprint("sorting all grs by position for safety...")
    t1 = timer()

    novel_exons = novel_exons.sort()
    novel_introns = novel_introns.sort()
    ref_exons = ref_exons.sort()
    ref_introns = ref_introns.sort()

    t2 = timer()
    eprint("took {} (s)".format(t2 - t1))

    eprint("adding intron_id column...")

    t3 = timer()
    novel_introns = intron_id(novel_introns)
    ref_introns = intron_id(ref_introns)
    t4 = timer()

    eprint("took {} s".format(t4 - t3))

    #2. Track number of introns in each novel transcript
    # novel_tx_intron_counts = (novel_introns.as_df()
                              # .groupby("transcript_id").size())


    # novel_introns, ref_introns

    # 3. Store intron_ids for each transcript, sorted by intron_number (where 1 = first intron regardless of strand) in a df/Series
    eprint("generating df of novel txipts sorted by intron number...")

    t5 = timer()

    novel_intron_ids_ordered = (novel_introns.as_df()
                                .groupby("transcript_id")
                                .apply(sort_introns_by_strand)
                                .reset_index(drop=True)
                                .rename({"index": "intron_number"}, axis="columns")
                               )
    novel_intron_ids_ordered["intron_number"] = novel_intron_ids_ordered["intron_number"].add(1)

    # df of txipt_id | intron_id | intron_number
    novel_intron_ids_ordered = novel_intron_ids_ordered.loc[:,["transcript_id","intron_id","intron_number"]]

    t6 = timer()
    eprint("took {} s".format(t6 - t5))
#     eprint(novel_intron_ids_ordered.dtypes)


    #4. Find novel introns with any overlap with reference introns
    # Inner join to add ref_rows to novel gr
    eprint("finding overlaps between novel and reference introns...")

    t7 = timer()

    # Have to convert Starts and Ends to np.int64 to prevent left-join error
    # https://github.com/biocore-ntnu/pyranges/issues/170
    # TO DO: re-report
    joined = pr.PyRanges(novel_introns.as_df(), int64=True).join(pr.PyRanges(ref_introns.as_df(), int64=True), strandedness="same", suffix ="_ref", nb_cpu=nb_cpu)

    t8 = timer()

    eprint("took {} s".format(t8 - t7))

    #5. Filter for overlaps that exactly match (or differ by given tolerance - for now not supported)
    eprint("filtering overlaps for exact matches...")

    t9 = timer()
    joined = joined.subset(lambda df: abs(df.Start - df.Start_ref) + abs(df.End - df.End_ref) <= 0, nb_cpu=nb_cpu)
    t10 = timer()

    eprint("took {} s".format(t10 - t9))

    # Minimal info needed on matches between novel and reference introns
    joined = joined.as_df()[["transcript_id","intron_id","transcript_id_ref","intron_id_ref"]]

#     eprint(joined.dtypes)

    #6. Join ordered novel introns with match info
    #7. Assign a simple tracker column 'match' of True (where intron is matched) and False (where intron is not matched)

    eprint("preparing for filtering intron matches...")
    t11 = timer()

    if match_type == "any":
        # Looking for intron to match any annotated intron, regardless of reference transcript
        novel_ref_match_info = novel_intron_ids_ordered.merge(joined,
                                                              how="left",
                                                              on="intron_id",
                                                              suffixes=["_novel","_match"]
                                                             )

        # Assign 'match' column for each intron.
        # Since we don't really care which intron it matches, & no matches will mean NaN
        novel_ref_match_info["match"] = novel_ref_match_info["transcript_id_ref"]
        novel_ref_match_info["match"] = novel_ref_match_info["match"].fillna(0)
        novel_ref_match_info["match"] = novel_ref_match_info["match"].replace("\w*", 1, regex=True)

        novel_ref_match_info = novel_ref_match_info.drop_duplicates(subset=["intron_id"])

        # Minimal informative info is novel tx, novel intron_id & number, match column
        novel_ref_match_info = novel_ref_match_info[["transcript_id_novel","intron_id","intron_number","match"]]


    elif match_type == "transcript":
        # Looking for introns (except last) to match the same reference transcript
        # merge_ordered can do a 'grouped merge' filling in empty rows (introns) for each transcript_id
        # This is especially useful if want to do transcript-specific intron matching
        # For each reference transcript, all novel introns will be filled with NaN if no overlap for given transcript_id
        # (i.e. novel txipt matches all but last intron of reference transcript)

        novel_ref_match_info = (pd.merge_ordered(novel_intron_ids_ordered,
                                    joined,
                                    how="left",
                                    on="intron_id",
                                    right_by="transcript_id_ref", # group matches by ref tx & join tx by tx
                                    suffixes=["_novel","_match"],
                                    fill_method=None)
                   .sort_values(by=["transcript_id_novel","intron_number"])
                               )

        # merge_ordered fills rows for each intron for each ref tx in df, regardless of whether any overlap
        # .dropna(axis="rows", subset=["intron_id_ref"])
        novel_ref_match_info = (novel_ref_match_info.groupby(["transcript_id_novel", "transcript_id_ref"])
                                .filter(lambda df: (df["intron_id_ref"].notna()).any()) # Retained if ref tx has >=1 matching introns
                                .reset_index(drop=True))

        # Make a match column where 1 = match, 0 = no match for each ref id and novel intron
        novel_ref_match_info["match"] = novel_ref_match_info["intron_id_ref"]
        novel_ref_match_info["match"] = novel_ref_match_info["match"].fillna(0)
        novel_ref_match_info["match"] = novel_ref_match_info["match"].replace("\w*", 1, regex=True)

    t12 = timer()
    eprint("took {} s".format(t12 - t11))


    # 8. Filter down matching transcripts to those that all ref introns except penultimate or all introns...
    eprint("filtering for valid intron chain matches...")
    t13 = timer()
    if match_type == "any":
        # Only need to check by novel transcript_id
        filt_novel_ref_match_info = (novel_ref_match_info.groupby("transcript_id_novel")
                                     .filter(lambda x: validate_matching_chain(x, max_terminal_non_match)
                                            )
                                    )

    elif match_type == "transcript":
        # Check novel tx vs each ref tx
        filt_novel_ref_match_info = (novel_ref_match_info.groupby(["transcript_id_novel","transcript_id_ref"])
                                     .filter(lambda x: validate_matching_chain(x, max_terminal_non_match)
                                            )
                                    )
    t14 = timer()
    eprint("took {} s".format(t14 - t13))


    # Return simplified df of novel transcript_id & matching transcript_id if applicable

    if match_type == "any":
        return filt_novel_ref_match_info["transcript_id_novel"].drop_duplicates()

    elif match_type == "transcript":
        return filt_novel_ref_match_info[["transcript_id_novel", "transcript_id_ref"]].drop_duplicates()


def _join_grs(gr_left, gr_right, strandedness=None, how=None, report_overlap=False, slack=0, suffix="_b", nb_cpu=1):
    '''
    '''

    return gr_left.join(gr_right, strandedness, how, report_overlap, slack, suffix, nb_cpu)


def check_three_end(gr):
    '''
    '''
    df = gr.as_df()

    return all(df.Start == df.End)


def filter_first_intron_tx(novel_exons, ref_exons, ref_introns, chain_match_info, nb_cpu=1):
    '''
    Function to return novel last exon transcripts occurring within annotated first exons,
    In which the 3'end boundary of the novel first exon exactly matches an annotated first exon
    These transcript types cannot have valid intron chain match but can still produce valid last exon isoforms
    A well known e.g. of this isoform type is TDP-43 depletion sensitive STMN2 isoform
    '''
    # start = timer()
    eprint("finding novel last exon isoforms contained within annotated first introns...")

    #1 - Pull out non-chain matched novel isoforms
    eprint("extracting non-chain matched novel isoforms...")
    s1 = timer()

    if isinstance(chain_match_info, pd.Series):
        # match type was any
        novel_exons_nm = novel_exons.subset(lambda df: ~df["transcript_id"].isin(set(chain_match_info.tolist())), nb_cpu=nb_cpu)

    elif isinstance(chain_match_info, pd.DataFrame):
        # match_by/match_type was transcript
        novel_exons_nm = novel_exons.subset(lambda df: ~df["transcript_id"].isin(set(chain_match_info["transcript_id_novel"].tolist())), nb_cpu=nb_cpu)

    e1 = timer()
    eprint("took {} s".format(e1 - s1))


    #2 -Extract last exons of non-chain matched novel isoforms
    eprint("extracting last exons of non-chain matched novel isoforms")
    s2 = timer()

    novel_nm_last_exons = get_terminal_regions(novel_exons_nm,
                                               region_number_col="exon_number",
                                               source="stringtie",
                                               filter_single=True,
                                               which_region="last",
                                               nb_cpu=nb_cpu)

    e2 = timer()

    eprint("took {} s".format(e2 - s2))

    #3 - Extract first introns from ref transcripts
    eprint("extracting first introns from reference transcripts...")

    s3 = timer()
    ref_first_introns = get_terminal_regions(ref_introns,
                                             feature_key="intron",
                                             region_number_col="intron_number",
                                             source=None,
                                             filter_single=False,
                                             which_region="first",
                                             nb_cpu=nb_cpu
                                             )

    e3 = timer()
    eprint("took {} s".format(e3 - s3))

    #2.3 - find last exons of non-matched txs completely contained within annotated first introns
    eprint("finding novel last exons completely contained within annotated first introns...")

    s4 = timer()
    novel_nm_fi = pr.PyRanges(novel_nm_last_exons.as_df(), int64=True).overlap(pr.PyRanges(ref_first_introns.as_df(), int64=True),
                                                                               how="containment",
                                                                               strandedness="same",
                                                                               #nb_cpu=nb_cpu
                                                                               )
    e4 = timer()
    eprint("took {} s".format(e4 - s4))


    tr_ids = set(novel_nm_fi.transcript_id.tolist())
    n_tr = len(tr_ids)
    eprint("number of novel tx with first intron contained annotated first introns - {}".format(n_tr))


    #2.5 - Get 3'ends of first exons of transcripts with first-intron contained last exons
    eprint("finding 3'ends of first exons of novel txipts with last exons fully contained within annotated introns...")
    s5 = timer()

    novel_nm_fi_fe_3p = (get_terminal_regions(novel_exons_nm.subset(lambda df: df["transcript_id"].isin(tr_ids),
                                                                    nb_cpu=nb_cpu
                                                                    ),
                                              which_region="first",
                                              source="stringtie",
                                              nb_cpu=nb_cpu
                                              )
                         .three_end()
                         )

    e5 = timer()
    eprint("took {} s".format(e5 - s5))


    #2.6 - Get 3'ends of first exons of reference transcripts
    eprint("finding 3'ends of first exons of reference transcripts...")
    s6 = timer()

    ref_first_exons_3p = (get_terminal_regions(ref_exons,
                                               region_number_col="exon_number",
                                               source=None,
                                               filter_single=True,
                                               which_region="first",
                                               nb_cpu=nb_cpu
                                               )
                          .three_end())

    e6 = timer()
    eprint("took {} s".format(e6 - s6))


    eprint("checking whether 3'end coordinates have been defined correctly...")

    eprint("checking reference first exon 3'ends...")
    if check_three_end(ref_first_exons_3p):
        eprint("Start and End values are equal (impossible to have overlap)")
        ref_first_exons_3p = ref_first_exons_3p.assign("End", lambda df: df.End + 1)

    eprint("checking novel first exon 3'ends...")
    if check_three_end(novel_nm_fi_fe_3p):
        eprint("Start and End values are equal (impossible to have overlap)")
        novel_nm_fi_fe_3p = novel_nm_fi_fe_3p.assign("End", lambda df: df.End + 1)


    # 2.7 - Find first-intron contained novel LE isoforms with the outgoing SJ of first exon matching a ref first exon
    eprint("finding first exon exact 3'end matches between novel and reference first exons...")
    s7 = timer()

    try:
        first_intron_contained_match = novel_nm_fi_fe_3p.join(ref_first_exons_3p,
                                                              strandedness="same",
                                                              suffix="_ref",
                                                              nb_cpu=nb_cpu
                                                              )
    except KeyError:
        # This is the 2nd join/introns in function call problem popping up again
        # how = kwargs["how"]; KeyError: 'how'
        eprint("multithreaded run failed... trying single-threaded")

        # first_intron_contained_match = novel_nm_fi_fe_3p.join(ref_first_exons_3p,
        #                                                       how=None,
        #                                                       strandedness="same",
        #                                                       suffix="_ref",
        #                                                       nb_cpu=1
        #                                                       )

        first_intron_contained_match = (_join_grs(novel_nm_fi_fe_3p,
                                                  ref_first_exons_3p,
                                                  how="left",
                                                  strandedness="same",
                                                  suffix="_ref",
                                                  nb_cpu=1)
                                        .subset(lambda df: df.Start_b != -1, nb_cpu=nb_cpu))

    e7 = timer()
    eprint("took {} s".format(e7 - s7))



    # eprint("first_intron_contained_match columns {}".format(first_intron_contained_match.columns))

    if len(first_intron_contained_match.as_df().index) == 0:
        eprint("No first intron contained novel last exon isoforms had eactly matching 3'boundaries of reference first exon. Returning chain_match_info with no additional IDs...")

        return chain_match_info

    else:

        first_intron_contained_match = first_intron_contained_match.as_df()[["transcript_id","transcript_id_ref"]].rename({"transcript_id": "transcript_id_novel"}, axis=1)

        if isinstance(chain_match_info, pd.Series):
            return pd.concat([pd.DataFrame(chain_match_info), first_intron_contained_match]).reset_index(drop=True)

        elif isinstance(chain_match_info, pd.DataFrame):
            return pd.concat([chain_match_info, first_intron_contained_match]).reset_index(drop=True)


def _df_add_intron_number(df, out_col):
    '''

    '''

    n_exons = len(df.index)

    # Note: (I think) dictionary unpacking is required so out_col can be a variable...

    if (df["Strand"] == "+").all():
        # first in order by Start position in each txipt = left-most start position (i.e. most 5')
        return df.assign(**{out_col: list(range(1, n_exons + 1))})

    elif (df["Strand"] == "-").all():
        # firs in order by Start position in each txipt = most 3' is left-most start position
        return df.assign(**{out_col: list(range(1, n_exons +1))[::-1]})


def add_intron_number(introns, id_col = "transcript_id", out_col="intron_number", nb_cpu=1):
    '''
    '''

    start = timer()

    assert len(set(introns.as_df().Feature.tolist())) == 1, "only one feature type (e.g. all introns, all exons) should be present in gr"

    # Sort by position (could add an nb_cpu here...)
    introns = introns.sort()

    introns_out = (introns.apply(lambda df:
                                 df.groupby(id_col)
                                 .apply(lambda x: _df_add_intron_number(x, out_col)),
                                 nb_cpu=nb_cpu
                                 ))


    end = timer()
    eprint("took {} s".format(end - start))

    return introns_out


def main(novel_path, ref_path, match_by, max_terminal_non_match, out_gtf, nb_cpu):
    '''
    '''
    start = timer()

    eprint("reading in input gtf files, this can take a while...")
    eprint("reading gtf containing novel assembled transcripts...")

    s1 = timer()
    novel = pr.read_gtf(novel_path)
    e1 = timer()

    eprint("took {} s".format(e1 - s1))

    eprint("reading gtf containing reference transcripts...")

    s2 = timer()
    ref = pr.read_gtf(ref_path)
    e2 = timer()

    eprint("took {} s".format(e2 - s2))

    eprint("extracting protein-coding & lncRNA gene types from reference annotation file...")

    start2 = timer()
    ref_pc = ref.subset(lambda df: df["gene_type"].isin(["protein_coding", "lncRNA"]), nb_cpu=nb_cpu)
    end2 = timer()

    eprint("took {} s".format(end2 - start2))

    eprint("extracting novel exons from input GTF file...")

    start4 = timer()
    novel_exons = novel.subset(lambda df: df["Feature"] == "exon", nb_cpu=nb_cpu)
    end4 = timer()

    eprint("took {} s".format(end4 - start4))

    eprint("extracting exons of protein_coding and lncRNA genes from reference annotation...")

    start5 = timer()
    ref_pc_exons = ref.subset(lambda df: df["Feature"] == "exon", nb_cpu=nb_cpu)
    end5 = timer()

    eprint("took {} s".format(end5 - start5))

    eprint("finding introns for each reference transcript...")

    start3 = timer()

    try:
        ref_pc_introns = ref_pc.features.introns(by="transcript", nb_cpu=nb_cpu)
    except KeyError:
        # Specific error with ray when execute introns func for the 2nd time in same script...
        # KeyError: 'by'
        # Whilst avoiding working out what's going on, I can run my super slow intron finding script...
        eprint("pr.features.introns returned KeyError ('by'), using my hacky intron finding workaround...")
        ref_pc_introns = introns_by_tx(ref_pc_exons, nb_cpu=nb_cpu)

    end3 = timer()

    eprint("took {} s".format(end3 - start3))

    eprint("finding introns for each novel transcript...")

    start1 = timer()

    try:
        novel_introns = novel.features.introns(by="transcript", nb_cpu=nb_cpu)
    except KeyError:
        # Specific error with ray when execute introns func for the 2nd time in same script...
        # KeyError: 'by'
        # Whilst avoiding working out what's going on, I can run my super slow intron finding script...
        eprint("pr.features.introns returned KeyError ('by'), using my hacky intron finding workaround...")
        novel_introns = introns_by_tx(novel_exons, nb_cpu=nb_cpu)

    end1 = timer()

    eprint("took {} s".format(end1 - start1))

    # Need/want to make sure each introns object has a strand-aware intron number by transcript_id (i.e. 1 = first intron
    try:
        assert "intron_number" in ref_pc_introns.as_df().columns.tolist()
    except AssertionError:
        eprint("adding intron_number column to reference introns object...")
        ref_pc_introns = add_intron_number(ref_pc_introns, nb_cpu=nb_cpu)


    try:
        assert "intron_number" in novel_introns.as_df().columns.tolist()

    except AssertionError:
        eprint("adding intron_number column to novel introns object...")
        novel_introns = add_intron_number(novel_introns, nb_cpu=nb_cpu)


    eprint("finding novel transcripts with valid matches in their intron chain to reference transcripts...")

    start6 = timer()
    valid_matches = filter_transcripts_by_chain(novel_exons,
                                                novel_introns,
                                                ref_pc_exons,
                                                ref_pc_introns,
                                                match_type=match_by,
                                                max_terminal_non_match=max_terminal_non_match,
                                                nb_cpu=nb_cpu)
    end6 = timer()

    eprint("took {} s".format(end6 - start6))

    eprint("finding novel last exons completely contained within annotated first introns...")
    s7 = timer()

    fi_valid_matches = filter_first_intron_tx(novel_exons,
                                              ref_pc_exons,
                                              ref_pc_introns,
                                              valid_matches,
                                              nb_cpu)

    e7 = timer()
    eprint("finding first intron novel isoforms - took {} s".format(e7 - s7))

    if isinstance(fi_valid_matches, pd.Series):
        # match type was any
        valid_novel = novel.subset(lambda df: df["transcript_id"].isin(set(fi_valid_matches.tolist())), nb_cpu=nb_cpu)
        valid_novel.to_gtf(out_prefix + ".gtf")

    elif isinstance(fi_valid_matches, pd.DataFrame):
        # match_by/match_type was transcript
        valid_novel = novel.subset(lambda df: df["transcript_id"].isin(set(fi_valid_matches["transcript_id_novel"].tolist())), nb_cpu=nb_cpu)
        valid_novel.to_gtf(out_prefix + ".gtf")

        valid_matches.to_csv(out_prefix + ".match_stats.tsv", sep="\t", header=True)


    end = timer()

    eprint("Completed: script took {} s / {} min (3dp) ".format(round(end - start, 3), round((end - start) / 60, 3)))





if __name__ == '__main__':

    descrpn = """Script to filter assembled transcripts for intron chain matches with reference transcripts
                 to identify novel 3'end transcript isoforms"""

    parser = argparse.ArgumentParser(description=descrpn)

    parser.add_argument("-i", "--input-transcripts", default='', dest="novel_gtf", help = "path to GTF file containing novel transcripts assembled by StringTie.", required=True)
    parser.add_argument("-r", "--reference-transcripts", default='', type=str, dest="ref_gtf", help="path to GTF file containing reference transcripts against which to match intron chains of novel transcripts. Should contain same chromosome naming scheme", required=True)
    parser.add_argument("-m", "--match-by", default="any", type=str, choices=["any", "transcript"], dest="match_by", help="Consider novel transcript a valid match if all but penultimate intron(s) match introns of any transcript ('any') or the same transcript ('transcript'). 'transcript' is CURRENTLY UNIMPLEMENTED (default: %(default)s)")
    parser.add_argument("-n", "--max-terminal-non-match", default=1, type=int, dest="max_terminal_non_match", help="Maximum number of uninterrupted reference-unmatched introns at 3'end of novel transcript for it to be considered a valid match (default: %(default)s)")
    parser.add_argument("-c", "--cores", default=1, type=int, help="number of cpus/threads for parallel processing (default: %(default)s)")
    parser.add_argument("-o", "--output-prefix", type=str, default="intron_chain_matched_transcripts", dest="output_prefix", help="Prefix for output files (GTF with valid matches, matching stats TSV etc.). '.<suffix>' added depending on output file type (default: %(default)s)")

    if len(sys.argv) == 1:
        parser.print_help()
        parser.exit()

    args = parser.parse_args()



    main(args.novel_gtf, args.ref_gtf, args.match_by, args.max_terminal_non_match, args.output_prefix, args.cores)
