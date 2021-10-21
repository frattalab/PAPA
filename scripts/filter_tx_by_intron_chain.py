#! /usr/bin/env python3
from __future__ import print_function
import pyranges as pr
import numpy as np
import pandas as pd
from papa_helpers import eprint, add_region_number, get_terminal_regions, get_internal_regions, introns_by_tx
import argparse
# import os
import sys
# import logging
from timeit import default_timer as timer

'''
Filter assembled transcripts for those that match reference transcripts in their intron chain up until their penultimate intron
This allows identification of novel last exons (conservatively), but limits addition of false positive transcripts
'''

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




def validate_matching_chain(df, max_terminal_non_match=1):
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


def _agg_validate_matching_chain(df, max_terminal_non_match=1, colnames=["match_class", "n_terminal_non_match", "runs", "starts", "vals"]):
    '''
    Internal function on grouped dataframe containing all matches for each of its introns
    Returns a summarised df of tx_id | match_class | n_terminal_non_match
    tx_id - transcript_id for novel tx
    match_class - str - 'valid' or 'not_valid' - does intron chain match until penultimate introns?
    max_terminal_non_match - int - number of consecutive introns at 3'end without a reference match

    Would be good to return the summarised rle
    '''

    runs, starts, vals = rle(df["match"])
    # eprint("runs[-1] - {}".format(runs[-1]))

    # Valid matches are:
    # All introns match (e.g. bleedthrough event)
    # All but last x introns match (usually one) (i.e. runs = 1,0)

    # str_match_runs = [",".join()]

    if np.all(vals):
        # All introns match (e.g. bleedthrough event, 3'UTR extension (and any reassembled reference transcripts))
        return pd.DataFrame({colnames[0]: ["valid"],
                             colnames[1]: [0],
                             colnames[2]: ":".join(runs.astype(str)),
                             colnames[3]: ":".join(starts.astype(str)),
                             colnames[4]: ":".join(vals.astype(str))
                             }
                            )

    elif np.array_equal(vals, [1,0]):
        # All but last x introns match
        if runs[-1] <= max_terminal_non_match:
            # n of 3'end non-matching introns is less than cut-off
            return pd.DataFrame({colnames[0]: ["valid"],
                                 colnames[1]: [runs[-1]],
                                 colnames[2]: ":".join(runs.astype(str)),
                                 colnames[3]: ":".join(starts.astype(str)),
                                 colnames[4]: ":".join(vals.astype(str))
                                 }
                                )

        else:
            # Too many unmatched at 3'end
            return pd.DataFrame({colnames[0]: ["not_valid"],
                                 colnames[1]: [runs[-1]],
                                 colnames[2]: ":".join(runs.astype(str)),
                                 colnames[3]: ":".join(starts.astype(str)),
                                 colnames[4]: ":".join(vals.astype(str))
                                 }
                                )

    else:
        # No exact match, no match starting at 5'end followed by only non-matched
        return pd.DataFrame({colnames[0]: ["not_valid"],
                             colnames[1]: [np.nan],
                             colnames[2]: ":".join(runs.astype(str)),
                             colnames[3]: ":".join(starts.astype(str)),
                             colnames[4]: ":".join(vals.astype(str))
                             }
                            )




def _join_grs(gr_left, gr_right, strandedness=None, how=None, report_overlap=False, slack=0, suffix="_b", nb_cpu=1):
    '''
    '''

    return gr_left.join(gr_right, strandedness, how, report_overlap, slack, suffix, nb_cpu)


def check_three_end(gr):
    '''
    '''
    df = gr.as_df()

    return all(df.Start == df.End)


def _check_int64(gr):
    '''
    Helper function to check if Start & End columns are int64, if not to convert
    Useful as a safety before pr.join as errors often arise if leave at int32
    '''

    if gr.dtypes.loc["Start"] != "int64":
        gr = gr.assign("Start", lambda df: df["Start"].astype("int64"))

    if gr.dtypes.loc["End"] != "int64":
        gr = gr.assign("End", lambda df: df["End"].astype("int64"))

    return gr


def _pd_merge_gr(df, df_to_merge, how, on, suffixes, to_merge_cols):
    '''
    Perform a pd.merge inside a pr.apply to add columns from gr_to_merge based on metadata in gr_df
    Here, will check the chromosome and strand of provided gr (gr_df)
    and subset gr_to_merge to chr/strand before converting to df
    This should cut down on memory requirements when convert to df (which requires pd.concat() x chrom & strand)
    For this kind of merge, only expect joins between regions on same chr/strand
    '''
    #chromsomes returns a list of chr names (always 1 val)
    assert isinstance(to_merge_cols, list)

    if df_to_merge.empty:
        eprint("df_to_merge for chr/strand pair {} is empty - returning to_merge cols filled with NaNs".format(",".join([df.Chromosome.iloc[0], df.Strand.iloc[0]])))

        df_cols = df.columns.tolist()
        # on won't have suffix added - need to remove as a target column
        to_merge_cols = [col for col in to_merge_cols if col != on]

        # eprint(df_cols)
        # eprint(to_merge_cols)

        # list of cols shared between dfs - need to add suffixes[1]
        # list of cols only in df_to_merge - these stay the same in a merge
        only_cols = [col for col in to_merge_cols if col not in df_cols]
        shared_cols = [col for col in to_merge_cols if col in df_cols]

        # eprint("only_cols - {}".format(only_cols))
        # eprint("shared_cols - {}".format(shared_cols))

        out_shared = [col + suffixes[1] for col in shared_cols]
        target_cols = out_shared + only_cols

        nrows = len(df.index)

        out_cols = {col: pd.Series([np.nan]*nrows) for col in target_cols}

        return df.assign(**out_cols)

    else:
        return df.merge(df_to_merge,
                        how=how,
                        on=on,
                        suffixes=suffixes)



def chain_match_by_tx(df, id_col="transcript_id", max_terminal_non_match=1):
    '''
    Applied to a pandas df (i.e. use within pyranges.apply())
    '''

    assert id_col in df.columns

    # Check that df has at least 1 not None id_col value
    if df.fillna(value=np.nan)[id_col].isna().all():
        # Will get key error if try to groupby column with all NA/None values
        return pd.DataFrame()

    else:
        # Try chain matching
        try:
            out_df = (df.groupby(id_col)
                        .apply(lambda x: _agg_validate_matching_chain(x,
                                                                      max_terminal_non_match))
                        # return to a column
                        .reset_index(id_col)
                        # drops weird 0 0 0 index made by my function... (To do: work out why)
                        .reset_index(drop=True)
                     )

        except KeyError:
            # sometimes get this... KeyError: 'Requested level (transcript_id_novel) does not match index name (None)'
            eprint("df returned 'transcript_id_novel' KeyError for chr/strand pair {} - returning empty df".format(",".join([df.Chromosome.iloc[0], df.Strand.iloc[0]])))
            out_df = pd.DataFrame()

    return out_df


def filter_transcripts_by_chain(novel_introns, ref_introns, match_type = "transcript", max_terminal_non_match=2, nb_cpu = 1):
    '''
    '''

    # Only keep essential columns (or those with potentially useful info) to save on memory
    # TO DO: this should really happen in main

    # novel_cols_to_keep = [col for col in ["Feature","transcript_id"] if col in novel_exons.columns.tolist()]
    # ref_cols_to_keep = [col for col in ["Feature", "transcript_id", "gene_id", "gene_name"] if col in ref_exons.columns.tolist()]

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

    # novel_exons = novel_exons.sort()
    novel_introns = novel_introns.sort()
    # ref_exons = ref_exons.sort()
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
    eprint("generating gr of novel txipts sorted by intron number...")

    t5 = timer()

    novel_intron_ids_ordered = novel_introns.apply(lambda df: df.sort_values(by=["transcript_id", "intron_number"], ascending=True))
    # gr of txipt_id | intron_id | intron_number
    novel_intron_ids_ordered = novel_intron_ids_ordered[["transcript_id","intron_id","intron_number"]]

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
    novel_introns = novel_introns.assign("Start", lambda df: df.Start.astype("int64")).assign("End", lambda df: df.End.astype("int64"))
    ref_introns = ref_introns.assign("Start", lambda df: df.Start.astype("int64")).assign("End", lambda df: df.End.astype("int64"))

    joined = novel_introns.join(ref_introns, strandedness="same", suffix ="_ref", nb_cpu=nb_cpu)

    t8 = timer()

    eprint("took {} s".format(t8 - t7))

    #5. Filter for overlaps that exactly match (or differ by given tolerance - for now not supported)
    eprint("filtering overlaps for exact matches...")

    t9 = timer()
    joined = joined.subset(lambda df: abs(df.Start - df.Start_ref) + abs(df.End - df.End_ref) <= 0, nb_cpu=nb_cpu)
    t10 = timer()

    eprint("took {} s".format(t10 - t9))

    # Minimal info needed on matches between novel and reference introns
    joined = joined[["transcript_id","intron_id","transcript_id_ref","intron_id_ref","intron_number_ref"]]

#     eprint(joined.dtypes)

    #6. Join ordered novel introns with match info
    #7. Assign a simple tracker column 'match' of True (where intron is matched) and False (where intron is not matched)

    eprint("preparing for filtering intron matches...")
    t11 = timer()

    if match_type == "any":
        # Looking for intron to match any annotated intron, regardless of reference transcript

        # i) Make a gr joining matched novel introns with all novel introns

        # If a chrom/strand pair has no matches (i.e. is an empty df)
        # pd.merge will not produce output,
        # and dfs in PyRanges object will be incompatible as don't have matching columns
        # Use this as set of cols to expect if matched_df (joined) is empty
        joined_cols = joined.columns.tolist()

        novel_ref_match_info = novel_intron_ids_ordered.apply_pair(joined,
                                                                   lambda all_novel, matched_novel:
                                                                   _pd_merge_gr(all_novel,
                                                                                matched_novel,
                                                                                how="left",
                                                                                on="intron_id",
                                                                                suffixes=[None, "_match"],
                                                                                to_merge_cols=joined_cols),
                                                                   nb_cpu=nb_cpu
                                                                   )

        # eprint(novel_ref_match_info)
        #
        # for chr, df in novel_ref_match_info.items():
        #     eprint(chr)
        #     eprint(df.columns)
        #     eprint(len(df.columns))

        # ii) Assign 'match' column for each intron_id
        # pd.merge puts an NaN in non-overlapping rows (i.e. intron not matched)
        novel_ref_match_info = novel_ref_match_info.assign("match",
                                                           lambda df:
                                                           pd.Series(np.where(np.isnan(df["Start_match"]),
                                                                              0,
                                                                              1)
                                                                     )
                                                           )

        # iii) check that first intron of each novel transcript matches first intron of a reference transcript
        # Replace 'match' with 0 if novel first intron doesn't match ref first intron

        novel_ref_match_info = novel_ref_match_info.assign("match",
                                                           lambda df: pd.Series(np.where(
                                                                                         (df["intron_number"] == 1) &
                                                                                         (df["intron_number_ref"] != 1),
                                                                                         0,
                                                                                         df["match"])))

        # iv) Make 'match' an ordered categorical column describing whether match is to first or other intron
        # With priority to a match over not match
        novel_ref_match_info = novel_ref_match_info.assign("match",
                                                           lambda df:
                                                           df["match"].astype("category")
                                                           .cat
                                                           .set_categories([1,0], ordered=True)
                                                           )

        # Now when drop duplicates a match is prioritised over no matches
        novel_ref_match_info = novel_ref_match_info.apply(lambda df:
                                                          df.sort_values(["transcript_id",
                                                                          "intron_number",
                                                                          "match"])
                                                          .drop_duplicates(subset=["intron_id","match"],
                                                                           keep="first"),
                                                          nb_cpu=nb_cpu
                                                          )

        # Subset to minimal required metadata for chain_match checking
        novel_ref_match_info = novel_ref_match_info.apply(lambda df: df.rename(columns={"transcript_id": "transcript_id_novel"}))
        novel_ref_match_info = novel_ref_match_info[["transcript_id_novel", "intron_id", "intron_number", "match"]]


    # elif match_type == "transcript":
        # Looking for introns (except last) to match the same reference transcript
        # merge_ordered can do a 'grouped merge' filling in empty rows (introns) for each transcript_id
        # This is especially useful if want to do transcript-specific intron matching
        # For each reference transcript, all novel introns will be filled with NaN if no overlap for given transcript_id
        # (i.e. novel txipt matches all but last intron of reference transcript)

        # novel_ref_match_info = (pd.merge_ordered(novel_intron_ids_ordered,
        #                             joined,
        #                             how="left",
        #                             on="intron_id",
        #                             right_by="transcript_id_ref", # group matches by ref tx & join tx by tx
        #                             suffixes=["_novel","_match"],
        #                             fill_method=None)
        #            .sort_values(by=["transcript_id_novel","intron_number"])
        #                        )
        #
        # # merge_ordered fills rows for each intron for each ref tx in df, regardless of whether any overlap
        # # .dropna(axis="rows", subset=["intron_id_ref"])
        # novel_ref_match_info = (novel_ref_match_info.groupby(["transcript_id_novel", "transcript_id_ref"])
        #                         .filter(lambda df: (df["intron_id_ref"].notna()).any()) # Retained if ref tx has >=1 matching introns
        #                         .reset_index(drop=True))
        #
        # # Make a match column where 1 = match, 0 = no match for each ref id and novel intron
        # novel_ref_match_info["match"] = novel_ref_match_info["intron_id_ref"]
        # novel_ref_match_info["match"] = novel_ref_match_info["match"].fillna(0)
        # novel_ref_match_info["match"] = novel_ref_match_info["match"].replace("\w*", 1, regex=True)

    t12 = timer()
    eprint("took {} s".format(t12 - t11))


    # 8. Filter down matching transcripts to those that all ref introns except penultimate or all introns...
    eprint("categorising intron chain matches as valid/invalid...")
    t13 = timer()

    if match_type == "any":
        # Only need to check by novel transcript_id
        # filt_novel_ref_match_info = (novel_ref_match_info.groupby("transcript_id_novel")
        #                              .filter(lambda x: validate_matching_chain(x, max_terminal_non_match)
        #                                     )
        #                             )
        # for chrom, df in novel_ref_match_info:
        #     eprint(chrom)
        #     eprint(df.columns)
        #     eprint("Empty df - {}".format(df.empty))
        #     eprint(df["transcript_id_novel"].value_counts(dropna=False))
        #     eprint(df["transcript_id_novel"].value_counts(dropna=False).loc[lambda x: x.index == None])

        novel_ref_match_info_agg = (novel_ref_match_info.apply(lambda df:
                                                               chain_match_by_tx(df,
                                                                                 "transcript_id_novel",
                                                                                 max_terminal_non_match),
                                                               as_pyranges=False, # Summarises by Tx and drops coord info
                                                               nb_cpu=nb_cpu)
                                    )


    # elif match_type == "transcript":
    #     # Check novel tx vs each ref tx
    #     filt_novel_ref_match_info = (novel_ref_match_info.groupby(["transcript_id_novel","transcript_id_ref"])
    #                                  .filter(lambda x: validate_matching_chain(x, max_terminal_non_match)
    #                                         )
    #                                 )

    t14 = timer()
    eprint("took {} s".format(t14 - t13))

    # Now turn into single pandas df (1 row per novel txipt)
    # Ignore_index to drop the PyRanges chrom/strand keys as not needed
    novel_ref_match_info_agg = pd.concat(novel_ref_match_info_agg,
                                         ignore_index=True)

    # eprint(novel_ref_match_info_agg)


    return novel_ref_match_info_agg


    # # Return simplified df of novel transcript_id & matching transcript_id if applicable
    # if match_type == "any":
    #     return filt_novel_ref_match_info["transcript_id_novel"].drop_duplicates()
    #
    # elif match_type == "transcript":
    #     return filt_novel_ref_match_info[["transcript_id_novel", "transcript_id_ref"]].drop_duplicates()





def filter_first_intron_tx(novel_first_exons,
                           novel_last_exons,
                           ref_first_exons,
                           ref_first_introns,
                           chain_match_info,
                           novel_source,
                           nb_cpu=1):
    '''
    Function to return novel last exon transcripts occurring within annotated first exons,
    In which the 3'end boundary of the novel first exon exactly matches an annotated first exon
    These transcript types cannot have valid intron chain match but can still produce valid last exon isoforms
    A well known e.g. of this isoform type is TDP-43 depletion sensitive STMN2 isoform
    '''

    assert isinstance(chain_match_info, pd.DataFrame)

    # start = timer()
    eprint("finding novel last exon isoforms contained within annotated first introns...")

    #1 - Pull out non-chain matched novel isoforms
    eprint("extracting non-chain matched novel isoforms...")
    s1 = timer()

    mask = chain_match_info["match_class"] == "valid"
    # matched_chain_match_info = chain_match_info[mask]
    nm_chain_match_info = chain_match_info[~mask]

    # eprint(nm_chain_match_info)


    e1 = timer()
    eprint("took {} s".format(e1 - s1))


    #2 -Extract last exons of non-chain matched novel isoforms
    eprint("extracting last exons of non-chain matched novel isoforms")

    s2 = timer()
    novel_last_exons_nm = (novel_last_exons.subset(lambda df:
                                                   df["transcript_id"]
                                                   .isin(set(nm_chain_match_info
                                                            ["transcript_id_novel"]
                                                            .tolist()
                                                             )
                                                         ),
                                                   nb_cpu=nb_cpu)
                           )
    e2 = timer()

    eprint("took {} s".format(e2 - s2))

    #3 - Extract first introns from ref transcripts
    # No longer needed...

    # eprint("extracting first introns from reference transcripts...")
    #
    # s3 = timer()
    # ref_first_introns = get_terminal_regions(ref_introns,
    #                                          feature_key="intron",
    #                                          region_number_col="intron_number",
    #                                          source=None,
    #                                          filter_single=False,
    #                                          which_region="first",
    #                                          nb_cpu=nb_cpu
    #                                          )
    #
    # e3 = timer()
    # eprint("took {} s".format(e3 - s3))

    #2.3 - find last exons of non-matched txs completely contained within annotated first introns
    eprint("finding novel last exons completely contained within annotated first introns...")
    novel_last_exons_nm = novel_last_exons_nm.assign("Start", lambda df: df.Start.astype("int64")).assign("End", lambda df: df.End.astype("int64"))
    ref_first_introns = ref_first_introns.assign("Start", lambda df: df.Start.astype("int64")).assign("End", lambda df: df.End.astype("int64"))

    s4 = timer()

    novel_nm_fi = novel_last_exons_nm.overlap(ref_first_introns,
                                              how="containment",
                                              strandedness="same")

    e4 = timer()
    eprint("took {} s".format(e4 - s4))


    try:
        # eprint(novel_nm_fi.columns)
        tr_ids = set(novel_nm_fi.transcript_id.tolist())
        n_tr = len(tr_ids)

    except AssertionError:
        # empty PyRanges so transcript_id column not present
        tr_ids = set()
        n_tr = 0


    if n_tr == 0:
        eprint("0 novel transcripts with last exons contained within " +
               "annotated exons - returning initial chain matching info df")
        return chain_match_info

    eprint("number of novel tx with first intron contained annotated first introns - {}".format(n_tr))

    #2.5 - Get 3'ends of first exons of transcripts with first-intron contained last exons
    eprint("finding 3'ends of first exons of novel txipts with last exons fully contained within annotated introns...")
    s5 = timer()

    novel_nm_fi_fe_3p = (novel_first_exons.subset(lambda df: df["transcript_id"].isin(tr_ids),
                                                  nb_cpu=nb_cpu)
                                          .three_end()
                         )

    e5 = timer()
    eprint("took {} s".format(e5 - s5))


    #2.6 - Get 3'ends of first exons of reference transcripts
    eprint("finding 3'ends of first exons of reference transcripts...")
    s6 = timer()

    ref_first_exons_3p = ref_first_exons.three_end()

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
        # Reformat first_intron_contained_match to match chain_match_info and update first intron matched transcripts to valid
        first_intron_contained_match = (first_intron_contained_match.as_df()
                                        [["transcript_id","transcript_id_ref"]]
                                        .rename({"transcript_id": "transcript_id_novel"}, axis=1)
                                        .assign(**{"match_class": "valid",
                                                   "n_terminal_non_match": 0,
                                                   "isoform_class": "first_intron_spliced"
                                                   }
                                                ))

        fi_ids = set(first_intron_contained_match["transcript_id_novel"].tolist())

        eprint("Number of filtered first intron novel last exon tx - {}".format(len(fi_ids)))

        # if isinstance(chain_match_info, pd.Series):
        #     return pd.concat([pd.DataFrame(chain_match_info), first_intron_contained_match]).reset_index(drop=True)
        # elif isinstance(chain_match_info, pd.DataFrame):

        # Update chain_match_info with previous first intron transcripts now classified as valid matches

        return (pd.concat([chain_match_info[~chain_match_info["transcript_id_novel"].isin(fi_ids)],
                          first_intron_contained_match
                           ]
                          )
                .reset_index(drop=True)
                )


def add_3p_extension_length(gr,
                            ref_gr,
                            id_col="transcript_id",
                            out_col="3p_extension_length",
                            nb_cpu=1):
    '''
    Add column '3p_extension_length' reporting 3'end extension of overlapping regions in gr & ref_gr (distance relative to gr)
    Note that for each unique ID (id_col) in gr, the smallest extension will be reported
    Avoids cases where overlaps with short & long isoforms, but tx is just reassembly of long isoform
    (so it's just reassembly of longer isoform, but is an extension relative to shorter ref isoform)
    '''

    # Find columns unique to ref_gr, so can drop at the end (no suffix added)
    not_cols = gr.columns.tolist()
    ref_to_drop = [col for col in ref_gr.columns if col not in not_cols]

    joined = gr.join(ref_gr,
                     strandedness="same",
                     how=None,
                     nb_cpu=nb_cpu)

    joined = joined.assign(out_col,
                           lambda df: df["End"] - df["End_b"] if (df["Strand"] == "+").all() else
                           df["Start_b"] - df["Start"], nb_cpu=nb_cpu)


    # To avoid capturing extensions of shorter isoforms (that really just are the longer known isoform)
    # Pick the smallest extension for each transcripts
    joined = joined.apply(lambda df: df.sort_values([id_col, out_col],
                                                    ascending=True).drop_duplicates(subset=[id_col],
                                                                                keep="first"),
                          nb_cpu=nb_cpu)


    return joined.drop(ref_to_drop).drop(like="_b$")


def filter_complete_match(novel_last_exons,
                          ref_first_exons,
                          ref_last_exons,
                          ref_exons,
                          ref_introns,
                          min_extension_length,
                          chain_match_info,
                          novel_source,
                          nb_cpu=1):
    '''
    Filter transcripts with complete intron chain matches for bleedthrough intronic last exons and 3'UTR extensions
    Reassembled reference transcripts will have a complete intron chain match to ref but not be meaningful novel isoforms
    This wrapper script will filter chain_match_info to remove reassembled ref transcripts
    And assign a 'isoform_class' column as well - 'exon_bleedthrough', 'utr_extension'
    '''

    assert isinstance(min_extension_length, int)

    #1. Extract complete match isoforms from chain_match_info
    exact_ids = set(chain_match_info.loc[(chain_match_info["match_class"] == "valid") &
                                         (chain_match_info["n_terminal_non_match"]
                                            .astype(float)
                                            .astype("Int64") == 0),
                                         "transcript_id_novel"].tolist())

    # eprint("ids with complete intron chain match - {}".format(",".join(exact_ids)))

    #2. Identify 'bleedthrough' intronic events, where:
    #### - 3'end lies within an annotated intron
    #### - Novel last exon overlaps with ref internal exon
    #### - 3'end of  novel last exon is downstream of overlapping ref internal intron

    # A - check three end  of last is within annotated intron.
    m_l_novel_exons = novel_last_exons.subset(lambda df: df["transcript_id"].isin(exact_ids), nb_cpu=nb_cpu)

    m_l_novel_exons_3p = m_l_novel_exons.three_end()

    if check_three_end(m_l_novel_exons_3p):
        m_l_novel_exons_3p = m_l_novel_exons_3p.assign("End", lambda df: df.End + 1)

    # eprint("\n3'ENDS OF NOVEL EXONS WITH EXACT CHAIN MATCHES")
    # eprint(m_l_novel_exons_3p)


    # eprint("\nTHIS IS REF INTRONS")
    # eprint(ref_introns)

    intron_cont_3p_ids = (set(m_l_novel_exons_3p.overlap(ref_introns,
                                                         strandedness="same",
                                                         how="containment")
                                                .as_df()
                                                ["transcript_id"].tolist()
                              )
                          )

    # eprint("\nTHIS IS IDS OF 3'ENDS COMPLETELY CONTAINED WITHIN ANNOTATED INTRONS")
    # eprint(intron_cont_3p_ids)

    m_l_novel_exons_bl = m_l_novel_exons.subset(lambda df: df["transcript_id"].isin(intron_cont_3p_ids), nb_cpu=nb_cpu)

    # B - Get reference internal exons
    eprint("finding reference internal exons...")
    s2 = timer()

    ref_exons_int = get_internal_regions(ref_exons)

    e2 = timer()
    eprint("took {} s".format(e2 - s2))

    # C - find last exons with 3'end downstream of annotated internal exon 3'end
    eprint("finding tx ids of novel isoforms with bleedthrough last exons...")

    # '3p_extension_length' column denoting length of 3'-end extension relative to olapping internal exons
    bleedthrough_gr_pre_l = add_3p_extension_length(m_l_novel_exons_bl, ref_exons_int, nb_cpu=nb_cpu)
    bleedthrough_gr_pre_l = bleedthrough_gr_pre_l.subset(lambda df: df["3p_extension_length"] > 0)

    bld_ids_pre_l = set(bleedthrough_gr_pre_l.transcript_id)
    bld_ext_len_dist = bleedthrough_gr_pre_l.as_df()["3p_extension_length"].describe(percentiles=[i * 0.1 for i in range(1,11,1)])

    eprint(f"Number of putative internal bleedthrough events - {len(bld_ids_pre_l)}")
    eprint(f"Internal bleedthrough extension length distribution\n{bld_ext_len_dist}")

    # Subset for extensions >= min_extension_length
    bleedthrough_gr_post_l = bleedthrough_gr_pre_l.subset(lambda df: df["3p_extension_length"] >= min_extension_length)
    bld_ids_post_l = set(bleedthrough_gr_post_l.transcript_id)
    eprint(f"After minimum length filter - {min_extension_length} - number of putative internal bleedthrough events - {len(bld_ids_post_l)}")

    # Set of IDs of putative events failing min length filter
    bld_ids_fail_l = bld_ids_pre_l.difference(bld_ids_post_l)

    eprint("checking if putative internal events also overlap with ref last exons...")
    # Find putative bleedthrough exons that do not overlap any known last exon
    # i.e. Overlaps an annotated 'hybrid' internal exon
    # and 3'end is contained within downstream corresponding intron
    # But exon can also act as terminal exon with downstream poly(A) site
    # the sometimes terminates further downstream
    # Putative bleedthrough events could just be reassembled hybrid terminal exons

    bleedthrough_ids_post_le = (set(m_l_novel_exons_bl.subset(lambda df: df["transcript_id"].isin(bld_ids_post_l))
                                                      .overlap(ref_last_exons,
                                                               strandedness="same",
                                                               invert=True)
                                                      .as_df()
                                                      ["transcript_id"]
                                                      .tolist()
                                    )
                                )

    eprint("number of internal bleedthrough events after filtering out those overlapping ref last exons - {}".format(len(bleedthrough_ids_post_le)))

    # Get a set of ids of valid bleedthrough events overlapping ref last exons
    bleedthrough_ids_le_olap = bld_ids_post_l.difference(bleedthrough_ids_post_le)

    # Get set of complete chain match IDs that aren't bleedthrough events
    exact_ids_n_bl = exact_ids.difference(bleedthrough_ids_post_le)


    ## Check for 3'UTR extensions

    # 3'end extension distance between olapping novel last exons ref last exons
    # Filter for minimum 3'end extension length
    # 3p end of novel last doesn't overlap with annotated first exon
    # (Coverage profile likely overlaps between genes, unlikely can predict correct termination event from overlapping profile...)
    m_l_novel_exons_n_bl = m_l_novel_exons.subset(lambda df: df["transcript_id"].isin(exact_ids_n_bl), nb_cpu=nb_cpu)

    m_l_novel_exons_n_bl_3p = m_l_novel_exons_n_bl.three_end()
    # eprint(m_l_novel_exons_n_bl_3p)

    if check_three_end(m_l_novel_exons_n_bl_3p):
        m_l_novel_exons_n_bl_3p = m_l_novel_exons_n_bl_3p.assign("End", lambda df: df.End + 1)


    not_fe_3p_ids = (set(m_l_novel_exons_n_bl_3p.overlap(ref_first_exons,
                                                         strandedness="same",
                                                         how="containment",
                                                         invert=True)
                                                         .as_df()
                                                         ["transcript_id"].tolist()
                         )
                     )

    m_l_novel_exons_n_bl = m_l_novel_exons_n_bl.subset(lambda df: df["transcript_id"].isin(not_fe_3p_ids), nb_cpu=nb_cpu)

    # Fin
    utr_extension_pre_l = add_3p_extension_length(m_l_novel_exons_n_bl, ref_last_exons, nb_cpu=nb_cpu)
    utr_extension_pre_l = utr_extension_pre_l.subset(lambda df: df["3p_extension_length"] > 0)
    utr_extension_ids_pre_l = set(utr_extension_pre_l.transcript_id)

    utr_ext_len_dist = utr_extension_pre_l.as_df()["3p_extension_length"].describe(percentiles=[i * 0.1 for i in range(1,11,1)])

    eprint(f"number of putative UTR/last exon extension events - {len(utr_extension_ids_pre_l)}")
    eprint(f"3'UTR extension length distribution\n{utr_ext_len_dist}")
    # Check for extensions >= min_extension_length
    utr_extension_ids_post_l = set(utr_extension_pre_l.subset(lambda df: df["3p_extension_length"] >= min_extension_length).transcript_id)

    # Set of IDs of putative events failing min length filter
    utr_extension_ids_fail_l = utr_extension_ids_pre_l.difference(utr_extension_ids_post_l)

    eprint(f"After minimum length filter - {min_extension_length} - number of putative 3'UTR extensions - {len(utr_extension_ids_post_l)}")

    # Check that 'putative 3'UTR extensions' do not overlap with any reference internal exons
    # Genes with proximal annotated last exons could have extensions that
    # correspond to intron retention with distal last exon
    # Or extensions that terminate within annotated internal exons
    # This is a heuristic to reduce the number of FP calls

    eprint("checking if 3'ends of putative extension events overlap with other exons...")

    # IDs that don't overlap with over ref exons
    utr_extension_ids_post_oe = (set(m_l_novel_exons_n_bl_3p.subset(lambda df: df["transcript_id"].isin(utr_extension_ids_post_l),
                                                                 nb_cpu=nb_cpu)
                                                         .overlap(ref_exons,
                                                                  strandedness="same",
                                                                  invert=True # return non-overlapping
                                                                  )
                                                         .as_df()["transcript_id"]
                                                         .tolist()
                                     )
                                 )

    eprint("number of UTR/last exon extension events after filtering out " +
           f"any overlapping ref exons - {len(utr_extension_ids_post_oe)}")

    # Get a set of IDs classified as extensions but olapping with other exons
    utr_extension_ids_e_olap = utr_extension_ids_post_l.difference(utr_extension_ids_post_oe)

    # both_classes = bleedthrough_ids.union(utr_extension_ids)

    def _temp_assign(df, bld_ids,
                     bld_ids_len_out, bld_ids_out,
                     ext_ids, ext_ids_len_out,
                     ext_ids_out):

        try:
            int(df["n_terminal_non_match"])

        except ValueError:
            return np.nan

        if int(df["n_terminal_non_match"]) == 0:

            if df["transcript_id_novel"] in bld_ids:
                return "internal_exon_bleedthrough"

            elif df["transcript_id_novel"] in bld_ids_len_out:
                return "internal_exon_bleedthrough_below_min_length"

            elif df["transcript_id_novel"] in bld_ids_out:
                return "internal_exon_bleedthrough_last_exon_overlap"

            elif df["transcript_id_novel"] in ext_ids:
                return "ds_3utr_extension"

            elif df["transcript_id_novel"] in ext_ids_len_out:
                return "ds_3utr_extension_below_min_length"

            elif df["transcript_id_novel"] in ext_ids_out:
                return "ds_3utr_extension_internal_exon_overlap"

            elif df["match_class"] == "valid" and not pd.isna(df["isoform_class"]):
                return df["isoform_class"]

            else:
                return "reassembled_reference"



    # Update chain_match_info with isoform_class column
    # conditions = [chain_match_info["n_terminal_non_match"].astype(float).astype("Int64").eq(0, fill_value=999) & chain_match_info["transcript_id_novel"].isin(bleedthrough_ids),
    #               chain_match_info["n_terminal_non_match"].astype(float).astype("Int64").eq(0, fill_value=999) & chain_match_info["transcript_id_novel"].isin(utr_extension_ids),
    #               chain_match_info["n_terminal_non_match"].astype(float).astype("Int64").eq(0, fill_value=999) & ~chain_match_info["transcript_id_novel"].isin(both_classes),
    #               ]
    #
    # choices = ["bleedthrough", "utr_extension", "reassembled_reference"]
    #
    # eprint("{}".format(chain_match_info["n_terminal_non_match"].astype(float).astype("Int64").eq(0, fill_value=9999) & chain_match_info["transcript_id_novel"].isin(bleedthrough_ids)))
    # chain_match_info["isoform_class"] = np.select(conditions, choices, default=np.nan)
    # n_loc = chain_match_info.columns.get_loc("n_terminal_non_match")
    # id_loc = chain_match_info.columns.get_loc("transcript_id_novel")

    chain_match_info["isoform_class"] = chain_match_info.apply(lambda df: _temp_assign(df,
                                                                                       bleedthrough_ids_post_le,
                                                                                       bld_ids_fail_l,
                                                                                       bleedthrough_ids_le_olap,
                                                                                       utr_extension_ids_post_oe,
                                                                                       utr_extension_ids_fail_l,
                                                                                       utr_extension_ids_e_olap),
                                                               axis="columns")

    # Update filtered out isoforms to not_valid so not included in filtered GTF

    chain_match_info["match_class"] = np.where(chain_match_info["isoform_class"].isin(["reassembled_reference",
                                                                                       "ds_3utr_extension_internal_exon_overlap",
                                                                                       "ds_3utr_extension_below_min_length",
                                                                                       "internal_exon_bleedthrough_last_exon_overlap",
                                                                                       "internal_exon_bleedthrough_below_min_length"]),
                                               "not_valid",
                                               chain_match_info["match_class"])

    #) chain_match_info.assign(isoform_class=lambda df: pd.Series([_temp_assign(row, bleedthrough_ids, utr_extension_ids, n_loc, id_loc)
    #                                                                               for row in df.itertuples(index = False)]))

    return chain_match_info


def annotate_3utr_introns(novel_last_introns=None,
                          ref_last_exons=None,
                          chain_match_info=None,
                          class_col="isoform_class",
                          class_key="3utr_intron_spliced",
                          nb_cpu=1):
    '''
    Identify intron chain matched isoforms with a novel last intron fully contained within an annotated 3'UTR/last exon
    '''

    assert novel_last_introns is not None
    assert ref_last_exons is not None
    assert chain_match_info is not None

    eprint("finding novel isoforms with 3'UTR introns...")

    #1. Extract currently unclassified valid isoforms for checking if 3'UTR introns

    valid_nc_ids = (set(chain_match_info.loc[lambda x: (x["match_class"] == "valid") &
                                                       (pd.isna(x["isoform_class"])),
                                             "transcript_id_novel"]
                                        .tolist())
                    )

    novel_last_introns_nc = novel_last_introns.subset(lambda df: df["transcript_id"].isin(valid_nc_ids), nb_cpu=nb_cpu)


    #2. Check if unclassified last exons are completely contained within annotated last exon
    novel_3utr_introns = novel_last_introns_nc.overlap(ref_last_exons,
                                                       how="containment",
                                                       strandedness="same")

    try:
        # eprint(novel_nm_fi.columns)
        utr3_tr_ids = set(novel_3utr_introns.as_df()["transcript_id"].tolist())
        n_tr = len(utr3_tr_ids)

    except AssertionError:
        utr3_tr_ids = set()
        n_tr = 0


    if n_tr == 0:
        eprint("0 novel transcripts with spliced 3'UTR intron" +
               "fully contained within annotated last exons found")
        return chain_match_info

    eprint("n of novel tx with spliced 3'UTR intron - {}".format(n_tr))


    #3. Update chain match info with reclassified 3'UTR intron transcripts
    chain_match_info[class_col] = np.where(chain_match_info["transcript_id_novel"].isin(utr3_tr_ids),
                                           class_key,
                                           chain_match_info[class_col])

    return chain_match_info


def annotate_distal_last_exons(novel_last_introns=None,
                               ref_last_exons=None,
                               chain_match_info=None,
                               class_col="isoform_class",
                               class_key="ds_alt_spliced",
                               nb_cpu=1):
    '''
    Annotate chain-matched novel isoforms as distal alternative last exons
    ref last exon must be completely contained within ref last intron
    '''

    #1. Extract currently unclassified valid isoforms for checking if alt distal LEs

    valid_nc_ids = (set(chain_match_info.loc[lambda x: (x["match_class"] == "valid") &
                                                       (pd.isna(x["isoform_class"])),
                                             "transcript_id_novel"]
                                        .tolist())
                    )

    novel_last_introns_nc = novel_last_introns.subset(lambda df: df["transcript_id"].isin(valid_nc_ids), nb_cpu=nb_cpu)


    #2. find ref last exons completely contained within novel last intron
    ref_last_exons_cnt = ref_last_exons.overlap(novel_last_introns_nc,
                                                strandedness="same",
                                                how="containment")

    #3. Get novel last introns overlaping with completely contained last exons
    ds_spliced = novel_last_introns_nc.overlap(ref_last_exons_cnt,
                                               strandedness="same")

    try:
        # eprint(novel_nm_fi.columns)
        ds_spliced_ids = set(ds_spliced.as_df()["transcript_id"].tolist())
        n_tr = len(ds_spliced_ids)

    except AssertionError:
        ds_spliced_ids = set()
        n_tr = 0


    if n_tr == 0:
        eprint("0 novel transcripts with distal last exons" +
               "found (where ref last exon is fully contained within novel last intron)")
        return chain_match_info

    eprint("n of novel tx with downstream spliced novel last exon - {}".format(n_tr))


    #3. Update chain match info with reclassified 3'UTR intron transcripts
    chain_match_info[class_col] = np.where(chain_match_info["transcript_id_novel"].isin(ds_spliced_ids),
                                           class_key,
                                           chain_match_info[class_col])

    return chain_match_info


def annotate_internal_spliced(novel_last_exons=None,
                              ref_introns=None,
                              ref_exons=None,
                              chain_match_info=None,
                              class_col="isoform_class",
                              class_key="internal_intron_spliced",
                              nb_cpu=1):
    '''
    Classify valid isoforms as 'internal_intron_spliced' if novel last exon is fully contained within annotated intron
    Also make a quick check that putative do not overlap with any known exons
    '''

    valid_nc_ids = (set(chain_match_info.loc[lambda x: (x["match_class"] == "valid") &
                                                       (pd.isna(x["isoform_class"])),
                                             "transcript_id_novel"]
                                        .tolist())
                    )

    novel_last_exons_nc = novel_last_exons.subset(lambda df: df["transcript_id"].isin(valid_nc_ids), nb_cpu=nb_cpu)

    # Find unclassified valid events completely contained within annotated introns
    ref_internal_introns = get_internal_regions(ref_introns,
                                                feature_col="Feature",
                                                feature_key="intron",
                                                id_col="transcript_id",
                                                region_number_col="intron_number")

    int_spliced = novel_last_exons_nc.overlap(ref_internal_introns,
                                              strandedness="same",
                                              how="containment")

    int_spliced_pre = set(int_spliced.as_df()["transcript_id"].tolist())

    eprint("before checking for overlap with ref exons" +
           ", number of putative novel internal spliced last exons" +
           " is {}".format(len(int_spliced_pre)))

    int_spliced_post = (set(int_spliced.overlap(ref_exons,
                                             strandedness="same",
                                             invert=True)
                                    .as_df()
                                    ["transcript_id"].tolist()
                            )
                        )

    eprint("after checking for overlap with ref exons" +
           ", number of putative novel internal spliced last exons" +
           " is {}".format(len(int_spliced_post))
           )

    # Get set of intron contained but exon overlapping tx ids
    int_spliced_ex_olap = int_spliced_pre.difference(int_spliced_post)

    # Update chain_match_info with isoform_class column
    conditions = [(chain_match_info["match_class"] == "valid") &
                  (chain_match_info["transcript_id_novel"].isin(int_spliced_post)),
                  (chain_match_info["match_class"] == "valid") &
                  (chain_match_info["transcript_id_novel"].isin(int_spliced_ex_olap))
                 ]

    choices = [class_key, class_key + "_exon_overlap"]

    chain_match_info[class_col] = np.select(conditions,
                                            choices,
                                            default=chain_match_info["isoform_class"])

    # exon_overlap events will no longer be considered valid
    # These could be last exons bleeding into internal exons,
    # where StringTie is unlikely to define 3'end robustly

    chain_match_info["match_class"] = np.where(chain_match_info[class_col] == class_key + "_exon_overlap",
                                               "not_valid",
                                               chain_match_info["match_class"]
                                               )

    return chain_match_info


def main(novel_path, ref_path,
         match_by, max_terminal_non_match,
         min_extension_length, out_prefix,
         novel_source, ref_anno_type,
         nb_cpu):
    '''
    '''
    start = timer()

    assert ref_anno_type in ["gencode", "ensembl"]

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

    # Define expected name of column/attribute defining the 'gene type'
    if ref_anno_type == "gencode":
        gene_type_col = "gene_type"

    elif ref_anno_type == "ensembl":
        gene_type_col = "gene_biotype"

    start2 = timer()

    if gene_type_col not in ref.columns.tolist():
        raise ValueError("Expected gene type attribute - {} - is not present in input reference GTF. Is reference source ('-s', '--reference-source') defined correctly?".format(gene_type_col))

    else:
        ref_pc = ref.subset(lambda df: df[gene_type_col].isin(["protein_coding", "lncRNA"]), nb_cpu=nb_cpu)

    end2 = timer()

    eprint("took {} s".format(end2 - start2))

    eprint("extracting novel exons from input GTF file...")

    start4 = timer()
    novel_exons = novel.subset(lambda df: df["Feature"] == "exon", nb_cpu=nb_cpu)
    end4 = timer()

    eprint("took {} s".format(end4 - start4))

    eprint("extracting first and last exons from input GTF file...")

    s3 = timer()

    if "exon_number" not in novel_exons.columns.tolist():
        eprint("Input GTF file does not contain 'exon' number column. Adding a 'strand-aware' exon number...")

        novel_exons = add_region_number(novel_exons,
                                        feature_key="exon",
                                        out_col="exon_number")

        # Now novel_source is irrelevant...
        novel_first_exons = get_terminal_regions(novel_exons,
                                                 region_number_col="exon_number",
                                                 source=None,
                                                 filter_single=True,
                                                 which_region="first"
                                                 )

        novel_last_exons = get_terminal_regions(novel_exons,
                                                region_number_col="exon_number",
                                                source=None,
                                                filter_single=True,
                                                which_region="last"
                                                )
    else:
        # Already has exon_number column no need to add it...
        novel_first_exons = get_terminal_regions(novel_exons,
                                                region_number_col="exon_number",
                                                source=novel_source,
                                                filter_single=True,
                                                which_region="first"
                                                 )

        novel_last_exons = get_terminal_regions(novel_exons,
                                                region_number_col="exon_number",
                                                source=novel_source,
                                                filter_single=True,
                                                which_region="last"
                                                )

    e3 = timer()
    eprint("extracting input first and last exons - took {} s".format(e3 - s3))


    eprint("extracting exons of protein_coding and lncRNA genes from reference annotation...")

    start5 = timer()
    ref_pc_exons = ref.subset(lambda df: df["Feature"] == "exon", nb_cpu=nb_cpu)
    end5 = timer()

    eprint("took {} s".format(end5 - start5))

    eprint("extracting first and last exons from reference annotation...")

    s4 = timer()
    ref_pc_first_exons = get_terminal_regions(ref_pc_exons,
                                              region_number_col="exon_number",
                                              source=None,
                                              filter_single=True,
                                              which_region="first"
                                              )

    ref_pc_last_exons = get_terminal_regions(ref_pc_exons,
                                             region_number_col="exon_number",
                                             source=None,
                                             filter_single=True,
                                             which_region="last"
                                             )

    e4 = timer()
    eprint("extracting reference first and last exons - took {} s".format(e4 - s4))

    eprint("finding introns for each reference transcript...")

    start3 = timer()

    try:
        ref_pc_introns = ref_pc.features.introns(by="transcript", nb_cpu=1)
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
        # Harcode nb_cpu to 1 to avoid 'KeyError' when invoke Ray for 2nd time
        # https://github.com/biocore-ntnu/pyranges/issues/201
        # https://github.com/biocore-ntnu/pyranges/issues/200
        # To do: try and downgrade Ray so this is avoided?
        novel_introns = novel.features.introns(by="transcript", nb_cpu=1)

    except KeyError:
        # Specific error with ray when execute introns func for the 2nd time in same script...
        # KeyError: 'by'
        # Whilst avoiding working out what's going on, I can run my super slow intron finding script...
        eprint("pr.features.introns returned KeyError ('by'), using my hacky intron finding workaround...")
        novel_introns = introns_by_tx(novel_exons, nb_cpu=nb_cpu)

    end1 = timer()

    eprint("took {} s".format(end1 - start1))

    # Need/want to make sure each introns object has a strand-aware intron number by transcript_id (i.e. 1 = first intron)
    s_in = timer()

    try:
        assert "intron_number" in ref_pc_introns.as_df().columns.tolist()
    except AssertionError:
        eprint("adding intron_number column to reference introns object...")
        ref_pc_introns = add_region_number(ref_pc_introns, nb_cpu=1)


    try:
        assert "intron_number" in novel_introns.as_df().columns.tolist()

    except AssertionError:
        eprint("adding intron_number column to novel introns object...")
        novel_introns = add_region_number(novel_introns, nb_cpu=1)

    e_in = timer()

    eprint("Adding intron number columns took {} s".format(e_in - s_in))

    eprint("Extracting reference first introns and novel last introns...")

    s5 = timer()
    ref_pc_first_introns = get_terminal_regions(ref_pc_introns,
                                                region_number_col="intron_number",
                                                feature_key="intron",
                                                source=None,
                                                filter_single=True,
                                                which_region="first"
                                                )

    novel_last_introns = get_terminal_regions(novel_introns,
                                              region_number_col="intron_number",
                                              feature_key="intron",
                                              source=None,
                                              filter_single=True,
                                              which_region="last"
                                              )

    e5 = timer()
    eprint("extracting reference first and novel last introns - took {} s".format(e5 - s5))

    eprint("finding novel transcripts with valid matches in their intron chain to reference transcripts...")

    start6 = timer()
    valid_matches = filter_transcripts_by_chain(novel_introns,
                                                ref_pc_introns,
                                                match_type=match_by,
                                                max_terminal_non_match=max_terminal_non_match,
                                                nb_cpu=nb_cpu)
    end6 = timer()

    eprint("took {} s".format(end6 - start6))

    # eprint("this is df with match info for each novel transcript...")
    # eprint(valid_matches)
    # eprint(valid_matches.dtypes)

    eprint("finding novel last exons completely contained within annotated first introns...")
    s7 = timer()


    fi_valid_matches = filter_first_intron_tx(novel_first_exons,
                                              novel_last_exons,
                                              ref_pc_first_exons,
                                              ref_pc_first_introns,
                                              valid_matches,
                                              novel_source,
                                              nb_cpu)

    e7 = timer()
    eprint("finding first intron novel isoforms - took {} s".format(e7 - s7))

    eprint("finding bleedthrough and 3'UTR extension events from transcripts with exact intron chain matches...")
    s8 = timer()

    # eprint(fi_valid_matches)

    bl_utr_valid_matches = filter_complete_match(novel_last_exons,
                                                 ref_pc_first_exons,
                                                 ref_pc_last_exons,
                                                 ref_pc_exons,
                                                 ref_pc_introns,
                                                 min_extension_length,
                                                 fi_valid_matches,
                                                 novel_source,
                                                 nb_cpu)

    e8 = timer()
    eprint("finding bleedthrough and 3'UTR extension events - took {} s".format(e8 - s8))


    eprint("classifying remaining valid isoforms...")
    # Will have some downstream alternatively spliced last exons, internal alt spliced and spliced 3'UTR introns flopping about
    ui3_bl_utr_valid_matches = annotate_3utr_introns(novel_last_introns,
                                                     ref_pc_last_exons,
                                                     bl_utr_valid_matches,
                                                     nb_cpu=nb_cpu)

    ui3_bl_utr_valid_matches = annotate_distal_last_exons(novel_last_introns,
                                                          ref_pc_last_exons,
                                                          ui3_bl_utr_valid_matches,
                                                          nb_cpu=nb_cpu)

    ui3_bl_utr_valid_matches = annotate_internal_spliced(novel_last_exons,
                                                         ref_pc_introns,
                                                         ref_pc_exons,
                                                         ui3_bl_utr_valid_matches,
                                                         nb_cpu=nb_cpu)

    # Have made all specific classifications and filtered where possible
    # Remaining valid isoforms likely have complex structures
    ui3_bl_utr_valid_matches["isoform_class"] = np.where((ui3_bl_utr_valid_matches["match_class"] == "valid") &
                                                       (pd.isna(ui3_bl_utr_valid_matches["isoform_class"])),
                                                        "other", ui3_bl_utr_valid_matches["isoform_class"])


    if isinstance(ui3_bl_utr_valid_matches, pd.Series):
        # match type was any
        valid_novel = novel.subset(lambda df: df["transcript_id"].isin(set(ui3_bl_utr_valid_matches.tolist())), nb_cpu=nb_cpu)
        valid_novel.to_gtf(out_prefix + ".gtf")

    elif isinstance(ui3_bl_utr_valid_matches, pd.DataFrame):
        # match_by/match_type was transcript
        valid_novel = novel.subset(lambda df: df["transcript_id"].isin(set(ui3_bl_utr_valid_matches.loc[ui3_bl_utr_valid_matches["match_class"] == "valid", "transcript_id_novel"].tolist())), nb_cpu=nb_cpu)
        valid_novel.to_gtf(out_prefix + ".gtf")

        summary_counts = (ui3_bl_utr_valid_matches.loc[lambda x: x["match_class"] == "valid", :]
                          .drop_duplicates(subset=["transcript_id_novel"])
                          ["isoform_class"]
                          .value_counts(dropna=False)
                          )

        summary_counts_nv = (ui3_bl_utr_valid_matches.loc[lambda x: x["match_class"] == "not_valid", :]
                          .drop_duplicates(subset=["transcript_id_novel"])
                          ["isoform_class"]
                          .value_counts(dropna=False)
                          )
        # eprint(summary_counts)

        pd.DataFrame(summary_counts).reset_index().to_csv(out_prefix + ".valid.class_summary_counts.tsv",
                                                          sep="\t",
                                                          header=["isoform_class","count"],
                                                          index=False,
                                                          na_rep="NA")

        pd.DataFrame(summary_counts_nv).reset_index().to_csv(out_prefix + ".not_valid.class_summary_counts.tsv",
                                                             sep="\t",
                                                             header=["isoform_class","count"],
                                                             index=False,
                                                             na_rep="NA")

        ui3_bl_utr_valid_matches.to_csv(out_prefix + ".match_stats.tsv", sep="\t", header=True, index=False,na_rep="NA")


    end = timer()

    eprint("Completed: script took {} s / {} min (3dp) ".format(round(end - start, 3), round((end - start) / 60, 3)))





if __name__ == '__main__':

    descrpn = """Script to filter assembled transcripts for intron chain matches with reference transcripts
                 to identify novel 3'end transcript isoforms"""

    parser = argparse.ArgumentParser(description=descrpn)

    parser.add_argument("-i", "--input-transcripts",
                        default='', dest="novel_gtf",
                        help = "path to GTF file containing novel transcripts assembled by StringTie.",
                        required=True)
    parser.add_argument("-r", "--reference-transcripts",
                        default='', type=str,
                        dest="ref_gtf",
                        help="path to GTF file containing reference transcripts against which to match intron chains of novel transcripts. Should contain same chromosome naming scheme",
                        required=True)
    parser.add_argument("-s", "--reference-source",
                        type=str, default="gencode",
                        choices=["gencode", "ensembl"], dest="ref_anno_type",
                        help="Whether reference annotation is sourced from Gencode 'gencode' or Ensembl. This is so can use correct attribute key to extract protein-coding and lncRNA genes (default: %(default)s)")
    parser.add_argument("--input-exon-number-format", default="stringtie",
                        choices=["stringtie","strand_aware"], dest="novel_exon_n_fmt",
                        help="Are 'exon numbers' in input transcripts assigned 1..n leftmost-rightmost ignoring strand (StringTie's convention, 'stringtie') or strand aware (Gencode & Ensembl annotation convention, 'strand_aware') (default: %(default)s)")
    parser.add_argument("-m", "--match-by",
                        default="any", type=str,
                        choices=["any", "transcript"], dest="match_by",
                        help="Consider novel transcript a valid match if all but penultimate intron(s) match introns of any transcript ('any') or the same transcript ('transcript'). 'transcript' is CURRENTLY UNIMPLEMENTED (default: %(default)s)")
    parser.add_argument("-n", "--max-terminal-non-match",
                        default=1, type=int,
                        dest="max_terminal_non_match",
                        help="Maximum number of uninterrupted reference-unmatched introns at 3'end of novel transcript for it to be considered a valid match (default: %(default)s)")
    parser.add_argument("-e", "--min-extension-length",
                        default=1, type=int,
                        dest="min_extension_length",
                        help="Minimum 3'end extension length relative to overlapping isoform for an 'extension' event to be retained (default: %(default)s)")
    parser.add_argument("-c", "--cores",
                        default=1, type=int,
                        help="number of cpus/threads for parallel processing (default: %(default)s)")
    parser.add_argument("-o", "--output-prefix",
                        type=str, default="intron_chain_matched_transcripts",
                        dest="output_prefix",
                        help="Prefix for output files (GTF with valid matches, matching stats TSV etc.). '.<suffix>' added depending on output file type (default: %(default)s)")

    if len(sys.argv) == 1:
        parser.print_help()
        parser.exit()

    args = parser.parse_args()

    if args.novel_exon_n_fmt == "strand_aware":
        exon_n_type = None
    else:
        exon_n_type = args.novel_exon_n_fmt

    main(args.novel_gtf, args.ref_gtf, args.match_by, args.max_terminal_non_match, args.min_extension_length, args.output_prefix, exon_n_type, args.ref_anno_type, args.cores)
