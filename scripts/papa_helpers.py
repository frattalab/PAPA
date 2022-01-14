##!/usr/bin/env python3

from __future__ import print_function
import pyranges as pr
import pandas as pd
import numpy as np
import os
import sys


'''
Collection of utility functions shared across pipeline scripts, mainly for working with PyRanges objects
'''

def eprint(*args, **kwargs):
    '''
    Nice lightweight function to print to STDERR (saves typing, I'm lazy)
    Credit: https://stackoverflow.com/questions/5574702/how-to-print-to-stderr-in-python (MarcH)
    '''
    print(*args, file=sys.stderr, **kwargs)


def _n_ids(gr, id_col):

    assert id_col in gr.columns

    return len(set(gr.as_df()[id_col]))


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


def sort_introns_by_strand(df):
    '''
    '''
    # first reset_index call removes the original index of the group (e.g. row 4005 in df)
    # second reset_index call adds the sorted index as a column to the dataframe (the order along exon in each transcript)
    if (df.Strand == '+').all():
        return df.sort_values(by=['End']).reset_index(drop=True).reset_index()
    elif (df.Strand == '-').all():
        return df.sort_values(by=['Start'], ascending=False).reset_index(drop=True).reset_index()


def get_terminal_regions(gr,
                         feature_col = "Feature",
                         feature_key = "exon",
                         id_col = "transcript_id",
                         region_number_col = "exon_number",
                         source = None,
                         which_region="last",
                         filter_single = False,
                         ):
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
    try:
        mod_gr = (gr.assign(region_number_col,
                            lambda df: df[region_number_col].astype(float).astype(int),
                            nb_cpu=1)
                  )
    except KeyError:
        # Currently getting weird KeyError with assign for certain chromosome
        # Mostly non-std chrom names
        # No error if do '.<exon_number>' to assign, but this makes inflexible to colname
        # Also no error if gr -> df assign -> gr
        eprint("pr.assign returned KeyError. Converting {} to int via pandas df conversion".format(region_number_col))

        mod_gr = gr.as_df()
        mod_gr[region_number_col] = mod_gr[region_number_col].astype(float).astype(int)
        mod_gr = pr.PyRanges(mod_gr)


    # Make sure gr is sorted by transcript_id & 'region number' (ascending order so 1..n)
    mod_gr = mod_gr.apply(lambda df: df.sort_values(by=[id_col, region_number_col], ascending=True),
                          nb_cpu=1)


    # Filter out single-exon transcripts
    if filter_single:
        eprint("Filtering for multi-exon transcripts...")
        eprint("Before: {}".format(len(set(mod_gr.as_df()[id_col].tolist()))))

        # Setting to 'False' marks all duplicates as True (so keep these)
        mod_gr = mod_gr.subset(lambda df: df.duplicated(subset=[id_col], keep=False), nb_cpu=1)

        eprint("After: {}".format(len(set(mod_gr.as_df()[id_col].tolist()))))


    if source is None:
        # source = None means that 1 = first region of group regardless of strand
        # Pick last region entry by max region number for each transcript (id_col)
        # Pick first region entry by min region number for each transcript (id_col)

        # keep="last" sets last in ID to 'False' and all others true (negate to keep last only)
        # keep="first" sets first in ID to 'False'

        out_gr = mod_gr.subset(lambda df: ~(df.duplicated(subset=[id_col], keep=which_region)),
                               nb_cpu=1
                              )


    elif source == "stringtie":
        # Numbering Doesn't respect strand
        # Need to flip selecting first/last in group depending on strand
        # minus strand - pick min if Minus strand, max if plus strand

        if which_region == "first":
            # + strand - pick first in group, - strand - pick last in group

            out_gr = (mod_gr.subset(lambda df:
                                    #1. plus strand & first in group/ID
                                    (df["Strand"] == "+") & ~(df.duplicated(subset=[id_col],
                                                                            keep="first")) |
                                    #2. minus strand & last in group/ID
                                    (df["Strand"] == "-") & ~(df.duplicated(subset=[id_col],
                                                                            keep="last")),
                                    nb_cpu=1)
                     )

        elif which_region == "last":
            # + strand - pick last in group/ID
            # - strand - pick first in group/ID
            out_gr = (mod_gr.subset(lambda df:
                                    #1. plus strand & last in group/ID
                                    (df["Strand"] == "+") & ~(df.duplicated(subset=[id_col],
                                                                            keep="last")) |
                                    #2. minus strand & first in group/ID
                                    (df["Strand"] == "-") & ~(df.duplicated(subset=[id_col],
                                                                            keep="first")),
                                    nb_cpu=1)
                     )


    return out_gr


def _df_add_region_number(df,id_col,sort_col="Start"):
    '''
    Return a Series of strand-aware region numbers (5'-3' in 1..n)
    Function to be used internally in a pr.assign (mainly by add_region_number)
    '''
    if id_col not in df.columns.tolist():
        raise KeyError(f"id_col - {id_col} - is not present in df for chr/strand pair {','.join([df.Chromosome.iloc[0], df.Strand.iloc[0]])}")

    elif (df.Strand == "+").all():
        # Start position smallest to largest = 5'-3'

        return df.groupby(id_col)[sort_col].rank(method="min", ascending=True)

    elif (df.Strand == "-").all():
        # Start position largest to smallest = 5'-3'

        return df.groupby(id_col)[sort_col].rank(method="min", ascending=False)

    elif df.empty:
        eprint("df is empty - returning empty pd.Series")
        return pd.Series()


def add_region_number(gr,
                      id_col="transcript_id",
                      feature_key="intron",
                      out_col="intron_number",
                      feature_col="Feature",
                      nb_cpu=1):
    '''
    Adds column to gr containing a strand aware region number column,
    ordered 5'-3' 1..n by a group of features (e.g. transcript)
    '''

    # Make sure only 'feature_key' rows are in the gr
    assert gr.as_df()[feature_col].drop_duplicates().tolist() == [feature_key], "only {} entries should be present in gr".format(feature_key)

    # Make sure sorted by position first.
    gr = gr.sort()

    # Add in region number column in strand aware manner, so 1 = most 5', n = most 3'

    gr = gr.assign(out_col, lambda df: _df_add_region_number(df, id_col), nb_cpu=nb_cpu)

    return gr


def get_internal_regions(gr,
                         feature_col="Feature",
                         feature_key="exon",
                         id_col="transcript_id",
                         region_number_col="exon_number",
                         ):
    '''
    Return gr of internal exons for each transcript_id
    In process, exon_number_col will be converted to type 'int'
    '''

    assert gr.as_df()[feature_col].drop_duplicates().tolist() == [feature_key], "only {} entries should be present in gr".format(feature_key)


    # Pull out exons, convert exon_number to int
    exons_gr = gr.assign(region_number_col,
                         lambda df: df[region_number_col].astype(float).astype("Int64"),
                         nb_cpu=1)

    # Make sure gr is sorted by transcript_id & 'region number' (ascending order so 1..n)
    exons_gr = exons_gr.apply(lambda df: df.sort_values(by=[id_col,
                                                            region_number_col
                                                            ],
                                                        ascending=True),
                              nb_cpu=1)

    # Filter out 1st + last exons for each ID
    # first exons for each transcript (.ne(1))
    # keep="last" sets last dup value to 'False' & all others True
    # This will filter out last exons

    out_gr = (exons_gr.subset(lambda df: (df[region_number_col].ne(1).astype(bool)) &
                     (df.duplicated(subset=["transcript_id"], keep="last")),
                     nb_cpu=1
                    )
             )

    return out_gr


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
