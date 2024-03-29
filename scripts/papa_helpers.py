#!/usr/bin/env python3

#     Collection of utility functions used by scripts in the PAPA pipeline
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

from __future__ import print_function
import pyranges as pr
import pandas as pd
import numpy as np
from pyranges.pyranges import PyRanges
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


def _n_ids(gr, id_col, is_df=False):

    assert id_col in gr.columns

    if not is_df:
        return len(set(gr.as_df()[id_col]))
    else:
        return gr[id_col].nunique()


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


def _check_to_extract_order(attribute, to_extract):
    '''
    Checks & re-orders list of keys present in a string
    attribute: str of the 'attribute' field from a GTF file
    to_extract: list of keys wish to extract from attribute string

    - returns to_extract values sorted from left-right occurence in attribute
    - if values are not found in attribute they are removed from the list

    '''

    # Get indexes of each
    idxs = [attribute.find(attr) for attr in to_extract]
    idx_ser = pd.Series(idxs, index=to_extract)

    # if not found (-1), need to remove from to_extract
    idx_ser = idx_ser[idx_ser != -1]

    return idx_ser.sort_values().index.tolist()


def _fetch_attributes(attribute,
                      to_extract=["gene_id",
                                  "transcript_id",
                                  "exon_number"
                                  "gene_name"],
                      ):
    '''
    Extracts specific key-value pairs from a GTF attribute string using a regex into individual columns
    attribute: Series of 'attribute' column string read in from GTF
    to_extract: list of keys to extract from the attribute string
    '''

    # Check order of extract cols matches the order they appear in atttribute column
    # If not present in 1st row then dropped
    to_extract_v = _check_to_extract_order(attribute.iloc[0], to_extract)


    # assert to_extract[0] == "gene_id", f"gene_id must be the first attribute to extract, following was found - {to_extract[0]}"

    # Remove quotes from attributes string
    no_quotes = attribute.str.replace('"', '').str.replace("'", "")

    # Assemble regex
    # Regex completely lifted from _fetch_gene_transcript_exon_id (pr.readers._fetch_gene_transcript_exon_id)
    # Specifically add the '\s' to properly parse cases where key is further present with a suffix (e.g. 'gene_id' & 'gene_id_b' keys)
    # https://github.com/biocore-ntnu/pyranges/blob/a36a2c7fac88f297bf41c529734c2cb6950bed3b/pyranges/readers.py#L210
    # "gene_id.?(.+?);(?:.*transcript_id.?(.+?);)?(?:.*exon_number.?(.+?);)?(?:.*exon_id.?(.+?);)?"

    first_attr = to_extract_v[0] + "\s(.+?);"

    others = ["(?:.*" + attr + "\s(.+?);)?" for attr in to_extract_v[1:]]

    attr_regex = "".join(list(first_attr) + others)

    # eprint(attr_regex)

    df = no_quotes.str.extract(attr_regex,
                               expand=True)  # .iloc[:, [1, 2, 3]]

    df.columns = to_extract_v

    # TODO: Why does _fetch_gene_transcript_exon_id need 'annotation' arg for Ensembl? What is it doing?

    return df


def read_gtf_specific(f,
                      attr_to_extract=["gene_id",
                                       "transcript_id",
                                       "exon_number",
                                       "gene_name"],
                      skiprows=None,
                      as_df=False,
                      nrows=None
                      ):
    '''
    GTF reader to extract a custom set of attributes. Almost entirely mimics read_gtf_restricted, except that can provide own list of attributes to extract
    No way near as flexible or feature-rich as normal function (yet)

    Notes:
    - Since the function uses a regex to extract key-value pairs, the order of keys in attr_to_extract is very important
    - The function splits by 'feature' (i.e. 'gene', 'exon' or 'transcript') to try to account for different keys between features (e.g. gene doesn't have 'exon_number')
    - The order is checked, but only against the FIRST value in each group (of a chunk of the input file). For reference GTFs this is probably sufficient, but if you have non-standard/rare keys in attr_to_extract they will most likely be parsed incorrectly
    '''

    dtypes = {
        "Chromosome": "category",
        "Feature": "category",
        "Strand": "category"
    }


    df_iter = pd.read_csv(
        f,
        sep="\t",
        comment="#",
        usecols=[0, 2, 3, 4, 5, 6, 8],
        header=None,
        names="Chromosome Feature Start End Score Strand Attribute".split(),
        dtype=dtypes,
        chunksize=int(1e5),
        skiprows=skiprows,
        nrows=nrows)

    dfs = []
    for df in df_iter:
        if sum(df.Score == ".") == len(df):
            cols_to_concat = "Chromosome Start End Strand Feature".split()
        else:
            cols_to_concat = "Chromosome Start End Strand Feature Score".split(
            )

        # Extract specific key-value pairs from the attribute string/column
        # Since different features (e.g. 'gene', 'exon') often have different attributes
        # Provided order in attr_to_extract may not be identical (would break the regex)
        # Extracts the provided attributes (if present) whilst checking the order they are observed (hackily in the first row of the group)
        extract = df.groupby("Feature")["Attribute"].apply(_fetch_attributes, attr_to_extract)

        if "exon_number" in extract.columns:
            extract.exon_number = extract.exon_number.astype(float)

        df = pd.concat([df[cols_to_concat], extract], axis=1, sort=False)

        dfs.append(df)

    df = pd.concat(dfs, sort=False)

    df.loc[:, "Start"] = df.Start - 1

    if not as_df:
        return pr.PyRanges(df)
    else:
        return df


def check_stranded(gr):
    '''
    Validate PyRanges object as stranded, removing offending rows if applicable
    '''

    if gr.stranded:
        return gr

    else:
        df = gr.as_df()
        eprint(f"Input gr is unstranded - 'Strand' values in input gr - {df['Strand'].drop_duplicates().tolist()}")

        # Get a list of Strand values that aren't + / -
        invalid_strand = list(set(df['Strand'].drop_duplicates()) - {"+", "-"})

        df = df.loc[df["Strand"].isin(["+", "-"]), :]

        # PyRanges requires Strand col to have '+' & '-' categories only
        df["Strand"] = df["Strand"].cat.remove_categories(invalid_strand)

        gr = pr.PyRanges(df)

        assert gr.stranded

        return gr


def check_concat(gr):
    '''
    Validate that underlying dfs of (concatenated) PyRanges object have same number of columns

    If dfs do not have same shape, methods like gr.assign will fail due to inconsistent index to add column
    This heuristic will find the first chr/str with the most columns
    and update all other chr/str dfs with missing cols from the max filled wih NaNs (making up to same number of cols)
    '''

    col_lens = {chr_str: len(df.columns) for chr_str, df in gr}

    if len(set(col_lens.values())) == 1:
        # All columns are of same length
        return gr

    else:
        # Some dfs have missing columns
        # Want to make sure each df has same number (and labels) of columns
        default_cols = ["Chromosome", "Start", "End", "Strand"]

        # Get all metadata cols from every df in gr
        # Will be lots of duplicates in here, but this way will avoid massively changing order of columns (from the first df in gr)
        extra_cols = [col for df in gr.values() for col in df.columns if col not in default_cols]

        # Generate list of target cols without duplicates (preserving order of appearance)
        all_cols = list(dict.fromkeys(default_cols + extra_cols))

        # Update all dfs so they have same columns (and order)
        # Col will be filled with NaNs if df didn't have before
        out_gr = gr.apply(lambda df: df.reindex(columns=all_cols))

        return out_gr


def _df_collapse_metadata(df, id_col, standard_cols, collapse_cols, collapse_uniq_cols, collapse_sep):
    '''
    Intended to be applied to internal dfs of PyRanges objects
    '''

    found_collapsed = [col for col in collapse_cols if col in df.columns]

    not_found_collapsed = set(collapse_cols) - set(found_collapsed)

    if len(not_found_collapsed) > 0:
        chr_strand = f"{df.Chromosome.drop_duplicates()[0]},{df.Strand.drop_duplicates()[0]}"
        eprint(f"following 'collapse_cols' columns not found in df (chr/strand) - {chr_strand} - {', '.join(not_found_collapsed)}")

    grouped = df.groupby(id_col)

    # Pick first entry for all standard_cols, these should be same for all rows of id_col
    # Leaves a df with id_col values as index
    std_collapsed = grouped[standard_cols].first()

    # For collapse cols, collapse to single row of delimited strings for each column
    # Again leave a df with id_col values as index labels
    clp_collapsed = grouped[found_collapsed].agg(lambda col: collapse_sep.join(col.astype(str)))

    if collapse_uniq_cols is not None:
        # Collapse these cols to single row of delimited strings whilst dropping duplicates
        # Again leave a df with id_col values as index labels
        clp_uniq_collapsed = grouped[collapse_uniq_cols].agg(lambda col: collapse_sep.join(list(dict.fromkeys(col.astype(str)))))

        int_collapsed = clp_collapsed.merge(clp_uniq_collapsed, left_index=True, right_index=True)

        collapsed = std_collapsed.merge(int_collapsed, left_index=True, right_index=True)

    else:
        # combine by id_col
        collapsed = std_collapsed.merge(clp_collapsed, left_index=True, right_index=True)

    return collapsed


def collapse_metadata(gr,
                      id_col="transcript_id",
                      standard_cols=["Chromosome", "Start", "End", "Strand"],
                      collapse_cols=None,
                      collapse_uniq_cols=None,
                      collapse_sep=","):
    '''
    Collapse to a single entry/row per ID entry whilst retaining/collapsing metadata on duplicate rows
    standard_cols: list of column labels that have the same value for all entries of id_col and do not need to be collapsed.
    This is essential for PyRanges standard columns, as you do not want to be changing their dtypes to string. All columns labels in this list retain their dtype, and the first value is retained

    collapse_cols: list of column labels containing metadata you'd like to collapse to a single row (separated by collapse_sep)
        If None, then all columns in gr except for standard_cols, id_col & collapse_uniq_cols will be collapsed
    collapse_uniq_cols: list of column labels containing metadata you'd like to collapse to a single row whilst dropping duplicate values. Values will maintain order of appearance in df
    '''

    assert all([True if col in gr.columns else False for col in standard_cols])

    if collapse_uniq_cols is not None:
        # Here just checking the columns are found in df
        assert all([True if col in gr.columns else False for col in collapse_uniq_cols])

    if collapse_cols is None:
        if collapse_uniq_cols is not None:
            def_cols = standard_cols + [id_col] + collapse_uniq_cols

        else:
            def_cols = standard_cols + [id_col]

        collapse_cols = [col for col in gr.columns if col not in def_cols]

    else:
        assert all([True if col in gr.columns else False for col in collapse_cols])


    return gr.apply(lambda df: _df_collapse_metadata(df,
                                                     id_col,
                                                     standard_cols,
                                                     collapse_cols,
                                                     collapse_uniq_cols,
                                                     collapse_sep
                                                     )
                    )
