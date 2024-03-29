#!/usr/bin/env python3

#     Combine reference last exons and novel last exons into a single Gene Transfer Format file (for quantification with Salmon)
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
from pyranges.readers import read_gtf_restricted
from papa_helpers import (
    eprint,
    get_terminal_regions,
    get_internal_regions,
    add_region_number,
    _df_add_region_number,
    _pd_merge_gr,
    read_gtf_specific,
    check_concat,
    collapse_metadata,
)
from get_novel_last_exons import (
    find_extension_events,
    find_spliced_events,
    _filter_gr_for_not_tx,
    add_region_rank, no_overlap_3p_ids
)
import sys
import argparse
from timeit import default_timer as timer


"""
Generate a GTF file of reference and novel last exons (or unique last exon segments for first/internal extensions)
Also annotate transcript IDs according to

"""


ref_attr_key_order = ["gene_id", "transcript_id", "gene_name", "exon_number"]
ref_attr_key_order_n_en = ["gene_id", "transcript_id", "gene_name"]


def cluster_to_region_number(
    gr, group_id_col, out_col="le_number", cluster_col="Cluster"
):
    """
    Returns gr with 'out_col' column added
    where out_col is leftmost to rightmost cluster_col converted to a
    strand-aware 1..n order by group_id_col
    1 = most 5' site in group_id_col
    n = most 3' site in group_id_col
    """

    # For groupby.rank to work appropriately, need a single row per 'Cluster'
    c2p = gr[[group_id_col, cluster_col]].apply(
        lambda df: df.drop_duplicates(subset=cluster_col)
    )

    # Add 1..n 5'-3' region number as a column
    c2p = c2p.assign(
        out_col, lambda df: _df_add_region_number(df, group_id_col, cluster_col)
    )

    # Return 'out_col' to original gr
    c2p_cols = c2p.columns.tolist()
    out_gr = gr.apply_pair(
        c2p,
        lambda df, df_to_merge: _pd_merge_gr(
            df,
            df_to_merge,
            how="left",
            on=cluster_col,
            suffixes=[None, "_match"],
            to_merge_cols=c2p_cols,
        ),
    )

    # avoid adding extra 'PyRanges' cols (Chromosome etc.) from c2p
    return out_gr.drop(like="_match$")


def _df_add_common_gene_id_col(df, gene_col, novel_id_str, novel_ref_gene_col):
    """
    Internal to add_common_gene_id_col
    returns pd.Series
    """

    if novel_ref_gene_col not in df.columns:
        # This means chr/strand pair df has no novel entries
        # gene_col = reference gene col
        return df[gene_col]

    else:
        return pd.Series(
            np.where(
                df[gene_col].str.contains(novel_id_str, regex=False),
                df[novel_ref_gene_col],  # olapping ref gene ID for novel
                df[gene_col],  # is a ref gene/last exon
            )
        )


def add_common_gene_id_col(
    gr,
    out_col="ref_gene_id",
    id_col1="gene_id",
    id_col1_str="PAPA",
    id_col2="gene_id_ref",
):

    """
    Add a common identifier column between two related columns based on value in one column.
    If values of the first column ('id_col1') contains a given string ('id_col1_str', exact match no regex), report the value of 'id_col2' for that row
    Otherwise use the value from 'id_col1'

    The typical use case is output of e.g. get_novel_last_exons.py / filter_tx_by_three_end.py, where a novel assembled event will have the overlapping reference gene ID/name in attribute column under 'gene_id_ref'
    An annotation GTF will typically have the same attribute stored under 'gene_id' key.
    This function will unify the values into a single new column, using the value of 'gene_id_ref' where a novel event is encountered, otherwise the value in 'gene_id'
    """

    return gr.assign(
        out_col,
        lambda df: _df_add_common_gene_id_col(df, id_col1, id_col1_str, id_col2),
    )


def _df_add_common_gene_name_col(
    df, col_to_check, col_to_check_str, id_col_t, id_col_f
):
    """
    Internal to add_common_gene_name_col
    returns pd.Series
    """

    if id_col_t not in df.columns:
        # This means chr/strand pair df has no novel entries
        # gene_col = reference gene col
        return df[id_col_f]

    else:
        return pd.Series(
            np.where(
                df[col_to_check].str.contains(col_to_check_str, regex=False),
                df[id_col_t],  # olapping ref gene name for novel
                df[id_col_f],  # is a ref gene/last exon (use ref gene_name)
            )
        )


def add_common_gene_name_col(
    gr,
    out_col="ref_gene_name",
    col_to_check="gene_id",
    col_to_check_str="PAPA",
    id_col_t="gene_name_ref",
    id_col_f="gene_name",
):

    """
    Add a common identifier column between two related columns based on value in one column.
    If values of the first column ('id_col1') contains a given string ('id_col1_str', exact match no regex), report the value of 'id_col2' for that row
    Otherwise use the value from 'id_col1'

    The typical use case is output of e.g. get_novel_last_exons.py / filter_tx_by_three_end.py, where a novel assembled event will have the overlapping reference gene ID/name in attribute column under 'gene_id_ref'
    An annotation GTF will typically have the same attribute stored under 'gene_id' key.
    This function will unify the values into a single new column, using the value of 'gene_id_ref' where a novel event is encountered, otherwise the value in 'gene_id'
    """

    return gr.assign(
        out_col,
        lambda df: _df_add_common_gene_name_col(
            df, col_to_check, col_to_check_str, id_col_t, id_col_f
        ),
    )


def _df_grp_update_le_number(
    df, le_id_col="le_id", le_num_col="le_number", is_ext_col="is_extension"
):
    """
    # TODO: update to cumulative sum based approach
    # https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.Series.cumsum.html?highlight=cumsum
    # https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.core.groupby.DataFrameGroupBy.cumsum.html?highlight=cumsum
    - Basically be, cumulative sum of 'is_extension' (modified so same le_id = 1 extension only (1st))
    - Should have same index, so add df.groupby(id_col).apply(lambda df: df[le_number].add(cumsum))
    https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.Series.add.html?highlight=add#pandas.Series.add

    Returns: pd.Series
    """

    # Track how many extension IDs for gene starting from 5'end
    # Once encounter an ext, all IDs downstream need to be shifted down by n of preceeding extensions
    ext_count = 0

    # Track le_id of most recently found ext event
    # If gene has alt extension isoforms of the same last exon
    # Will group together for simplification's sake
    prev_ext_le_id = ""

    out = []

    for _, row in df.iterrows():
        if row[is_ext_col] == 1:
            # Is an extension
            # Check if same LE as before
            if row[le_id_col] == prev_ext_le_id:
                # Tx is a diff extension of same le - group together
                out.append(row[le_num_col] + ext_count)

            else:
                # New distinct extension event
                ext_count += 1
                out.append(row[le_num_col] + ext_count)

                # Track le_id in case another for the same last exon
                prev_ext_le_id = row[le_id_col]

        elif row[is_ext_col] == 0:
            # Not an extension
            # If no previous extension, le number is unchaged
            # Otherwise shifted by n extensions
            out.append(row[le_num_col] + ext_count)

    return pd.Series(out, index=df.index)


def update_extension_le_number(df, id_col="ref_gene_id", out_col="le_number_ext"):
    """ """

    # assert

    df = df.groupby(id_col).apply(
        lambda grp: grp.assign(**{out_col: _df_grp_update_le_number(grp)})
    )

    return df


def update_ext_le_ids(
    gr,
    type_col="event_type",
    le_ext_key="last_exon_extension",
    id_col="ref_gene_id",
    number_col="le_number",
    le_id_col="le_id",
):
    """ """

    gr = gr.assign(
        "is_extension",
        lambda df: pd.Series(
            np.where(df[type_col].str.contains(le_ext_key, na=False, regex=False), 1, 0)
        ),
    )

    # Extra sort so extension events come after non-extensions of the same le_id
    gr = gr.apply(lambda df: df.sort_values(by=[id_col, number_col, "is_extension"]))

    # Adds column 'le_number_ext'
    gr = gr.apply(lambda df: update_extension_le_number(df, id_col=id_col))

    # eprint(combined_ext[["ref_gene_id", "le_id", "le_number", "event_type", "is_extension", "le_number_ext"]].print(n=50))

    # Reassign le_id with updated number, return to same shape as no extensions gr
    gr = (
        gr.drop([le_id_col, number_col, "is_extension"])
        .apply(lambda df: df.rename(columns={"le_number_ext": number_col}))
        .assign(
            le_id_col,
            lambda df: df[id_col] + "_" + df[number_col].astype(int).astype(str),
        )
    )

    return gr


def annotate_le_ids(
    novel_le,
    ref_le,
    novel_id_col="gene_id_ref",
    ref_id_col="gene_id",
    le_id_outcol="le_id",
    ref_extensions=False,
):
    """ """

    # Identify last exon Ids containing novel extensions - these need to be handled separately
    # set of cluster IDs for each chr, strand tuple
    # {(chr, strand): {cluster_IDs}}
    d_ext_gene_ids = novel_le.apply(
        lambda df: set(
            df[df["event_type"].str.contains("extension", regex=False)][novel_id_col]
        ),  # should be ref gene ID in novel events (foun)
        as_pyranges=False,
    )

    # Combine into a single set
    # https://blog.finxter.com/union-multiple-sets-in-python/
    ext_gene_ids = set().union(*d_ext_gene_ids.values())

    if ref_extensions:
        # Want to also consider specific reference events as if novel extensions
        # {(chr, strand): {gene_id1, gene_id2}}
        ref_ext_gene_ids = ref_le.apply(
            lambda df: set(
                df[df["event_type"].str.contains("extension", regex=False)][ref_id_col]
            ),
            as_pyranges=False,
        )

        n_ref_ext_ids = sum(len(ids_set) for ids_set in ref_ext_gene_ids.values())
        eprint(
            f"Number of genes containing extension events sourced from labelled annotated/reference transcripts - {n_ref_ext_ids}"
        )

        # Add gene_ids to existing list
        ext_gene_ids = ext_gene_ids.union(*ref_ext_gene_ids.values())

    # Separate novel & ref grs into genes with novel extensions ('<ref/novel>_le_ext') & those without ('<ref/novel>_le_n_ext')

    ref_le_ext = ref_le.subset(lambda df: df[ref_id_col].isin(ext_gene_ids))
    ref_le_n_ext = ref_le.subset(lambda df: ~(df[ref_id_col].isin(ext_gene_ids)))

    novel_le_ext = novel_le.subset(lambda df: df[novel_id_col].isin(ext_gene_ids))
    novel_le_n_ext = novel_le.subset(lambda df: ~(df[novel_id_col].isin(ext_gene_ids)))

    # Make combined GTF of reference and novel last exons
    # Keeping extensions separate
    combined_ext = pr.concat([ref_le_ext, novel_le_ext])
    combined_n_ext = pr.concat([novel_le_n_ext, ref_le_n_ext])

    # Make sure all dfs of gr have same columns (number & labels)
    # This *should* be the case, but if chr/strand df is unique to one of the concatenated grs then can see different num of cols
    combined_ext = check_concat(combined_ext)
    combined_n_ext = check_concat(combined_n_ext)

    # make sure gene_id corresponds to reference gene ID
    # Otherwise novel + ref of same gene will be considered separately (annotated in diff column)
    # New common col = "ref_gene_id"

    if len(combined_ext) > 0:

        combined_ext = add_common_gene_id_col(combined_ext)
        combined_ext = add_common_gene_name_col(combined_ext)
        combined_ext = combined_ext.cluster(by="ref_gene_id")

        # Assign 5'-3' 1..n 'last_exon number' for each gene
        # Group together overlapping exons with a common identifier
        # .cluster(strand=None) groups as ('I') expect i.e. only overlapping intervals on the same strand can be merged
        combined_ext = cluster_to_region_number(combined_ext, "ref_gene_id")

        eprint("assigning 'le_id' (last exon ID) for each gene...")

        combined_ext = combined_ext.assign(
            le_id_outcol,
            lambda df: df["ref_gene_id"]
            + "_"
            + df["le_number"].astype(int).astype(str),
        )

        # No extensions last exon number is ready to go
        # Extension events - need to update le number so extension = own event (for last_extensions)
        combined_ext = update_ext_le_ids(
            combined_ext,
            type_col="event_type",
            le_ext_key="last_exon_extension",
            id_col="ref_gene_id",
            le_id_col=le_id_outcol,
        )

    combined_n_ext = add_common_gene_id_col(combined_n_ext)
    combined_n_ext = add_common_gene_name_col(combined_n_ext)

    combined_n_ext = combined_n_ext.cluster(by="ref_gene_id")

    combined_n_ext = cluster_to_region_number(combined_n_ext, "ref_gene_id")

    eprint("assigning 'le_id' (last exon ID) for each gene...")

    combined_n_ext = combined_n_ext.assign(
        "le_id",
        lambda df: df["ref_gene_id"] + "_" + df["le_number"].astype(int).astype(str),
    )

    if len(combined_ext) > 0:
        # GTF containing defined last exons (with le_id etc. defined)
        combined = pr.concat([combined_ext, combined_n_ext])

    else:
        combined = combined_n_ext

    return combined


def annotate_ref_event_types(
    le: pr.PyRanges,
    li: pr.PyRanges,
    ref_nl_exons: pr.PyRanges,
    ref_nl_introns: pr.PyRanges,
    suffix="_ref"
):
    '''_summary_

    _extended_summary_

    Parameters
    ----------
    le : pr.PyRanges
        _description_
    li : pr.PyRanges
        _description_
    ref_nl_exons : pr.PyRanges
        _description_
    ref_nl_introns : pr.PyRanges
        _description_
    suffix : str, optional
        _description_, by default "_ref"

    Returns
    -------
    _type_
        _description_
    '''

    def_rank_col = "region_rank"

    # must be present otherwise event type definition fails
    assert def_rank_col in ref_nl_exons.columns
    if def_rank_col in le.columns:
        rank_col_exons = def_rank_col + suffix # will be suffixed when join rows
    else:
        rank_col_exons = def_rank_col

    assert def_rank_col in ref_nl_introns.columns
    if def_rank_col in li.columns:
        rank_col_introns = def_rank_col + suffix # will be suffixed when join rows
    else:
        rank_col_introns = def_rank_col

    # First, check that 3'ends of last exons do not overlap with first/internal exons
    # (Note: these would be removed downstream anyway)
    no_nl_exon_overlap_ids = no_overlap_3p_ids(le, ref_nl_exons)
    eprint(f"Number of last exons without 3'ends internal to first/internal exons - {len(no_nl_exon_overlap_ids)}")

    le = le.subset(lambda df: df["transcript_id"].isin(no_nl_exon_overlap_ids))
    li = li.subset(lambda df: df["transcript_id"].isin(no_nl_exon_overlap_ids))

    # First identify extensions - for ref events these will just be first/internal
    # Don't care about checking 5' ends of last exons as trusting reference
    extensions = find_extension_events(
        le,
        ref_nl_exons,
        min_extension_length=1,  # min_extension_length - if completely contained will be removed later
        first_5p_tolerance=0,  # dummy value  as not checking
        other_5p_tolerance=0,  # dummy value as not checking
        rank_col=rank_col_exons, # le & ref_le contain this col so will be suffixed when join
        tolerance_filter=False,
        return_filtered_ids=False,
        suffix=suffix
    )

    extensions = extensions.drop(rank_col_exons)

    # Identify spliced events (in theory this should be all other events, just want to annotate their relative location in gene)
    # First filter novel events for non-extensions

    li_n_ext = li.apply_pair(
        extensions, lambda df, df2: _filter_gr_for_not_tx(df, df2), as_pyranges=True
    )

    # finds and annotates LIs that are internal to the gene (i.e. 5'ss matches a first/internal SJ of another tx)
    spliced_nl = find_spliced_events(
        li_n_ext, ref_nl_introns, tolerance_filter=True, rank_col=rank_col_introns
    )  # will check for 5'ss matches as just reporting alt SJs for ref le (prevents overloading)

    # Everything else remaining will be 'last exon spliced'
    # Want to keep same metadata as first/internal, so will re-run against self
    li_n_ext_les = li_n_ext.subset(lambda df: ~df["transcript_id"].isin(set(spliced_nl.transcript_id)))

    spliced_l = find_spliced_events(
        li_n_ext_les, li_n_ext_les, tolerance_filter=True, rank_col=rank_col_introns
    )

    spliced = pr.concat([spliced_nl, spliced_l])

    # Spliced is last introns, want corresponding last exons in output
    spliced_le = le.subset(
        lambda df: df["transcript_id"].isin(set(spliced.transcript_id))
    )

    # Want to retain metadata from matching reference regions
    # These coordinates/metadata correspond to matching overlapping ref last intron/SJ
    spliced = spliced[
        [
            "transcript_id",
            "gene_id_ref",
            "transcript_id_ref",
            "gene_name_ref",
            "Start_ref",
            "End_ref",
            "event_type",
        ]
    ]

    spliced_cols = spliced.columns.tolist()

    spliced_le = spliced_le.apply_pair(
        spliced,
        lambda df, df2: _pd_merge_gr(
            df,
            df2.drop(columns=["Chromosome", "Start", "End", "Strand"]),
            how="left",
            on="transcript_id",
            suffixes=[None, "_spl"],
            to_merge_cols=spliced_cols,
        ),
    )

    combined = pr.concat([extensions, spliced_le])

    # Finally, collapse metadata/duplicate attribute values for each last exon (transcript ID)
    # This can occur if same last exon matches to multiple reference transcripts/junctions
    # Cols added in script that always have 1 value per last exon (& always assigned a non-NaN value)
    assigned_cols = ["event_type"]

    # columns in original input file
    query_cols = le.columns.tolist()

    cols_to_collapse = [
        col
        for col in combined.columns
        if col not in query_cols and col not in assigned_cols
    ]

    # dict of {<collapse_col>: np.nan} - nan values in these cols only are replaced (to prevent empty string "" in output GTF)
    # Doing it this way means cols with all nan values will be dropped in final GTF
    replace_dict = {col: np.nan for col in cols_to_collapse}

    combined = pr.PyRanges(combined.as_df().replace(replace_dict, "NULL"))

    combined = collapse_metadata(
        combined,
        standard_cols=query_cols,
        collapse_cols=cols_to_collapse,
        collapse_uniq_cols=assigned_cols,
    )

    return combined


def main(
    novel_le_path,
    ref_gtf_path,
    ref_attr_key_order,
    trust_exon_number_ref,
    ref_extensions_string,
    output_prefix,
):
    """ """

    eprint("Reading in reference GTF, this can take a while...")

    start = timer()

    if trust_exon_number_ref:
        ref = read_gtf_specific(ref_gtf_path, attr_to_extract=ref_attr_key_order)

    else:
        ref = read_gtf_specific(ref_gtf_path, attr_to_extract=ref_attr_key_order_n_en)

    end = timer()

    eprint(f"Reading reference GTF complete - took {end - start} s")

    eprint("Extracting last exons for each transcript in reference GTF...")

    if trust_exon_number_ref:
        ref_e = ref.subset(lambda df: df["Feature"] == "exon")
        # Annotate as first, internal or last relative to tx
        ref_e = add_region_rank(ref_e)
        ref_le = get_terminal_regions(ref)

    else:
        ref_e = ref.subset(lambda df: df["Feature"] == "exon")
        ref_e = add_region_number(ref_e, feature_key="exon", out_col="exon_number")
        # Annotate as first, internal or last relative to tx
        ref_e = add_region_rank(ref_e)
        ref_le = get_terminal_regions(ref_e)

    # Also get a gr of non-last reference exons
    ref_e_nl = pr.concat(
        [get_terminal_regions(ref_e, which_region="first"), get_internal_regions(ref_e)]
    )

    eprint("Extracting introns for each transcript in reference GTF...")

    ref_i = ref.features.introns(by="transcript")
    # Add region number and annotate as first, internal or last
    ref_i = add_region_number(ref_i, feature_key="intron", out_col="intron_number")
    ref_i = add_region_rank(ref_i, region_number_col="intron_number")

    ref_li = get_terminal_regions(
        ref_i, feature_key="intron", region_number_col="intron_number"
    )

    ref_i_nl = pr.concat(
        [get_terminal_regions(ref_i, feature_key="intron", region_number_col="intron_number", which_region="first"),
         get_internal_regions(ref_i, feature_key="intron", region_number_col="intron_number")]
        )

    # Assign event types for reference LEs
    if ref_extensions_string is not None:
        # Want to consider specific reference events as if they were extensions
        # annotate 'event_type' as 'last_exon_extension' (same value as novel events) if their 'transcript_id' value contains ref_extensions_string
        eprint(
            f"Identifying reference last exons annotated as extensions i.e. transcript_id containing {ref_extensions_string}"
        )
        ref_le_ext = ref_le.subset(
            lambda df: df["transcript_id"].str.contains(
                ref_extensions_string, regex=False
            )
        )
        ref_le_ext.event_type = "last_exon_extension"

        if len(ref_le_ext) == 0:
            raise Exception(
                f"No reference transcript_id containing provided string - {ref_extensions_string} - were found. Correct or do not pass --ref-extensions-string to avoid this error message"
            )

        n_ref_ext = ref_le_ext.event_type.loc[lambda x: x == "last_exon_extension"].size

        eprint(
            f"Number of reference transcript ids containing extensions string - {n_ref_ext}"
        )

        ref_ext = True

        ref_le_n_ext = ref_le.subset(
            lambda df: ~df["transcript_id"].str.contains(
                ref_extensions_string, regex=False
            )
        )

        # Annotate as first/internal extensions if applicable, otherwise a distinct spliced last exon
        ref_le_n_ext = annotate_ref_event_types(
            ref_le_n_ext, ref_li, ref_e_nl, ref_i_nl
        )

        ref_le = pr.concat([ref_le_ext, ref_le_n_ext])

    else:
        eprint(
            "--ref-extensions-string not provided - all reference transcripts will be 'collapsed' to same last exon ID given overlap"
        )
        ref_ext = False

        ref_le = annotate_ref_event_types(
            ref_le, ref_li, ref_e_nl, ref_i_nl
        )

    eprint("Reading in input GTF of novel last exons...")

    # novel_le = pr.read_gtf(novel_le_path)
    novel_le = read_gtf_specific(
        novel_le_path,
        [
            "gene_id",
            "transcript_id",
            "exon_number",
            "gene_id_ref",
            "gene_name_ref",
            "Start_ref",
            "End_ref",
            "event_type",
        ],
    )

    # Make sure gene_id_ref has a single value string if only 1 distinct gene ID (Otherwise collapse to non-redundant comma-separated string of gene IDs)
    # TODO: looping over Series is a dirty way to do this, find a more canonical approach
    # TODO: why do some last exons have NaNs for gene_id_ref?

    novel_le = novel_le.assign(
        "gene_id_ref",
        lambda df: pd.Series(
            [
                ",".join(list(dict.fromkeys(str(ref_id).split(",")))) # convert to string in case of NaNs
                for ref_id in df["gene_id_ref"]
            ]
        ),
    )

    # These events need to be handled separately - will create a 'metagene' combining
    n_mult = novel_le.as_df()["gene_id_ref"].str.contains(",", regex=False).sum()
    eprint(f"Number of novel events matching multiple reference genes - {n_mult}")

    if n_mult == 0:
        combined_mult = pr.PyRanges()
        combined_single = annotate_le_ids(
            novel_le, ref_le, le_id_outcol="le_id", ref_extensions=ref_ext
        )

    else:
        # Split novel_le into novel_le_single (overlaps 1 ref gene ID) & novel_le_mult (overlaps > 1 ref gene ID)
        # These will be annotated separately then combined in the final step
        novel_le_single = novel_le.subset(
            lambda df: ~df["gene_id_ref"].str.contains(",", regex=False)
        )
        novel_le_mult = novel_le.subset(
            lambda df: df["gene_id_ref"].str.contains(",", regex=False)
        )

        # Need a set of ref gene_ids from novel_le_mult to subset reference GTF
        # Present in novel_le_mult like <gene_id>,<gene_id_b> - need to split by 'comma' and expand into sinlge list
        mult_ref_ids = set(
            novel_le_mult.as_df()["gene_id_ref"].str.split(",").explode()
        )

        # Split ref le into _single & _mult
        ref_le_single = ref_le.subset(lambda df: ~df["gene_id"].isin(mult_ref_ids))
        ref_le_mult = ref_le.subset(lambda df: df["gene_id"].isin(mult_ref_ids))

        # Need to update the reference gene_id for _mult events
        # in novel - <gene_id_a>,<gene_id_b>
        # in ref - <gene_id_a>
        # create dict of {<ref_id_1>: <ref_id_1>,<ref_id_2>, <ref_id_2>: <ref_id_1>,<ref_id_2>}
        mult_ids_dict = {}
        for mult_id in novel_le_mult.as_df()["gene_id_ref"]:
            splt = mult_id.split(",")
            for id in splt:
                mult_ids_dict[id] = mult_id

        ref_le_mult = ref_le_mult.assign(
            "gene_id", lambda df: df["gene_id"].apply(lambda x: mult_ids_dict[x])
        )

        # eprint(ref_le_mult[["gene_id"]])

        # Combine across ref & novel, annotate le_ids
        # Note that novel extension events are annotated separately, as these will overlap with their reference counterpart (so need to update assigned le_number)
        eprint(
            "Combining ref & novel last exons objects and grouping last exons based on overlap..."
        )
        combined_single = annotate_le_ids(
            novel_le_single, ref_le_single, le_id_outcol="le_id", ref_extensions=ref_ext
        )
        combined_mult = annotate_le_ids(
            novel_le_mult, ref_le_mult, le_id_outcol="le_id", ref_extensions=ref_ext
        )

    # GTF containing defined last exons (with le_id etc. defined)
    combined = pr.concat([combined_single, combined_mult])

    # Want a GTF containing 'unique regions' (relative to first/internal exons) for each last exon
    # These regions will be used for quantification (don't want to assign expression in shared region only to the last exon (internal exons aren't included))

    eprint(
        "Extracting unique regions for last exons overlapping reference first/internal exons"
    )

    # Need to manually set strandedness to compare on same strand whilst awaiting clarification on behaviour
    # https://github.com/biocore-ntnu/pyranges/issues/255
    quant_combined = combined.subtract(ref_e_nl, strandedness="same")

    # Some le_ids can be dropped if they are completely contained within non-last exons
    le_ids_dropped = set(combined.le_id) - set(quant_combined.le_id)
    eprint(
        f"Number of last exon IDs dropped due to complete containment inside ref overlapping exons - {len(le_ids_dropped)}"
    )

    # eprint(combined)
    # eprint(quant_combined)

    eprint("Generating tx2le, le2gene assignment tables...")

    eprint(
        f"Writing 'tx2gene' (transcript_id | gene_id) to TSV... - {output_prefix + '.tx2gene.tsv'}"
    )

    (
        quant_combined.subset(
            lambda df: df.duplicated(subset=["ref_gene_id"], keep=False)
        )  # remove single isoform genes (keep='False' marks all duplicates as True (so keep these))
        .as_df()[["transcript_id", "ref_gene_id"]]
        .drop_duplicates()
        .rename(columns={"ref_gene_id": "gene_id"})
        .to_csv(output_prefix + ".tx2gene.tsv", sep="\t", index=False, header=True)
    )

    eprint(
        f"Writing 'tx2le' (transcript_id | le_id) to TSV... - {output_prefix + '.tx2le.tsv'}"
    )

    (
        quant_combined.subset(
            lambda df: df.duplicated(subset=["ref_gene_id"], keep=False)
        )
        .as_df()[["transcript_id", "le_id"]]
        .drop_duplicates()
        .sort_values(by="le_id")
        .to_csv(output_prefix + ".tx2le.tsv", sep="\t", index=False, header=True)
    )

    eprint(
        f"Writing 'le2gene' (le_id | gene_id) to TSV... - {output_prefix + '.le2gene.tsv'}"
    )
    (
        quant_combined.subset(
            lambda df: df.duplicated(subset=["ref_gene_id"], keep=False)
        )
        .as_df()[["le_id", "ref_gene_id"]]
        .drop_duplicates()
        .sort_values(by="le_id")
        .rename(columns={"ref_gene_id": "gene_id"})
        .to_csv(output_prefix + ".le2gene.tsv", sep="\t", index=False, header=True)
    )

    eprint(
        f"Writing 'le2genename' (le_id | gene_name) to TSV... - {output_prefix + '.le2genename.tsv'}"
    )
    (
        quant_combined.subset(
            lambda df: df.duplicated(subset=["ref_gene_id"], keep=False)
        )
        .as_df()[["le_id", "ref_gene_name"]]
        .drop_duplicates()
        .sort_values(by="le_id")
        .rename(columns={"ref_gene_name": "gene_name"})
        .to_csv(output_prefix + ".le2genename.tsv", sep="\t", index=False, header=True)
    )

    eprint(
        f"writing 'info' table (tx_id | le_id | gene_id | gene_name | event_type | coords | annot_status) to file - {output_prefix + '.info.tsv'}"
    )
    (
        combined.subset(lambda df: df.duplicated(subset=["ref_gene_id"], keep=False))
        .subset(
            lambda df: ~df["le_id"].isin(le_ids_dropped)
        )  # remove LEs completely contained within known exons
        .as_df()[
            [
                "transcript_id",
                "le_id",
                "ref_gene_id",
                "ref_gene_name",
                "event_type",
                "Chromosome",
                "Start",
                "End",
                "Strand",
            ]
        ]
        .rename(columns={"ref_gene_id": "gene_id", "ref_gene_name": "gene_name"})
        .drop_duplicates()
        .assign(
            **{
                "annot_status": lambda df: np.where(
                    df["transcript_id"].str.contains("PAPA", regex=False),
                    "novel",
                    "annotated",
                ),
            }
        )
        .sort_values(by="le_id")
        .to_csv(output_prefix + ".info.tsv", sep="\t", index=False, header=True)
    )

    eprint("Writing last exon GTFs to file...")

    eprint(
        f"Writing quantification-ready last exons GTF to file - {output_prefix + '.quant.last_exons.gtf'}"
    )
    quant_combined.drop(["gene_id_ref", "gene_name_ref", "Cluster"]).to_gtf(
        output_prefix + ".quant.last_exons.gtf"
    )

    eprint(f"Writing last exons GTF to file - {output_prefix + '.last_exons.gtf'}")
    (
        combined.drop(["gene_id_ref", "gene_name_ref", "Cluster"])
        .subset(
            lambda df: ~df["le_id"].isin(le_ids_dropped)
        )  # remove LEs completely contained within known exons
        .to_gtf(output_prefix + ".last_exons.gtf")
    )


if __name__ == "__main__":

    start = timer()

    descrpn = """Generate quantification ready GTF of last exons & group transcripts according to shared last exon"""

    parser = argparse.ArgumentParser(
        description=descrpn,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,  # Add defaults to end of help strings
    )

    parser.add_argument(
        "-i",
        "--input-gtf",
        type=str,
        default=argparse.SUPPRESS,
        required=True,
        help="Path to input GTF file containing last exons",
    )

    parser.add_argument(
        "-r",
        "--reference-gtf",
        required=True,
        type=str,
        default=argparse.SUPPRESS,
        help="Path to GTF file containing reference transcripts from which last exons will be quantified",
    )

    parser.add_argument(
        "-o",
        "--output-prefix",
        dest="output_prefix",
        type=str,
        default="novel_ref_combined",
        help="""path to/prefix for output files.
                                '.quant.last_exons.gtf' is suffixed for GTF of unique last exons regions for quantification,
                                '.last_exons.gtf' for GTF of last exons,
                                '.tx2le.tsv' suffixed for (transcript_id | le_id) TSV,
                                '.tx2gene.tsv' for (transcript_id | gene_id) TSV,
                                '.le2gene.tsv' for (le_id | gene_id) TSV,
                                """,
    )

    parser.add_argument(
        "--trust-ref-exon-number",
        action="store_true",
        default=False,
        help="Whether to 'trust' the exon number attribute in reference GTF as being strand-aware",
    )

    parser.add_argument(
        "--ref-extensions-string",
        type=str,
        default=None,
        help="Treat 'transcript_id' values in reference GTF containing this string as if novel extension events (i.e. they will be considered a distinct 'last exon isoform' to shorter isoform)",
    )

    if len(sys.argv) == 1:
        parser.print_help()
        parser.exit()

    args = parser.parse_args()

    main(
        args.input_gtf,
        args.reference_gtf,
        ref_attr_key_order,
        args.trust_ref_exon_number,
        args.ref_extensions_string,
        args.output_prefix,
    )

    end = timer()

    eprint(
        f"Script complete: took {round(end - start, 3)} s / {round((end - start) / 60, 3)} min (3 dp)"
    )
