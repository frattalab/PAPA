#!/usr/bin/env python3

import pyranges as pr
import pandas as pd
import numpy as np
from pyranges.readers import read_gtf_restricted
from papa_helpers import eprint, get_terminal_regions, get_internal_regions, add_region_number, _df_add_region_number, _pd_merge_gr, read_gtf_specific, check_concat
import sys
import argparse
from timeit import default_timer as timer


'''
Generate a GTF file of reference and novel last exons (or unique last exon segments for first/internal extensions)
Also annotate transcript IDs according to

'''


ref_attr_key_order = ["gene_id", "transcript_id", "gene_name", "exon_number"]
ref_attr_key_order_n_en = ["gene_id", "transcript_id", "gene_name"]


def cluster_to_region_number(gr, group_id_col, out_col="le_number", cluster_col="Cluster"):
    '''
    Returns gr with 'out_col' column added
    where out_col is leftmost to rightmost cluster_col converted to a
    strand-aware 1..n order by group_id_col
    1 = most 5' site in group_id_col
    n = most 3' site in group_id_col
    '''

    # For groupby.rank to work appropriately, need a single row per 'Cluster'
    c2p = (gr[[group_id_col, cluster_col]]
           .apply(lambda df: df.drop_duplicates(subset=cluster_col)
                  )
           )

    # Add 1..n 5'-3' region number as a column
    c2p = c2p.assign(out_col,
                     lambda df: _df_add_region_number(df, group_id_col, cluster_col))

    # Return 'out_col' to original gr
    c2p_cols = c2p.columns.tolist()
    out_gr = gr.apply_pair(c2p, lambda df, df_to_merge: _pd_merge_gr(df,
                                                                     df_to_merge,how="left",
                                                                     on=cluster_col,
                                                                     suffixes=[None, "_match"],
                                                                     to_merge_cols=c2p_cols)
                           )

    # avoid adding extra 'PyRanges' cols (Chromosome etc.) from c2p
    return out_gr.drop(like="_match$")


def _df_add_common_gene_id_col(df, gene_col, novel_id_str, novel_ref_gene_col):
    '''
    Internal to add_common_gene_id_col
    returns pd.Series
    '''

    if novel_ref_gene_col not in df.columns:
        # This means chr/strand pair df has no novel entries
        # gene_col = reference gene col
        return df[gene_col]

    else:
        return pd.Series(np.where(df[gene_col].str.contains(novel_id_str, regex=False),
                                  df[novel_ref_gene_col], # olapping ref gene ID for novel
                                  df[gene_col] # is a ref gene/last exon
                                  )
                         )


def add_common_gene_id_col(gr,
                           out_col="ref_gene_id",
                           gene_col="gene_id",
                           novel_id_str="PAPA",
                           novel_ref_gene_col="gene_id_ref",
                           ):

    '''
    '''

    return gr.assign(out_col,
                     lambda df: _df_add_common_gene_id_col(df, gene_col, novel_id_str, novel_ref_gene_col)
                     )


def _df_grp_update_le_number(df,
                             le_id_col="le_id",
                             le_num_col="le_number",
                             is_ext_col="is_extension"):
    '''
    # TODO: update to cumulative sum based approach
    # https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.Series.cumsum.html?highlight=cumsum
    # https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.core.groupby.DataFrameGroupBy.cumsum.html?highlight=cumsum
    - Basically be, cumulative sum of 'is_extension' (modified so same le_id = 1 extension only (1st))
    - Should have same index, so add df.groupby(id_col).apply(lambda df: df[le_number].add(cumsum))
    https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.Series.add.html?highlight=add#pandas.Series.add

    Returns: pd.Series
    '''

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
    '''
    '''

    # assert

    df = (df.groupby(id_col)
          .apply(lambda grp: grp.assign(**{out_col: _df_grp_update_le_number(grp)
                                           }
                                        )
                 )
          )


    return df


def update_ext_le_ids(gr,
                      type_col="event_type",
                      le_ext_key="last_exon_extension",
                      id_col="ref_gene_id",
                      number_col="le_number",
                      le_id_col="le_id"):
    '''
    '''

    gr = gr.assign("is_extension",
                   lambda df: pd.Series(np.where(df[type_col].str.contains(le_ext_key,
                                                                           na=False,
                                                                           regex=False),
                                                 1,
                                                 0
                                                 )
                                        )
                   )


    # Extra sort so extension events come after non-extensions of the same le_id
    gr = gr.apply(lambda df:
                  df.sort_values(by=[id_col, number_col, "is_extension"])
                  )

    # Adds column 'le_number_ext'
    gr = gr.apply(lambda df: update_extension_le_number(df))

    # eprint(combined_ext[["ref_gene_id", "le_id", "le_number", "event_type", "is_extension", "le_number_ext"]].print(n=50))

    # Reassign le_id with updated number, return to same shape as no extensions gr
    gr = (gr.drop([le_id_col, number_col, "is_extension"])
            .apply(lambda df: df.rename(columns={"le_number_ext": number_col}))
            .assign(le_id_col,
                    lambda df: df[id_col] + "_" + df[number_col].astype(int).astype(str)
                    )
          )

    return gr


def annotate_le_ids(novel_le, ref_le, novel_id_col="gene_id_ref", ref_id_col="gene_id", le_id_outcol="le_id"):
    '''
    '''

    # Identify last exon Ids containing novel extensions - these need to be handled separately
    # set of cluster IDs for each chr, strand tuple
    # {(chr, strand): {cluster_IDs}}
    d_ext_gene_ids = novel_le.apply(lambda df: set(df[df["event_type"].str.contains("extension", regex=False)][novel_id_col]) # should be ref gene ID in novel events (foun)
                                  ,
                                  as_pyranges=False)

    # Combine into a single set
    # https://blog.finxter.com/union-multiple-sets-in-python/
    ext_gene_ids = set().union(*d_ext_gene_ids.values())

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
        combined_ext = combined_ext.cluster(by="ref_gene_id")

        # Assign 5'-3' 1..n 'last_exon number' for each gene
        # Group together overlapping exons with a common identifier
        # .cluster(strand=None) groups as ('I') expect i.e. only overlapping intervals on the same strand can be merged
        combined_ext = cluster_to_region_number(combined_ext, "ref_gene_id")

        eprint("assigning 'le_id' (last exon ID) for each gene...")

        combined_ext = combined_ext.assign(le_id_outcol,
                                           lambda df: df["ref_gene_id"] + "_" + df["le_number"].astype(int).astype(str)
                                           )

        # No extensions last exon number is ready to go
        # Extension events - need to update le number so extension = own event (for last_extensions)
        combined_ext = update_ext_le_ids(combined_ext,
                                         type_col="event_type",
                                         le_ext_key="last_exon_extension",
                                         id_col="ref_gene_id",
                                         le_id_col=le_id_outcol)


    combined_n_ext = add_common_gene_id_col(combined_n_ext)

    combined_n_ext = combined_n_ext.cluster(by="ref_gene_id")

    combined_n_ext = cluster_to_region_number(combined_n_ext, "ref_gene_id")

    eprint("assigning 'le_id' (last exon ID) for each gene...")

    combined_n_ext = combined_n_ext.assign("le_id",
                                           lambda df: df["ref_gene_id"] + "_" + df["le_number"].astype(int).astype(str)
                                           )

    if len(combined_ext) > 0:
        # GTF containing defined last exons (with le_id etc. defined)
        combined = pr.concat([combined_ext, combined_n_ext])

    else:
        combined = combined_n_ext

    return combined




def main(novel_le_path,
         ref_gtf_path,
         ref_attr_key_order,
         trust_exon_number_ref,
         output_prefix
         ):
    '''
    '''

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
        ref_le = get_terminal_regions(ref)

    else:
        ref_e = ref.subset(lambda df: df["Feature"] == "exon")
        ref_e = add_region_number(ref_e, feature_key="exon", out_col="exon_number")
        ref_le = get_terminal_regions(ref_e)


    eprint("Reading in input GTF of novel last exons...")

    # novel_le = pr.read_gtf(novel_le_path)
    novel_le = read_gtf_specific(novel_le_path, ["gene_id", "transcript_id", "exon_number", "gene_id_ref", "Start_ref", "End_ref", "event_type"])


    # Make sure gene_id_ref has a single value string if only 1 distinct gene ID (Otherwise collapse to non-redundant comma-separated string of gene IDs)
    # TODO: looping over Series is a dirty way to do this, find a more canonical approach
    novel_le = novel_le.assign("gene_id_ref",
                               lambda df: pd.Series([",".join(list(dict.fromkeys(ref_id.split(","))))
                                                     for ref_id in df["gene_id_ref"]])
                               )

    # These events need to be handled separately - will create a 'metagene' combining
    eprint(f"Number of novel events matching multiple reference genes - {novel_le.as_df()['gene_id_ref'].str.contains(',', regex=False).sum()}")

    # Split novel_le into novel_le_single (overlaps 1 ref gene ID) & novel_le_mult (overlaps > 1 ref gene ID)
    # These will be annotated separately then combined in the final step
    novel_le_single = novel_le.subset(lambda df: ~df['gene_id_ref'].str.contains(',', regex=False))
    novel_le_mult = novel_le.subset(lambda df: df['gene_id_ref'].str.contains(',', regex=False))

    # Need a set of ref gene_ids from novel_le_mult to subset reference GTF
    # Present in novel_le_mult like <gene_id>,<gene_id_b> - need to split by 'comma' and expand into sinlge list
    mult_ref_ids = set(novel_le_mult.as_df()['gene_id_ref'].str.split(',').explode())

    # Split ref le into _single & _mult
    ref_le_single = ref_le.subset(lambda df: ~df['gene_id'].isin(mult_ref_ids))
    ref_le_mult = ref_le.subset(lambda df: df['gene_id'].isin(mult_ref_ids))

    # Need to update the reference gene_id for _mult events
    # in novel - <gene_id_a>,<gene_id_b>
    # in ref - <gene_id_a>
    # create dict of {<ref_id_1>: <ref_id_1>,<ref_id_2>, <ref_id_2>: <ref_id_1>,<ref_id_2>}
    mult_ids_dict = {}
    for mult_id in novel_le_mult.as_df()['gene_id_ref']:
        splt = mult_id.split(",")
        for id in splt:
            mult_ids_dict[id] = mult_id

    ref_le_mult = ref_le_mult.assign("gene_id",
                                     lambda df: df["gene_id"].apply(lambda x: mult_ids_dict[x]))

    # eprint(ref_le_mult[["gene_id"]])

    # Combine across ref & novel, annotate le_ids
    # Note that novel extension events are annotated separately, as these will overlap with their reference counterpart (so need to update assigned le_number)
    eprint("Combining ref & novel last exons objects and grouping last exons based on overlap...")
    combined_single = annotate_le_ids(novel_le_single, ref_le_single, le_id_outcol="le_id")
    combined_mult = annotate_le_ids(novel_le_mult, ref_le_mult, le_id_outcol="le_id")

    # GTF containing defined last exons (with le_id etc. defined)
    combined = pr.concat([combined_single, combined_mult])

    # Want a GTF containing 'unique regions' (relative to first/internal exons) for each last exon
    # These regions will be used for quantification (don't want to assign expression in shared region only to the last exon (internal exons aren't included))

    eprint("Extracting unique regions for last exons overlapping reference first/internal exons")
    ref_e_nl = pr.concat([get_terminal_regions(ref_e, which_region="first"),
                          get_internal_regions(ref_e)]
                         )

    eprint("Generating 'unique regions' for last exons overlapping non-last reference exons...")

    # Need to manually set strandedness to compare on same strand whilst awaiting clarification on behaviour
    # https://github.com/biocore-ntnu/pyranges/issues/255
    quant_combined = combined.subtract(ref_e_nl, strandedness="same")

    # Some le_ids can be dropped if they are completely contained within non-last exons
    le_ids_dropped = set(combined.le_id) - set(quant_combined.le_id)
    eprint(f"Number of last exon IDs dropped due to complete containment inside ref overlapping exons - {len(le_ids_dropped)}")

    # eprint(combined)
    # eprint(quant_combined)

    eprint("Generating tx2le, le2gene assignment tables...")

    eprint(f"Writing 'tx2gene' (transcript_id | gene_id) to TSV... - {output_prefix + '.tx2gene.tsv'}")

    (quant_combined.subset(lambda df: df.duplicated(subset=["gene_id"], keep=False)) # remove single isoform genes (keep='False' marks all duplicates as True (so keep these))
     .as_df()
                   [["transcript_id", "ref_gene_id"]]
                   .drop_duplicates()
                   .rename(columns={"ref_gene_id": "gene_id"})
                   .to_csv(output_prefix + ".tx2gene.tsv",
                           sep="\t",
                           index=False,
                           header=True)
     )


    eprint(f"Writing 'tx2le' (transcript_id | le_id) to TSV... - {output_prefix + '.tx2le.tsv'}")

    (quant_combined.subset(lambda df: df.duplicated(subset=["gene_id"], keep=False))
     .as_df()
     [["transcript_id", "le_id"]]
     .drop_duplicates()
     .sort_values(by="le_id")
     .to_csv(output_prefix + ".tx2le.tsv",
             sep="\t",
             index=False,
             header=True)
     )

    eprint(f"Writing 'le2gene' (le_id | gene_id) to TSV... - {output_prefix + '.le2gene.tsv'}")
    (quant_combined.subset(lambda df: df.duplicated(subset=["gene_id"], keep=False))
     .as_df()
     [["le_id", "ref_gene_id"]]
     .drop_duplicates()
     .sort_values(by="le_id")
     .rename(columns={"ref_gene_id": "gene_id"})
     .to_csv(output_prefix + ".le2gene.tsv",
             sep="\t",
             index=False,
             header=True)
     )

    eprint(f"writing 'info' table (tx_id | le_id | gene_id | gene_name | event_type | coords | annot_status) to file - {output_prefix + '.info.tsv'}")
    (quant_combined.subset(lambda df: df.duplicated(subset=["gene_id"], keep=False))
     .as_df()
     [["transcript_id", "le_id", "ref_gene_id", "gene_name", "event_type", "Chromosome", "Start", "End", "Strand"]]
     .rename(columns={"ref_gene_id": "gene_id"})
     .drop_duplicates()
     .assign(**{"annot_status": lambda df: np.where(df["transcript_id"].str.contains("PAPA", regex=False),
                                                    "novel",
                                                    "annotated"),
                }
             )
     .sort_values(by="le_id")
     .to_csv(output_prefix + ".info.tsv",
             sep="\t",
             index=False,
             header=True)

    )




    eprint("Writing last exon GTFs to file...")

    eprint(f"Writing quantification-ready last exons GTF to file - {output_prefix + '.quant.last_exons.gtf'}")
    quant_combined.to_gtf(output_prefix + ".quant.last_exons.gtf")

    eprint(f"Writing last exons GTF to file - {output_prefix + '.last_exons.gtf'}")
    combined.to_gtf(output_prefix + ".last_exons.gtf")



if __name__ == '__main__':

    start = timer()

    descrpn = """Generate quantification ready GTF of last exons & group transcripts according to shared last exon"""

    parser = argparse.ArgumentParser(description=descrpn,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter, # Add defaults to end of help strings
                                     )

    parser.add_argument("-i",
                        "--input-gtf",
                        type=str,
                        default=argparse.SUPPRESS,
                        required=True,
                        help="Path to input GTF file containing last exons")

    parser.add_argument("-r",
                        "--reference-gtf",
                        required=True,
                        type=str,
                        default=argparse.SUPPRESS,
                        help="Path to GTF file containing reference transcripts from which last exons will be quantified")

    parser.add_argument("-o","--output-prefix",
                        dest="output_prefix",
                        type=str,
                        default="novel_ref_combined",
                        help="""path to/prefix for output files.
                                '.quant.last_exons.gtf' is suffixed for GTF of unique last exons regions for quantification,
                                '.last_exons.gtf' for GTF of last exons,
                                '.tx2le.tsv' suffixed for (transcript_id | le_id) TSV,
                                '.tx2gene.tsv' for (transcript_id | gene_id) TSV,
                                '.le2gene.tsv' for (le_id | gene_id) TSV,
                                """)

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
         ref_attr_key_order,
         args.trust_ref_exon_number,
         args.output_prefix)

    end = timer()

    eprint(f"Script complete: took {round(end - start, 3)} s / {round((end - start) / 60, 3)} min (3 dp)")
