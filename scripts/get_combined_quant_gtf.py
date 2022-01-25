#!/usr/bin/env python3

import pyranges as pr
import pandas as pd
import numpy as np
from pyranges.readers import read_gtf_restricted
from papa_helpers import eprint, get_terminal_regions, add_region_number, _df_add_region_number, _pd_merge_gr, read_gtf_specific
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


def add_common_gene_id_col(gr,
                           out_col="ref_gene_id",
                           gene_col="gene_id",
                           novel_id_str="PAPA",
                           novel_ref_gene_col="gene_id_b",
                           ):

    '''
    '''

    return gr.assign(out_col,
                     lambda df: pd.Series(np.where(df[gene_col].str.contains(novel_id_str, regex=False),
                                                   df[novel_ref_gene_col], # olapping ref gene ID for novel
                                                   df[gene_col] # is a ref gene/last exon
                                                   )
                                          )
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
    - Should have same index, so add df.group_by(id_col).apply(lambda df: df[le_number].add(cumsum))
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
    novel_le = read_gtf_specific(novel_le_path, ["gene_id", "transcript_id", "exon_number", "gene_id_b", "Start_b", "End_b", "event_type"])

    eprint("Combining ref & novel last exons objects and grouping last exons based on overlap...")

    # Identify last exon Ids containing novel extensions - these need to be handled separately
    # set of cluster IDs for each chr, strand tuple
    # {(chr, strand): {cluster_IDs}}
    d_ext_gene_ids = novel_le.apply(lambda df: set(df[df["event_type"].str.contains("extension", regex=False)]["gene_id_b"]) # should be ref gene ID in novel events (foun)
                                                ,
                                    as_pyranges=False)

    # eprint(d_ext_gene_ids)
    #
    # for k,v in d_ext_gene_ids.items():
    #     eprint(k)
    #     eprint(d_ext_gene_ids[k])

    # Combine into a single set
    # https://blog.finxter.com/union-multiple-sets-in-python/
    ext_gene_ids = set().union(*d_ext_gene_ids.values())

    ref_le_ext = ref_le.subset(lambda df: df["gene_id"].isin(ext_gene_ids))
    ref_le_n_ext = ref_le.subset(lambda df: ~(df["gene_id"].isin(ext_gene_ids)))

    # eprint("Reference genes containing novel extensions")
    # eprint(ref_le_ext)
    #
    # eprint("Reference genes without novel extensions")
    # eprint(ref_le_n_ext)

    novel_le_ext = novel_le.subset(lambda df: df["gene_id_b"].isin(ext_gene_ids))

    # eprint("Genes with novel last exon (extensions)")
    # eprint(novel_le_ext)

    # eprint(novel_le.columns)

    novel_le_n_ext = novel_le.subset(lambda df: ~(df["gene_id_b"].isin(ext_gene_ids)))

    # eprint("Genes with novel last exon (no extensions)")
    # eprint(novel_le_n_ext)

    # eprint(novel_le_n_ext.columns)
    # eprint(novel_le_n_ext)


    # Make combined GTF of reference and novel last exons
    # Keeping extensions separate
    combined_ext = pr.concat([ref_le_ext, novel_le_ext])
    combined_n_ext = pr.concat([ref_le_n_ext, novel_le_n_ext])

    # Group together overlapping exons with a common identifier
    combined_ext = combined_ext.cluster()
    combined_n_ext = combined_n_ext.cluster()

    # eprint(combined_ext[["gene_id", "Cluster", "event_type"]])
    # eprint(combined_ext.columns)
    #
    # eprint(combined_n_ext[["gene_id", "Cluster", "event_type"]])
    # eprint(combined_n_ext.columns)

    # make sure gene_id corresponds to reference gene ID
    # Otherwise novel + ref of same gene will be considered separately (annotated in diff column)
    combined_ext = add_common_gene_id_col(combined_ext)
    combined_n_ext = add_common_gene_id_col(combined_n_ext)

    # Assign 5'-3' 1..n 'last_exon number' for each gene
    combined_ext = cluster_to_region_number(combined_ext, "ref_gene_id")

    # DF blows up if I don't do this
    # eprint(len(combined_ext))
    # combined_ext = combined_ext.apply(lambda df: df.drop_duplicates(subset="transcript_id"))
    # eprint(len(combined_ext))

    combined_n_ext = cluster_to_region_number(combined_n_ext, "ref_gene_id")

    eprint("assigning 'le_id' (last exon ID) for each gene...")

    combined_ext = combined_ext.assign("le_id",
                                       lambda df: df["ref_gene_id"] + "_" + df["le_number"].astype(int).astype(str)
                                       )

    combined_n_ext = combined_n_ext.assign("le_id",
                                           lambda df: df["ref_gene_id"] + "_" + df["le_number"].astype(int).astype(str)
                                           )

    # eprint(len(combined_n_ext))
    # combined_n_ext = combined_n_ext.apply(lambda df: df.drop_duplicates(subset="transcript_id"))
    # eprint(len(combined_n_ext))

    # No extensions last exon number is ready to go
    # Extension events - need to update le number so extension = own event
    combined_ext = combined_ext.assign("is_extension",
                                       lambda df: pd.Series(np.where(df["event_type"].str.contains("last_exon_extension",
                                                                                                   na=False,
                                                                                                   regex=False),
                                                                     1,
                                                                     0
                                                                     )
                                                            )
                                       )


    # Extra sort so extension events come after non-extensions of the same le_id
    combined_ext = combined_ext.apply(lambda df:
                                      df.sort_values(by=["ref_gene_id", "le_number", "is_extension"])
                                      )

    combined_ext = combined_ext.apply(lambda df: update_extension_le_number(df))

    eprint(combined_ext[["ref_gene_id", "le_id", "le_number", "event_type", "is_extension", "le_number_ext"]].print(n=50))

    # Reassign le_id with updated number, return to same shape as no extensions gr
    combined_ext = (combined_ext.drop(["le_id", "le_number", "is_extension"])
                    .apply(lambda df: df.rename(columns={"le_number_ext": "le_number"}))
                    .assign("le_id",
                            lambda df: df["ref_gene_id"] + "_" + df["le_number"].astype(int).astype(str)
                            )
                    )

    eprint(combined_ext)

    # Get set of regions from which need to extract unique regions
    int_mask = (combined_ext.event_type.isin(["first_exon_extension", "internal_exon_extension"]))

    combined_ext_int = combined_ext[int_mask]

    eprint(combined_ext_int[["gene_id", "event_type", "Start_b", "End_b"]])

    olap_exons = combined_ext_int.apply(lambda df:
                                        (df[["Chromosome", "Start_b", "End_b", "Strand"]]
                                         .rename(columns={"Start_b": "Start", "End_b": "End"})
                                         .astype({"Start": "int32", "End": "int32"})))

    eprint(ref_e[ref_e.gene_id == "ENSG00000000457.14"])
    eprint(combined_ext_int[combined_ext_int.ref_gene_id == "ENSG00000000457.14"])
    eprint(combined_ext_int)

    eprint("Getting extension-specific regions for first & internal events")
    # Substract overlapping exon region from first & internal extensions

    combined_ext_int = combined_ext_int.subtract(ref_e)

    eprint(combined_ext_int)

    # combined_ext_int.to_gtf("wtf_internal.gtf")

    quant_combined_ext = pr.concat([combined_ext[~int_mask], combined_ext_int])

    # GTF containing last exon regions for quantification purposes
    quant_combined = pr.concat([quant_combined_ext, combined_n_ext])

    # GTF containing defined last exons
    combined = pr.concat([combined_ext, combined_n_ext])

    eprint("Generating tx2le, le2gene assignment tables...")

    eprint(f"Writing 'tx2gene' (transcript_id | gene_id) to TSV... - {output_prefix + '.tx2gene.tsv'}")
    quant_combined.as_df()[["transcript_id",
                            "ref_gene_id"]].drop_duplicates().to_csv(output_prefix + ".tx2gene.tsv",
                                                                     sep="\t",
                                                                     index=False,
                                                                     header=True)


    eprint(f"Writing 'le2pas' (transcript_id | pas_id) to TSV... - {output_prefix + '.tx2le.tsv'}")

    quant_combined.as_df()[["transcript_id",
                            "le_id"]].drop_duplicates().to_csv(output_prefix + ".tx2le.tsv",
                                                               sep="\t",
                                                               index=False,
                                                               header=True)


    eprint(f"Writing 'le2gene' (le_id | gene_id) to TSV... - {output_prefix + '.le2gene.tsv'}")
    quant_combined.as_df()[["le_id",
                            "ref_gene_id"]].drop_duplicates().to_csv(output_prefix + ".le2gene.tsv",
                                                                     sep="\t",
                                                                     index=False,
                                                                     header=True)




    eprint("Writing last exon GTFs to file...")

    eprint(f"Writing quantification-ready last exons GTF to file - {output_prefix + '.quant.last_exons.gtf'}")
    quant_combined.to_gtf(output_prefix + ".quant.last_exons.gtf")

    eprint(f"Writing last exons GTF to file - {output_prefix + '.last_exons.gtf'}")
    combined.to_gtf(output_prefix + ".last_exons.gtf")



if __name__ == '__main__':

    # eprint(sys.argv[1])
    # eprint(sys.argv[2])
    main(sys.argv[1], sys.argv[2], ref_attr_key_order, False, "its_ready")
