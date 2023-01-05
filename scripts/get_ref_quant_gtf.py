#!/usr/bin/env python3

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
    add_region_rank,
)
from get_combined_quant_gtf import cluster_to_region_number, update_ext_le_ids, annotate_ref_event_types
import sys
import argparse
from timeit import default_timer as timer


"""
Generate a GTF file of reference last exons (or unique last exon segments for first/internal extensions) for quantification with Salmon
Note: this follows the workflow of get_combined_quant_gtf.py, the main distinction being that this script will not include a novel GTF
"""


ref_attr_key_order = ["gene_id", "transcript_id", "gene_name", "exon_number"]
ref_attr_key_order_n_en = ["gene_id", "transcript_id", "gene_name"]


def annotate_le_ids_ref(ref_le: pr.PyRanges,ref_id_col="gene_id",le_id_outcol="le_id",ref_extensions=False):
    '''_summary_

    _extended_summary_
    Mimics annotate_le_ids fro get_combined_quant_gtf.py, just allows for not having novel LEs
    TODO: make a general function that makes novel_le optional

    Parameters
    ----------
    ref_le : pr.PyRanges
        _description_
    ref_id_col : str, optional
        _description_, by default "gene_id"
    le_id_outcol : str, optional
        _description_, by default "le_id"
    ref_extensions : bool, optional
        _description_, by default False
    '''
    
    if ref_extensions:
        # Want to also consider specific reference events as if novel extensions
        # Will handle these genes separately
        # {(chr, strand): {gene_id1, gene_id2}}
        d_ref_ext_gene_ids = ref_le.apply(
            lambda df: set(
                df[df["event_type"].str.contains("extension", regex=False)][ref_id_col]
            ),
            as_pyranges=False,
        )

        n_ref_ext_ids = sum(len(ids_set) for ids_set in ref_ext_gene_ids.values())
        eprint(
            f"Number of genes containing extension events sourced from labelled annotated/reference transcripts - {n_ref_ext_ids}"
        )
        
        # Combine into a single set
        # https://blog.finxter.com/union-multiple-sets-in-python/
        ref_ext_gene_ids = set().union(*d_ref_ext_gene_ids.values())
        
    else:
        ref_ext_gene_ids = set()
    
    if len(ref_ext_gene_ids) == 0:
        # no extensions in ref GTF to treat as separate events
        # Assign 5'-3' 1..n 'last_exon number' for each gene
        # Group together overlapping exons with a common identifier
        # .cluster(strand=None) groups as ('I') expect i.e. only overlapping intervals on the same strand can be merged
        ref_le_ext = ref_le
        ref_le_n_ext = pr.PyRanges()
        ref_le = ref_le.cluster(by=ref_id_col, strand=None)
        ref_le = cluster_to_region_number(ref_le, ref_id_col)
        
        eprint("assigning 'le_id' (last exon ID) for each gene...")

        ref_le = ref_le.assign(le_id_outcol,
                               lambda df: df[ref_id_col] + "_" + df["le_number"].astype(int).astype(str),
                    )
        
        out = ref_le
    
    else:
        # Need to assign le_ids for genes with/without extensions separately
        ref_le_ext = ref_le.subset(lambda df: df[ref_id_col].isin(ref_ext_gene_ids))
        ref_le_n_ext = ref_le.subset(lambda df: ~(df[ref_id_col].isin(ref_ext_gene_ids)))
        
    
    ref_le_n_ext = ref_le_n_ext.cluster(by=ref_id_col, strand=None)
    ref_le_n_ext = cluster_to_region_number(ref_le_n_ext, ref_id_col)
    
    ref_le_ext = ref_le_ext.cluster(by=ref_id_col, strand=None)
    ref_le_ext = cluster_to_region_number(ref_le_ext, ref_id_col)
        
    eprint("assigning 'le_id' (last exon ID) for each gene...")
    ref_le_n_ext = ref_le_n_ext.assign(le_id_outcol,
                           lambda df: df[ref_id_col] + "_" + df["le_number"].astype(int).astype(str),
                        )
    
    ref_le_ext = ref_le_ext.assign(le_id_outcol,
                           lambda df: df[ref_id_col] + "_" + df["le_number"].astype(int).astype(str)
    )
    
    # Extension events - need to update le number so extension = own event (for last_extensions)
    ref_le_ext = update_ext_le_ids(ref_le_ext,
                                       type_col="event_type",
                                       le_ext_key="last_exon_extension",
                                       id_col=ref_id_col,
                                       le_id_col=le_id_outcol
                                       )
    
    out = pr.concat([ref_le_ext, ref_le_n_ext])
    
    return out
    

def main(ref_gtf_path,
         ref_attr_key_order,
         trust_exon_number_ref,
         ref_extensions_string,
         output_prefix,
):

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

    eprint("Extracting introns for each transcript in reference GTF...")

    ref_i = ref.features.introns(by="transcript")
    # Add region number and annotate as first, internal or last
    ref_i = add_region_number(ref_i, feature_key="intron", out_col="intron_number")
    ref_i = add_region_rank(ref_i, region_number_col="intron_number")

    ref_li = get_terminal_regions(
        ref_i, feature_key="intron", region_number_col="intron_number"
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

        n_ref_ext = ref_le.event_type.loc[lambda x: x == "last_exon_extension"].sum()

        if n_ref_ext == 0:
            raise Exception(
                f"No reference transcript_id containing provided string - {ref_extensions_string} - were found. Correct or do not pass --ref-extensions-string to avoid this error message"
            )

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
            ref_le_n_ext, ref_li, ref_e, ref_i, remove_last_exon_extensions=True
        )

        ref_le = pr.concat([ref_le_ext, ref_le_n_ext])

    else:
        eprint(
            "--ref-extensions-string not provided - all reference transcripts will be 'collapsed' to same last exon ID given overlap"
        )
        ref_ext = False

        ref_le = annotate_ref_event_types(
            ref_le, ref_li, ref_e, ref_i, remove_last_exon_extensions=True
        )
        
    # Add an identifier grouping overlapping last exons together
    ref_le = annotate_le_ids_ref(ref_le, ref_extensions=ref_ext)
    
    # Want a GTF containing 'unique regions' of last exons 
    # These regions will be used for quantification (don't want to assign expression in shared region only to the last exon (internal exons aren't included))
    
    eprint(
        "Extracting unique regions for last exons overlapping reference first/internal exons"
    )
    ref_e_nl = pr.concat(
        [get_terminal_regions(ref_e, which_region="first"), get_internal_regions(ref_e)]
    )

    eprint(
        "Generating 'unique regions' for last exons overlapping non-last reference exons..."
    )

    # Need to manually set strandedness to compare on same strand whilst awaiting clarification on behaviour
    # https://github.com/biocore-ntnu/pyranges/issues/255
    quant_ref_le = ref_le.subtract(ref_e_nl, strandedness="same")

    # Some le_ids can be dropped if they are completely contained within non-last exons
    le_ids_dropped = set(ref_le.le_id) - set(quant_ref_le.le_id)
    eprint(
        f"Number of last exon IDs dropped due to complete containment inside ref overlapping exons - {len(le_ids_dropped)}"
    )
    
    eprint("Generating tx2le, le2gene assignment tables...")

    eprint(
        f"Writing 'tx2gene' (transcript_id | gene_id) to TSV... - {output_prefix + '.tx2gene.tsv'}"
    )

    (
        quant_ref_le.subset(
            lambda df: df.duplicated(subset=["gene_id"], keep=False)
        )  # remove single isoform genes (keep='False' marks all duplicates as True (so keep these))
        .as_df()[["transcript_id", "gene_id"]]
        .drop_duplicates()
        .to_csv(output_prefix + ".tx2gene.tsv", sep="\t", index=False, header=True)
    )

    eprint(
        f"Writing 'tx2le' (transcript_id | le_id) to TSV... - {output_prefix + '.tx2le.tsv'}"
    )

    (
        quant_ref_le.subset(
            lambda df: df.duplicated(subset=["gene_id"], keep=False)
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
        quant_ref_le.subset(
            lambda df: df.duplicated(subset=["gene_id"], keep=False)
        )
        .as_df()[["le_id", "gene_id"]]
        .drop_duplicates()
        .sort_values(by="le_id")
        .to_csv(output_prefix + ".le2gene.tsv", sep="\t", index=False, header=True)
    )

    eprint(
        f"Writing 'le2genename' (le_id | gene_name) to TSV... - {output_prefix + '.le2genename.tsv'}"
    )
    (
        quant_ref_le.subset(
            lambda df: df.duplicated(subset=["gene_id"], keep=False)
        )
        .as_df()[["le_id", "gene_name"]]
        .drop_duplicates()
        .sort_values(by="le_id")
        .to_csv(output_prefix + ".le2genename.tsv", sep="\t", index=False, header=True)
    )

    eprint(
        f"writing 'info' table (tx_id | le_id | gene_id | gene_name | event_type | coords | annot_status) to file - {output_prefix + '.info.tsv'}"
    )
    (
        ref_le.subset(lambda df: df.duplicated(subset=["gene_id"], keep=False))
        .subset(
            lambda df: ~df["le_id"].isin(le_ids_dropped)
        )  # remove LEs completely contained within known exons
        .as_df()[
            [
                "transcript_id",
                "le_id",
                "gene_id",
                "gene_name",
                "event_type",
                "Chromosome",
                "Start",
                "End",
                "Strand",
            ]
        ]
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
    quant_ref_le.drop(["Cluster"]).to_gtf(
        output_prefix + ".quant.last_exons.gtf"
    )

    eprint(f"Writing last exons GTF to file - {output_prefix + '.last_exons.gtf'}")
    (
        ref_le.drop(["Cluster"])
        .subset(
            lambda df: ~df["le_id"].isin(le_ids_dropped)
        )  # remove LEs completely contained within known exons
        .to_gtf(output_prefix + ".last_exons.gtf")
    )
    
if __name__ == "__main__":

    start = timer()

    descrpn = """Generate quantification ready GTF of reference last exons & group transcripts according to shared last exon"""

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
        ref_attr_key_order,
        args.trust_ref_exon_number,
        args.ref_extensions_string,
        args.output_prefix,
    )

    end = timer()

    eprint(
        f"Script complete: took {round(end - start, 3)} s / {round((end - start) / 60, 3)} min (3 dp)"
    )
