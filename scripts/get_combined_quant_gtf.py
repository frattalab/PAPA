#!/usr/bin/env python3

import pyranges as pr
import pandas as pd
import numpy as np
from pyranges.readers import read_gtf_restricted
from papa_helpers import eprint, get_terminal_regions, add_region_number, read_gtf_specific
import sys
import argparse
from timeit import default_timer as timer


'''
Generate a GTF file of reference and novel last exons (or unique last exon segments for first/internal extensions)
Also annotate transcript IDs according to

'''

ref_attr_key_order = ["gene_id", "transcript_id", "gene_name", "exon_number"]
ref_attr_key_order_n_en = ["gene_id", "transcript_id", "gene_name"]

def main(novel_le_path,
         ref_gtf_path,
         ref_attr_key_order,
         trust_exon_number_ref,
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
        ref = add_region_number(ref_e, feature_key="exon", out_col="exon_number")
        ref_le = get_terminal_regions(ref)

    eprint("Reading in input GTF of novel last exons...")

    novel_le = pr.read_gtf(novel_le_path)

    eprint("Combining ref & novel last exons objects and grouping last exons based on overlap...")

    combined = pr.concat([ref_le, novel_le])

    combined = combined.cluster()

    # Identify last exon Ids containing novel extensions - these need to be handled separately
    # set of cluster IDs for each chr, strand tuple
    # {(chr, strand): {cluster_IDs}}
    cluster_ext = combined.apply(lambda df: set(df.loc[df["event_type"].str.contains("extension"),
                                                       "gene_id_b"] # should be ref gene ID in novel events (foun)
                                                ),
                               as_pyranges=False)

    # # Combine into a single set
    # le_id_ext = set().union(*le_id_ext)

    # Get separate grs for genes containing novel extensions and those that do not
    combined_ext = combined.df(lambda df: (df["gene_id"].isin(le_id_ext[(df["Chromosome"].iloc[0],
                                                                         df["Strand"].iloc[0])
                                                                        ]
                                                              )
                                           ) |
                                          (df["gene_id_b"].isin(le_id_ext[(df["Chromosome"].iloc[0],
                                                                           df["Strand"].iloc[0]
                                                                           )
                                                                          ]
                                                                )
                                           )
                                )

    combined_n_ext = combined.df(lambda df: (~df["gene_id"].isin(le_id_ext[(df["Chromosome"].iloc[0],
                                                                            df["Strand"].iloc[0])
                                                                           ]
                                                              )
                                             ) |
                                 (~df["gene_id_b"].isin(le_id_ext[(df["Chromosome"].iloc[0],
                                                                   df["Strand"].iloc[0]
                                                                   )
                                                                  ]
                                                        )
                                  )
                                 )

    eprint(combined_ext[["event_type"]])

    eprint(combined_n_ext[["event_type"]])


if __name__ == '__main__':


    main(sys.argv[1], sys.argv[2], ref_attr_key_order, False)
