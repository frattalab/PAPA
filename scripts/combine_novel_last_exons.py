#!/usr/bin/env python3

import pyranges as pr
import pandas as pd
import os
import sys
import argparse

'''
Merge GTFs of putative novel last exons into a single GTF, grouping last exons with a common ID if they have any overlap
'''

def read_process_le_gtf(path, sample_id_col, fname_suffix):
    '''
    '''

    sample_id = os.path.basename(path).replace(fname_suffix, "")

    eprint(f"Reading in last exons at {path}")
    gtf = pr.read_gtf(path)

    gtf = gtf.assign(sample_id_col,
                     lambda df: pd.Series([sample_id]*len(df), index=df.index())
                     )

    return gtf


def main(input_gtf_list,
         sample_id_col,
         fname_suffix,
         le_id_col,
         out_gtf):
    '''
    '''

    le_gtfs = [read_process_le_gtf(pth, sample_id_col) for pth in input_gtf_list]

    les = pr.concat(le_gtfs)

    # Assign a 'last exon ID' based on whether last exons overlap
    les = (pr.cluster()
             .apply(lambda df: df.rename({"Cluster": le_id_col},
                                         axis="columns")
                    )
           )

    eprint(les.columns)

    les.to_gtf(out_gtf)
