#!/usr/bin/env python3

import pyranges as pr
import pandas as pd
from papa_helpers import eprint, _n_ids
import os
import sys
import argparse
from timeit import default_timer as timer

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
                     lambda df: pd.Series([sample_id]*len(df), index=df.index)
                     )

    return gtf


def main(input_gtf_list,
         sample_id_col,
         fname_suffix,
         le_id_col,
         out_gtf):
    '''
    '''

    le_gtfs = [read_process_le_gtf(pth, sample_id_col, fname_suffix) for pth in input_gtf_list]

    les = pr.concat(le_gtfs)

    # TODO: consider using 5'ends (w/ some tolerance) to group together LEs of the same gene

    # Assign a 'last exon ID' based on whether last exons overlap
    les = (les.cluster()
             .apply(lambda df: df.rename({"Cluster": le_id_col},
                                         axis="columns")
                    )
           )

    # eprint(les.columns)
    eprint(f"Number of last exons after grouping by overlap - {_n_ids(les, le_id_col)}")

    les.to_gtf(out_gtf)


if __name__ == '__main__':

    start = timer()

    descrpn = """Combine multiple GTFs containing last exons into a single GTF, grouping last exons across GTFs according to overlap and annotating according to source GTF"""

    parser = argparse.ArgumentParser(description=descrpn,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter, # Add defaults to end of help strings
                                     )

    parser.add_argument("-i",
                        "--input-gtfs",
                        required=True,
                        nargs="+",
                        default=argparse.SUPPRESS,
                        type=str,
                        help="Path to input GTF files of last exons to merge. If multiple, pass paths consecutively and space-separated")

    parser.add_argument("-s",
                        "--gtf-suffix",
                        default="_last_exons.gtf",
                        type=str,
                        help="Suffix string to remove from GTF filenames when assigning source sample ID")

    parser.add_argument("-o",
                        "--output-gtf",
                        type=str,
                        default="merged_last_exons.gtf",
                        help="Name of/path to output GTF file containing combined last exons across input GTFs")

    parser.add_argument("--sample-id-attribute-name",
                        type=str,
                        default="sample_id",
                        dest="sample_id_col",
                        help="Name of attribute storing sample information/source GTF file for last exon")

    parser.add_argument("--last-exon-id-attribute-name",
                        type=str,
                        default="last_exon_id",
                        dest="le_id_col",
                        help="Name of attribute storing last exon grouping identifier")

    if len(sys.argv) == 1:
        parser.print_help()
        parser.exit()

    args = parser.parse_args()

    main(args.input_gtfs,
         args.sample_id_col,
         args.gtf_suffix,
         args.le_id_col,
         args.output_gtf)

    end = timer()

    eprint(f"Script complete: took {round(end - start, 3)} s / {round((end - start) / 60, 3)} min (3 dp)")
