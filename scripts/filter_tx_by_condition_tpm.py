##!/usr/bin/env python3

import pyranges as pr
import pandas as pd
from papa_helpers import eprint
import os
import sys
import argparse
from timeit import default_timer as timer

'''
'''


def main(gtf_path,tracking_path,gtf_list_path,min_mean_tpm,output_prefix):
    '''
    '''


    #1. Read list of input GTF files - these fill out columns post the 4 default cols in tracking file
    eprint("Collecting input GTF file names and reading in '.tracking' TSV...")
    with open(gtf_list_path) as infile:
        samples_list = [os.path.basename(line.rstrip("\n")) for line in infile]


    # tx_id_in_merged | super locus/gene ID | reference gene ID | overlap class code
    tracking_default_cols = ["transcript_id",
                             "gene_id",
                             "ref_gene_id",
                             "overlap_class_code"]

    tracking_cols = tracking_default_cols + samples_list

    tracking = pd.read_csv(tracking_path,
                           sep="\t",
                           names=tracking_cols)

    # eprint(tracking)

    #2. Extract TPM from each 'sample' string
    # 5th 'column' is the TPM value
    # q1:PAPA.2|PAPA.2.2|2|0.106867|0.277613|0.795843|535
    eprint("Extracting TPM value from matching transcripts strings...")
    tracking[samples_list] = tracking[samples_list].apply(lambda x: x.str.split("\|", expand=True)[4], axis=0)

    #3. Txipts not always found in all samples - replace extracted TPM (None) with 0
    eprint("Replacing missing TPM values (tx not expressed in sample) with 0...")
    tracking = tracking.fillna(0)

    #4. Calculate mean across samples
    #Convert all sample_cols to numeric
    eprint("Calculating mean TPM per transcript...")
    tracking = tracking.astype({sample: "float" for sample in samples_list})

    tracking["mean_tpm"] = tracking[samples_list].mean(axis=1)

    #5. Get set of transcripts passing min mean TPM filter
    eprint(f"Filtering transcripts based on min mean TPM - {min_mean_tpm}")
    eprint(f"Number of transcripts pre mean TPM >= {min_mean_tpm} filter - {len(set(tracking.transcript_id))}")

    valid_tx_ids = set(tracking.loc[tracking["mean_tpm"] >= min_mean_tpm,
                                    "transcript_id"])

    eprint(f"Number of transcripts post mean TPM >= {min_mean_tpm} filter - {len(valid_tx_ids)}")

    #6. Filter GTF for valid_tx_ids
    eprint("Reading in input GTF...")
    gtf = pr.read_gtf(gtf_path)

    eprint("Subsetting input GTF for transcripts passing filter...")
    gtf = gtf.subset(lambda df: df.transcript_id.isin(valid_tx_ids))

    eprint("Writing filtered GTF to file...")
    gtf.to_gtf(output_prefix + ".gtf")

    #7. Generate TPM matrix & output to file
    # transcript_id | gene_id | mean_tpm | sample1..n
    eprint("Generating TPM matrix...")
    tpm_cols = ["transcript_id", "gene_id", "mean_tpm"] + samples_list
    tracking = tracking[tpm_cols]

    tracking.to_csv(output_prefix + ".tpm_matrix.tsv",
                    sep="\t",
                    index=False,
                    header=True)





if __name__ == '__main__':

    start = timer()

    descrpn = """Script to filter assembled transcripts based on a minimum mean TPM expression across samples of the same condition"""

    parser = argparse.ArgumentParser(description=descrpn,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter, # Add defaults to end of help strings
                                     )

    parser.add_argument("-i",
                        "--input-gtf",
                        dest="input_gtf",
                        type=str,
                        default=argparse.SUPPRESS,
                        required=True,
                        help="Path to merged GTF file for samples of the same condition output by GFFcompare")

    parser.add_argument("-t",
                        "--tracking",
                        type=str,
                        default=argparse.SUPPRESS,
                        required=True,
                        help="Path to '<prefix>.tracking' TSV file produce by GFFcompare for input GTF's merge")

    parser.add_argument("-l",
                        "--gtf-list",
                        type=str,
                        default=argparse.SUPPRESS,
                        required=True,
                        help="Path to list of GTF files used as input to GFFcompare merge call. Used to replace 'query IDs' with more identifiable sample names in output TPM matrix.")

    parser.add_argument("-m",
                        "--min-tpm",
                        type=float,
                        default=1.0,
                        help="Minimum mean TPM (across all samples) for a transcript to be retained in output GTF"
                        )

    parser.add_argument("-o","--output-prefix",
                        dest="output_prefix",
                        type=str,
                        default="mean_tpm_filtered_transcripts",
                        help="Prefix for output files. '.<suffix>' added depending on output file type ('.gtf' for valid transcripts GTF, '.tpm_matrix.tsv' for TSV with TPM expression across samples)")

    if len(sys.argv) == 1:
        parser.print_help()
        parser.exit()

    args = parser.parse_args()

    main(args.input_gtf, args.tracking,args.gtf_list,args.min_tpm,args.output_prefix)

    end = timer()

    eprint(f"Complete: took {round(end - start, 3)} s / {round((end - start) / 60, 3)} min (3 dp)")
