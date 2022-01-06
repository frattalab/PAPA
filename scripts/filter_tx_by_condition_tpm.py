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


def _filter_gtf(col, output_suffix, tx_id_col="transcript_id", mask_ids=[None]):
    '''
    Output a GTF file filtered for transcript_ids in a series

    col: Named Series of transcript ID strings (e.g. a column from a df). Name should be the path to the input GTF file
    output_suffix: Suffix string to add to output filtered GTF file (must end in '.gtf')
    tx_id_col: Name of column/attribute containing transcript IDs
    mask_ids: list of IDs to remove from set of IDs used for filtering

    returns col unmodified (but w)

    Intended to be applied internally to col-wise pd.apply with GFFcompare 'tracking' file
    '''

    assert output_suffix.endswith(".gtf")

    eprint(f"Filtering following file for valid transcript IDs - {col.name}")
    gtf = pr.read_gtf(col.name)


    ids = set(col) - set(mask_ids)

    eprint(f"Number of tx ids to retain - {len(ids)}")

    eprint("Filtering for valid transcript IDs...")

    gtf = gtf.filter(lambda df: df[tx_id_col].isin(ids))

    out_gtf = col.name.rstrip(".gtf") + output_suffix

    gtf.to_gtf(out_gtf)

    return col




def main(tracking_path,
         gtf_list_path,
         min_mean_tpm,
         gtf_output_suffix):
    '''
    '''


    #1. Read list of input GTF files - these fill out columns post the 4 default cols in tracking file
    eprint("Collecting input GTF file names and reading in '.tracking' TSV...")
    with open(gtf_list_path) as infile:
        samples_list = [line.rstrip("\n") for line in infile]


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

    eprint(f"Number of transcripts post mean TPM >= {min_mean_tpm} filter - {len(set(valid_tx_ids.transcript_id))}")

    # Extract transcript id from the 'tracking' string for each sample
    # 2nd 'column' is the transcript_id in that sample1
    # q1:PAPA.2|PAPA.2.2|2|0.106867|0.277613|0.795843|535
    # Note: if no tx in that sample column contains '-' (split returns None)
    valid_tx_ids[samples_list] = (valid_tx_ids[samples_list]
                                  .apply(lambda col: col.str.split("\|", expand=True)[1],
                                         axis="index")
                                  )

    #6. Filter individual GTFs for valid_tx_ids
    # Each column, read in GTF, filter & output
    # Note: this function returns the series unmodified
    valid_tx_ids[samples_list] = valid_tx_ids[samples_list].apply(lambda col: _filter_gtf(col, gtf_output_suffix), axis="index")


    # #7. Generate TPM matrix & output to file
    # # transcript_id | gene_id | mean_tpm | sample1..n
    # eprint("Generating TPM matrix...")
    # tpm_cols = ["transcript_id", "gene_id", "mean_tpm"] + samples_list
    # tracking = tracking[tpm_cols]
    #
    # tracking.to_csv(output_prefix + ".tpm_matrix.tsv",
    #                 sep="\t",
    #                 index=False,
    #                 header=True)





if __name__ == '__main__':

    start = timer()

    descrpn = """Script to filter assembled transcripts based on a minimum mean TPM expression across samples of the same condition. Corresponding individual GTF assemblies are filtered for transcripts passing this threshold and output to file"""

    parser = argparse.ArgumentParser(description=descrpn,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter, # Add defaults to end of help strings
                                     )

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
                        help="Path to list of GTF files used as input to GFFcompare merge call. This must be the list used to generate '.tracking' file and is the set of individual GTFs that will be filtered and output separately by the script")

    parser.add_argument("-m",
                        "--min-tpm",
                        type=float,
                        default=1.0,
                        help="Minimum mean TPM (across all samples) for a transcript to be retained in output GTF"
                        )

    parser.add_argument("-o","--output-suffix",
                        dest="output_suffix",
                        type=str,
                        default=".mean_tpm_filtered.gtf",
                        help="Suffix for individual filtered GTF files containing transcripts passing mean expression threshold (MUST end with '.gtf')")

    if len(sys.argv) == 1:
        parser.print_help()
        parser.exit()

    args = parser.parse_args()

    main(args.tracking, args.gtf_list, args.min_tpm, args.output_suffix)

    end = timer()

    eprint(f"Complete: took {round(end - start, 3)} s / {round((end - start) / 60, 3)} min (3 dp)")
