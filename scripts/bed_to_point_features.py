#!/usr/bin/env python3


import pyranges as pr
import pandas as pd
import argparse
import os
import sys

'''
Convert PolyASite atlas BED file to a 'point features' TSV file ready for use with StrigTie -ptf



PolyASite (v2) file format (BED/TSV)
chromosome | cluster start | cluster end | cluster ID <chr:representative coord:strand> | Average tags per million |
| (col 6) strand | % samples supporting | n protocols supporting | ave expression all samples | cluster annotation |
| (col 11) polya signals

Since polyASite atlas reports poly(A) site clusters, will take the representative coordinate from cluster ID (Name column of BED file)
as the single nucleotide position

'Point-features' file format (TSV):
Column 1 - chromosome | Column 2 - coordinate | Column 3 - Strand | Column 4 - Feature type <CPAS/TSS>


'''



def main(in_bed, add_chr, out_tsv):

    bed = pr.read_bed(in_bed)

    #1. Select required columns (Chromosome, Strand, Name)
    # Chromosome & Strand are 'inbuilt', subset to Name
    # print(bed.columns)

    bed = bed[["Name"]]

    #2. Get representative coordinate from 'Name'/cluster ID
    # 1:1234:+
    bed.coordinate = bed.Name.str.split(":", expand=True)[1]

    # print(bed)

    #3. Add feature type column - all are poly(A) sites so 'CPAS'
    bed.feature = "CPAS"

    #4. add 'chr' prefix to chromosome if wanted (so matches GTF annotation)
    # PolyASite atlas uses Ensembl chromosome names ()
    bed = bed.as_df()

    if add_chr:
        bed["Chromosome"] = "chr" + bed["Chromosome"].astype(str)

    else: # no need to add prefix
        pass

    bed[["Chromosome", "coordinate", "Strand", "feature"]].to_csv(out_tsv, sep="\t", header=False, index=False)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Script to convert PolyASite v2 BED file to StringTie 'point features' TSV file")

    parser.add_argument("-i","--input-bed", type=str, default="", dest="input_bed", help="path to PolyASite atlas BED file (must be uncompressed)")
    parser.add_argument("--prefix_chr", default=False, dest="add_chr", action='store_true', help="should 'chr' prefix be added to chromosome names reported in PolyASite BED file?")
    parser.add_argument("-o", type=str, default="pas_point_features.tsv", dest="output_tsv", help="path/name of output 'point features' TSV file (default: %(default)s)")
    parser.add_argument("--overwrite", default=False, dest="overwrite", action="store_true", help="overwrite file specified by -o if it exists")

    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help()
        parser.exit()

    # check whether specified output file already exists (and clashes with --overwrite option)
    if os.path.exists(args.output_tsv):
        if not args.overwrite:
            sys.exit("specified output file - {0} - already exists. Pass --overwrite to overwrite existing file")

    else: # output file doesn't already exist (doesn't matter if overwrite set or not)
        pass

    main(args.input_bed, args.add_chr, args.output_tsv)
