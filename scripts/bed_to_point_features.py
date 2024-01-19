#!/usr/bin/env python3

#     Convert PolyASite atlas BED file to a 'point features' TSV file ready for use with StrigTie
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



def main(in_bed, add_chr, standard_bed6, out_tsv):

    bed = pr.read_bed(in_bed)

    #3. Add feature type column - all are poly(A) sites so 'CPAS'
    bed.feature = "CPAS"

    if not standard_bed6:
        #1. Select required columns (Chromosome, Strand, Name)
        # Chromosome & Strand are 'inbuilt', subset to Name
        # print(bed.columns)

        bed = bed[["Name"]]

        #2. Get representative coordinate from 'Name'/cluster ID
        # 1:1234:+
        bed.coordinate = bed.Name.str.split(":", expand=True)[1]

        # print(bed)

        bed = bed.as_df()

        #4. add 'chr' prefix to chromosome if wanted (so matches GTF/BAM file seqnamses)
        # PolyASite atlas uses Ensembl chromosome names
        if add_chr:
            bed["Chromosome"] = "chr" + bed["Chromosome"].astype(str)

        else: # no need to add prefix
            pass

        bed[["Chromosome", "coordinate", "Strand", "feature"]].to_csv(out_tsv, sep="\t", header=False, index=False)

    else:
        # Just taking start/end coords in BED col 2 & 3 as coordinates to report
        # Assuming wants 'position' format style coordinates, so coord is 1-based
        # So to convert BED start to 1-based, need to add 1
        bed.coordinate = bed.Start + 1

        bed = bed.as_df()

        # add 'chr' prefix to chromosome if wanted
        if add_chr:
            bed["Chromosome"] = "chr" + bed["Chromosome"].astype(str)

        bed[["Chromosome", "coordinate", "Strand", "feature"]].to_csv(out_tsv, sep="\t", header=False, index=False)



if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Script to convert PolyASite v2 BED file to StringTie 'point features' TSV file")

    parser.add_argument("-i","--input-bed", type=str, default="", dest="input_bed", help="path to PolyASite atlas BED file (must be uncompressed)")
    parser.add_argument("--prefix_chr", default=False, dest="add_chr", action='store_true', help="should 'chr' prefix be added to chromosome names reported in PolyASite BED file?")
    parser.add_argument("--standard-bed6", default=False, action='store_true', help="Whether to consider input BED file as a 'conventional' BED6 file and report the Start column as the PAS coordinate (in 1-based)")
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

    main(args.input_bed, args.add_chr, args.standard_bed6, args.output_tsv)
