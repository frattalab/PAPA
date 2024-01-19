#     Generate simulated FASTQ files for testing the pipeline
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

if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install(version = "3.16")
}

if (!require("polyester",quietly = T)) {
  BiocManager::install("polyester")
}

if (!require("pacman", quietly = TRUE)) {
  install.packages("pacman")
}

library(pacman)

p_load(polyester, Biostrings)


path_fasta <- "tests/test_annotations/test_transcripts.fa"
path_gtf <- "tests/test_annotations/test_transcripts_combined.gtf"

fa <- readDNAStringSet(path_fasta)
tx_names <- names(fa)

tx_fcs <- data.frame(ctl = rep(1,times = length(tx_names)),
                     kd = rep(1,times = length(tx_names)),
                     row.names = tx_names)

tx_fcs

# set the reference txipts to have no change in KD
tx_fcs[grepl("^ref_g", row.names(tx_fcs)), "kd"] <- 1
tx_fcs

# Set the first and last novel txipts to have a FC of 4
# Set the  transcript of ref_g3 to have a FC of 4
# Set the first ref & '_ext' transcript of ref_g3 FC of 5 & 10 in KD respectively
tx_p <- grep("^nov_g.*_p$", row.names(tx_fcs))
tx_up1 <- row.names(tx_fcs)[tx_p[1]]
tx_up2 <- row.names(tx_fcs)[tx_p[length(tx_p)]]
tx_up1
tx_up2

ref_trs_up <- c("ref_g3_tr_1", "ref_g3_tr_5_ext")

tx_fcs[c(tx_up1, tx_up2), "kd"] <- 4
tx_fcs

tx_fcs[ref_trs_up[1], "kd" ] <- 5
tx_fcs[ref_trs_up[2], "kd"] <- 10
tx_fcs

# 100 bp reads, PE, stranded, 20x coverage
# ~20x coverage ----> reads per transcript = transcriptlength/readlength * 20
reads_per_tx <- round(20 * width(fa) / 100)
reads_per_tx

simulate_experiment(fasta = path_fasta,
                    # gtf = path_gtf,
                    reads_per_transcript = reads_per_tx,
                    num_reps = c(3,3),
                    fold_changes = as.matrix(tx_fcs),
                    paired = TRUE,
                    strand_specific = TRUE,
                    outdir = "tests/test_sequence")

# want to make it easier to know which samples are which - will add prefix of CTL / KD for group1/group2 respecvitley
rep2group <- read.table(file = "tests/test_sequence/sim_rep_info.txt", sep = "\t", header = T)


# Want to convert simulated read FASTAs to FASTQs
sim_fastas <- list.files(path = "tests/test_sequence",
                         pattern = ".fasta$",
                         full.names = T)

for (f in sim_fastas) {
  
  outdir <- dirname(f)
  bname <- basename(f)
  # sample _01 _1 fasta
  bname_noext <- strsplit(bname, ".", fixed = TRUE)[[1]][1]
  sname <- paste(strsplit(bname, "_", fixed = TRUE)[[1]][1:2], collapse = "_")
  # print(sname)
  grpname <- paste("group", rep2group[rep2group$rep_id == sname, "group"], sep = "")
  # print(grpname)
  outpath <- file.path(outdir, paste(paste(grpname,
                                     bname_noext,
                                     sep = "_"
                                     ),                   
                                     ".fastq.gz", sep = "")
                       )
  print(outpath)
  tmp_fa <- readDNAStringSet(filepath = f)
  
  # Set dummy quality scores (high)
  # All reads are 100bp so width is same every time, but diff sampels may have diff n reads (length)
  # Setting phred score to 1e-4 (99.99% basecall accuracy ) just to have high dummy value
  tmp_quals <- rep(PhredQuality(rep(1e-4, 100)), # quals for 100bp read
                   length(tmp_fa) # quals for every read in sample
                   )
  
  tmp_fa_qual <- QualityScaledDNAStringSet(tmp_fa, quality = tmp_quals) 
  
  writeQualityScaledXStringSet(x = tmp_fa_qual,
                               filepath = outpath,
                               compress = T)
  
}


# tmp_fq <- readDNAStringSet(filepath = "tests/test_sequence/sample_01_1.fasta")
# # All reads are 100bp, need quality assignments for every position, and every read in FASTA
# # Setting phred score to 1e-4 (99.99% basecall accuracy ) just to have dummy value
# quals <- rep(PhredQuality(rep(1e-4, 100)), length(tmp_fq))
# length(tmp_fq)
# 
# 40 * width(tmp_fq)
# 
# tmp_fq_qual <- QualityScaledDNAStringSet(tmp_fq,quality = quals
#                                          )
# #quality = 1e-4*width(tmp_fq)
# writeQualityScaledXStringSet(x = tmp_fq_qual,compress = T, filepath = "tests/test_sequence/sample_01_1.fastq.gz")




