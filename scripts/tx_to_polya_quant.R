library(optparse)
library(tximport)

option_list <- list(make_option(c("-s", "--sample-table"),
                                type="character",
                                dest="sample_table",
                                help="Path to sample table CSV file used as input to PAPA pipeline"),
                    make_option(c("-d", "--salmon-dir"),
                                dest = "salmon_dir",
                                type = "character",
                                help = "Path to top-level directory under which per-sample Salmon quantification outputs are stored"),
                    make_option(c("-t","--tx2le"),
                                type="character",
                                help = "Path to <prefix>'.tx2le.tsv file storing transcript ID to last exon ID assignment (output by get_combined_quant_gtf.py)"),
                    make_option(c("-g","--tx2gene"),
                                type="character",
                                help = "Path to <prefix>'.tx2le.tsv file storing transcript ID to Gene ID assignment (output by get_combined_quant_gtf.py)"),
                    make_option(c("-o", "--output-prefix"),
                                dest = "output_prefix",
                                default = "summarised_pas_quantification",
                                help = "Prefix to names of output files storing per-poly(A) site summarised counts (<output_prefix>.counts.tsv) or TPM values (<output_prefix>.tpm.tsv) ([default= %default])")
                    )

opt_parser <- OptionParser(option_list = option_list)

if (length(commandArgs(trailingOnly = TRUE)) == 0) {
  print_help(opt_parser)
  stop()
}

opt <- parse_args(opt_parser)

#1. Read in sample table & construct paths to all Salmon quantification files

# opt_parser$sample_table
sample_tbl <- read.table(opt$sample_table,
                         header = TRUE,
                         sep = ",",
                         stringsAsFactors = FALSE)

sample_names <- sample_tbl$sample_name

quant_paths <- file.path(opt$salmon_dir, sample_names, "quant.sf")
names(quant_paths) <- sample_names

if (isFALSE(all(file.exists(quant_paths)))) {
  stop(paste("Not all provided paths to quant.sf files exist, has -d/--salmon-dir been provided correctly?",
             paste(quant_paths, file.exists(quant_paths), sep = " - "), sep = " ")
       )
}

#2 read in tx2pas
tx2pas <- read.table(opt$tx2pas,
                     header = T,
                     sep = "\t",
                     stringsAsFactors = FALSE)

tx2gene <- read.table(opt$tx2gene,
                      header = T,
                      sep = "\t",
                      stringsAsFactors = FALSE)

txi.tx <- tximport(quant_paths,
                   type = "salmon",
                   txOut = TRUE,
                   tx2gene = tx2gene,
                   countsFromAbundance = "lengthScaledTPM")

#Sum abundance for transcripts sharing the same polyA site
txi.pas <- summarizeToGene(txi.tx, tx2pas, countsFromAbundance = "lengthScaledTPM")

# print(head(txi.pas$abundance))
# print(head(txi.pas$counts))


# Add 'poly(A) site IDs to counts df
# could be pas_id or le_id depending on input file type
isoform_id_col <- colnames(tx2pas)[-1]
message(paste("isoform ID column:", isoform_id_col, sep=" "))

pas_counts <- cbind(rownames(txi.pas$counts), as.data.frame(txi.pas$counts))
pas_tpm <- cbind(rownames(txi.pas$abundance), as.data.frame(txi.pas$abundance))
# print(head(pas_counts))

colnames(pas_counts) <- c(isoform_id_col, colnames(as.data.frame(txi.pas$counts)))
colnames(pas_tpm) <- c(isoform_id_col, colnames(as.data.frame(txi.pas$abundance)))

# print(head(pas_counts))
# Return gene ID to to counts df (so have polyASite ID & gene IDs to group events)
# gsub("_[0-9]+$", "", head(pas_counts$pas_id))

# print("tx2gene")
# print(head(tx2gene))
# print("tx2pas")
# print(head(tx2pas))

message("Adding gene ID to counts & TPM df...")

iso2gene <- merge(tx2gene, tx2pas, by = "transcript_id")[,c("gene_id", isoform_id_col)]
# print(head(iso2gene))

# print(head(pas_counts))
# print(colnames(pas_counts))
pas_counts <- merge(pas_counts,
                    iso2gene,
                    by = isoform_id_col)

pas_tpm <- merge(pas_tpm,
                 iso2gene,
                 by = isoform_id_col)

# print(head(pas_counts))

# reorder so pas_id | gene_id are the first two columns
pas_counts <- pas_counts[, c(1, ncol(pas_counts), seq(2, ncol(pas_counts) - 1, 1))]
pas_tpm <- pas_tpm[, c(1, ncol(pas_tpm), seq(2, ncol(pas_tpm) - 1, 1))]

# Just in case, remove duplicate pas_id rows from matrix
pas_counts <- pas_counts[!duplicated(pas_counts[, isoform_id_col]), ]
pas_tpm <- pas_tpm[!duplicated(pas_tpm[, isoform_id_col]), ]

# "data/tximport_pas_quantification.tsv"
message(paste("Writing poly(A) isoform counts file to - ", opt$output_prefix, ".counts.tsv", sep = ""))

write.table(pas_counts,
            file = paste(opt$output_prefix, ".counts.tsv", sep=""),
            sep = "\t",
            quote = FALSE,
            col.names = TRUE,
            row.names = FALSE)

message(paste("Writing poly(A) isoform tpm file to - ", opt$output_prefix, ".tpm.tsv", sep = ""))

write.table(pas_tpm,
            file = paste(opt$output_prefix, ".tpm.tsv", sep = ""),
            sep = "\t",
            quote = FALSE,
            col.names = TRUE,
            row.names = FALSE)
