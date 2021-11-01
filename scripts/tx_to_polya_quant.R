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
                    make_option(c("-t","--tx2pas"),
                                type="character",
                                help = "Path to <prefix>'.tx2pas.tsv file storing transcript ID to PolyASite ID assignment output by assign_tx_to_pas.py"),
                    make_option(c("-g","--tx2gene"),
                                type="character",
                                help = "Path to <prefix>'.tx2pas.tsv file storing transcript ID to Gene ID assignment output by assign_tx_to_pas.py"),
                    make_option(c("-o", "--output-file"),
                                dest = "output_file", 
                                default = "summarised_pas_quantification.tsv",
                                help = "Name of output TSV file storing summarise per-poly(A) site quantification ([default= %default])")
                    )

opt_parser <- OptionParser(option_list = option_list)
opt <-parse_args(opt_parser)


if (length(commandArgs(trailingOnly = TRUE)) == 0) {
  print_help(opt_parser)
  stop()
}

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

txi.pas <- summarizeToGene(txi.tx, tx2pas, countsFromAbundance = "lengthScaledTPM")

# Add 'poly(A) site IDs to counts df
pas_counts <- cbind(pas_id = rownames(txi.pas$counts), as.data.frame(txi.pas$counts))

# Return gene ID to to counts df (so have polyASite ID & gene IDs to group events)
# gsub("_[0-9]+$", "", head(pas_counts$pas_id))

pas_counts <- merge(pas_counts,
                    merge(tx2gene, tx2pas, by = "transcript_id")[,c("pas_id", "gene_id")],
                    by = "pas_id")

# reorder so pas_id | gene_id are the first two columns
pas_counts <- pas_counts[, c(1, ncol(pas_counts), seq(2, ncol(pas_counts) - 1, 1))]

# "data/tximport_pas_quantification.tsv"
write.table(pas_counts,
            file = opt$output_file,
            sep = "\t",
            quote = FALSE,
            col.names = TRUE,
            row.names = FALSE)




# input:
# tx2pas dataframe output by assign_tx_to_pas.py
# Directory storing Salmon quantification
# Input sample table from pipeline