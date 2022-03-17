suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tximport))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(glue))

option_list <- list(make_option(c("-s", "--sample-table"),
                                type="character",
                                dest="sample_table",
                                help="Path to sample table CSV file used as input to PAPA pipeline. By default the first value in the 'condition' column is taken as the 'base key'"),
                    make_option(c("-d", "--salmon-dir"),
                                dest = "salmon_dir",
                                type = "character",
                                help = "Path to top-level directory under which per-sample Salmon quantification outputs are stored"),
                    make_option(c("-t","--tx2le"),
                                type="character",
                                help = "Path to <prefix>'.tx2le.tsv file storing transcript ID to last exon ID assignment (output by get_combined_quant_gtf.py)"),
                    make_option(c("-g","--le2gene"),
                                type="character",
                                help = "Path to <prefix>.le2gene.tsv file storing last exon ID to Gene ID assignment (output by get_combined_quant_gtf.py)"),
                    make_option(c("-o", "--output-prefix"),
                                dest = "output_prefix",
                                default = "summarised_pas_quantification",
                                help = "Prefix to names of output files storing per-event summarised counts (<output_prefix>.counts.tsv), per-event summarised TPM values (<output_prefix>.tpm.tsv), per-event PPAU values (<output_prefix>.ppau.tsv) and per-gene summarised TPM values (<output_prefix>.gene_tpm.tsv) ([default= %default])")
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

n_cond <- length(unique(sample_tbl$condition))

if ( n_cond != 2) {

  stop(paste("Sample table must only contain 2 distinct conditions -",
             n_cond,
             "were found"))
}

output_prefix <- opt$output_prefix

sample_names <- sample_tbl$sample_name

base_key <- sample_tbl$condition[1]
treat_key <- unique(sample_tbl$condition[sample_tbl$condition != base_key])

base_sample_names <- sample_tbl[sample_tbl$condition == base_key, "sample_name"]
treat_sample_names <- sample_tbl[sample_tbl$condition != base_key, "sample_name"]

quant_paths <- file.path(opt$salmon_dir, sample_names, "quant.sf")
names(quant_paths) <- sample_names

if (isFALSE(all(file.exists(quant_paths)))) {
  stop(paste("Not all provided paths to quant.sf files exist, has -d/--salmon-dir been provided correctly?",
             paste(quant_paths, file.exists(quant_paths), sep = " - "), sep = " ")
       )
}

#2 read in tx2pas
tx2pas <- read.table(opt$tx2le,
                     header = T,
                     sep = "\t",
                     stringsAsFactors = FALSE)

tx2gene <- read.table(opt$le2gene,
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
message(paste("Writing poly(A) isoform counts file to - ", output_prefix, ".counts.tsv", sep = ""))

write.table(pas_counts,
            file = paste(output_prefix, ".counts.tsv", sep=""),
            sep = "\t",
            quote = FALSE,
            col.names = TRUE,
            row.names = FALSE)

message(paste("Writing poly(A) isoform tpm file to - ", output_prefix, ".tpm.tsv", sep = ""))

write.table(pas_tpm,
            file = paste(output_prefix, ".tpm.tsv", sep = ""),
            sep = "\t",
            quote = FALSE,
            col.names = TRUE,
            row.names = FALSE)

# Generate PPAU matrices
# Rows = last exon ID,
# columns = PPAU values for each sample (sample_name as col names),
# mean_{condition} for every provided condition &
# delta_PPAU_{condition}_{base} for every provided contrast
# Note: EVERY PROVIDED CONSTRAST IS CURRENTLY JUST A SINGLE ONE.

message("Calculating per-sample PPAU for each event...")

# Generate per-sample PPAU for each last exon/event
## PPAU = isoform TPM / sum(isoform TPM for gene)
pas_ppau <- pas_tpm %>%
  # Put sample columns of TPM values into rows
  pivot_longer(all_of(sample_names),
               names_to = "sample_name",
               values_to = "TPM") %>%
  # Calculate total TPM for each gene in each sample
  group_by(gene_id, sample_name) %>%
  mutate(total_TPM = sum(TPM)) %>%
  ungroup() %>%
  # Calculate PPAU (isoform expression fraction of total expression)
  mutate(PPAU = TPM / total_TPM)

# Get a matrix/df of event ID as rows and PPAU for each sample as columns
ppau <- pas_ppau %>%
  select(-c(TPM, total_TPM)) %>%
  pivot_wider(names_from = sample_name, values_from = PPAU)

# Also get a matrix/df of gene_id as rows and total TPM for each sample as columns
ppau_tot <- pas_ppau %>%
  select(-c(!!sym(isoform_id_col), TPM, PPAU)) %>%
  # drop duplicate rows by gene and sample (total TPM for sample same value in all rows)
  distinct(gene_id, sample_name, .keep_all = TRUE) %>%
  # Convert sample_name gene TPMs to columns of sample name + gene TPM
  pivot_wider(names_from = sample_name, values_from = total_TPM)

# Calculate mean PPAU per-condition & delta PPAU between treatment and base conditions
message("Calculate condition-means & delta PPAU for specified contrasts...")
ppau <- ppau %>%
  # 0 / 0 returns NaN (i.e. gene has no expression)
  mutate(across(all_of(sample_names), ~ replace_na(., 0))) %>%
  rowwise() %>%
  mutate("mean_PPAU_{base_key}" := mean(c_across(all_of(base_sample_names))),
         "mean_PPAU_{treat_key}" := mean(c_across(all_of(treat_sample_names))),
         "delta_PPAU_{treat_key}_{base_key}" := !!sym(glue("mean_PPAU_{treat_key}")) - !!sym(glue("mean_PPAU_{base_key}"))) %>%
  as_tibble()

# head(ppau)

# Calculate mean & median gene TPM per condition
message("Calculating mean & median gene TPMs for each condition...")
ppau_tot <- ppau_tot %>%
  rowwise() %>%
  mutate("mean_gene_TPM_{base_key}" := mean(c_across(all_of(base_sample_names))),
         "mean_gene_TPM_{treat_key}" := mean(c_across(all_of(treat_sample_names))),
         "median_gene_TPM_{base_key}" := median(c_across(all_of(base_sample_names))),
         "median_gene_TPM_{treat_key}" := median(c_across(all_of(treat_sample_names)))
         ) %>%
  as_tibble()

# head(ppau_tot)

message("Writing per-event PPAU matrix to {output_prefix}.ppau.tsv")
write.table(ppau,
            file = glue("{output_prefix}.ppau.tsv"),
            #paste(output_prefix, ".tpm.tsv", sep = ""),
            sep = "\t",
            quote = FALSE,
            col.names = TRUE,
            row.names = FALSE)

message("Writing per-gene TPM matrix to {output_prefix}.gene_tpm.tsv")
write.table(ppau_tot,
            file = glue("{output_prefix}.gene_tpm.tsv"),
            #paste(output_prefix, ".tpm.tsv", sep = ""),
            sep = "\t",
            quote = FALSE,
            col.names = TRUE,
            row.names = FALSE)
