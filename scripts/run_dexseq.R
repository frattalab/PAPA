library(optparse)


##args

option_list <- list(make_option(c("-i", "--input-counts"),
                                type="character",
                                dest = "counts",
                                help="Path to <prefix>.tx_counts.tsv file storing count matrix for each APA site (output by qapa_tximport.R)"),
                    make_option(c("-g","--tx2gene"),
                                type="character",
                                help = "Path to <prefix>.tx2gene.tsv file storing transcript ID to Gene ID assignment (output by generate_tx2apaid_tbls.py)"),
                    make_option(c("-s", "--sample-table"),
                                type="character",
                                dest="sample_table",
                                help="Path to sample table CSV file used as input to QAPA pipeline"),
                    make_option(c("-f", "--formulas"),
                                type = "character",
                                dest = "formulas",
                                help = "Path to 2-line TXT file denoting full and reduced model formulas for differential usage analysis"
                                ),
                    make_option(c("--base-key"),
                                type="character",
                                dest="base_key",
                                help="Value in 'condition' column of sample table which denotes samples belonging to the 'base' group for contrasts."),
                    make_option(c("--contrast-key"),
                                type="character",
                                dest="contrast_key",
                                help="Value in 'condition' column of sample table which denotes samples belonging to the 'contrast' group for contrasts."),
                    make_option(c("-n", "--constrast-name"),
                                type = "character",
                                dest = "contrast_name",
                                default = "contrast_vs_base",
                                help = "Name of contrast to identify comparison of interest tested with current script to be added as column to end of results table (default = %default)"),
                    make_option(c("--condition-col"),
                                type = "character",
                                dest = "condition_col",
                                default = "condition",
                                help = "Name of contrast to identify comparison of interest tested with current script to be added as column to end of results table (default = %default)"),
                    make_option(c("-c", "--cores"),
                                type="integer",
                                default = 1,
                                help = "Number of cores to use for parallel processing (default = %default)"),
                    make_option(c("-m", "--min-mean-count"),
                                type = "integer",
                                dest = "min_mean_count",
                                default = 10,
                                help = "Min mean count in any condition for an isoform to be retained for differential usage analysis (default = %default)"),
                    make_option(c("-r", "--min-relative-usage"),
                                type = "numeric",
                                dest = "min_rel_usage",
                                default = 0.05,
                                help = "Min mean relative usage in any condition for an isoform to be retained for differential usage analysis (default = %default)"),
                    make_option(c("-o", "--output-prefix"),
                                type = "character",
                                dest = "output_prefix",
                                default="differential_usage",
                                help = "Prefix to names of output files - <prefix>.image.RData for the Rdata object, <prefix>.results.tsv for SatuRn output dataframe (default = %default)")
)

opt_parser <- OptionParser(option_list = option_list)

if (length(commandArgs(trailingOnly = TRUE)) == 0) {
  print_help(opt_parser)
  stop()
}

opt <- parse_args(opt_parser)


#
counts_path <- opt$counts
sample_tbl_path <- opt$sample_table
tx2gene_path <- opt$tx2gene
formulas_path <- opt$formulas
condition_col <- opt$condition_col
base_key <- opt$base_key
contrast_key <- opt$contrast_key
contrast_name <- opt$contrast_name
min_rel_usage <- opt$min_rel_usage
min_mean_reads <- opt$min_mean_count
cores <- opt$cores
output_prefix <- opt$output_prefix

library(SummarizedExperiment)
library(DEXSeq)
library(readr)
library(dplyr)
library(tidyr)
library(tibble)
library(glue)
library(stringr)

### FUNCTIONS


#' filter isoforms/events with for minimum mean count (after adjusting for library size composition) in either base/contrast samples
min_mean_count_filter <- function(counts, base_samples, treat_samples, min_mean_reads = 10) {
  #1. First perform low count filter to remove isoforms with low expressiom
  # Just pick isoforms that have an library size adjusted mean count of at least n reads in either of the conditions
  # Using DESeq2's normalisation method (so count = count / size factor for sample)

  message(glue::glue("Filtering isoforms for minimum mean count in any condition of at least - {min_mean_reads}"))

  # Named vector of size factors for each sample in matrix
  sfs <- DESeq2::estimateSizeFactorsForMatrix(counts)


  # Divide each count column (sample) by its corresponding size factor
  norm_counts <- sweep(x = counts,
                       MARGIN = 2, # operate on columns
                       STATS = sfs,
                       FUN = '/')

  # Get list of samples (column names) for each condition
  conds <- list(base_samples, treat_samples)

  names(conds) <- c("base", "contrast")


  # Calculate means for each condition
  # Matrix with cols base_cond | treat_cond (mean count)
  means_norm_counts <- sapply(X = conds, FUN = function(x) {rowMeans(norm_counts[, x])})

  # Get the max condition mean count for each sample
  # Just a vector of max counts
  maxs <- rowMaxs(means_norm_counts)

  # Filter original count matrices for rows with a condition mean > min_mean_reads in either condition
  counts <- counts[maxs > min_mean_reads, ]

  message(glue::glue("Number of isoforms after filtering of min mean reads in either condition - {nrow(counts)}"))

  return(counts)
}

#' Filter count matrix for isoforms with at least a minimum mean relative usage in either base or contrast condition
min_rel_filter <- function(counts, tx2gene, base_samples, treat_samples, min_rel_usage = 0.05) {

  if (!("le_id" %in% colnames(counts))) {
    counts_df <- rownames_to_column(counts, "le_id")
  } else {
    counts_df <- counts
  }

  message(glue::glue("Filtering isoforms for minimum fractional mean relative usage in any condition of at least - {min_rel_usage}"))

  # join in gene ID
  counts_df <- dplyr::left_join(counts_df, tx2gene, by = "le_id")

  # convert to long format, add condition
  counts_long <- tidyr::pivot_longer(counts_df,
                              cols = all_of(c(base_samples, treat_samples)),
                              names_to = "sample_name",
                              values_to = "count") %>%
    dplyr::mutate(condition = if_else(sample_name %in% base_samples, "base", "contrast"))

  # calc per-sample rel usage for each tx
  rel_use <- counts_long %>%
    dplyr::group_by(gene_id, sample_name) %>%
    dplyr::mutate(rel_usage = count / sum(count)) %>%
    dplyr::ungroup()

  # calculate condition-wise mean relative usage for each tx
  mean_rel_use <- rel_use %>%
    dplyr::group_by(le_id, condition) %>%
    dplyr::summarise(mean_rel_usage = mean(rel_usage)) %>%
    dplyr::ungroup()

  # calculate max mean relative usage
  max_mean_rel_use <- mean_rel_use %>%
    dplyr::group_by(le_id) %>%
    dplyr::summarise(max_mean_rel_usage = max(mean_rel_usage))

  pass_ids <- dplyr::filter(max_mean_rel_use, max_mean_rel_usage > min_rel_usage) %>% dplyr::pull(le_id)



  if (!("le_id" %in% colnames(counts))) {
    out_counts <- counts[which(rownames(counts) %in% pass_ids), ]
  } else {
    out_counts <- dplyr::filter(counts, le_id %in% pass_ids)
  }

  message(glue::glue("Number of isoforms after filtering of min mean reads in either condition - {nrow(out_counts)}"))

  return(out_counts)

}

#' Remove genes that only have a single isoform present in the count matrix
#' (Removes rows in which a column value (e.g. gene_id) only appears once in a dataframe)
rm_single_iso_genes <- function(df, id_col = "gene_id") {

  out_df <- df %>%
    dplyr::group_by(!!sym(id_col)) %>%
    dplyr::filter(dplyr::n_distinct(le_id) > 1) %>%
    dplyr::ungroup()

  message(glue::glue("Number of isoforms after filtering out single isoform genes - {nrow(out_df)}"))
  message(glue::glue("Number of genes with multiple isoforms - {n_distinct(out_df[[id_col]])}"))

  return(out_df)
}



### Main script logic

counts <- read_tsv(counts_path)
sample_tbl <- read_csv(sample_tbl_path)
tx2gene <- read_tsv(tx2gene_path)
formulas <- read_lines(formulas_path)

# first filter sample_tbl for specified keys
sample_tbl <- dplyr::filter(sample_tbl, !!sym(condition_col) %in% c(base_key, contrast_key))

# convert condition col to factor, setting base_key as reference level
sample_tbl <- mutate(sample_tbl, !!sym(condition_col) := factor(!!sym(condition_col), levels = c(base_key, contrast_key)))

# Check if additional co-variate columns in reduced formula - these should also be converted to factor
# 2nd element in formulas is reduced model
red_form_spl <- unlist(str_split(formulas[2], " ", simplify = F))
red_form_cov <- red_form_spl[str_detect(red_form_spl, ":")]

if (length(red_form_cov) > 0) {
  # extract covariate column names from formula
  # e.g. pair:exon -> pair
  covs <- str_remove_all(red_form_cov, ":exon")
  glue::glue("Extracted covariate column names (space separated) - {paste(covs, collapse = ' ')}")

  # convert covariate cols to factors
  sample_tbl <- mutate(sample_tbl, across(all_of(covs), as.factor))
}


# Filter counts matrices for lowly expressed isoforms and those with low relative usage
base_samples <- filter(sample_tbl, !!sym(condition_col) == base_key) %>% pull(sample_name)
contrast_samples <- filter(sample_tbl, !!sym(condition_col) == contrast_key) %>% pull(sample_name)

# filter counts matrix for base/contrast samples
# also remove gene_id column which is also present in count matrix output by tx_to_polya_quant.R - will otherwise break min_mean_count_filter()
counts_filt <- counts[, colnames(counts) %in% c("le_id", base_samples, contrast_samples)]

counts_filt <- min_mean_count_filter(column_to_rownames(counts_filt, "le_id"),
                      base_samples,
                      contrast_samples,
                      min_mean_reads = min_mean_reads) %>%
  min_rel_filter(., tx2gene, base_samples, contrast_samples, min_rel_usage = min_rel_usage) %>%
  rownames_to_column("le_id") %>%
  left_join(tx2gene, by = "le_id") %>%
  rm_single_iso_genes()

# Now filter tx2gene so only retain transcripts found in count matrix
tx2gene_filt <- filter(tx2gene, le_id %in% pull(counts_filt, le_id)) %>%
  # https://stackoverflow.com/questions/67684167/how-to-use-match-in-dplyr-to-order-a-column-based-on-an-external-vector
  arrange(ordered(le_id, counts_filt$le_id))

# Generate summarised experiment object
message("Generating DEXSeq object...")
dexseq_se <- SummarizedExperiment(assays = list( counts = round(column_to_rownames(dplyr::select(counts_filt, - gene_id), "le_id"))),
                                           colData = sample_tbl,
                                           rowData = rename(tx2gene_filt, tx_name = le_id)
)


# Generate DEXSeq dataset object ready for analysis
dxd <- DEXSeqDataSetFromSE(dexseq_se,
                           design = as.formula(formulas[1]) # full model
                           )

# DEXSeq analysis
message("Estimating size factors...")
system.time({dxd <- estimateSizeFactors( dxd ) })

message("Estimating dispersions..")
system.time({dxd <- estimateDispersions( dxd, BPPARAM = MulticoreParam(cores))})

message("Testing for differential usage...")
system.time({dxd <- testForDEU(dxd,
                               reducedModel = as.formula(formulas[2]),
                               BPPARAM=MulticoreParam(cores))
}
)

message("Estimating exon usage fold changes...")
system.time({ dxd <- estimateExonFoldChanges( dxd,
                                              fitExpToVar= condition_col,
                                              BPPARAM = MulticoreParam(cores),
                                              denominator = base_key)
}
)

message("Extracting results table...")
dxr = DEXSeqResults( dxd )

message("Calculating per-gene q-values")
# named vector of geneID - qval
qvals <- perGeneQValue(dxr)
qvals_df <- enframe(qvals, "groupID", "gene.qvalue")

message("Reformatting results table...")
# generate a tidier output table
dxr_out <- as_tibble(dxr, rownames = "binID") %>%
  mutate(le_id = unlist(transcripts)) %>% # QAPA only produces 1 ID per bin, so simple unlist produces 1 element vector for each row
  select(- transcripts, -genomicData) # don't provide ranges as input so column unncecessary

dxr_out <- left_join(dxr_out, qvals_df, by = "groupID")

# separate count data from DEXSeq results - too much extra info
counts_out <- dplyr::select(dxr_out, binID, groupID, featureID, le_id, starts_with("countData."))
dxr_out <- dplyr::select(dxr_out, -starts_with("countData."))
dxr_out <- mutate(dxr_out, contrast_name = contrast_name)

# output
save.image(paste(output_prefix, ".Rdata", sep = ""))
write_tsv(counts_out, paste(output_prefix, ".filtered_countData.tsv", sep = ""), col_names = T)
write_tsv(dxr_out, paste(output_prefix, ".results.tsv", sep = ""), col_names = T)
