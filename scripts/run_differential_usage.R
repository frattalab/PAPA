
start_time <- Sys.time()
suppressPackageStartupMessages(library(optparse))

option_list <- list(make_option(c("-i", "--input-counts"),
                                type="character",
                                dest = "counts",
                                help="Path to <prefix>.counts.tsv file storing count matrix for each last exon isoform (output by tx_to_polya_quant.R)"),
                    make_option(c("-g","--le2gene"),
                                type="character",
                                help = "Path to <prefix>.le2gene.tsv file storing last exon ID to Gene ID assignment (output by get_combined_quant_gtf.py)"),
                    make_option(c("-s", "--sample-table"),
                                type="character",
                                dest="sample_table",
                                help="Path to sample table CSV file used as input to PAPA pipeline. By default the first value in the 'condition' column is taken as the 'base key'"),
                    make_option(c("-c", "--cores"),
                                type="integer",
                                default = 1,
                                help = "Number of cores to use for parallel processing (default = %default)"
                                ),
                    make_option(c("-m", "--min-mean-count"),
                                type = "integer",
                                dest = "min_mean_count",
                                default = 10,
                                help = "Min mean count in any condition for an isoform to be retained for differential usage analysis (default = %default)"),
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

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(glue))
suppressPackageStartupMessages(library(satuRn))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(BiocParallel))

## Define set of functions to run diff usage analysis

#' Removes rows in which a column value (e.g. gene_id) only appears once in a dataframe
rm_single_iso_genes <- function(df, id_col = "gene_id") {
  
  subset(df, duplicated(df[, id_col]) | duplicated(df[ , id_col], fromLast = TRUE))
  
}

#' Match up order of rows between two dataframes
#' df - dataframe which want to re-order rows relative to match_to
#' match_to = dataframe with rownames wish to match order of df
#' id_col = name of column in df which contains values found in rownames(match_to)
match_col_order <- function(df, match_to, id_col = "isoform_id") {
  
  df[match(rownames(match_to), df[, id_col]), ]
  
}

#' Make a SummarizedExperiment object 
make_summarized_exp <- function(in_counts, in_colData, in_rowData, cond_col = "condition") {
  
  se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = in_counts),
                                                   colData = in_colData,
                                                   rowData = in_rowData)
  
  # Specify design formula - extract variable of interest (condition) from colData
  metadata(se)$formula <- ~ 0 + as.factor(colData(se)[, cond_col])
  
  return(se)
}

#' Return a design matrix for simple 2-group comparison
dtu_design <- function(condition, base_cond, treat_cond) {
  
  desgn <- model.matrix(~ 0 + as.factor(condition))
  colnames(desgn) <- levels(condition)
  
  return(desgn)
  
}

#' Return a contrasts matrix for comparison of interest
dtu_contrasts <- function(desgn, base_cond, treat_cond) {
  
  x <- paste(c(treat_cond, base_cond), collapse = "-")
  
  contrsts <- limma::makeContrasts(contrasts = x,
                                   levels = desgn)
  
  return(contrsts)
  
}


#' Wrapper function to run different
dtu_wrapper <- function(counts, iso2gene, sample_tbl, min_mean_reads = 10, cores = 1, id_col = "gene_id", iso_id_col = "isoform_id", cond_col = "condition", base_cond = "CTL", treat_cond = "KD") {
  
  
  # First double check if any single-isoform genes in input isoform_id | gene_id df
  message("Filtering single-isoform genes...")
  print(glue::glue("Number of isoforms pre-filtering - {nrow(iso2gene)}"))
  iso2gene <- rm_single_iso_genes(iso2gene, id_col)
  
  print(glue::glue("Number of isoforms post-filtering - {nrow(iso2gene)}"))
  
  counts <- counts[which(rownames(counts) %in% iso2gene[, iso_id_col]), ]

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
  conds <- list(sample_tbl[sample_tbl[, cond_col] == base_cond, "sample_name"],
                sample_tbl[sample_tbl[, cond_col] == treat_cond, "sample_name"]
      )
  
  names(conds) <- c(base_cond, treat_cond)

  
  # Calculate means for each condition      
  # Matrix with cols base_cond | treat_cond (mean count)
  means_norm_counts <- sapply(X = conds, FUN = function(x) {rowMeans(norm_counts[, x])})

  # Get the max condition mean count for each sample
  # Just a vector of max counts
  maxs <- rowMaxs(means_norm_counts)

  # Filter original count matrices for rows with a condition mean > 10
  counts <- counts[maxs > min_mean_reads, ]
  
  if (length(counts) == 0) {
    stop(glue::glue("Filtered counts matrix is empty - consider if the specified min count filter - {min_mean_reads} - is too stringent"))
  }
  
  #2. Match up iso2gene to filtered counts matrix
  iso2gene <- iso2gene[which(iso2gene[, iso_id_col] %in% rownames(counts)), ]
  
  message(glue::glue("Number of isoforms after filtering for min expression - {nrow(iso2gene)}"))
  
  # Filter out single isoform genes again.
  iso2gene <- rm_single_iso_genes(iso2gene, id_col)
  
  message(glue::glue("Number of isoforms after removing single-isoform genes - {nrow(iso2gene)}"))
  message(glue::glue("Number of multi-isoform genes to be tested - {length(unique(iso2gene$gene_id))}"))

  # Remove single isoform genes from count matrix
  counts <- counts[which(rownames(counts) %in% iso2gene[, iso_id_col]), ]
  
  # Now make sure isoforms are in same order as count matrix
  iso2gene <- iso2gene[match(rownames(counts), iso2gene[, iso_id_col]), ]
  
  
  #3. Make Summarised experiment object
  
  se <- make_summarized_exp(counts,
                            in_colData = sample_tbl,
                            in_rowData = iso2gene,
                            cond_col = cond_col)
  
  # 4. fit DTU parameters across condition
  form_dtu <- as.formula(paste("~ 0", cond_col, sep = " + "))

  message(glue::glue("DTU formula - {paste(as.character(form_dtu), collapse = ' ')}"))

  message("Estimating parameters for differential usage...")
  
  if (cores == 1) {
    
    system.time({ 
      se <- satuRn::fitDTU(object = se,
                           formula = form_dtu,
                           parallel = FALSE,
                           verbose= TRUE)
      }
      )
  
  } else {
    
    system.time({ 
      se <- satuRn::fitDTU(object = se,
                           formula = form_dtu,
                           parallel = TRUE,
                           BPPARAM = BiocParallel::MulticoreParam(cores),
                           verbose = TRUE)
    })
    
  }
  
  # 5. Set up test for DTU 
  message("Testing for differential usage across specified contrast...")

  dsgn <- dtu_design(sample_tbl[, cond_col], base_cond, treat_cond)
  
  cntrsts <- dtu_contrasts(dsgn, base_cond, treat_cond)
  
  se <- testDTU(object = se,
                contrasts = cntrsts
  )
  
  
  out_list <- list("se" = se,
                   "filtered_counts" = counts,
                   "filtered_iso2gene" = iso2gene,
                   "design" = dsgn,
                   "contrasts" = cntrsts
                   )
  
  return(out_list)

}

# Read in and process data

sample_tbl <- read.table(opt$sample_table,
                         header = TRUE,
                         sep = ",",
                         colClasses = c("condition"="factor"))

counts <- read.table(opt$counts,
                     sep = "\t",
                     header = TRUE,
                     row.names = "le_id")

# count matrix as output by tx_to_polya_quant.R contains gene_id column
# Remove it if it's there...

if ("gene_id" %in% colnames(counts)) {
  counts <- counts[, setdiff(colnames(counts), c("gene_id"))]
}

# print(colnames(counts))
# print(sample_tbl$sample_name)

le2gene <- read.table(opt$le2gene,
                      header = TRUE,
                      sep = "\t",
                      stringsAsFactors = FALSE
                      )

# Must be isoform_id | gene_id for SatuRn
colnames(le2gene) <- c("isoform_id", "gene_id")

# Get base condition from 1st row in sample table
base_key <- sample_tbl[, "condition"][1]
treat_key <- sample_tbl[sample_tbl[, "condition"] != base_key, "condition"][1]

message(glue("Base condition key - {base_key}"))
message(glue("Treatment condition key - {treat_key}"))
        
message("Performing differential usage analysis...")

dtu <- dtu_wrapper(counts = counts,
                   iso2gene = le2gene,
                   sample_tbl = sample_tbl,
                   min_mean_reads = opt$min_mean_count,
                   cores = opt$cores,
                   id_col = "gene_id",
                   iso_id_col = "isoform_id",
                   cond_col = "condition",
                   base_cond = base_key,
                   treat_cond = treat_key
                   )


message("Differential usage analysis complete, writing output table and Rdata objects...")

out_res <- paste(opt$output_prefix, ".results.tsv", sep = "")
out_rda <- paste(opt$output_prefix, ".image.RData", sep = "")

message(glue("Writing results dataframe to {out_res}"))

res <- rowData(dtu$se)[[paste("fitDTUResult_", treat_key, "-", base_key, sep = "")]] %>%
                        rownames_to_column("isoform_id")
write.table(res,
            file = out_res,
            quote = FALSE,
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE
            )

message(glue("Writing workspace image to {out_rda}"))

save.image(out_rda)

message("Script complete - took...")
comp <- Sys.time() - start_time
print(comp)

