suppressPackageStartupMessages(library(optparse))

option_list <- list(make_option(c("-i", "--input-tsv"),
                                type="character",
                                dest = "results",
                                help="Path to <prefix>.results.tsv file storing results table from SatuRn differential usage analysis (output by run_differential_usage.R)"),
                    make_option(c("-a","--annot-info"),
                                type="character",
                                dest="annot_info",
                                help = "Path to <prefix>.info.tsv file storing last exon ID annotation information (output by get_combined_quant_gtf.py)"),
                    make_option(c("-p", "--ppau"),
                                type="character",
                                help = "path to <prefix>.ppau.tsv file storing matrix of sample and condition-wise PPAU values for each last exon isoform"),
                    make_option(c("-o", "--output-prefix"),
                                type = "character",
                                dest = "output_prefix",
                                default="differential_usage.results",
                                help = "Prefix to name of augmented SatuRn results output table <prefix>.processed.tsv (default = %default)")
)

opt_parser <- OptionParser(option_list = option_list)

if (length(commandArgs(trailingOnly = TRUE)) == 0) {
  print_help(opt_parser)
  stop()
}

opt <- parse_args(opt_parser)

res_path <- opt$results
info_path <- opt$annot_info
ppau_path <- opt$ppau
output_prefix <- opt$output_prefix

suppressPackageStartupMessages(library(dplyr))

# read in input files

# satuRn results table
res <- read.table(file = res_path,
                  header = TRUE,
                  sep = "\t",
                  stringsAsFactors = F
                  )

info <- read.table(file = info_path,
                  header = TRUE,
                  sep = "\t",
                  stringsAsFactors = F
                  )

ppau <- read.table(file = ppau_path,
                   header = TRUE,
                   sep = "\t",
                   stringsAsFactors = F
)

# isoform_id is equivalent to le_id - rename to match with info
res <- rename(res, le_id = isoform_id)

# info_path has mixed case column names - convert all to lower case to make tidier
info <- rename_with(info, tolower)

# reorder cols so coordinate info comes at the end (just move annot_status before chromosome)
info <- relocate(info, annot_status, .before = chromosome)

# The same le_id can have multiple transcripts (& multiple coords) contributing to it
# Collapse these to a single row for each LE (saves the results table getting over-duplicated)

# vector of colnames from info that expect to have a single value for each le_id
exp_uniq_cols <- c("gene_id", "gene_name", "event_type", "annot_status", "chromosome", "strand")
# vector of colnames from info where expect le_id to have multiple distinct values across rows
exp_collapse_cols <- c("transcript_id", "start", "end")

# TODO: implement check for whether le_id does not have 1 value for exp_uniq_cols (if so separate, parse separately)

# n_le_ids <- n_distinct(info$le_id)
# info %>%
#   group_by(le_id) %>%
#   # for each le, count number of unique values in each col (expect = 1)
#   summarise(across(all_of(exp_uniq_cols),
#             n_distinct)
#             )

info_clpsd <- info %>%
  group_by(le_id) %>%
  summarise(across(all_of(exp_uniq_cols),
                               unique),
            across(all_of(exp_collapse_cols),
                   ~ paste(.x, collapse = ",")
                   ),

            )

# move chromosome & strand back next to start & end coords
info_clpsd <- relocate(info_clpsd, chromosome, strand, .before = start)


# now join info with the results tbl
res_jned <- left_join(res, info_clpsd, by = "le_id")
# res_jned

# now join with minimal PPAU information (means & deltas for each group)
res_jned <- left_join(res_jned,
                      select(ppau, le_id, starts_with("mean"), starts_with("delta")),
                      by = "le_id")

write.table(x = res_jned,
            file = paste(output_prefix, ".processed.tsv", sep = ""),
            sep = "\t",
            col.names = T,
            row.names = F, 
            quote = F
            )



