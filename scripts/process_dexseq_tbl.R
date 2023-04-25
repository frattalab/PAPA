suppressPackageStartupMessages(library(optparse))

option_list <- list(make_option(c("-i", "--input-tsv"),
                                type="character",
                                dest = "results",
                                help="Path to <prefix>.results.tsv file storing results table from DEXSeq differential usage analysis (output by run_dexseq.R)"),
                    make_option(c("-a","--annot-info"),
                                type = "character",
                                dest = "annot_info",
                                help = "Path to <prefix>.info.tsv file storing last exon ID annotation information (output by get_combined_quant_gtf.py)"),
                    make_option(c("-p","--ppau"),
                                type="character",
                                help = "path to <prefix>.ppau.tsv file storing matrix of sample and condition-wise PPAU values for each last exon isoform"),
                    make_option(c("-o", "--output-prefix"),
                                type = "character",
                                dest = "output_prefix",
                                default="differential_usage.results",
                                help = "Prefix to name of augmented DEXSeq results output table <prefix>.processed.tsv (default = %default)")
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


suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))


###

# satuRn results table
res <- read_tsv(res_path)
info <- read_tsv(info_path)
ppau <- read_tsv(ppau_path)

# info_path has mixed case column names - convert all to lower case to make tidier
info <- rename_with(info, tolower)

# reorder cols so coordinate info comes at the end (just move annot_status before chromosome)
info <- relocate(info, annot_status, .before = chromosome)

# TODO: only interested in le_ids analysed by DEXSeq - should subset for these only (e.g. sufficient expression)
info <- filter(info, le_id %in% pull(res, le_id))

# vector of colnames from info that expect to have a single value for each le_id
exp_uniq_cols <- c("gene_id", "gene_name", "event_type", "annot_status", "chromosome", "strand")
# vector of colnames from info where expect le_id to have multiple distinct values across rows
exp_collapse_cols <- c("transcript_id", "start", "end")

# Check whether le_id does not have 1 value for exp_uniq_cols (if so separate, parse separately)
# summarise(across(x, ~ f)) expects every element in x to produce an identical length/size output from the function f (I think)
# If expected unique cols have > 2 different counts of unique values then will get an error with mismatched sizes
# these le_ids will not have their values collapsed for now...

# NB: this is often caused by 'event_type', any cols where column values are collapsed if multiple annotated partner txs (e.g. 'gene_name').
# TODO: consider updating how retain annotation information for novel events of output GTFs... (separate df?)
# TODO: really expect event_type to be unique for a le_id (which is just combination of overlapping isoforms)?

# For these le_ids, info rows will not be collapsed to a single row
n_exp_uniq <- info %>%
  group_by(le_id) %>%
  # for each le, count number of unique values in each col (expect = 1)
  summarise(across(all_of(exp_uniq_cols),
            n_distinct),
            ) %>%
  # le_id | exp_uniq_cols ...
  # Keep le_ids with > 1 value in uniq cols (potential troublemakers...)
  filter(if_any(all_of(exp_uniq_cols), ~ .x != 1)) %>%
  pivot_longer(cols = -le_id,
               names_to = "col_name",
               values_to = "n_uniq_vals") %>%
  group_by(le_id) %>%
  # Number of unique per-column unique value counts 
  mutate(n = n_distinct(n_uniq_vals)) %>%
  # If > 2 then parsing below will break
  filter(n > 2)

if (length(n_exp_uniq) != 0) {
  problem_le_ids <- unique(n_exp_uniq$le_id)
} else {
  print("All expected columns have a single unique value for all le_ids")
  problem_le_ids <- c()
}

info_gd <- filter(info, !le_id %in% problem_le_ids)
info_bd <- filter(info,  le_id %in% problem_le_ids)

info_clpsd_gd <- info_gd %>%
  group_by(le_id) %>%
  summarise(across(all_of(exp_uniq_cols),
                               unique),
            across(all_of(exp_collapse_cols),
            ~ paste(.x, collapse = ","))
            )

# move chromosome & strand back next to start & end coords
# line up colnames before concatenating 
# TODO: is lining up necessary?
info_clpsd_gd <- relocate(info_clpsd_gd, chromosome, strand, .before = start)
info_bd <- relocate(info_bd, transcript_id,  chromosome, strand, .before = start)

# change start & end to characters 
info_bd <- mutate(info_bd, across(.cols = all_of(c("start", "end")), as.character))

info_clpsd <- bind_rows(info_clpsd_gd, info_bd)

# now join info with the results tbl
res_jned <- left_join(res, info_clpsd, by = "le_id")
# res_jned

# now join with minimal PPAU information (means & deltas for each group)
res_jned <- left_join(res_jned,
                      select(ppau, le_id, starts_with("mean"), starts_with("delta")),
                      by = "le_id")



# Sort so proximal site always comes first for each gene. Also put significant genes at top of table for quick easy browsing
res_jned <- arrange(res_jned, gene.qvalue, gene_id, le_id)

write_tsv(res_jned, file = paste(output_prefix, ".processed.tsv", sep = ""), col_names = T)
