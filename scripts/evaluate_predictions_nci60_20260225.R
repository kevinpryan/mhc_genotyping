module_path <- "/hlamajority-paper/external/nf-hlamajority/bin/"
options(box.path = module_path)
box::use(lib/df_to_list[...])
box::use(lib/are_vectors_identical[...])
box::use(lib/majority_voting[...])
setwd("/hlamajority-paper/external/mhc_genotyping/")
gold.standard.nci60 <- readRDS("data/gold_standard_nci60.rds")
hlamajority.in <- read.table("../../data/cell-lines/benchmark-cell-lines-all-kourami-3-63-0-majority-vote/combined_results/nf_hlamajority_votes_combined_sorted.tsv", sep = "\t", header = T)
all.in <- read.table("../../data/cell-lines/benchmark-cell-lines-all-kourami-3-63-0-majority-vote/combined_results/nf_hlamajority_all_calls_sorted.tsv", sep = "\t", header = T)
#all.in.coverage <- read.table("../../data/cell-lines/benchmark-cell-lines-all/combined_results/nf_hlamajority_mean_depth_exons2_3_hla_classI_sorted.tsv", sep = "\t", header = T)
colnames(gold.standard.nci60)[1] <- "sample"
colnames(gold.standard.nci60)[2:7] <- gsub(pattern = "\\.", replacement = "", colnames(gold.standard.nci60[2:7]))
gold.standard.nci60 <- as.data.frame(gold.standard.nci60)
#na_samples_gs <- gold.standard.nci60 %>% 
#  dplyr::filter(if_any(c(A1, A2, B1, B2, C1, C2), is.na))
#na_samples_gs

source("scripts/functions/ggroup_mapper.R")
source("scripts/functions/evaluate_predictions_functions.R")

hlamajority_long <- hlamajority.in %>%
  mutate(tool = "hlamajority") %>%
  select(sample, tool, everything()) # Ensure column order matches

cell.line.in.wide <- hlamajority.in %>%
  pivot_wider(
    id_cols = sample,             # keep 'sample' as the identifier
    names_from = gene,            # use 'gene' values to create new column groups
    values_from = c(allele1, allele2), # these are the values to spread
    names_glue = "{gene}{.value}" # creates columns like Aallele1, Aallele2
  ) %>%
  dplyr::rename(
    A1 = `HLA-Aallele1`,
    A2 = `HLA-Aallele2`,
    B1 = `HLA-Ballele1`,
    B2 = `HLA-Ballele2`,
    C1 = `HLA-Callele1`,
    C2 = `HLA-Callele2`
  ) %>%
  mutate(tool = "hlamajority")
  

cell.line.in.wide

# 2. Combine into one Master Dataframe
# This allows us to loop through "unique(master_df$tool)" 
master_df <- bind_rows(all.in, cell.line.in.wide)

# 3. List of tools to analyze
tools_to_analyze <- unique(master_df$tool)
genes_to_analyze <- c("A", "B", "C")

# Map SRR to cell line names
sample_names_nci60 <- readRDS("data/sample_names_nci60.rds")
map_sample_name_nci60 <- function(res_df) {
  res_df <- res_df %>%
    left_join(sample_names_nci60, by = c("sample" = "Experiment")) %>%
    mutate(sample_id = coalesce(sample_id.new, sample)) %>%
    select(-sample_id.new, -sample_id.srr) %>% 
    dplyr::mutate(sample = sample_id.gs) %>%
    dplyr::select(-sample_id, -sample_id.old, -sample_id.gs)
}

map_sample_name_nci60_full <- function(res_df) {
  res_df <- res_df %>%
    left_join(sample_names_nci60, by = c("sample" = "Experiment")) %>%
    mutate(sample_id = coalesce(sample_id.new, sample)) %>%
    select(-sample_id.new, -sample_id.srr) %>% 
    dplyr::mutate(sample = sample_id.gs) %>%
    dplyr::select(-sample_id.old, -sample_id.gs)
}

master_df_mapped <- map_sample_name_nci60(master_df)                   
master_df_mapped_full <- map_sample_name_nci60_full(master_df)                   
saveRDS(master_df_mapped_full, file = "data/results/hlamajority/nci-map.Rds")
# run_full_benchmark <- function(master_df, gold_standard, genes, tools) {
#   
#   # Initialize storage
#   summary_stats <- data.frame() # For Accuracy/Call Rate table
#   detailed_artifacts <- list()  # For specific error tables/heatmaps
#   
#   for (gene in genes) {
#     print(paste("=== Processing HLA-", gene, " ===", sep=""))
#     
#     # 1. Prepare Gold Standard for this gene
#     #    (Extract cols A1/A2, convert to list, Clean G-groups)
#     cols <- c(paste0(gene, "1"), paste0(gene, "2"))
#     gs_list_raw <- df_to_list_fix(gold_standard, cols = cols)
#     gs_list_clean <- clean_list(gs_list_raw, gene)
#     # Gold standard callable samples (non-NA alleles)
#     gs_callable <- sapply(gs_list_clean, function(x) !any(is.na(x)))
#     gs_callable_ids <- names(gs_callable)[gs_callable]
#     
#     n_gs_called <- length(gs_callable_ids)
#     for (tool in tools) {
#       # 2. Prepare Tool Data
#       tool_data <- master_df %>% filter(tool == !!tool)
#       
#       # Skip if tool didn't call this gene (safety check)
#       if (nrow(tool_data) == 0) next
#       
#       # Convert to list and clean
#       tool_list_raw <- df_to_list_fix(tool_data, cols = cols)
#       tool_list_clean <- clean_list(tool_list_raw, gene)
#       
#       # 3. Run Accuracy Comparison (excludes NAs)
#       #    Note: ensure run_comparison returns the list we defined previously
#       metrics <- run_comparison(gs_list_clean, tool_list_clean, tool_name = tool)
#       
#       # 4. Run Error Typing (Zygosity checks)
#       #    Note: We use the cleaned lists here
#       error_types <- analyze_error_types(gs_list_clean, tool_list_clean, tool)
#       
#       # 5. Run Difficult Allele Analysis
#       #    (Extracts the 'incorrect_df' from the metrics object)
#       diff_alleles <- analyze_difficult_alleles(gs_list_clean, metrics$gold_standard_vs_tool_incorrect_calls)
#       
#       # 6. Store Summary Stats
#       # Calculate Call Rate
#       # common_ids count vs valid_indices count from run_comparison logic
#       # (We can reconstruct it here roughly or extract if run_comparison returns it)
#       # Simpler approach: 
#       # n_total <- length(intersect(names(gs_list_clean), names(tool_list_clean)))
#       # # Count how many in tool list are NOT NA
#       # n_called <- sum(!sapply(tool_list_clean, function(x) any(is.na(x))))
#       # call_rate <- n_called / n_total
#       # Restrict to samples where GS made a call
#       tool_sub <- tool_list_clean[gs_callable_ids]
#       
#       # Tool callable among GS-callable samples
#       tool_callable <- sapply(tool_sub, function(x) !any(is.na(x)))
#       
#       n_total <- length(gs_callable_ids)        # gold-standard calls
#       n_called <- sum(tool_callable)             # tool calls where GS is callable
#       n_excluded_na <- n_total - n_called        # guaranteed ≥ 0
#       
#       call_rate <- n_called / n_total
#       
#       summary_stats <- rbind(summary_stats, data.frame(
#         Gene = gene,
#         Tool = tool,
#         Accuracy = metrics$final_accuracy,
#         Call_Rate = call_rate,
#         Num_Samples = n_total,              # GS callable
#         Num_GS_Called = n_gs_called,         # same as Num_Samples (explicit)
#         Num_Tool_Called = n_called,
#         Num_Excluded_NA = n_excluded_na
#       ))
#       
#       # 7. Store Detailed Artifacts (nested list)
#       # Structure: List$Gene$Tool$Object
#       detailed_artifacts[[gene]][[tool]] <- list(
#         metrics = metrics,
#         error_types = error_types,
#         difficult_alleles = diff_alleles
#       )
#     }
#   }
#   
#   return(list(summary = summary_stats, details = detailed_artifacts))
# }

# Run the benchmark
benchmark_results <- run_full_benchmark(
  master_df = master_df_mapped,
  gold_standard = gold.standard.nci60,
  genes = genes_to_analyze,
  tools = tools_to_analyze
)
#benchmark_results$gold_standard_missing %>% group_by(Gene) %>% summarise(n = n())
saveRDS(benchmark_results, file = "data/results/hlamajority/nci-full-results-hlamajority-majority-vote.Rds")
full_stats <- calculate_overall_stats(benchmark_results$summary)
write.csv(full_stats, file = "data/results/hlamajority/nci-full-stats-hlamajority-majority-vote.csv", row.names = F, quote = F)
# 3. Create the Pretty Table
clean_table <- format_publication_table(full_stats)

# Print it nicely
library(knitr)
kable(clean_table, caption = "Accuracy by Gene and Overall (Excluding NAs)")

all_results_hlamajority <- vroom("../../data/cell-lines/benchmark-cell-lines-all-kourami-3-63-0-majority-vote/combined_results/nf_hlamajority_votes_combined_sorted.tsv")
all_results_hlamajority$support
all_results_hlamajority
benchmark_results
sample_names_nci60
#left_join(all_results_hlamajority, sample_names_nci60)
all_results_hlamajority_cell_line_id <- map_sample_name_nci60(all_results_hlamajority)
# extract_scores <- function(results, gene) {
#   x <- results$details[[gene]]$hlamajority$metrics$scores_per_person
# 
#   #tab <- table(x)
#   
#   data.frame(
#     gene  = paste("HLA-",gene, sep = ""),
#     sample  = names(x),
#     Score = as.integer(x),
#     row.names = NULL
#   )
# }
# 
# genes <- c("A", "B", "C")
# final_df <- do.call(
#   rbind,
#   lapply(genes, function(g) {
#     extract_scores(benchmark_results, g)
#       })
# )
# final_df
# head(master_df_mapped_full)
# master_df_mapped_full_sample_sampleid <- master_df_mapped_full %>% dplyr::select(sample, sample_id) %>% distinct(sample, sample_id)
# score_support <- left_join(final_df, all_results_hlamajority_cell_line_id, by = c("sample", "gene"))
# score_support_sampleid <- left_join(score_support, master_df_mapped_full_sample_sampleid)

extract_scores_tool <- function(results, gene, tool) {
  x <- results$details[[gene]][[tool]]$metrics$scores_per_person
  
  #tab <- table(x)
  
  data.frame(
    gene  = paste("HLA-",gene, sep = ""),
    sample  = names(x),
    Score = as.integer(x),
    tool = tool,
    row.names = NULL
  )
}

extract_scores_tool(benchmark_results, "A", "kourami")

tools <-  c("kourami", "hlala", "polysolver", "optitype", "hlamajority")
genes <- c("A", "B", "C")

final_df_score <- do.call(
  rbind,
  lapply(genes, function(g) {
    do.call(
      rbind,
      lapply(tools, function(t) {
        extract_scores_tool(benchmark_results, g, t)
      })
    )
  })
)

# tools <-  c("kourami", "hlala", "polysolver", "optitype", "hlamajority")
# 
# final_df <- do.call(
#   rbind,
#   lapply(genes, function(g) {
#     do.call(
#       rbind,
#       lapply(tools, function(t) {
#         summarise_gene_tool(results, g, t)
#       })
#     )
#   })
# )

#master_df_mapped_full_sample_sampleid <- master_df_mapped_full %>% dplyr::select(sample, sample_id) %>% distinct(sample, sample_id)
all_results_hlamajority_cell_line_id_for_join <- all_results_hlamajority_cell_line_id %>% dplyr::select(-allele1, -allele2, -matching_tools)
score_support <- left_join(final_df_score, all_results_hlamajority_cell_line_id_for_join, by = c("sample", "gene")) 
# score_support_sampleid <- left_join(score_support, master_df_mapped_full_sample_sampleid)
write.csv(score_support, file = "data/results/hlamajority/nci60-score-depth-per-sample-per-tool.csv", row.names = F, quote= F)

cor(score_support$support, score_support$Score, method = "spearman")

plot(score_support$support, score_support$Score)
