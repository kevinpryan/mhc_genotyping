module_path <- "/hlamajority-paper/external/nf-hlamajoirty/bin/"
options(box.path = module_path)
box::use(lib/df_to_list[...])
box::use(lib/are_vectors_identical[...])
box::use(lib/majority_voting[...])
setwd("/hlamajority-paper/external/mhc_genotyping/")
# read in 1kg gold standard
gold.standard.nci60 <- readRDS("data/gold_standard_nci60.rds")
hlamajority.in <- read.table("../../data/cell-lines/benchmark-cell-lines-all/combined_results/nf_hlamajority_results_majority_vote_combined_sorted.tsv", sep = "\t", header = T)
all.in <- read.table("../../data/cell-lines/benchmark-cell-lines-all/combined_results/nf_hlamajority_hlatyping_results_all_calls_sorted.tsv", sep = "\t", header = T)
all.in.coverage <- read.table("../../data/cell-lines/benchmark-cell-lines-all/combined_results/nf_hlamajority_mean_depth_exons2_3_hla_classI_sorted.tsv", sep = "\t", header = T)
colnames(gold.standard.nci60)[1] <- "sample"
colnames(gold.standard.nci60)[2:7] <- gsub(pattern = "\\.", replacement = "", colnames(gold.standard.nci60[2:7]))
gold.standard.nci60 <- as.data.frame(gold.standard.nci60)
source("scripts/functions/ggroup_mapper.R")
source("scripts/functions/evaluate_predictions_functions.R")

hlamajority_long <- hlamajority.in %>%
  mutate(tool = "hlamajority") %>%
  select(sample, tool, everything()) # Ensure column order matches

# 2. Combine into one Master Dataframe
# This allows us to loop through "unique(master_df$tool)" 
master_df <- bind_rows(all.in, hlamajority_long)

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

master_df_mapped <- map_sample_name_nci60(master_df)                   
run_full_benchmark <- function(master_df, gold_standard, genes, tools) {
  
  # Initialize storage
  summary_stats <- data.frame() # For Accuracy/Call Rate table
  detailed_artifacts <- list()  # For specific error tables/heatmaps
  
  for (gene in genes) {
    print(paste("=== Processing HLA-", gene, " ===", sep=""))
    
    # 1. Prepare Gold Standard for this gene
    #    (Extract cols A1/A2, convert to list, Clean G-groups)
    cols <- c(paste0(gene, "1"), paste0(gene, "2"))
    gs_list_raw <- df_to_list_fix(gold_standard, cols = cols)
    gs_list_clean <- clean_list(gs_list_raw, gene)
    
    for (tool in tools) {
      # 2. Prepare Tool Data
      tool_data <- master_df %>% filter(tool == !!tool)
      
      # Skip if tool didn't call this gene (safety check)
      if (nrow(tool_data) == 0) next
      
      # Convert to list and clean
      tool_list_raw <- df_to_list_fix(tool_data, cols = cols)
      tool_list_clean <- clean_list(tool_list_raw, gene)
      
      # 3. Run Accuracy Comparison (excludes NAs)
      #    Note: ensure run_comparison returns the list we defined previously
      metrics <- run_comparison(gs_list_clean, tool_list_clean, tool_name = tool)
      
      # 4. Run Error Typing (Zygosity checks)
      #    Note: We use the cleaned lists here
      error_types <- analyze_error_types(gs_list_clean, tool_list_clean, tool)
      
      # 5. Run Difficult Allele Analysis
      #    (Extracts the 'incorrect_df' from the metrics object)
      diff_alleles <- analyze_difficult_alleles(gs_list_clean, metrics$gold_standard_vs_tool_incorrect_calls)
      
      # 6. Store Summary Stats
      # Calculate Call Rate
      # common_ids count vs valid_indices count from run_comparison logic
      # (We can reconstruct it here roughly or extract if run_comparison returns it)
      # Simpler approach: 
      n_total <- length(intersect(names(gs_list_clean), names(tool_list_clean)))
      # Count how many in tool list are NOT NA
      n_called <- sum(!sapply(tool_list_clean, function(x) any(is.na(x))))
      call_rate <- n_called / n_total
      
      summary_stats <- rbind(summary_stats, data.frame(
        Gene = gene,
        Tool = tool,
        Accuracy = metrics$final_accuracy,
        Call_Rate = call_rate,
        Num_Samples = n_total,
        Num_Excluded_NA = n_total - n_called
      ))
      
      # 7. Store Detailed Artifacts (nested list)
      # Structure: List$Gene$Tool$Object
      detailed_artifacts[[gene]][[tool]] <- list(
        metrics = metrics,
        error_types = error_types,
        difficult_alleles = diff_alleles
      )
    }
  }
  
  return(list(summary = summary_stats, details = detailed_artifacts))
}

# Run the benchmark
benchmark_results <- run_full_benchmark(
  master_df = master_df_mapped,
  gold_standard = gold.standard.nci60,
  genes = genes_to_analyze,
  tools = tools_to_analyze
)

full_stats <- calculate_overall_stats(benchmark_results$summary)

# 3. Create the Pretty Table
clean_table <- format_publication_table(full_stats)

# Print it nicely
library(knitr)
kable(clean_table, caption = "Accuracy by Gene and Overall (Excluding NAs)")
