module_path <- "/hlamajority-paper/external/nf-hlamajoirty/bin/"
options(box.path = module_path)
box::use(lib/df_to_list[...])
box::use(lib/are_vectors_identical[...])
box::use(lib/majority_voting[...])
setwd("/hlamajority-paper/external/mhc_genotyping/")
# read in 1kg gold standard
gold.standard.1kg <- readRDS("data/gold_standard_1kg.rds")
hlamajority.in <- read.table("../../data/1000-genomes/benchmark-1000genomes-nfhlamajority-local/combined_results/nf_hlamajority_results_majority_vote_combined_sorted.tsv", sep = "\t", header = T)
all.in <- read.table("../../data/1000-genomes/benchmark-1000genomes-nfhlamajority-local/combined_results/nf_hlamajority_hlatyping_results_all_calls_sorted.tsv", sep = "\t", header = T)
all.in.coverage <- read.table("../../data/1000-genomes/benchmark-1000genomes-nfhlamajority-local/combined_results/nf_hlamajority_mean_depth_exons2_3_hla_classI_sorted.tsv", sep = "\t", header = T)

colnames(gold.standard.1kg)[1] <- "sample"
colnames(gold.standard.1kg)[2:7] <- gsub(pattern = "\\.", replacement = "", colnames(gold.standard.1kg[2:7]))
gold.standard.1kg <- as.data.frame(gold.standard.1kg)
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
  master_df = master_df,
  gold_standard = gold.standard.1kg,
  genes = genes_to_analyze,
  tools = tools_to_analyze
)

# --- VIEW RESULTS ---

# 1. The Summary Table (for your App Note Table 1)
print(benchmark_results$summary)

# Example formatting for easier reading:
library(knitr)
benchmark_results$summary %>% 
  mutate(Accuracy = paste0(round(Accuracy * 100, 2), "%"),
         Call_Rate = paste0(round(Call_Rate * 100, 2), "%")) %>%
  arrange(Gene, desc(Accuracy)) %>%
  kable()

benchmark_results_detailed <- benchmark_results$summary %>% 
  mutate(Accuracy = paste0(round(Accuracy * 100, 2), "%"),
         Call_Rate = paste0(round(Call_Rate * 100, 2), "%")) %>%
  arrange(Gene, desc(Accuracy)) 

benchmark_results_detailed %>%
  kable()
write.csv(benchmark_results_detailed, file = "/hlamajority-paper/data/1000-genomes/benchmark-1000genomes-nfhlamajority-local/hlamajority-1000-genomes-benchmark-results-detailed.csv", quote = F, row.names = F)

# Combine error dataframes for all tools for Gene B
errors_B <- bind_rows(lapply(benchmark_results$details$B, function(x) x$error_types))

ggplot(errors_B %>% filter(Type != "Correct"), aes(x = Tool, fill = Type)) +
  geom_bar(position = "stack") +
  theme_minimal() +
  labs(title = "HLA-B Error Profiles", y = "Count of Errors")

errors_A <- bind_rows(lapply(benchmark_results$details$A, function(x) x$error_types))

ggplot(errors_A %>% filter(Type != "Correct"), aes(x = Tool, fill = Type)) +
  geom_bar(position = "stack") +
  theme_minimal() +
  labs(title = "HLA-A Error Profiles", y = "Count of Errors")

errors_C <- bind_rows(lapply(benchmark_results$details$C, function(x) x$error_types))

ggplot(errors_C %>% filter(Type != "Correct"), aes(x = Tool, fill = Type)) +
  geom_bar(position = "stack") +
  theme_minimal() +
  labs(title = "HLA-C Error Profiles", y = "Count of Errors")

kourami_C_diff <- benchmark_results$details$C$kourami$difficult_alleles

# Show top 5 hardest alleles
head(kourami_C_diff, 5)

benchmark_results$summary %>% 
  filter(Gene == "A") %>% 
  arrange(desc(Call_Rate))

polysolver_B_diff <- benchmark_results$details$B$polysolver$difficult_alleles

benchmark_results$details$A$polysolver$metrics$proportion_excluded

# Show top 5 hardest alleles
head(polysolver_B_diff, 5)
polysolver_B_diff
full_stats <- calculate_overall_stats(benchmark_results$summary)

# 3. Create the Pretty Table
clean_table <- format_publication_table(full_stats)

# Print it nicely
library(knitr)
kable(clean_table, caption = "Accuracy by Gene and Overall (Excluding NAs)")
write.csv(clean_table, file = "/hlamajority-paper/data/1000-genomes/benchmark-1000genomes-nfhlamajority-local/hlamajority-1000-genomes-benchmark-results-summary.csv", quote = F, row.names = F)
analyze_difficult_alleles(full_gold_standard = benchmark_results$details$B$polysolver$metrics$gold_standard_list,incorrect_df = benchmark_results$details$B$polysolver$metrics$gold_standard_vs_tool_incorrect_calls
                          )

optitype_A_errors <- benchmark_results$details$A$optitype$metrics$gold_standard_vs_tool_incorrect_calls

analyze_difficult_alleles(full_gold_standard = benchmark_results$details$B$hlamajority$metrics$gold_standard_list,incorrect_df = benchmark_results$details$B$hlamajority$metrics$gold_standard_vs_tool_incorrect_calls)
plot_confusion_heatmap(benchmark_results$details$B$polysolver$metrics$gold_standard_vs_tool_incorrect_calls, min_count = 3)
