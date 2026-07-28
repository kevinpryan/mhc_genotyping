library(vroom)
library(knitr)
module_path <- "/hlamajority-paper/external/nf-hlamajority/bin/"
options(box.path = module_path)
box::use(lib/df_to_list[...])
box::use(lib/are_vectors_identical[...])
box::use(lib/majority_voting[...])
setwd("/hlamajority-paper/external/mhc_genotyping/")
gold.standard.nci60 <- readRDS("data/gold_standard_nci60.rds")
#hlamajority.in <- read.table("../../data/cell-lines/benchmark-cell-lines-all-kourami-3-63-0-majority-vote/combined_results/nf_hlamajority_votes_combined_sorted.tsv", sep = "\t", header = T)
#hlamajority.in <- read.table("../../data/raw/cell-lines/benchmark-cell-lines-all-kourami-3-63-0-majority-vote/combined_results/nf_hlamajority_votes_combined_sorted.tsv", sep = "\t", header = T)
hlamajority.in <- read.table("../../data/raw/cell-lines-before-polysolver-change/majority/combined_results/nf_hlamajority_votes_combined_sorted.tsv", sep = "\t", header = T)

#all.in <- read.table("../../data/raw/cell-lines/benchmark-cell-lines-all-kourami-3-63-0-majority-vote/combined_results/nf_hlamajority_all_calls_sorted.tsv", sep = "\t", header = T)
all.in <- read.table("../../data/raw/cell-lines-before-polysolver-change/majority/combined_results/nf_hlamajority_all_calls_sorted.tsv", sep = "\t", header = T)
#all.in.coverage <- read.table("../../data/cell-lines/benchmark-cell-lines-all/combined_results/nf_hlamajority_mean_depth_exons2_3_hla_classI_sorted.tsv", sep = "\t", header = T)
colnames(gold.standard.nci60)[1] <- "sample"
colnames(gold.standard.nci60)[2:7] <- gsub(pattern = "\\.", replacement = "", colnames(gold.standard.nci60[2:7]))
gold.standard.nci60 <- as.data.frame(gold.standard.nci60)
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
master_df %>% dplyr::filter(tool == "polysolver")
master_df %>% dplyr::filter(tool == "nf-hlamajority" & sample == "SRX4239529")
gold.standard.nci60$sample
# Map SRR to cell line names
#sample_names_nci60 <- readRDS("data/sample_names_nci60.rds")
sample_names_nci60 <- readRDS("data/sample_names_nci60_srx.rds") %>% rename(sample = "sample_id.srx")
sample_names_nci60 %>% dplyr::filter(sample == "SRX4239529")
sample_names_nci60 %>% dplyr::filter(sample == "SRX4239554")

map_sample_name_nci60 <- function(res_df) {
  res_df <- res_df %>%
    #left_join(sample_names_nci60, by = c("sample" = "Experiment")) %>%
    left_join(sample_names_nci60, by = "sample") %>%
    mutate(sample_id = coalesce(sample_id.new, sample)) %>%
    dplyr::select(-sample_id.new
                  #, -sample_id.srr
                  ) %>% 
    dplyr::mutate(sample = sample_id) 
    #dplyr::select(-sample_id, -sample_id.old, -sample_id.gs)
}

# map_sample_name_nci60_full <- function(res_df) {
#   res_df <- res_df %>%
#     #left_join(sample_names_nci60, by = c("sample" = "Experiment")) %>%
#     left_join(sample_names_nci60, by = "sample") %>%
#     mutate(sample_id = coalesce(sample_id.new, sample)) %>%
#     dplyr::mutate(sample = sample_id)
#     #select(-sample_id.new, -sample_id.srr) %>% 
#     #dplyr::select(-sample_id.old, -sample_id.gs)
# }
map_sample_name_nci60_full <- function(res_df) {
  res_df <- res_df %>%
    #left_join(sample_names_nci60, by = c("sample" = "Experiment")) %>%
    left_join(sample_names_nci60, by = "sample") %>%
    mutate(sample_id = coalesce(sample_id.new, sample),
           sample_id.srx = sample) %>%
    dplyr::mutate(sample = sample_id)
  #select(-sample_id.new, -sample_id.srr) %>% 
  #dplyr::select(-sample_id.old, -sample_id.gs)
}

print("master_df...")
head(master_df)
master_df_mapped <- map_sample_name_nci60(master_df)                   
master_df_mapped_full <- map_sample_name_nci60_full(master_df)                   
#saveRDS(master_df_mapped_full, file = "data/results/hlamajority/nci-map.Rds")
#outdir <- c("../../data/processed/results/hlamajority/cell-lines/")
outdir <- c("../../data/processed/cell-lines-before-polysolver-change/majority/")

dir.create(file.path(outdir), showWarnings = FALSE, recursive = TRUE)
saveRDS(master_df_mapped_full, file = paste(outdir, "nci-map.Rds", sep = ""))
#"../../data/processed/results/hlamajority/cell-lines/nci-map.Rds")

# Run the benchmark
benchmark_results <- run_full_benchmark(
  master_df = master_df_mapped,
  gold_standard = gold.standard.nci60,
  genes = genes_to_analyze,
  tools = tools_to_analyze
)
# optitype_calls_before_polysolver_change <- master_df_mapped %>% dplyr::filter(tool == "optitype")
# write.csv(optitype_calls_before_polysolver_change, file = "../optitype-calls-before-polysolver-change.csv", quote = F, row.names = F)
#benchmark_results$gold_standard_missing %>% group_by(Gene) %>% summarise(n = n())
#saveRDS(benchmark_results, file = "data/results/hlamajority/nci-full-results-hlamajority-majority-vote.Rds")
#saveRDS(benchmark_results, file = "../../data/processed/results/hlamajority/cell-lines/nci-full-results-hlamajority-majority-vote.Rds")
saveRDS(benchmark_results, file = paste(outdir, "nci-full-results-hlamajority-majority-vote.Rds", sep = ""))

full_stats <- calculate_overall_stats(benchmark_results$summary)
#write.csv(full_stats, file = "data/results/hlamajority/nci-full-stats-hlamajority-majority-vote.csv", row.names = F, quote = F)
#write.csv(full_stats, file = "../../data/processed/results/hlamajority/cell-lines/nci-full-stats-hlamajority-majority-vote.csv", row.names = F, quote = F)
write.csv(full_stats, file = paste(outdir, "nci-full-stats-hlamajority-majority-vote.csv", sep = ""), row.names = F, quote = F)
# 3. Create the Pretty Table
clean_table <- format_publication_table(full_stats)

# Print it nicely
kable(clean_table, caption = "Accuracy by Gene and Overall (Excluding NAs)")

#all_results_hlamajority <- vroom("../../data/cell-lines/benchmark-cell-lines-all-kourami-3-63-0-majority-vote/combined_results/nf_hlamajority_votes_combined_sorted.tsv")
all_results_hlamajority <- vroom("../../data/raw/cell-lines-before-polysolver-change/majority/combined_results/nf_hlamajority_votes_combined_sorted.tsv")

all_results_hlamajority_cell_line_id <- map_sample_name_nci60(all_results_hlamajority)

extract_scores_tool <- function(results, gene, tool) {
  x <- results$details[[gene]][[tool]]$metrics$scores_per_person
  
  data.frame(
    gene  = paste("HLA-",gene, sep = ""),
    sample  = names(x),
    Score = as.integer(x),
    tool = tool,
    row.names = NULL
  )
}

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

all_results_hlamajority_cell_line_id_for_join <- all_results_hlamajority_cell_line_id %>% dplyr::select(-allele1, -allele2, -matching_tools)
score_support <- left_join(final_df_score, all_results_hlamajority_cell_line_id_for_join, by = c("sample", "gene")) 
#write.csv(score_support, file = "data/results/hlamajority/nci60-score-depth-per-sample-per-tool.csv", row.names = F, quote= F)
#write.csv(score_support, file = "../../data/processed/results/hlamajority/cell-lines/nci60-score-depth-per-sample-per-tool.csv", row.names = F, quote= F)
write.csv(score_support, file = paste(outdir, "nci60-score-depth-per-sample-per-tool.csv", sep = ""), row.names = F, quote= F) #"../../data/processed/results/hlamajority/cell-lines/nci60-score-depth-per-sample-per-tool.csv", row.names = F, quote= F)

#cor(score_support$support, score_support$Score, method = "spearman")
#plot(score_support$support, score_support$Score)
