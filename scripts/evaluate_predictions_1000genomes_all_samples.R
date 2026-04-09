library(dplyr)
library(knitr)
module_path <- "/hlamajority-paper/external/nf-hlamajority/bin/"
options(box.path = module_path)
box::use(lib/df_to_list[...])
box::use(lib/are_vectors_identical[...])
box::use(lib/majority_voting[...])
setwd("/hlamajority-paper/external/mhc_genotyping/")
source("scripts/functions/ggroup_mapper.R")
source("scripts/functions/evaluate_predictions_functions.R")
gold.standard.1kg <- readRDS("data/gold_standard_1kg.rds")
#hlamajority.in <- read.table("../../data/raw/1000-genomes/benchmark-1000genomes-nfhlamajority-local-update-db-exclude-trim-majority-all-samples/benchmark-1000genomes-nfhlamajority-all-20260309-majority-handle-error-kourami-hlala/combined_results/nf_hlamajority_votes_combined_sorted.tsv", sep = "\t", header = T)
hlamajority.in <- read.table("../../data/raw/1000-genomes/majority/all_samples/combined_results/nf_hlamajority_votes_combined_sorted.tsv", sep = "\t", header = T)
colnames(gold.standard.1kg)[1] <- "sample"
colnames(gold.standard.1kg)[2:7] <- gsub(pattern = "\\.", replacement = "", colnames(gold.standard.1kg[2:7]))
gold.standard.1kg <- as.data.frame(gold.standard.1kg)
na_samples_gs <- gold.standard.1kg %>% 
  dplyr::filter(if_any(c(B1, B2, C1, C2), is.na))
#outdir <- c("../../data/processed/results/hlamajority/1000genomes-all-samples/")
outdir <- c("../../data/processed/1000-genomes/majority/")
dir.create(file.path(outdir), showWarnings = FALSE, recursive = TRUE)
#dir.create(file.path("../../data/processed/", "results/hlamajority/1000genomes-all-samples/"), showWarnings = FALSE, recursive = TRUE)
write.csv(na_samples_gs, file = paste(outdir, "1000-genomes-gs-na-samples.csv", sep = ""), row.names = F, quote = F)
#all.in <- read.table("../../data/1000-genomes/benchmark-1000genomes-nfhlamajority-local-update-db-exclude-trim-majority-all-samples/benchmark-1000genomes-nfhlamajority-all-20260309-majority-handle-error-kourami-hlala/combined_results/nf_hlamajority_all_calls_sorted.tsv", sep = "\t", header = T)
#all.in.coverage <- read.table("../../data/1000-genomes/benchmark-1000genomes-nfhlamajority-local-update-db-exclude-trim-majority-all-samples/benchmark-1000genomes-nfhlamajority-all-20260309-majority-handle-error-kourami-hlala/combined_results/nf_hlamajority_depth_sorted.tsv", sep = "\t", header = T)
all.in <- read.table("../../data/raw/1000-genomes/majority/all_samples/combined_results/nf_hlamajority_all_calls_sorted.tsv", sep = "\t", header = T)
all.in.coverage <- read.table("../../data/raw/1000-genomes/majority/all_samples/combined_results/nf_hlamajority_depth_sorted.tsv", sep = "\t", header = T)
hlamajority_long <- hlamajority.in %>%
  mutate(tool = "hlamajority") %>%
  dplyr::select(sample, tool, everything()) # Ensure column order matches

hlamajority.in.wide <- hlamajority.in %>%
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

master_df <- bind_rows(all.in, hlamajority.in.wide)
tools_to_analyze <- unique(master_df$tool)
genes_to_analyze <- c("A", "B", "C")

benchmark_samples <- master_df$sample
gold.standard.1kg.samples <- gold.standard.1kg$sample
all.samples <- intersect(benchmark_samples, gold.standard.1kg.samples)
master_df <- master_df %>% 
                      dplyr::filter(sample %in% all.samples)
gold.standard.1kg <- gold.standard.1kg %>% 
                             dplyr::filter(sample %in% all.samples)
source("scripts/functions/evaluate_predictions_functions.R")
# Run the benchmark
benchmark_results <- run_full_benchmark(
  master_df = master_df,
  gold_standard = gold.standard.1kg,
  genes = genes_to_analyze,
  tools = tools_to_analyze
)

full_stats <- calculate_overall_stats(benchmark_results$summary)
write.csv(full_stats, file = paste(outdir, "1000-genomes-full-stats-hlamajority-majority-vote.csv", sep = ""), row.names = F, quote = F)
saveRDS(benchmark_results, file = paste(outdir, "1000-genomes-full-results-hlamajority-majority-vote.Rds", sep = ""))
#benchmark_results$details$A$kourami
# 3. Create the Pretty Table
clean_table <- format_publication_table(full_stats)

# Print it nicely
kable(clean_table, caption = "Accuracy by Gene and Overall (Excluding NAs)")


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

#extract_scores_tool(benchmark_results, "A", "kourami")

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


write.csv(final_df_score, file = paste(outdir, "1000genomes-score-per-sample.csv", sep = ""), row.names = F, quote= F)
#final_df_score
