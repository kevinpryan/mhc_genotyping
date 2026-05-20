library(vroom)
library(knitr)
file.exists("/hlamajority-paper/external/mhc_genotyping/scripts/functions/Conversion/arcasHLA_conversion_LENS.R")
module_path <- "/hlamajority-paper/external/mhc_genotyping/scripts/functions"
# options(box.path = module_path)
# box::use(Conversion/arcasHLA_conversion_LENS[
#   toolOutputToR.arcasHLA_conversion_LENS
# ])
setwd("/hlamajority-paper/external/mhc_genotyping/")
source("scripts/functions/ggroup_mapper.R")
source("scripts/functions/evaluate_predictions_functions.R")
gold.standard.nci60 <- readRDS("data/gold_standard_nci60.rds")
colnames(gold.standard.nci60)[1] <- "sample"
colnames(gold.standard.nci60)[2:7] <- gsub(pattern = "\\.", replacement = "", colnames(gold.standard.nci60[2:7]))
gold.standard.nci60 <- as.data.frame(gold.standard.nci60)

sample_names_nci60 <- readRDS("data/sample_names_nci60_srx.rds") %>% rename(sample = "sample_id.srx")

# read in results
lens_results <- read.csv("../../data/processed/cell-lines-after-polysolver-change/lens-majority-benchmark/combined-results/lens-v1.2-dev-v1.8-results-standardised.csv") %>%
                rename("sample_name" = "sample") %>% 
                mutate(sample = str_split_fixed(srx_srr, pattern = "_", n = 2)[,1])

sample_names_nci60 <- readRDS("data/sample_names_nci60_srx.rds") %>% rename(sample = "sample_id.srx")

map_sample_name_nci60_full <- function(res_df) {
  res_df <- res_df %>%
    left_join(sample_names_nci60, by = "sample") %>%
    mutate(sample_id = coalesce(sample_id.new, sample)) %>%
    dplyr::mutate(sample = sample_id)
}
master_df_mapped_full_lens <- map_sample_name_nci60_full(lens_results)        
master_df_mapped_full_lens <- master_df_mapped_full_lens %>% dplyr::select(c("sample", "tool", "A1", "A2", "B1", "B2", "C1", "C2"))
gold.standard.nci60 <- gold.standard.nci60 %>% dplyr::filter(sample != "UO-31")

# benchmark_results_lens <- run_full_benchmark(
#   master_df = master_df_mapped_full_lens,
#   gold_standard = gold.standard.nci60,
#   genes = genes_to_analyze,
#   tools = tools_to_analyze
# )
# for HLAmajority/LENS comparison, exclude UO-31 (single-end sample) - rerun evaluation without this sample
hlamajority.in <- read.table("../../data/raw/cell-lines-after-polysolver-change/majority/combined_results/nf_hlamajority_votes_combined_sorted.tsv", sep = "\t", header = T)
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
all.in <- read.table("../../data/raw/cell-lines-after-polysolver-change/majority/combined_results/nf_hlamajority_all_calls_sorted.tsv", sep = "\t", header = T)
master_df <- bind_rows(all.in, cell.line.in.wide)
genes_to_analyze <- c("A", "B", "C")
sample_names_nci60 <- readRDS("data/sample_names_nci60_srx.rds") %>% rename(sample = "sample_id.srx")
map_sample_name_nci60 <- function(res_df) {
  res_df <- res_df %>%
    #left_join(sample_names_nci60, by = c("sample" = "Experiment")) %>%
    left_join(sample_names_nci60, by = "sample") %>%
    mutate(sample_id = coalesce(sample_id.new, sample)) %>%
    dplyr::select(-sample_id.new) %>% 
    dplyr::mutate(sample = sample_id) 
  #dplyr::select(-sample_id, -sample_id.old, -sample_id.gs)
}
master_df_mapped_hlamajority <- map_sample_name_nci60(master_df)  
master_df_mapped_hlamajority <- master_df_mapped_hlamajority %>% dplyr::select(c("sample", "tool", "A1", "A2", "B1", "B2", "C1", "C2"))

all_data_combined <- rbind(master_df_mapped_full_lens, master_df_mapped_hlamajority)
all_data_combined <- all_data_combined %>% dplyr::filter(sample != "UO-31")
tools_to_analyze <- unique(all_data_combined$tool)
genes_to_analyze <- c("A", "B", "C")
# Run the benchmark
benchmark_results <- run_full_benchmark_concordance(
  master_df = all_data_combined,
  gold_standard = gold.standard.nci60,
  genes = genes_to_analyze,
  tools = tools_to_analyze
)
full_stats <- calculate_overall_stats(benchmark_results$summary)
df <- benchmark_results$tool_gene_discordance
outdir <- c("../../data/processed/cell-lines-after-polysolver-change/lens-majority-benchmark/")
saveRDS(benchmark_results, file = paste(outdir, "nci-full-results-lens-compare-hlamajority.Rds", sep = ""))
write.csv(full_stats, file = paste(outdir, "nci-full-stats-lens-compare-hlamajority.csv", sep = ""), row.names = F, quote = F)
