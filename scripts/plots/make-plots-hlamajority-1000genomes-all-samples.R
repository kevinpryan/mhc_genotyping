library(ggplot2)
library(ggthemes)
library(dplyr)
library(vroom)
library(forcats)
if (!require("Biostrings")) install.packages("BiocManager"); BiocManager::install("Biostrings")
BiocManager::install("pwalign")
install.packages("svglite")
library(Biostrings)
library(pwalign)
#library(ordinal)
library(tidyr)
setwd("/hlamajority-paper/external/mhc_genotyping/")
source("scripts/functions/evaluate_predictions_functions.R")
depth <- vroom("../../data/1000-genomes/benchmark-1000genomes-nfhlamajority-local-update-db-exclude-trim-majority-all-samples/benchmark-1000genomes-nfhlamajority-all-20260309-majority-handle-error-kourami-hlala/combined_results/nf_hlamajority_depth_sorted.tsv")
#depth <- depth %>% distinct(mean_depth_hla_exons_2_3_classI, .keep_all = T) 
all.in <- read.table("../../data/1000-genomes/benchmark-1000genomes-nfhlamajority-local-update-db-exclude-trim-majority-all-samples/benchmark-1000genomes-nfhlamajority-all-20260309-majority-handle-error-kourami-hlala/combined_results/nf_hlamajority_all_calls_sorted.tsv", sep = "\t", header = T)
results <- readRDS("data/results/hlamajority/1000genomes-all-samples/1000-genomes-full-results-hlamajority-majority-vote.Rds")
results$summary$Tool <- factor(
  results$summary$Tool,
  levels = c("kourami", "hlala", "polysolver", "optitype", "hlamajority"),
  labels = c("Kourami", "HLA*LA", "Polysolver", "Optitype", "nf-hlamajority")
)
df <- read.csv("data/results/hlamajority/1000genomes-all-samples/1000-genomes-full-stats-hlamajority-majority-vote.csv")

df <- df %>% dplyr::mutate(Accuracy = 100*Accuracy)
df$Gene <- factor(df$Gene, levels = c("Overall", "A", "B", "C"))
# df$Tool <- factor(df$Tool, 
#                   levels = c("kourami", "hlala", "polysolver", "optitype", "hlamajority"),
#                   labels = c("Kourami", "HLA-LA", "Polysolver", "Optitype", "nf-hlamajority")
# )
df$Tool <- factor(
  df$Tool,
  levels = c("kourami", "hlala", "polysolver", "optitype", "hlamajority"),
  labels = c("Kourami", "HLA*LA", "Polysolver", "Optitype", "nf-hlamajority")
)

my_gene_labels <- c(
  "Overall" = "Overall Accuracy",
  "A"     = "HLA-A",
  "B"     = "HLA-B",
  "C"     = "HLA-C"
) 

# palette_mistake_types <- c(
#                            Correct = "#BDD5EA",
#                            Dropout = "#577399",
#                            Hallucination = "#495867",
#                            "Partial Mismatch" = "#F7F7FF",
#                            "Complete Mismatch" = "#FE5F55"
#                            )

palette_mistake_types <- c(
  Correct = "#016FB9",
  Dropout = "#22AED1",
  Hallucination = "#6D8EA0",
  "Partial Mismatch" = "#AFA98D",
  "Complete Mismatch" = "#182825"
)

p <- ggplot(df, aes(x = Tool, y = Accuracy, fill = Tool)) +
  
  # Create bars
  geom_col(position = position_dodge(), width = 0.7, color = "black", size = 0.2) +
  
  # Facet by Gene
  facet_wrap(~Gene, scales = "fixed", 
             ncol = 2,
             labeller = as_labeller(my_gene_labels)
  ) +
  
  # Add text labels on top of bars (rounded to 1 decimal)
  geom_text(#aes(label = sprintf("%.1f", accuracy)), 
    aes(label = paste(round(Accuracy, 1), "%", sep = "")),
    position = position_dodge(width = 0.9), 
    vjust = -0.5, 
    size = 6) +
  
  # Colors: Highlight Hlamajority (assuming it's the last factor level)
  # You can customize these colors. 
  # Here: Greys for others, Red/Blue for Hlamajority
  scale_fill_manual(values = c("#999999", "#999999", "#999999", "#999999", "#E69F00")) +
  
  # Scales
  scale_y_continuous(limits = c(0, 120), breaks = seq(0, 100, 25), expand = c(0,0)) +
  
  # Labels
  labs(
    #title = "Genotyping Accuracy by Tool and Gene",
    #subtitle = "Comparison with NCI-60 WES Dataset",
    y = "Accuracy (%)",
    x = "Tool",
    fill = "Tool") +
  
  # Theme customization
  theme_bw() +
  coord_cartesian(ylim = c(0, 115)) +     # Set the visible limits here
  theme(
    #axis.text.x = element_text(angle = 45, hjust = 1), # Rotate x labels
    strip.background = element_rect(fill = "#f0f0f0"), # Facet header background
    strip.text = element_text(face = "bold", size = 18),
    panel.grid.major.x = element_blank(),
    legend.position = "none", # Hide legend since x-axis has labels
    axis.title = element_text(size = 18), 
    axis.text = element_text(size = 16), 
    axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"),
    axis.text.y = element_text(colour = "black")
  )

# 4. Display Plot
print(p)
ggsave(plot = p, filename = "data/results/hlamajority/plots/results-all-samples-1000genomes/hlamajority-accuracy-per-gene-1000genomes-all-samples.png", width = 10, height = 10)
ggsave(plot = p, filename = "data/results/hlamajority/plots/results-all-samples-1000genomes/hlamajority-accuracy-per-gene-1000genomes-all-samples.svg", width = 10, height = 10)

p_1row <- ggplot(df, aes(x = Tool, y = Accuracy, fill = Tool)) +
  
  # Create bars
  geom_col(position = position_dodge(), width = 0.7, color = "black", size = 0.2) +
  
  # Facet by Gene
  facet_wrap(~Gene, scales = "fixed", 
             ncol = 4,
             labeller = as_labeller(my_gene_labels)
  ) +
  
  # Add text labels on top of bars (rounded to 1 decimal)
  geom_text(#aes(label = sprintf("%.1f", accuracy)), 
    aes(label = paste(round(Accuracy, 1), "%", sep = "")),
    position = position_dodge(width = 0.9), 
    vjust = -0.5, 
    size = 6) +
  
  # Colors: Highlight Hlamajority (assuming it's the last factor level)
  # You can customize these colors. 
  # Here: Greys for others, Red/Blue for Hlamajority
  scale_fill_manual(values = c("#999999", "#999999", "#999999", "#999999", "#E69F00")) +
  
  # Scales
  scale_y_continuous(limits = c(0, 108), breaks = seq(0, 100, 25), expand = c(0,0)) +
  
  # Labels
  labs(
    #title = "Genotyping Accuracy by Tool and Gene",
    #subtitle = "Comparison with NCI-60 WES Dataset",
    y = "Accuracy (%)",
    x = "",
    fill = "Tool") +
  
  # Theme customization
  theme_bw() +
  coord_cartesian(ylim = c(0, 108)) +     # Set the visible limits here
  theme(
    #axis.text.x = element_text(angle = 45, hjust = 1), # Rotate x labels
    strip.background = element_rect(fill = "#f0f0f0"), # Facet header background
    strip.text = element_text(face = "bold", size = 24),
    panel.grid.major.x = element_blank(),
    legend.position = "none", # Hide legend since x-axis has labels
    axis.title = element_text(size = 20, face = "bold"), 
    axis.text = element_text(size = 18), 
    axis.text.x = element_text(angle = 45, hjust = 1, colour = "black", size = 20),
    axis.text.y = element_text(colour = "black"),
  )

# 4. Display Plot
print(p_1row)
ggsave(plot = p_1row, filename = "data/results/hlamajority/plots/results-all-samples-1000genomes/hlamajority-accuracy-per-gene-1000genomes-all-samples-1row.png", width = 15, height = 7)
ggsave(plot = p_1row, filename = "data/results/hlamajority/plots/results-all-samples-1000genomes/hlamajority-accuracy-per-gene-1000genomes-all-samples-1row.svg", width = 15, height = 7)

data_in_orig <- read.csv("../../data/claeys-et-al/claeys-et-al-benchmarking-results.csv")
data_in <- data_in_orig[1:2,1:4]
data_in_long <- pivot_longer(data = data_in, cols = c("A", "B", "C"), names_to = "HLA_Allele")
mhci <- ggplot(data_in_long, aes(fill=Feature, y=as.numeric(value), x=HLA_Allele)) +
  geom_col(position="dodge") +
  scale_y_continuous(expand = c(0, 0)) + # Keep this to remove bottom padding
  coord_cartesian(ylim=c(0,104)) +
  ylab("Accuracy (%)") +
  xlab("MHC Class I Gene") +
  scale_fill_tableau() +
  theme_bw() +
  scale_fill_discrete(c(""), labels = c("Best Individual Tool", "Metaclassifier")) +
  theme(axis.title=element_text(size=20), axis.text = element_text(size=16), legend.text = element_text(size = 16)) +
  geom_text( 
    position = position_dodge(width = 0.9), # Matches the bar dodge width
    aes(
      label = paste(round(as.numeric(value), 1), "%", sep = ""),
      group = Feature
    ),
    vjust = -0.2, # Position text just above the bar
    size = 4
  )
data_in_long$Method <- "Original_benchmark"
colnames(data_in_long) <- c("Tool", "Gene", "Accuracy", "Method")
#ylim(95.0,100.0) 
mhci

best_wes <- df %>% dplyr::filter(Tool != "nf-hlamajority" & Gene != "Overall") %>% group_by(Gene) %>% dplyr::filter(Accuracy == max(Accuracy)) %>% dplyr::select(Gene, Tool, Accuracy) %>% mutate(Tool = "best_wes") 
meta_wes <- df %>% dplyr::filter(Tool == "nf-hlamajority" & Gene != "Overall") %>%  dplyr::select(Gene, Tool, Accuracy)

best_wes$Tool <- "best_wes"
meta_wes$Tool <- "meta_wes"
full_results <- rbind.data.frame(best_wes, meta_wes)

full_results$Method <- "nf-hlamajority"
mhci_hlamajority <- ggplot(as.data.frame(full_results), aes(fill=Tool, y=as.numeric(Accuracy), x=Gene)) +
  geom_col(position="dodge") +
  scale_y_continuous(expand = c(0, 0)) + # Keep this to remove bottom padding
  coord_cartesian(ylim=c(0,104)) +
  ylab("Accuracy (%)") +
  xlab("MHC Class I Gene") +
  scale_fill_tableau() +
  theme_bw() +
  scale_fill_discrete(c(""), labels = c("Best Individual Tool", "Metaclassifier")) +
  theme(axis.title=element_text(size=20), axis.text = element_text(size=16), legend.text = element_text(size = 16)) +
  geom_text( 
    position = position_dodge(width = 0.9), # Matches the bar dodge width
    aes(
      label = paste(round(as.numeric(Accuracy), 1), "%", sep = ""),
      group = Tool
    ),
    vjust = -0.2, # Position text just above the bar
    size = 4
  )
mhci_hlamajority

# combine results
all_results <- rbind.data.frame(data_in_long, full_results)
all_results$Accuracy <- as.numeric(all_results$Accuracy)
all_results$Method <- factor(
  all_results$Method,
  levels = c("Original_benchmark", "nf-hlamajority"),
  labels = c("Original benchmark", "nf-hlamajority")
)

all_results$Gene <- factor(
  all_results$Gene,
  levels = c("A", "B", "C"),
  labels = c("HLA-A", "HLA-B", "HLA-C")
)

all_results$Tool <- factor(
  all_results$Tool,
  levels = c("best_wes", "meta_wes"),
  labels = c("Best WES", "Meta WES")
)

compare_hlamajority_claeys <- ggplot(all_results, aes(x = Gene, y = Accuracy, fill = Tool)) +
  
  # Create bars
  geom_col(position = position_dodge(), width = 0.7, color = "black", size = 0.2) +
  
  facet_wrap(~Method, scales = "fixed", 
             ncol = 2#,
             #labeller = as_labeller(my_gene_labels)
  ) +
  
  # Add text labels on top of bars (rounded to 1 decimal)
  geom_text(#aes(label = sprintf("%.1f", accuracy)), 
    aes(label = paste(round(as.numeric(Accuracy), 1), "%", sep = "")),
    position = position_dodge(width = 0.9), 
    vjust = -0.5, 
    size = 6) +
  
  # Colors: Highlight Hlamajority (assuming it's the last factor level)
  # You can customize these colors. 
  # Here: Greys for others, Red/Blue for Hlamajority
  #scale_fill_manual(values = c("#999999", "#999999", "#999999", "#999999", "#E69F00")) +
  
  # Scales
  scale_y_continuous(limits = c(0, 120), breaks = seq(0, 100, 25), expand = c(0,0)) +
  
  # Labels
  labs(
    #title = "Genotyping Accuracy by Tool and Gene",
    #subtitle = "Comparison with NCI-60 WES Dataset",
    y = "Accuracy (%)",
    x = "Gene",
    fill = "Tool") +
  
  # Theme customization
  theme_bw() +
  coord_cartesian(ylim = c(0, 115)) +     # Set the visible limits here
  # theme(
  #   #axis.text.x = element_text(angle = 45, hjust = 1), # Rotate x labels
  #   strip.background = element_rect(fill = "#f0f0f0"), # Facet header background
  #   strip.text = element_text(face = "bold", size = 18),
  #   panel.grid.major.x = element_blank(),
  #   #legend.position = "none", # Hide legend since x-axis has labels
  #   axis.title = element_text(size = 18), 
  #   axis.text = element_text(size = 16), 
  #   axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"),
  #   axis.text.y = element_text(colour = "black")
  # )
theme(
  #axis.text.x = element_text(angle = 45, hjust = 1), # Rotate x labels
  strip.background = element_rect(fill = "#f0f0f0"), # Facet header background
  strip.text = element_text(face = "bold", size = 24),
  panel.grid.major.x = element_blank(),
  #legend.position = "none", # Hide legend since x-axis has labels
  legend.text = element_text(size = 18),
  legend.title = element_text(size = 20, face = "bold"),
  axis.title = element_text(size = 20, face = "bold"), 
  axis.text = element_text(size = 18), 
  axis.text.x = element_text(angle = 45, hjust = 1, colour = "black", size = 20),
  axis.text.y = element_text(colour = "black"),
)
compare_hlamajority_claeys
ggsave(plot = compare_hlamajority_claeys, filename = "data/results/hlamajority/plots/results-all-samples-1000genomes/compare-claeys-hlamajority-1000genomes-all-samples.png", width = 15, height = 10)
ggsave(plot = compare_hlamajority_claeys, filename = "data/results/hlamajority/plots/results-all-samples-1000genomes/compare-claeys-hlamajority-1000genomes-all-samples.svg", width = 15, height = 10)

# redo this plot, but now compare the claeys et al benchmark with nf-hlamajority for each tool, for each gene
claeys.results.full <- read.csv("../../data/claeys-et-al/benchmarking_results_claeys_cleaned.csv")
df.for.comparison <- df %>% 
                     dplyr::select(Gene, Tool, Accuracy) %>% 
                     dplyr::filter(Gene != "Overall")
df.for.comparison$Tool <- gsub("nf-hlamajority", "Metaclassifier", df.for.comparison$Tool)
df.for.comparison$Study <- "nf-hlamajority"
df.for.comparison$Accuracy <- round(df.for.comparison$Accuracy, 1)
claeys.results.full$tool <- gsub(pattern = "hlala", replacement = "HLA*LA", x = claeys.results.full$tool)
claeys.results.full$tool <- gsub(pattern = "kourami", replacement = "Kourami", x = claeys.results.full$tool)
claeys.results.full$tool <- gsub(pattern = "optitype", replacement = "Optitype", x = claeys.results.full$tool)
claeys.results.full$tool <- gsub(pattern = "polysolver", replacement = "Polysolver", x = claeys.results.full$tool)
claeys.results.full.long <- tidyr::pivot_longer(data = claeys.results.full, cols = c("A", "B", "C"), names_to = "Gene", values_to = "Accuracy") %>% 
                            mutate(Accuracy = 100*Accuracy) %>% 
                            rename(Tool = "tool")
# now add metaclassifier results
# new_name = old_name
data_in_orig_meta_wes <- data_in_orig %>% dplyr::filter(Feature == "meta_wes") %>% 
                                          dplyr::select(Feature, A, B, C) %>% 
                                          pivot_longer(cols = c("A", "B", "C"), names_to = "Gene", values_to = "Accuracy") %>% 
                                          rename(`Tool` = "Feature") 
data_in_orig_meta_wes$Tool <- gsub("meta_wes", "Metaclassifier", data_in_orig_meta_wes$Tool)         
data_in_orig_meta_wes$Study <- "Original_benchmark"
claeys.results.full.long$Study <- "Original_benchmark"
claeys.results.full.long.combined <- rbind.data.frame(data_in_orig_meta_wes, claeys.results.full.long)
data_for_plotting <- rbind.data.frame(claeys.results.full.long.combined, df.for.comparison)
data_for_plotting$Accuracy <- as.numeric(data_for_plotting$Accuracy)

data_for_plotting$Study <- factor(
  data_for_plotting$Study,
  levels = c("Original_benchmark", "nf-hlamajority"),
  labels = c("Original benchmark", "nf-hlamajority")
)

data_for_plotting$Tool <- factor(
  data_for_plotting$Tool,
  levels = c("Kourami", "HLA*LA", "Polysolver", "Optitype", "Metaclassifier")#,
  #labels = c("Kourami", "HLA*LA", "Polysolver", "Optitype", "nf-hlamajority")
)

data_for_plotting$Gene <- factor(
  data_for_plotting$Gene,
  levels = c("A", "B", "C"),
  labels = c("HLA-A", "HLA-B", "HLA-C")
)

compare_hlamajority_claeys_per_gene_per_tool <- ggplot(data_for_plotting, aes(x = Study, y = Accuracy, fill = Study)) +
  # Create bars
  geom_col(position = position_dodge(), width = 0.7, color = "black", size = 0.2) +
 # facet_wrap(~Tool, scales = "free_y") + 
  facet_grid(
    Gene ~ Tool,
    scales = "free_y"#,
    #labeller = labeller(gene = my_gene_labels)
  ) +
  # Add text labels on top of bars (rounded to 1 decimal)
  geom_text(#aes(label = sprintf("%.1f", accuracy)), 
    aes(label = paste(round(as.numeric(Accuracy), 1), "%", sep = "")),
    position = position_dodge(width = 0.9), 
    vjust = -0.5, 
    size = 6) +
  # Scales
  scale_y_continuous(limits = c(0, 120), breaks = seq(0, 100, 25), expand = c(0,0)) +
  # Labels
  labs(
    y = "Accuracy (%)",
    x = "Tool",
    fill = "Study") +
  
  # Theme customization
  theme_bw() +
  coord_cartesian(ylim = c(0, 125)) +     # Set the visible limits here
  theme(
    #axis.text.x = element_text(angle = 45, hjust = 1), # Rotate x labels
    strip.background = element_rect(fill = "#f0f0f0"), # Facet header background
    strip.text = element_text(face = "bold", size = 20),
    panel.grid.major.x = element_blank(),
    legend.position = "bottom", # Hide legend since x-axis has labels
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 20, face = "bold"), 
    axis.text = element_text(size = 18), 
    axis.text.x = element_blank(),
    axis.text.y = element_text(colour = "black"),
  )
compare_hlamajority_claeys_per_gene_per_tool
ggsave(plot = compare_hlamajority_claeys_per_gene_per_tool, filename = "data/results/hlamajority/plots/results-all-samples-1000genomes/compare-hlamajority-claeys-per-gene-tool.png", width = 15, height = 7)
ggsave(plot = compare_hlamajority_claeys_per_gene_per_tool, filename = "data/results/hlamajority/plots/results-all-samples-1000genomes/compare-hlamajority-claeys-per-gene-tool.svg", width = 15, height = 7)
ggsave(plot = compare_hlamajority_claeys_per_gene_per_tool, filename = "data/results/hlamajority/plots/results-all-samples-1000genomes/compare-hlamajority-claeys-per-gene-tool.pdf", width = 15, height = 7)
ggsave(plot = compare_hlamajority_claeys_per_gene_per_tool, filename = "data/results/hlamajority/plots/results-all-samples-1000genomes/compare-hlamajority-claeys-per-gene-tool.pdf", width = 15, height = 7)

df_results <- results$summary %>% dplyr::mutate(Call_Rate = 100*Call_Rate)
df_results$Gene <- factor(df_results$Gene, levels = c("Overall", "A", "B", "C"))
# df_results$Tool <- factor(
#   df_results$Tool,
#   levels = c("kourami", "hlala", "polysolver", "optitype", "hlamajority"),
#   labels = c("Kourami", "HLA*LA", "Polysolver", "Optitype", "nf-hlamajority")
# )

p_call <- ggplot(df_results, aes(x = Tool, y = Call_Rate, fill = Tool)) +
  
  # Create bars
  geom_col(position = position_dodge(), width = 0.7, color = "black", size = 0.2) +
  
  # Facet by Gene
  facet_wrap(~Gene, scales = "fixed", 
             ncol = 3,
             labeller = as_labeller(my_gene_labels)
  ) +
  
  # Add text labels on top of bars (rounded to 1 decimal)
  geom_text(#aes(label = sprintf("%.1f", accuracy)), 
    aes(label = paste(round(Call_Rate, 1), "%", sep = "")),
    position = position_dodge(width = 0.9), 
    vjust = -0.5, 
    size = 7) +
  
  # Colors: Highlight Hlamajority (assuming it's the last factor level)
  # You can customize these colors. 
  # Here: Greys for others, Red/Blue for Hlamajority
  scale_fill_manual(values = c("#999999", "#999999", "#999999", "#999999", "#E69F00")) +
  
  # Scales
  scale_y_continuous(limits = c(0, 120), breaks = seq(0, 100, 25), expand = c(0,0)) +
  
  # Labels
  labs(
    y = "Call Rate (%)",
    x = "Tool",
    fill = "Tool") +
  
  # Theme customization
  theme_bw() +
  coord_cartesian(ylim = c(0, 112)) +     # Set the visible limits here
  theme(
    strip.background = element_rect(fill = "#f0f0f0"), # Facet header background
    strip.text = element_text(face = "bold", size = 24),
    panel.grid.major.x = element_blank(),
    legend.position = "none", # Hide legend since x-axis has labels
    axis.title = element_text(size = 20, face = "bold"), 
    axis.text = element_text(size = 18), 
    axis.text.x = element_text(angle = 45, hjust = 1, colour = "black", size = 20),
    axis.text.y = element_text(colour = "black"),
  )

# 4. Display Plot
print(p_call)

ggsave(plot = p_call, filename = "data/results/hlamajority/plots/results-all-samples-1000genomes/hlamajority-call-rate-per-gene-1000genomes-all-samples.png", width = 15, height = 7)
ggsave(plot = p_call, filename = "data/results/hlamajority/plots/results-all-samples-1000genomes/hlamajority-call-rate-per-gene-1000genomes-all-samples.svg", width = 15, height = 7)

# look at relationship between coverage and ability to call Kourami
#all.in.kourami <- 
all.in.kourami <- all.in %>% 
                  dplyr::filter(tool == "kourami") %>% 
                  mutate(`HLA-A` = paste(A1, A2, sep = ","),
                         `HLA-B` = paste(B1, B2, sep = ","),
                         `HLA-C` = paste(C1, C2, sep = ",")
                         )
all.in.kourami.long <- all.in.kourami %>%  pivot_longer(cols = c("HLA-A", "HLA-B", "HLA-C"), 
                                                    names_to = "gene",
                                                    values_to = "Genotype"
                                                    ) %>% 
                                          dplyr::select(sample, gene, tool, Genotype) %>% 
                                          mutate(is_called = ifelse(grepl("NA", Genotype), 0, 1))

all.in.kourami.depth <- depth %>% dplyr::select(sample, gene, mean_depth_hla_exons_2_3_gene) %>% right_join(all.in.kourami.long, by = c("sample", "gene"))
all.in.kourami.depth
head(all.in.kourami.depth)
# add depth to scores
# read in scores
scores <- read.csv("data/results/hlamajority/1000genomes-all-samples/1000genomes-score-per-sample.csv")
depth_scores <- scores %>%
  complete(gene, sample, tool) %>%   # create missing combinations
  left_join(depth, by = c("gene", "sample"))
# remove samples that are NA in the gold standard
gs.na <- read.csv("data/results/hlamajority/1000genomes-all-samples/1000-genomes-gs-na-samples.csv")
depth_scores_rm_na <- depth_scores %>% dplyr::filter(
  !(sample == "NA12234" & gene == "HLA-C") &
  !(sample == "NA12249" & gene == "HLA-B") &
  !(sample == "NA18548" & gene == "HLA-C") 
)
depth_scores_rm_na$correct_flag <- ifelse(depth_scores_rm_na$Score != 2 | is.na(depth_scores_rm_na$Score), "Incorrect", "Correct")
depth_scores_rm_na$correct_flag <- factor(depth_scores_rm_na$correct_flag, levels = c("Incorrect", "Correct"))

depth_scores_rm_na$tool <- factor(
  depth_scores_rm_na$tool,
  levels = c("kourami", "hlala", "polysolver", "optitype", "hlamajority"),
  labels = c("Kourami", "HLA*LA", "Polysolver", "Optitype", "nf-hlamajority")
)
summary_stats <- depth_scores_rm_na %>%
  group_by(tool, gene, correct_flag) %>%
  summarise(
    n = n(),
    mean_depth = mean(mean_depth_hla_exons_2_3_gene),
    median_depth = median(mean_depth_hla_exons_2_3_gene),
    sd_depth = sd(mean_depth_hla_exons_2_3_gene)
  ) %>%
  ungroup()

library(purrr)

depth_scores_rm_na %>%
  group_by(tool, gene) %>%
  summarise(p_value = wilcox.test(mean_depth_hla_exons_2_3_gene ~ correct_flag)$p.value)

ggplot(depth_scores_rm_na, aes(x = gene, y = mean_depth_hla_exons_2_3_gene, fill = correct_flag)) +
  geom_boxplot(outlier.alpha = 0.2, position = position_dodge(0.8)) +
  labs(
    x = "HLA gene",
    y = "log10(Mean HLA exon depth)",
    fill = "Call correctness"
  ) +
  # Facet by Tool with 2 columns to reduce awkward spacing
  facet_wrap(~tool, ncol = 2) +
  scale_fill_manual(values = c("Incorrect" = "#E69F00", "Correct" = "#56B4E9")) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),  # no rotation needed for HLA-A/B/C
    strip.background = element_rect(fill = "gray90", color = NA), # subtle facet background
    strip.text = element_text(face = "bold"),
    legend.position = "top"
  ) +
  scale_y_log10()

# library(lme4)
# 
# fit_lmm <- lmer(mean_depth_hla_exons_2_3_gene ~ correct_flag + (1|sample), data = depth_scores)
# summary(fit_lmm)

p_gene <- ggplot(depth_scores_rm_na, aes(x = correct_flag, y = mean_depth_hla_exons_2_3_gene, fill = correct_flag)) +
  geom_boxplot(outlier.shape = NA) + # Hides outlier points for a cleaner look
  geom_jitter(width = 0.2, alpha = 0.6) + # Adds individual data points
  facet_wrap(~tool, scales = "free_y") + # Creates a separate plot for each tool
  labs(title = "Gene-Specific Coverage by Tool and Call Accuracy",
       x = "Call Type",
       y = "Mean Depth of HLA Exons 2 & 3 (Gene)") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

print(p_gene)
library(ggpubr)
library(rstatix) 
stat_test_gene <- depth_scores_rm_na %>%
  group_by(tool) %>%
  wilcox_test(mean_depth_hla_exons_2_3_gene ~ correct_flag) %>%
  ungroup() %>%
  # 2. Adjust the p-values across all tests
  adjust_pvalue(method = "holm") %>%
  # 3. Add significance stars (optional, but nice for plots)
  add_significance("p.adj") %>%
  # 4. Get y-position for plotting the labels on the graph
  add_xy_position(x = "correct_flag", fun = "max", data = depth_scores_rm_na)
print(stat_test_gene)

p_gene_adjusted <- ggplot(depth_scores_rm_na, aes(x = correct_flag, y = mean_depth_hla_exons_2_3_gene, fill = correct_flag)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6) +
  facet_wrap(~tool, scales = "free_y") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  labs(title = "Gene-Specific Coverage by Tool and Call Accuracy",
       subtitle = "Holm test corrected p-values",
       x = "Call Type",
       y = "Mean Depth of HLA Exons 2 & 3 (Gene)") +
  theme_bw() +
  theme(legend.position = "none",
        axis.title = element_text(size = 18), 
        axis.text.x = element_text(size = 15, colour = "black"),
        axis.text.y = element_text(size = 13, colour = "black"),
        
        strip.text = element_text(size = 18, colour = "black")
        ) +
  # Add the adjusted p-values and significance bars from our table
  stat_pvalue_manual(
    stat_test_gene,
    label = "p.adj = {p.adj}, {p.adj.signif}", # Custom label
    tip.length = 0.01,
    bracket.nudge.y = 0.05,
    inherit.aes = FALSE
  )

print(p_gene_adjusted)
ggsave(plot = p_gene_adjusted, filename = "data/results/hlamajority/plots/compare-hlamajority-accuracy-per-gene-correct-incorrect-1000genomes-all-samples.png", width = 10, height = 7)

stat_test_gene_tool <- depth_scores_rm_na %>%
  group_by(gene, tool) %>%
  wilcox_test(mean_depth_hla_exons_2_3_gene ~ correct_flag) %>%
  ungroup() %>%
  adjust_pvalue(method = "holm") %>%      # adjust across ALL tests
  add_significance("p.adj") %>%
  add_xy_position(
    x = "correct_flag",
    fun = "max",
    data = depth_scores_rm_na
  )



p_gene_adjusted_gene_tool <- ggplot(depth_scores_rm_na, aes(x = correct_flag, y = mean_depth_hla_exons_2_3_gene, fill = correct_flag)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6) +
  #facet_wrap(~tool, scales = "free_y") +
  facet_grid(
    gene ~ tool,
    scales = "fixed"#,
    #labeller = labeller(gene = my_gene_labels)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  labs(title = "Gene-Specific Coverage by Tool and Call Accuracy",
       subtitle = "Holm test corrected p-values",
       x = "Call Type",
       y = "Mean Depth of HLA Exons 2 & 3") +
  theme_bw() +
  theme(legend.position = "none",
        axis.title = element_text(size = 18), 
        axis.text.x = element_text(size = 15, colour = "black"),
        axis.text.y = element_text(size = 13, colour = "black"),
        
        strip.text = element_text(size = 18, colour = "black")
  ) +
  # Add the adjusted p-values and significance bars from our table
  stat_pvalue_manual(
    stat_test_gene_tool,
    label = "p.adj = {p.adj}, {p.adj.signif}", # Custom label
    tip.length = 0.01,
    bracket.nudge.y = 0.05,
    inherit.aes = FALSE
  )
p_gene_adjusted_gene_tool
ggsave(plot = p_gene_adjusted_gene_tool, filename = "data/results/hlamajority/plots/compare-hlamajority-accuracy-per-gene-correct-incorrect-1000genomes-all-samples-facet-gene-tool.png", width = 10, height = 7)
ggsave(plot = p_gene_adjusted_gene_tool, filename = "data/results/hlamajority/plots/compare-hlamajority-accuracy-per-gene-correct-incorrect-1000genomes-all-samples-facet-gene-tool.svg", width = 10, height = 7)

# table of medians per group for reference
median_depth_correct_incorrect <- depth_scores_rm_na %>%
  group_by(gene, tool, correct_flag) %>% 
  summarise(median_coverage = median(mean_depth_hla_exons_2_3_gene))
write.csv(median_depth_correct_incorrect, file = "data/results/hlamajority/median-depth-incorrect-correct-per-tool-gene-1000genomes-1002samples.csv", row.names = F)

depth_scores_rm_na_all <- depth_scores_rm_na %>% group_by(gene) %>% distinct(sample, .keep_all = T) %>% dplyr::select(gene, sample, mean_depth_hla_exons_2_3_gene)
depth_scores_rm_na_all %>% group_by(gene) %>% summarise(median_coverage = median(mean_depth_hla_exons_2_3_gene))



p_gene_adjusted_gene_tool
summarise_gene_tool <- function(results, gene, tool) {
  
  x <- results$details[[gene]][[tool]]$error_types$Type
  
  tab <- table(x)
  
  data.frame(
    Gene  = gene,
    Tool  = tool,
    Type  = names(tab),
    Count = as.integer(tab),
    row.names = NULL
  )
}
tools <-  c("kourami", "hlala", "polysolver", "optitype", "hlamajority")
genes <- c("A", "B", "C")
final_df <- do.call(
  rbind,
  lapply(genes, function(g) {
    do.call(
      rbind,
      lapply(tools, function(t) {
        summarise_gene_tool(results, g, t)
      })
    )
  })
)

final_df_mismatches <- final_df |>
  dplyr::filter(Type != "Correct") |>
  group_by(Gene, Tool) |>
  mutate(
    Percent = Count / sum(Count) * 100
  ) |>
  ungroup()

final_df_mismatches$Type <- factor(final_df_mismatches$Type, 
                        levels = c("Dropout (Hetero -> Homo)", "Hallucination (Homo -> Hetero)", "Partial Mismatch", "Complete Mismatch"),
                        labels = c("Dropout", "Hallucination", "Partial Mismatch", "Complete Mismatch")
)

final_df_mismatches$Tool <- factor(final_df_mismatches$Tool, 
                        levels = c("kourami", "hlala", "polysolver", "optitype", "hlamajority"),
                        labels = c("Kourami", "HLA*LA", "Polysolver", "Optitype", "nf-hlamajority")
)

p_types_percent <- ggplot(
  final_df_mismatches,
  aes(x = Type, y = Percent, fill = Type)
) +
  geom_col(
    width = 0.7,
    color = "black",
    linewidth = 0.2
  ) +
  facet_grid(
    Gene ~ Tool,
    scales = "fixed",
    labeller = labeller(Gene = my_gene_labels)
  ) +
  scale_y_continuous(
    expand = c(0, 0)
  ) +
  labs(
    x = "Outcome type",
    y = "Percent of errors",
    fill = "Type"
  ) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "#f0f0f0"),
    strip.text = element_text(face = "bold", size = 18),
    panel.grid.major.x = element_blank(),
    # axis.text.x = element_text(
    #   angle = 45,
    #   hjust = 1,
    #   colour = "black"
    # ),
    axis.text.x = element_blank(),
    axis.text.y = element_text(colour = "black"),
    axis.title = element_text(size = 18),
    legend.position = "bottom",
    legend.title = element_text(size =18),
    legend.text = element_text(size = 16)
  ) +
  coord_cartesian(ylim = c(0, 100)) +
  scale_fill_manual(values = palette_mistake_types)

p_types_percent

# depth

ggsave(plot = p_types_percent, filename = "data/results/hlamajority/plots/call-types-1000genomes-1012samples-percent-errors-only.png", width = 10, height = 7)
ggsave(plot = p_types_percent, filename = "data/results/hlamajority/plots/call-types-1000genomes-1012samples-percent-errors-only.svg", width = 10, height = 7)

# Run plot
#plot_confusion_heatmap(incorrect_df, min_count = 3)
results$details$A$hlamajority$difficult_alleles
results$details$A$optitype$difficult_alleles

confusion_matrix_A <- plot_confusion_heatmap(results$details$A$hlamajority$metrics$gold_standard_vs_tool_incorrect_calls, min_count = 2, gene = "HLA-A", tool = "nf-hlamajority")
confusion_matrix_A
ggsave(filename = "data/results/hlamajority/plots/confusion-matrix/confusion-matrix-1000genomes-1012samples-hla-a.png", plot = confusion_matrix_A, width = 5)
results$details$A$hlamajority$metrics$gold_standard_vs_tool_incorrect
confusion_matrix_B <- plot_confusion_heatmap(results$details$B$hlamajority$metrics$gold_standard_vs_tool_incorrect_calls, min_count = 2, gene = "HLA-B", tool = "nf-hlamajority")
ggsave(filename = "data/results/hlamajority/plots/confusion-matrix/confusion-matrix-1000genomes-1012samples-hla-b.png", plot = confusion_matrix_B, width = 5)
confusion_matrix_C <- plot_confusion_heatmap(results$details$C$hlamajority$metrics$gold_standard_vs_tool_incorrect_calls, min_count = 2, gene = "HLA-C", tool = "nf-hlamajority")
ggsave(filename = "data/results/hlamajority/plots/confusion-matrix/confusion-matrix-1000genomes-1012samples-hla-c.png", plot = confusion_matrix_C, width = 5)

# confusion matrix optitype
confusion_matrix_A_optitype <- plot_confusion_heatmap(results$details$A$optitype$metrics$gold_standard_vs_tool_incorrect_calls, min_count = 2, gene = "HLA-A", tool = "Optitype")
confusion_matrix_A_optitype
ggsave(filename = "data/results/hlamajority/plots/confusion-matrix/optitype-hla-a-confusion-matrix-1000genomes-1012samples.png", plot = confusion_matrix_A_optitype, width = 5)
confusion_matrix_B_optitype <- plot_confusion_heatmap(results$details$B$optitype$metrics$gold_standard_vs_tool_incorrect_calls, min_count = 2, gene = "HLA-B", tool = "Optitype")
print(confusion_matrix_B_optitype)
ggsave(filename = "data/results/hlamajority/plots/confusion-matrix/optitype-hla-b-confusion-matrix-1000genomes-1012samples.png", plot = confusion_matrix_B_optitype, width = 5)
confusion_matrix_C_optitype <- plot_confusion_heatmap(results$details$C$optitype$metrics$gold_standard_vs_tool_incorrect_calls, min_count = 2, gene = "HLA-C", tool = "Optitype")
ggsave(filename = "data/results/hlamajority/plots/confusion-matrix/optitype-hla-c-confusion-matrix-1000genomes-1012samples.png", plot = confusion_matrix_C_optitype, width = 5)

# confusion matrix polysolver
confusion_matrix_A_polysolver <- plot_confusion_heatmap(results$details$A$polysolver$metrics$gold_standard_vs_tool_incorrect_calls, min_count = 2, gene = "HLA-A", tool = "Polysolver")
confusion_matrix_A_polysolver
ggsave(filename = "data/results/hlamajority/plots/confusion-matrix/polysolver-hla-a-confusion-matrix-1000genomes-1012samples.png", plot = confusion_matrix_A_polysolver, width = 5)
confusion_matrix_B_polysolver <- plot_confusion_heatmap(results$details$B$polysolver$metrics$gold_standard_vs_tool_incorrect_calls, min_count = 2, gene = "HLA-B", tool = "Polysolver")
print(confusion_matrix_B_polysolver)
ggsave(filename = "data/results/hlamajority/plots/confusion-matrix/polysolver-hla-b-confusion-matrix-1000genomes-1012samples.png", plot = confusion_matrix_B_polysolver, width = 5)
confusion_matrix_C_polysolver <- plot_confusion_heatmap(results$details$C$polysolver$metrics$gold_standard_vs_tool_incorrect_calls, min_count = 2, gene = "HLA-C", tool = "Polysolver")
ggsave(filename = "data/results/hlamajority/plots/confusion-matrix/polysolver-hla-c-confusion-matrix-1000genomes-1012samples.png", plot = confusion_matrix_C_polysolver, width = 5)

# confusion matrix HLA*LA
confusion_matrix_A_hlala <- plot_confusion_heatmap(results$details$A$hlala$metrics$gold_standard_vs_tool_incorrect_calls, min_count = 2, gene = "HLA-A", tool = "HLA*LA")
confusion_matrix_A_hlala
ggsave(filename = "data/results/hlamajority/plots/confusion-matrix/hlala-hla-a-confusion-matrix-1000genomes-1012samples.png", plot = confusion_matrix_A_hlala, width = 5)
confusion_matrix_B_hlala <- plot_confusion_heatmap(results$details$B$hlala$metrics$gold_standard_vs_tool_incorrect_calls, min_count = 2, gene = "HLA-B", tool = "HLA*LA")
print(confusion_matrix_B_hlala)
ggsave(filename = "data/results/hlamajority/plots/confusion-matrix/hlala-hla-b-confusion-matrix-1000genomes-1012samples.png", plot = confusion_matrix_B_hlala, width = 5)
confusion_matrix_C_hlala <- plot_confusion_heatmap(results$details$C$hlala$metrics$gold_standard_vs_tool_incorrect_calls, min_count = 2, gene = "HLA-C", tool = "HLA*LA")
ggsave(filename = "data/results/hlamajority/plots/confusion-matrix/hlala-hla-c-confusion-matrix-1000genomes-1012samples.png", plot = confusion_matrix_C_hlala, width = 5)

# confusion matrix Kourami
confusion_matrix_A_kourami <- plot_confusion_heatmap(results$details$A$kourami$metrics$gold_standard_vs_tool_incorrect_calls, min_count = 2, gene = "HLA-A", tool = "Kourami")
confusion_matrix_A_kourami
ggsave(filename = "data/results/hlamajority/plots/confusion-matrix/kourami-hla-a-confusion-matrix-1000genomes-1012samples.png", plot = confusion_matrix_A_kourami, width = 5)
confusion_matrix_B_kourami <- plot_confusion_heatmap(results$details$B$kourami$metrics$gold_standard_vs_tool_incorrect_calls, min_count = 2, gene = "HLA-B", tool = "Kourami")
print(confusion_matrix_B_kourami)
ggsave(filename = "data/results/hlamajority/plots/confusion-matrix/kourami-hla-b-confusion-matrix-1000genomes-1012samples.png", plot = confusion_matrix_B_kourami, width = 5)
confusion_matrix_C_kourami <- plot_confusion_heatmap(results$details$C$kourami$metrics$gold_standard_vs_tool_incorrect_calls, min_count = 2, gene = "HLA-C", tool = "Kourami")
ggsave(filename = "data/results/hlamajority/plots/confusion-matrix/kourami-hla-c-confusion-matrix-1000genomes-1012samples.png", plot = confusion_matrix_C_kourami, width = 5)


fasta_path <- "https://raw.githubusercontent.com/ANHIG/IMGTHLA/refs/heads/3630/hla_nuc.fasta"

hla_sequences <- readDNAStringSet(fasta_path)
hla_sequences
# 2. Helper Function to Clean Names
# The FASTA headers look like: "HLA:HLA00001 A*01:01:01:01 1098 bp"
# We want to extract just "A*01:01:01:01"
names(hla_sequences) <- sapply(strsplit(names(hla_sequences), " "), `[`, 2)

get_allele_distance("02:07", "02:01", "A", hla_sequences = hla_sequences)

hlaA <- hla_sequences[grepl("^A\\*", names(hla_sequences))]

# Extract allele name (e.g. A*02:01:01:01)
allele_names <- str_extract(names(hlaA), "A\\*[^ ]+")

# Extract exon 2+3 (positions 74–619)
exon23 <- subseq(hlaA, start = 74, end = 619)

# Create dataframe
hla_df <- tibble(
  allele_full = allele_names,
  allele_2field = str_extract(allele_full, "A\\*\\d+:\\d+"),
  seq = as.character(exon23)
)

bw <- 2.5
p_depth <- ggplot(
  depth,
  aes(mean_depth_hla_exons_2_3_classI)
) +
  geom_histogram(
    binwidth = bw,
    #fill = "grey70",
    color = "white",      # creates the gap
    linewidth = 0.6
  ) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 120)) +
  xlab("Mean sequencing depth at exons 2 and 3\nacross MHC Class I genes") +
  ylab("Number of samples") +
  theme_bw() +
  theme(axis.title = element_text(size= 15),
        axis.text = element_text(size = 15)) 
print(p_depth)

detailed.hlamajority.depth$mean_depth_hla_exons_2_3_classI

t.test(detailed.hlamajority.depth$mean_depth_hla_exons_2_3_classI, detailed.hlamajority$mean_depth_hla_exons_2_3)
mean(detailed.hlamajority.depth$mean_depth_hla_exons_2_3_classI)
mean(detailed.hlamajority$mean_depth_hla_exons_2_3)
dim(detailed.hlamajority)


# difference in depth correct not correct

