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
# read in data
setwd("/hlamajority-paper/external/mhc_genotyping/")
source("scripts/functions/ggroup_mapper.R")
source("scripts/functions/evaluate_predictions_functions.R")

results <- readRDS("data/results/hlamajority/nci-full-results-hlamajority-majority-vote.Rds")
results$summary$Tool <- factor(
  results$summary$Tool,
  levels = c("kourami", "hlala", "polysolver", "optitype", "hlamajority"),
  labels = c("Kourami", "HLA*LA", "Polysolver", "Optitype", "nf-hlamajority")
)
df <- read.csv("data/results/hlamajority/nci-full-stats-hlamajority-majority-vote.csv")
#df$Accuracy <- 100*(results_summary$Accuracy)
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
  coord_cartesian(ylim = c(0, 112)) +     # Set the visible limits here
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
    x = "Tool",
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
#ggsave(plot = p_1row, filename = "data/results/hlamajority/plots/hlamajority-accuracy-per-gene-nci60-1row.png", width = 15, height = 10)
#ggsave(plot = p_1row, filename = "data/results/hlamajority/plots/hlamajority-accuracy-per-gene-nci60-1row.svg", width = 15, height = 10)

ggsave(plot = p_1row, filename = "data/results/hlamajority/plots/hlamajority-accuracy-per-gene-nci60-1row.png", width = 15, height = 7)
ggsave(plot = p_1row, filename = "data/results/hlamajority/plots/hlamajority-accuracy-per-gene-nci60-1row.svg", width = 15, height = 7)

ggsave(plot = p_1row, filename = "../../figures/figures-thesis/figure-accuracy-cell-lines/hlamajority-accuracy-per-gene-nci60-1row.png", width = 15, height = 7)
ggsave(plot = p_1row, filename = "../../figures/figures-thesis/figure-accuracy-cell-lines/hlamajority-accuracy-per-gene-nci60-1row.svg", width = 15, height = 7)
ggsave(plot = p_1row, filename = "../../figures/figures-thesis/figure-accuracy-cell-lines/hlamajority-accuracy-per-gene-nci60-1row.pdf", width = 15, height = 7)

# call rate
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
  # theme(
  #   #axis.text.x = element_text(angle = 45, hjust = 1), # Rotate x labels
  #   strip.background = element_rect(fill = "#f0f0f0"), # Facet header background
  #   strip.text = element_text(face = "bold", size = 18),
  #   panel.grid.major.x = element_blank(),
  #   legend.position = "none", # Hide legend since x-axis has labels
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
    legend.position = "none", # Hide legend since x-axis has labels
    axis.title = element_text(size = 20, face = "bold"), 
    axis.text = element_text(size = 18), 
    axis.text.x = element_text(angle = 45, hjust = 1, colour = "black", size = 20),
    axis.text.y = element_text(colour = "black"),
  )

# 4. Display Plot
print(p_call)

ggsave(plot = p_call, filename = "data/results/hlamajority/plots/hlamajority-call-rate-per-gene-nci60.png", width = 15, height = 7)
ggsave(plot = p_call, filename = "data/results/hlamajority/plots/hlamajority-call-rate-per-gene-nci60.svg", width = 15, height = 7)

# coverage
# read in detailed table
#detailed.hlamajority <- read.table("../../data/cell-lines/benchmark-cell-lines-all-kourami-3-63-0-majority-vote/combined_results/nf_hlamajority_stats_combined_sorted.tsv", 
#                                   sep = "\t", 
#                                   #col.names = T, 
 #                                  )
detailed.hlamajority <- vroom("../../data/cell-lines/benchmark-cell-lines-all-kourami-3-63-0-majority-vote/combined_results/nf_hlamajority_stats_combined_sorted.tsv")
detailed.hlamajority.depth <- detailed.hlamajority %>% 
                              distinct(sample, mean_depth_hla_exons_2_3_classI) 
hist(detailed.hlamajority.depth$mean_depth_hla_exons_2_3_classI)

bw <- 2.5
p_depth <- ggplot(
  detailed.hlamajority.depth,
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
  coord_cartesian(ylim = c(0, 14)) +
  xlab("Mean sequencing depth at exons 2 and 3\nacross MHC Class I genes") +
  ylab("Number of samples") +
  theme_bw() +
  theme(axis.title = element_text(size= 15),
        axis.text = element_text(size = 15)) 
print(p_depth)

ggsave(plot = p_depth, filename = "data/results/hlamajority/plots/hlamajority-depth-nci60.png", width = 7, height = 7)
ggsave(plot = p_depth, filename = "data/results/hlamajority/plots/hlamajority-depth-nci60.svg", width = 7, height = 7)

results$details$A$hlamajority$difficult_alleles
results$details$B$hlamajority$difficult_alleles
results$details$C$hlamajority$difficult_alleles

A_results_summary <- table(results$details$A$hlamajority$error_types$Type)
B_results_summary <- table(results$details$B$hlamajority$error_types$Type)
C_results_summary <- table(results$details$C$hlamajority$error_types$Type)

# results_to_table <- function(results, gene){
#   results_summary <- table(results$details[[gene]]$hlamajority$error_types$Type)
#   results_table <- data.frame(Gene = `gene`,
#                               Type = names(results_summary),
#                               Count = as.vector(results_summary)
#   )
#   return(results_table)
# }

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

# final_df <- do.call(
#   rbind,
#   lapply(genes, function(g) results_to_table(results, g))
# )

final_df$Type <- factor(final_df$Type, 
                  levels = c("Correct", "Dropout (Hetero -> Homo)", "Hallucination (Homo -> Hetero)", "Partial Mismatch", "Complete Mismatch"),
                  labels = c("Correct", "Dropout", "Hallucination", "Partial Mismatch", "Complete Mismatch")
)

final_df$Tool <- factor(final_df$Tool, 
                  levels = c("kourami", "hlala", "polysolver", "optitype", "hlamajority"),
                  labels = c("Kourami", "HLA*LA", "Polysolver", "Optitype", "nf-hlamajority")
)

p_types <- ggplot(
  final_df,
  aes(x = Type, y = Count, fill = Type)
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
    y = "Number of samples",
    fill = "Type"
  ) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "#f0f0f0"),
    strip.text = element_text(face = "bold", size = 14),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(colour = "black"),
    axis.title = element_text(size = 16),
    legend.position = "bottom"
  ) +
  coord_cartesian(ylim = c(0, 35)) +
  scale_fill_manual(values = palette_mistake_types)
  
p_types
ggsave(plot = p_types, filename = "data/results/hlamajority/plots/call-types-nci60.png", width = 10, height = 7)
#install.packages("svglite")
ggsave(plot = p_types, filename = "data/results/hlamajority/plots/call-types-nci60.svg", width = 10, height = 7)

final_df_mismatches <- final_df |>
  dplyr::filter(Type != "Correct") |>
  group_by(Gene, Tool) |>
  mutate(
    Percent = Count / sum(Count) * 100
  ) |>
  ungroup()

# final_df_mismatches$Type <- factor(final_df_mismatches$Type, 
#                                    levels = c("Dropout (Hetero -> Homo)", "Hallucination (Homo -> Hetero)", "Partial Mismatch", "Complete Mismatch"),
#                                    labels = c("Dropout", "Hallucination", "Partial Mismatch", "Complete Mismatch")
# )
# 
# final_df_mismatches$Tool <- factor(final_df_mismatches$Tool, 
#                                    #levels = c("kourami", "hlala", "polysolver", "optitype", "hlamajority"),
#                                    levels = c("Kourami", "HLA*LA", "Polysolver", "Optitype", "nf-hlamajority")
# )

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

final_df <- final_df |>
  group_by(Gene, Tool) |>
  mutate(
    Percent = Count / sum(Count) * 100
  ) |>
  ungroup()

# p_types_percent <- ggplot(
#   final_df,
#   aes(x = Type, y = Percent, fill = Type)
# ) +
#   geom_col(
#     width = 0.7,
#     color = "black",
#     linewidth = 0.2
#   ) +
#   facet_grid(
#     Gene ~ Tool,
#     scales = "fixed",
#     labeller = labeller(Gene = my_gene_labels)
#   ) +
#   scale_y_continuous(
#     expand = c(0, 0)
#   ) +
#   labs(
#     x = "Outcome type",
#     y = "Percent of samples",
#     fill = "Type"
#   ) +
#   theme_bw() +
#   theme(
#     strip.background = element_rect(fill = "#f0f0f0"),
#     strip.text = element_text(face = "bold", size = 14),
#     panel.grid.major.x = element_blank(),
#     # axis.text.x = element_text(
#     #   angle = 45,
#     #   hjust = 1,
#     #   colour = "black"
#     # ),
#     axis.text.x = element_blank(),
#     axis.text.y = element_text(colour = "black"),
#     axis.title = element_text(size = 16),
#     legend.position = "bottom"
#   ) +
#   coord_cartesian(ylim = c(0, 100)) +
#   scale_fill_manual(values = palette_mistake_types)
# 
# p_types_percent

ggsave(plot = p_types_percent, filename = "data/results/hlamajority/plots/call-types-nci60-percent-errors-only.png", width = 10, height = 7)
#install.packages("svglite")
ggsave(plot = p_types_percent, filename = "data/results/hlamajority/plots/call-types-nci60-percent-errors-only.svg", width = 10, height = 7)

gold.standard.calls <- results$summary %>% dplyr::distinct(Gene, Num_GS_Called)
gold.standard.calls$Gene <- factor(gold.standard.calls$Gene, 
                  levels = c("A", "B", "C"),
                  labels = c("HLA-A", "HLA-B", "HLA-C")
)
p_gold_standard_calls <- ggplot(gold.standard.calls, aes(Num_GS_Called, fct_reorder(Gene, Num_GS_Called))) +
  geom_col() +
  xlab("Number of unambiguous gold standard calls") +
  ylab(NULL) +
  theme_minimal() +
  scale_x_continuous(expand = c(0, 0)) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    #axis.text.y = element_text(colour = "black", face = "bold"),
    axis.text = element_text(colour = "black", size = 12),
    axis.title.x = element_text(colour = "black", size = 12)
    
  ) +
  geom_text(aes(label = Num_GS_Called),
            position = position_dodge(width = 0.9),
            hjust = -0.5,
            size = 4) +
    coord_cartesian(xlim = c(0, 42))
p_gold_standard_calls

ggsave(plot = p_gold_standard_calls, filename = "data/results/hlamajority/plots/count-calls-gold-standard-nci60.png", width = 10, height = 7)
#install.packages("svglite")
ggsave(plot = p_gold_standard_calls, filename = "data/results/hlamajority/plots/count-calls-gold-standard-nci60.svg", width = 10, height = 7)
ggsave(plot = p_gold_standard_calls, filename = "data/results/hlamajority/plots/count-calls-gold-standard-nci60.pdf", width = 10, height = 7)

results$details$A$hlamajority$metrics$gold_standard_vs_tool_incorrect_calls
results$details$A$hlamajority$metrics$scores_per_person

results$details$A$hlamajority$difficult_alleles

confusion_matrix_A <- plot_confusion_heatmap(results$details$A$hlamajority$metrics$gold_standard_vs_tool_incorrect_calls, min_count = 1, gene = "HLA-A")
confusion_matrix_A
ggsave(filename = "data/results/hlamajority/plots/confusion-matrix-nci60-hla-a.png", plot = confusion_matrix_A, width = 5)
confusion_matrix_B <- plot_confusion_heatmap(results$details$B$hlamajority$metrics$gold_standard_vs_tool_incorrect_calls, min_count = 1, gene = "HLA-B")
ggsave(filename = "data/results/hlamajority/plots/confusion-matrix-nci60-hla-b.png", plot = confusion_matrix_B, width = 5)
confusion_matrix_C <- plot_confusion_heatmap(results$details$C$hlamajority$metrics$gold_standard_vs_tool_incorrect_calls, min_count = 1, gene = "HLA-C")
ggsave(filename = "data/results/hlamajority/plots/confusion-matrix-nci60-hla-c.png", plot = confusion_matrix_C, width = 5)

# Sequence similarity

fasta_path <- "https://raw.githubusercontent.com/ANHIG/IMGTHLA/refs/heads/3630/hla_nuc.fasta"

hla_sequences <- readDNAStringSet(fasta_path)

# 2. Helper Function to Clean Names
# The FASTA headers look like: "HLA:HLA00001 A*01:01:01:01 1098 bp"
# We want to extract just "A*01:01:01:01"
names(hla_sequences) <- sapply(strsplit(names(hla_sequences), " "), `[`, 2)

get_allele_distance("33:01", "25:01", "A", hla_sequences = hla_sequences)
get_allele_distance("33:01", "29:01", "A", hla_sequences = hla_sequences)
get_allele_distance("33:01", "30:03", "A", hla_sequences = hla_sequences)

get_allele_distance("15:05", "57:01", "B")

get_allele_distance("07:02", "07:17", "C")
get_allele_distance("12:02", "12:03", "C")

results$details$A$hlamajority$metrics$gold_standard_vs_tool_incorrect_calls
cell.lines.3301 <- c("T47D", "IGROV1", "TK-10")
cell.lines.3301.SRX <- 
all.calls <- vroom("../../data/cell-lines/benchmark-cell-lines-all-kourami-3-63-0-majority-vote/combined_results/nf_hlamajority_all_calls_sorted.tsv")
cell.line.map <- readRDS("data/results/hlamajority/nci-map.Rds")
cell.lines.3301.SRX <- cell.line.map %>% dplyr::filter(sample %in% cell.lines.3301) %>% dplyr::pull(sample_id) %>% unique() 
all.calls %>% dplyr::filter(sample %in% cell.lines.3301.SRX)

fasta_path_prot <- "https://raw.githubusercontent.com/ANHIG/IMGTHLA/refs/heads/3630/hla_prot.fasta"
hla_sequences_prot <- readDNAStringSet(fasta_path_prot)
names(hla_sequences_prot) <- sapply(strsplit(names(hla_sequences_prot), " "), `[`, 2)


# depth
hlamajority.depth <- detailed.hlamajority %>% dplyr::select(sample, gene, mean_depth_hla_exons_2_3_gene)
depth.scores <- read.csv("data/results/hlamajority/nci60-score-depth-per-sample-per-tool.csv", header = T)
master_df_mapped_full <- readRDS("data/results/hlamajority/nci-map.Rds") %>% dplyr::select(sample, sample_id) %>% distinct(sample, sample_id)

# read in all tool results
cell.lines.all.results <- read.table("../../data/cell-lines/benchmark-cell-lines-all-kourami-3-63-0-majority-vote/combined_results/nf_hlamajority_all_calls_sorted.tsv", sep = "\t", header = T)
cell.lines.all.results.samples <- cell.lines.all.results %>% dplyr::select(sample, tool)
# read in all hlamajority results 
hlamajority.in <- read.table("../../data/cell-lines/benchmark-cell-lines-all-kourami-3-63-0-majority-vote/combined_results/nf_hlamajority_votes_combined_sorted.tsv", sep = "\t", header = T)
samples <- cell.lines.all.results.samples %>% distinct(sample)
samples$tool <- "hlamajority"
samples.tools.complete <- rbind.data.frame(samples, cell.lines.all.results.samples)
names(samples.tools.complete) <- c("sample_id", "tool")
samples.tools.complete <- full_join(samples.tools.complete, master_df_mapped_full)
Gene <- c("A", "B", "C")
samples.tools.complete.gene <- samples.tools.complete
samples.tools.complete.gene$Gene <- Gene
samples.tools.complete.gene <- samples.tools.complete.gene %>% complete(tool, Gene, sample) 
# now only keep those not missing in gold standard
samples.missing.gs <- results$gold_standard_missing
samples.tools.complete.rm.na.gs <- samples.tools.complete %>% dplyr::filter(
  !(sample %in% "NA12234" & gene == "HLA-C") &
    !(sample == "NA12249" & gene == "HLA-B") &
    !(sample == "NA18548" & gene == "HLA-C") 
)

samples.tools.complete.gene.filtered <- samples.tools.complete.gene %>%
  dplyr::anti_join(
    samples.missing.gs,
    by = c("sample" = "Sample", "Gene" = "Gene")
  )
samples.tools.complete.gene.filtered$Gene <- paste("HLA-", samples.tools.complete.gene.filtered$Gene, sep = "")
# left join these with scores - NA score will be due to missing score
left_join(samples.tools.complete.gene.filtered, depth.scores, by = c("Gene" = "gene", "sample" = "sample", "tool" = "tool")) %>% View() #%>% nrow()
# filter to only ones with not missing gold standard


#scores <- read.csv("data/results/hlamajority/nci60-score-per-sample-per-tool.csv", header = T)
detailed.hlamajority.depth <- hlamajority.depth %>% rename(`sample_id` = "sample") %>% left_join(master_df_mapped_full) #%>% dplyr::select(-all_of(c(A1, A2, B1, B2, C1, C2))) #rename(sample_id = sample)
#depth_scores <- scores %>%
#  complete(gene, sample, tool) %>%   # create missing combinations
#  left_join(depth, by = c("gene", "sample"))
depth_scores_complete <- depth.scores %>%
  complete(gene, sample, tool) %>%   # create missing combinations
  left_join(detailed.hlamajority.depth, by = c("gene", "sample"))

depth_scores <- left_join(depth.scores, detailed.hlamajority.depth)
depth_scores
depth_scores$tool <- factor(
  depth_scores$tool,
  levels = c("kourami", "hlala", "polysolver", "optitype", "hlamajority"),
  labels = c("Kourami", "HLA*LA", "Polysolver", "Optitype", "nf-hlamajority")
)
depth_scores$correct_flag <- ifelse(depth_scores$Score == 2, "Correct", "Incorrect")
depth_scores$correct_flag <- factor(depth_scores$correct_flag, levels = c("Incorrect", "Correct"))

library(ggpubr)
library(rstatix) 
library(broom)

stat_test_gene <- depth_scores %>%
  group_by(tool) %>%
  wilcox_test(mean_depth_hla_exons_2_3_gene ~ correct_flag) %>%
  ungroup() %>%
  # 2. Adjust the p-values across all tests
  adjust_pvalue(method = "holm") %>%
  # 3. Add significance stars (optional, but nice for plots)
  add_significance("p.adj") %>%
  # 4. Get y-position for plotting the labels on the graph
  add_xy_position(x = "correct_flag", fun = "max", data = depth_scores)
print(stat_test_gene)

# wilcox_results <- depth_scores %>%
#   group_by(gene, tool) %>%
#   filter(n_distinct(correct_flag) == 2) %>%  # ensure both groups exist
#   summarise(
#     test = list(wilcox.test(mean_depth_hla_exons_2_3_gene ~ correct_flag,
#                             exact = FALSE)),
#     .groups = "drop"
#   ) %>%
#   mutate(tidy = purrr::map(test, broom::tidy)) %>%
#   tidyr::unnest(tidy) %>%
#   select(gene, tool, statistic, p.value)

stat_test_gene_tool <- depth_scores %>%
  group_by(gene, tool) %>%
  wilcox_test(mean_depth_hla_exons_2_3_gene ~ correct_flag) %>%
  ungroup() %>%
  adjust_pvalue(method = "holm") %>%      # adjust across ALL tests
  add_significance("p.adj") %>%
  add_xy_position(
    x = "correct_flag",
    fun = "max",
    data = depth_scores
  )

print(stat_test_gene_tool)


p_gene_adjusted <- ggplot(depth_scores, aes(x = correct_flag, y = mean_depth_hla_exons_2_3_gene, fill = correct_flag)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6) +
  facet_wrap(~tool, scales = "free_y") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  labs(title = "Gene-Specific Coverage by Tool and Call Accuracy",
       subtitle = "Benjamini-Hochberg corrected p-values",
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

stat_test_gene_tool <- depth_scores %>%
  group_by(gene, tool) %>%
  wilcox_test(mean_depth_hla_exons_2_3_gene ~ correct_flag) %>%
  ungroup() %>%
  adjust_pvalue(method = "holm") %>%      # adjust across ALL tests
  add_significance("p.adj") %>%
  add_xy_position(
    x = "correct_flag",
    fun = "max",
    data = depth_scores
  )



p_gene_adjusted_gene_tool <- ggplot(depth_scores, aes(x = correct_flag, y = mean_depth_hla_exons_2_3_gene, fill = correct_flag)) +
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
ggsave(plot = p_gene_adjusted_gene_tool, filename = "data/results/hlamajority/plots/compare-hlamajority-accuracy-per-gene-correct-incorrect-nci60-facet-gene-tool.png", width = 10, height = 7)
ggsave(plot = p_gene_adjusted_gene_tool, filename = "data/results/hlamajority/plots/compare-hlamajority-accuracy-per-gene-correct-incorrect-nci60-facet-gene-tool.svg", width = 10, height = 7)
