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
colnames(gold.standard.1kg)[1] <- "sample"
colnames(gold.standard.1kg)[2:7] <- gsub(pattern = "\\.", replacement = "", colnames(gold.standard.1kg[2:7]))
gold.standard.1kg <- as.data.frame(gold.standard.1kg)
source("scripts/functions/ggroup_mapper.R")
# 
# # 2. Define the Mapping Function
# # This takes a specific allele (e.g. "07:01") and converts it to the G-Group representative
get_g_group <- function(allele_in, gene) {
  # Clean input (remove 'A*', handle NAs)
  if(is.na(allele_in)) return(NA)

  # Ensure 4 digit format for lookup if input is 2 digit
  # (This is a simplified logic, you might need the more complex 'ga_map' from Claeys)

  # Check if this allele is part of a group in the g_groups table
  # You must filter g_groups by the specific gene (A, B, or C)
  match <- g_groups %>%
    filter(gene == !!gene, allele == paste0(gene, "*", allele_in))

  if(nrow(match) > 0) {
    # Return the group name (trimmed to 2-field)
    return(str_remove(match$key[1], "^[ABC]\\*"))
  } else {
    # If no G-group, just return the allele (trimmed to 2-field)
    return(allele_in)
  }
}
# 
# # Wrapper to apply to a whole list of dataframes
normalize_predictions <- function(prediction_list, gene) {
  lapply(prediction_list, function(pred_vec) {
    # pred_vec is c("01:01", "02:01")
    mapped <- sapply(pred_vec, get_g_group, gene = gene)
    return(unname(mapped))
  })
}

clean_list <- function(input_list, gene_name) {
  lapply(input_list, function(x) {
    # 1. Map to G-group (using Claeys' ga_map function)
    # ga_map usually expects "A", "B", "C" as gene name
    mapped_x <- ga_map(gene_name, x)
    
    # 2. If mapped_x returns multiple options (ambiguity), 
    # Claeys took the first one. You should too for consistency.
    final_x <- sapply(mapped_x, function(y) head(y, 1))
    
    return(final_x)
  })
}

gold.standard.1kg <- gold.standard.1kg %>% 
                     dplyr::filter(sample %in% hlamajority.in$sample)

gold.standard.1kg <- gold.standard.1kg[match(hlamajority.in$sample, gold.standard.1kg$sample),]

# fix 
df_to_list_fix <- function(df, cols) {
  # 1. Identify Sample IDs to use as list keys
  # If a 'sample' column exists, use it; otherwise use rownames
  if ("sample" %in% colnames(df)) {
    sample_ids <- df$sample
  } else {
    sample_ids <- rownames(df)
  }
  
  # 2. Extract only the relevant columns safely
  # as.matrix converts tibbles/dataframes to a base matrix
  # as.character ensures Factors are converted to Strings (prevents sort errors)
  sub_matrix <- as.matrix(df[, cols, drop = FALSE])
  mode(sub_matrix) <- "character" 
  
  # 3. Apply logic row by row using lapply (faster than for loop)
  output <- lapply(1:nrow(sub_matrix), function(i) {
    vec <- sub_matrix[i, ]
    
    # Check if row is entirely empty/NA
    if (all(is.na(vec) | vec == "")) {
      vec_sorted <- rep(NA, length(cols))
    } else {
      # 4. Sort safely
      # na.last = TRUE ensures NAs are kept at the end, preserving vector length
      vec_sorted <- sort(vec, na.last = TRUE)
    }
    
    # 5. Assign names (length is guaranteed to match now)
    names(vec_sorted) <- cols
    return(vec_sorted)
  })
  
  # 6. Assign sample IDs to the list
  names(output) <- sample_ids
  
  return(output)
}

# A_list_gold_standard_notna <- A_list_gold_standard[not_na(A_list_gold_standard)]
# A_list_hlamajority_for_comparison <- A_list_hlamajority[not_na(A_list_gold_standard)]
# all(identical(names(A_list_hlamajority_for_comparison), names(A_list_gold_standard_notna)))

# Function to count correct alleles (0, 1, or 2)
count_correct_alleles <- function(gold, call) {
  # 1. Check for NAs
  if (any(is.na(gold)) || any(is.na(call))) return(0)
  
  # 2. Make copies to manipulate
  g <- as.character(gold)
  c <- as.character(call)
  
  score <- 0
  
  # 3. Check first called allele
  if (c[1] %in% g) {
    score <- score + 1
    # Remove the matched allele from gold so it isn't matched twice
    # (This handles the homozygous vs heterozygous distinction)
    match_idx <- match(c[1], g)
    g <- g[-match_idx]
  }
  
  # 4. Check second called allele against the REMAINING gold allele
  if (length(c) > 1 && c[2] %in% g) {
    score <- score + 1
  } 
  return(score)
}

run_comparison <- function(gold.standard, calls, tool_name) {
  common_ids <- intersect(names(gold.standard), names(calls))
  
  # Initial alignment
  gold_aligned <- gold.standard[common_ids]
  call_aligned <- calls[common_ids]
  
  # --- NEW: FILTERING STEP ---
  # Identify samples where the CALLS have any NAs
  # (We assume if Gold Standard has NA, we should also exclude it)
  valid_indices <- mapply(function(g, c) {
    !any(is.na(c)) && !any(is.na(g))
  }, gold_aligned, call_aligned)
  
  # Filter lists to keep ONLY valid samples
  gold_aligned <- gold_aligned[valid_indices]
  call_aligned <- call_aligned[valid_indices]
  
  # Check if we removed everyone (edge case)
  if (length(gold_aligned) == 0) {
    print("Warning: No valid calls remaining after NA exclusion.")
    return(list(final_accuracy = 0, scores_per_person = NULL))
  }
  
  # Calculate scores for every individual at once
  # This returns a vector of 0, 1, or 2 for each person
  scores_per_person <- mapply(count_correct_alleles, gold_aligned, call_aligned)
  scores_incorrect_calls <- scores_per_person[which(scores_per_person < 2)]
  ids_incorrect_calls <- names(scores_incorrect_calls)
  gold_standard_incorrect_calls <- gold_aligned[ids_incorrect_calls]
  incorrect_calls <- call_aligned[ids_incorrect_calls]
  # Calculate final accuracy
  total_correct_alleles <- sum(scores_per_person)
  max_possible_alleles <- 2 * length(scores_per_person)
  
  # make gold standard of incorrect calls into a table
  gold_standard_calls_from_incorrect_tool_calls <- do.call(rbind.data.frame, gold_standard_incorrect_calls)
  colnames(gold_standard_calls_from_incorrect_tool_calls) <- c("allele1_gold_standard", "allele2_gold_standard")
  rownames(gold_standard_calls_from_incorrect_tool_calls) <- names(gold_standard_incorrect_calls)

  # make incorrect calls into a table
  incorrect_tool_calls <- do.call(rbind.data.frame, incorrect_calls)
  colnames(incorrect_tool_calls) <- c("allele1_tool", "allele2_tool")
  rownames(incorrect_tool_calls) <- names(incorrect_calls)
  
  # join tables
  gold_standard_vs_tool_incorrect_calls <- cbind.data.frame(gold_standard_calls_from_incorrect_tool_calls,
                                                            incorrect_tool_calls
                                                            )
  
  gold_standard_vs_tool_incorrect_calls$score <- scores_incorrect_calls
  
  final_accuracy <- total_correct_alleles / max_possible_alleles
  print(paste("Tool name: ", tool_name, sep = ""))
  print(paste("Accuracy:", round(final_accuracy * 100, 2), "%"))
  num_excluded <- length(common_ids) - length(valid_indices[valid_indices == TRUE])
  proportion_excluded <- num_excluded/length(common_ids)
  print(paste("Excluded ", num_excluded, " samples (", round(100*proportion_excluded, 2), "%) due to NA calls.", sep = ""))
  print(paste("Accuracy (Callers only):", round(final_accuracy * 100, 2), "%"))
  
  outs <- list(final_accuracy = final_accuracy, 
               scores_per_person = scores_per_person, 
               ids_incorrect_calls = ids_incorrect_calls, 
               gold_standard_vs_tool_incorrect_calls = gold_standard_vs_tool_incorrect_calls,
               proportion_excluded = proportion_excluded
               )
  return(outs)
}

analyze_error_types <- function(gold_std, tool_calls, tool_name) {
  
  common_ids <- intersect(names(gold_std), names(tool_calls))
  
  # Initialize results dataframe
  results <- data.frame(
    Sample = character(),
    Type = character(), # Homo, Hetero
    Outcome = character(), # Correct, Dropout, Hallucination, Partial, Mismatch
    stringsAsFactors = FALSE
  )
  
  for (id in common_ids) {
    g <- as.character(gold_std[[id]])
    c <- as.character(tool_calls[[id]])
    
    # Skip NAs
    if (any(is.na(g)) || any(is.na(c))) next
    
    # 1. Determine Zygosity States
    g_is_homo <- g[1] == g[2]
    c_is_homo <- c[1] == c[2]
    
    # 2. Count Correct Alleles (Using your logic)
    # We use a temp version of g to handle the counting
    g_temp <- g
    score <- 0
    if (c[1] %in% g_temp) {
      score <- score + 1
      g_temp <- g_temp[-match(c[1], g_temp)]
    }
    if (length(c) > 1 && c[2] %in% g_temp) {
      score <- score + 1
    }
    
    # 3. Categorize
    outcome <- ""
    
    if (score == 2) {
      outcome <- "Correct"
    } else if (score == 0) {
      outcome <- "Complete Mismatch"
    } else if (score == 1) {
      # This is the interesting part: Why did it get 1 wrong?
      
      if (g_is_homo && !c_is_homo) {
        # Gold: A, A | Tool: A, B
        outcome <- "Hallucination (Homo -> Hetero)"
      } else if (!g_is_homo && c_is_homo) {
        # Gold: A, B | Tool: A, A
        outcome <- "Dropout (Hetero -> Homo)"
      } else {
        # Gold: A, B | Tool: A, C
        # Both are hetero, but tool got one wrong
        outcome <- "Partial Mismatch"
      }
    }
    
    results <- rbind(results, data.frame(Sample = id, Type = outcome, Tool = tool_name))
  }
  
  return(results)
}

A_list_gold_standard <- df_to_list_fix(gold.standard.1kg, cols = c("A1", "A2"))
B_list_gold_standard <- df_to_list_fix(gold.standard.1kg, cols = c("B1", "B2"))
C_list_gold_standard <- df_to_list_fix(df = gold.standard.1kg, cols = c("C1", "C2"))

tools <- c("optitype", "polysolver", "hlala", "kourami")

# get accuracy for all tools A
accuracy_A <- c()
A_results <- list()
for (i in 1:length(tools)){
  A_list_tool <- all.in %>% dplyr::filter(tool == tools[i]) %>% df_to_list_fix(., cols = c("A1", "A2"))
  A_list_tool_cleaned <- clean_list(A_list_tool, "A")
  all_results <- run_comparison(A_list_gold_standard, A_list_tool_cleaned, tool_name = tools[i])
  accuracy_A[i] <- all_results$final_accuracy
  names(accuracy_A)[i] <- tools[i]
  A_results[[length(A_results) + 1]] <- all_results
  names(A_results)[i] <- tools[i]
  df <- analyze_error_types(A_list_gold_standard, A_list_tool_cleaned, tool_name = tools[i])
  
}

A_list_hlamajority <- df_to_list(hlamajority.in, cols = c("A1", "A2"))
names(A_list_hlamajority) <- hlamajority.in$sample
A_list_majority_clean <- clean_list(A_list_hlamajority, "A")
all_results <- run_comparison(A_list_gold_standard, A_list_majority_clean, tool_name = "hlamajority")
accuracy_A[length(tools)+1] <- all_results$final_accuracy
A_results[[length(A_results) + 1]] <- all_results
names(A_results)[length(A_results)] <- "hlamajority"

table(df$Tool, df$Type)

# Filter out "Correct" to zoom in on the errors
errors_only <- df %>% filter(Type != "Correct")

# Calculate proportions for the plot
errors_summary <- errors_only %>%
  group_by(Tool, Type) %>%
  summarise(Count = n()) %>%
  mutate(Proportion = Count / sum(Count))

ggplot(errors_summary, aes(x = Tool, y = Count, fill = Type)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_minimal() +
  labs(title = "Error Profiles by Tool (HLA-A)",
       subtitle = "breakdown of incorrect calls only",
       y = "Number of Errors",
       fill = "Error Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_brewer(palette = "Set1")



tab_incorrect <- A_results$kourami$gold_standard_vs_tool_incorrect_calls
tab_incorrect_hlamajority <- A_results$hlamajority$gold_standard_vs_tool_incorrect_calls

optitype_incorrect_calls_A <- A_results$optitype$incorrect_tool_calls
optitype_incorrect_calls_A_gold_standard <- A_results$optitype$gold_standard_calls_from_incorrect_tool_calls

accuracy_B <- c()
B_results <- list()

for (i in 1:length(tools)){
  B_list_tool <- all.in %>% dplyr::filter(tool == tools[i]) %>% df_to_list_fix(., cols = c("B1", "B2"))
  B_list_tool_cleaned <- clean_list(B_list_tool, "B")
  all_results <- run_comparison(B_list_gold_standard, B_list_tool_cleaned, tool_name = tools[i])
  accuracy_B[i] <- all_results$final_accuracy
  names(accuracy_B)[i] <- tools[i]
  B_results[[length(B_results) + 1]] <- all_results
  names(B_results)[i] <- tools[i]
}

B_list_hlamajority <- df_to_list_fix(hlamajority.in, cols = c("B1", "B2"))
names(B_list_hlamajority) <- hlamajority.in$sample
B_list_majority_clean <- clean_list(B_list_hlamajority, "B")
all_results <- run_comparison(B_list_gold_standard, B_list_majority_clean, tool_name = "hlamajority")
accuracy_B[length(tools)+1] <- all_results$final_accuracy
B_results[[length(B_results) + 1]] <- all_results
names(B_results)[length(B_results)] <- "hlamajority"

accuracy_C <- c()
C_results <- list()
for (i in 1:length(tools)){
  C_list_tool <- all.in %>% dplyr::filter(tool == tools[i]) %>% df_to_list_fix(., cols = c("C1", "C2"))
  C_list_tool_cleaned <- clean_list(C_list_tool, "C")
  all_results <- run_comparison(C_list_gold_standard, C_list_tool_cleaned, tool_name = tools[i])
  accuracy_C[i] <- all_results$final_accuracy
  names(accuracy_C)[i] <- tools[i]
  C_results[[length(C_results) + 1]] <- all_results
  names(C_results)[i] <- tools[i]
}
C_list_hlamajority <- df_to_list_fix(hlamajority.in, cols = c("C1", "C2"))
names(C_list_hlamajority) <- hlamajority.in$sample
C_list_majority_clean <- clean_list(C_list_hlamajority, "C")
all_results <- run_comparison(C_list_gold_standard, C_list_majority_clean, tool_name = "hlamajority")
accuracy_C[length(tools)+1] <- all_results$final_accuracy
C_results[[length(C_results) + 1]] <- all_results
names(C_results)[length(C_results)] <- "hlamajority"

results <- cbind.data.frame(accuracy_A, accuracy_B, accuracy_C)
colnames(results) <- c("A", "B", "C")
rownames(results)[nrow(results)] <- "hlamajority"
results

A_results$optitype$ids_incorrect_calls
names(A_results$optitype$incorrect_calls)
names(A_results$optitype$gold_standard_incorrect_calls)
test <- do.call(rbind.data.frame, A_results$optitype$gold_standard_incorrect_calls)
colnames(test) <- c("allele1", "allele2")
rownames(test) <- names(A_results$optitype$gold_standard_incorrect_calls)

# extract_mismatches <- function(incorrect_df) {
#   # Initialize empty lists to store results
#   missed_alleles <- c()
#   wrong_calls <- c()
#   sample_ids <- c()
#   
#   # Loop through each row of the incorrect dataframe
#   for (i in 1:nrow(incorrect_df)) {
#     # Get the row
#     row <- incorrect_df[i, ]
#     id <- rownames(incorrect_df)[i]
#     
#     # Create vectors for GS and Tool
#     gs <- c(row$allele1_gold_standard, row$allele2_gold_standard)
#     tool <- c(row$allele1_tool, row$allele2_tool)
#     
#     # Remove NAs just in case
#     gs <- gs[!is.na(gs)]
#     tool <- tool[!is.na(tool)]
#     
#     # LOGIC: 
#     # 1. Find matches
#     # 2. Remove matches from both lists
#     # 3. What remains in GS was "Missed"
#     # 4. What remains in Tool was the "Confusion"
#     
#     # We use a while loop to handle homozygous matching logic correctly
#     # (e.g. GS: 01, 01 vs Tool: 01, 02 -> Match one 01, Miss one 01)
#     
#     gs_remaining <- gs
#     tool_remaining <- tool
#     
#     for (g in gs) {
#       if (g %in% tool_remaining) {
#         # It's a match, remove it from tool_remaining so we don't match it twice
#         match_idx <- match(g, tool_remaining)
#         tool_remaining <- tool_remaining[-match_idx]
#         
#         # Remove from gs_remaining
#         match_gs_idx <- match(g, gs_remaining)
#         gs_remaining <- gs_remaining[-match_gs_idx]
#       }
#     }
#     
#     # Now gs_remaining contains alleles the tool missed
#     # tool_remaining contains what the tool guessed instead
#     
#     if (length(gs_remaining) > 0) {
#       missed_alleles <- c(missed_alleles, gs_remaining)
#       # If tool made a call but it was wrong, record it. 
#       # If tool returned fewer alleles (unlikely in your setup), fill NA
#       fill_tool <- c(tool_remaining, rep(NA, length(gs_remaining) - length(tool_remaining)))
#       wrong_calls <- c(wrong_calls, fill_tool)
#       sample_ids <- c(sample_ids, rep(id, length(gs_remaining)))
#     }
#   }
#   
#   return(data.frame(
#     Sample = sample_ids,
#     True_Allele = missed_alleles,
#     Predicted_Allele = wrong_calls,
#     stringsAsFactors = FALSE
#   ))
# }
# 
# analyze_difficult_alleles <- function(full_gold_standard, incorrect_df) {
#   
#   # 1. Get the specific mismatch pairs
#   mismatches <- extract_mismatches(incorrect_df)
#   
#   # 2. Calculate "Missed Counts" (Numerator)
#   missed_counts <- mismatches %>%
#     group_by(True_Allele) %>%
#     summarise(Errors = n())
#   
#   # 3. Calculate "Total Counts" in the Gold Standard (Denominator)
#   # We flatten the entire Gold Standard list into one big vector of alleles
#   all_gs_alleles <- unlist(full_gold_standard)
#   # Remove NAs
#   all_gs_alleles <- all_gs_alleles[!is.na(all_gs_alleles)]
#   
#   total_counts <- as.data.frame(table(all_gs_alleles))
#   colnames(total_counts) <- c("True_Allele", "Total_Occurrences")
#   
#   # 4. Merge
#   stats <- merge(total_counts, missed_counts, by = "True_Allele", all.x = TRUE)
#   
#   # Fill NAs with 0 (alleles that were never missed)
#   stats$Errors[is.na(stats$Errors)] <- 0
#   
#   # 5. Calculate Metrics
#   stats$Accuracy <- 1 - (stats$Errors / stats$Total_Occurrences)
#   stats$Error_Rate <- (stats$Errors / stats$Total_Occurrences)
#   
#   # 6. Sort by Error Rate (Desc) then Total Occurrences (Desc)
#   # We usually care about alleles that appear at least a few times (e.g. > 5)
#   stats <- stats %>%
#     arrange(desc(Error_Rate), desc(Total_Occurrences))
#   
#   return(stats)
# }

tool_name <- "kourami"
# Retrieve the object returned by run_comparison for Kourami
kourami_results <- A_results[[tool_name]] 

incorrect_df <- kourami_results$gold_standard_vs_tool_incorrect_calls

# 2. Analyze Difficulty
diff_stats <- analyze_difficult_alleles(A_list_gold_standard, incorrect_df)

# 3. View the "Most Difficult" alleles (min 5 occurrences to avoid noise)
top_difficult <- diff_stats %>% 
  filter(Total_Occurrences >= 5) %>%
  head(10)

print(top_difficult)

plot_confusion_heatmap <- function(incorrect_df, min_count = 2) {
  
  # Get specific pairs
  mismatches <- extract_mismatches(incorrect_df)
  
  # Count pairs
  pair_counts <- mismatches %>%
    group_by(True_Allele, Predicted_Allele) %>%
    summarise(Count = n(), .groups = "drop") %>%
    filter(Count >= min_count) # Filter out one-off errors to clean the plot
  
  ggplot(pair_counts, aes(x = Predicted_Allele, y = True_Allele, fill = Count)) +
    geom_tile() +
    scale_fill_gradient(low = "yellow", high = "red") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    labs(title = "Allele Confusion Matrix",
         subtitle = paste("Mistakes occuring >=", min_count, "times"),
         y = "True Allele (Gold Standard)",
         x = "Tool Prediction")
}

# Run plot
plot_confusion_heatmap(incorrect_df, min_count = 3)

if (!require("Biostrings")) install.packages("BiocManager"); BiocManager::install("Biostrings")
BiocManager::install("pwalign")

fasta_path <- "https://raw.githubusercontent.com/ANHIG/IMGTHLA/refs/heads/3630/hla_nuc.fasta"

hla_sequences <- readDNAStringSet(fasta_path)

# 2. Helper Function to Clean Names
# The FASTA headers look like: "HLA:HLA00001 A*01:01:01:01 1098 bp"
# We want to extract just "A*01:01:01:01"
names(hla_sequences) <- sapply(strsplit(names(hla_sequences), " "), `[`, 2)

get_allele_distance <- function(allele1, allele2, gene) {
  # Format input to match FASTA headers (e.g., "02:01" -> "A*02:01")
  # We assume inputs are 2-field, but FASTA is full resolution.
  # We pick the FIRST match (usually the reference/lowest number) for that 2-field group.
  
  a1_query <- paste0(gene, "*", allele1)
  a2_query <- paste0(gene, "*", allele2)
  # Find matching sequences (prefix matching)
  # grep finds "A*02:01:01:01" when searching for "A*02:01"
  idx1 <- grep(a1_query, names(hla_sequences))[1]
  idx2 <- grep(a2_query, names(hla_sequences))[1]
  if (is.na(idx1) || is.na(idx2)) return(NA)
  
  seq1 <- hla_sequences[[idx1]]
  seq2 <- hla_sequences[[idx2]]
  
  # Perform Pairwise Alignment
  # type="global" ensures we align the whole available sequence
  alignment <- pairwiseAlignment(seq1, seq2, type="global")
  
  # Calculate Edit Distance (Levenshtein)
  # This counts substitutions, insertions, and deletions
  dist <- nedit(alignment)
  
  # Calculate Percent Identity
  pid <- pid(alignment, type="PID1")
  
  return(list(distance = dist, percent_identity = pid))
}
result <- get_allele_distance("02:06", "02:07", "A")
print(paste("A*02:06 vs A*02:07 Distance:", result$distance, "bases"))

add_distance_to_errors <- function(incorrect_df, gene_name) {
  
  # 1. Handle input format
  # If the input is the raw matrix, convert it. If it's already converted, skip.
  if(!"True_Allele" %in% colnames(incorrect_df)) {
    incorrect_df <- extract_mismatches(incorrect_df)
  }
  
  distances <- numeric(nrow(incorrect_df))
  
  for(i in 1:nrow(incorrect_df)) {
    true_a <- incorrect_df$True_Allele[i]
    pred_a <- incorrect_df$Predicted_Allele[i]
    
    # 2. Skip if prediction or truth is NA
    if (is.na(pred_a) || is.na(true_a)) {
      distances[i] <- NA
      next
    }
    
    # 3. Calculate distance
    res <- get_allele_distance(true_a, pred_a, gene_name)
    
    # 4. SAFETY CHECK: Did we get a valid list back?
    if (is.list(res)) {
      distances[i] <- res$distance
      # Optional: Print progress for long loops
      # print(paste(true_a, "vs", pred_a, ":", res$distance))
    } else {
      # This handles the case where get_allele_distance returned NA
      # because the allele wasn't found in the FASTA
      warning(paste("Could not calculate distance for:", true_a, "vs", pred_a))
      distances[i] <- NA
    }
  }
  
  # 5. Assign and Return
  # Note: I fixed the variable name 'incorrect_df_cp' to 'incorrect_df'
  incorrect_df$Nucleotide_Distance <- distances
  return(incorrect_df)
}

# Run it
# (Assuming your incorrect_df is for HLA-A)
analyzed_errors <- add_distance_to_errors(incorrect_df, "A")

# View the results
head(analyzed_errors)

grep("A*02:01", names(hla_sequences))



# Define your tools and lists as before
tools_map <- list(
  "Optitype" = A_list_optitype_clean,
  "Polysolver" = A_list_polysolver_clean,
  "HLA-LA" = A_list_hlala_clean,
  "Kourami" = A_list_kourami_clean,
  "Majority" = A_list_majority_clean
)

all_errors_df <- data.frame()

for (t_name in names(tools_map)) {
  # Run the analysis
  df <- analyze_error_types(A_list_gold_standard, tools_map[[t_name]], t_name)
  all_errors_df <- rbind(all_errors_df, df)
}

# View summary table
table(all_errors_df$Tool, all_errors_df$Type)
