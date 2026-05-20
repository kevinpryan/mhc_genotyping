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
    # Claeys took the first one.
    final_x <- sapply(mapped_x, function(y) head(y, 1))
    
    return(final_x)
  })
}

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
  print(head(df))
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
  
  # FILTERING STEP
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
               proportion_excluded = proportion_excluded,
               gold_standard_list = gold_aligned
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

extract_mismatches <- function(incorrect_df) {
  # Initialize empty lists to store results
  missed_alleles <- c()
  wrong_calls <- c()
  sample_ids <- c()
  
  # Loop through each row of the incorrect dataframe
  for (i in 1:nrow(incorrect_df)) {
    # Get the row
    row <- incorrect_df[i, ]
    id <- rownames(incorrect_df)[i]
    
    # Create vectors for GS and Tool
    gs <- c(row$allele1_gold_standard, row$allele2_gold_standard)
    tool <- c(row$allele1_tool, row$allele2_tool)
    
    # Remove NAs just in case
    gs <- gs[!is.na(gs)]
    tool <- tool[!is.na(tool)]
    
    # LOGIC: 
    # 1. Find matches
    # 2. Remove matches from both lists
    # 3. What remains in GS was "Missed"
    # 4. What remains in Tool was the "Confusion"
    
    # We use a while loop to handle homozygous matching logic correctly
    # (e.g. GS: 01, 01 vs Tool: 01, 02 -> Match one 01, Miss one 01)
    
    gs_remaining <- gs
    tool_remaining <- tool
    
    for (g in gs) {
      if (g %in% tool_remaining) {
        # It's a match, remove it from tool_remaining so we don't match it twice
        match_idx <- match(g, tool_remaining)
        tool_remaining <- tool_remaining[-match_idx]
        
        # Remove from gs_remaining
        match_gs_idx <- match(g, gs_remaining)
        gs_remaining <- gs_remaining[-match_gs_idx]
      }
    }
    
    # Now gs_remaining contains alleles the tool missed
    # tool_remaining contains what the tool guessed instead
    
    if (length(gs_remaining) > 0) {
      missed_alleles <- c(missed_alleles, gs_remaining)
      # If tool made a call but it was wrong, record it. 
      # If tool returned fewer alleles (unlikely in your setup), fill NA
      fill_tool <- c(tool_remaining, rep(NA, length(gs_remaining) - length(tool_remaining)))
      wrong_calls <- c(wrong_calls, fill_tool)
      sample_ids <- c(sample_ids, rep(id, length(gs_remaining)))
    }
  }
  
  return(data.frame(
    Sample = sample_ids,
    True_Allele = missed_alleles,
    Predicted_Allele = wrong_calls,
    stringsAsFactors = FALSE
  ))
}

analyze_difficult_alleles <- function(full_gold_standard, incorrect_df) {
  
  # 1. Get the specific mismatch pairs
  mismatches <- extract_mismatches(incorrect_df)
  
  # 2. Calculate "Missed Counts" (Numerator)
  missed_counts <- mismatches %>%
    group_by(True_Allele) %>%
    summarise(Errors = n())
  
  # 3. Calculate "Total Counts" in the Gold Standard (Denominator)
  # We flatten the entire Gold Standard list into one big vector of alleles
  all_gs_alleles <- unlist(full_gold_standard)
  # Remove NAs
  all_gs_alleles <- all_gs_alleles[!is.na(all_gs_alleles)]
  
  total_counts <- as.data.frame(table(all_gs_alleles))
  colnames(total_counts) <- c("True_Allele", "Total_Occurrences")
  
  # 4. Merge
  stats <- merge(total_counts, missed_counts, by = "True_Allele", all.x = TRUE)
  
  # Fill NAs with 0 (alleles that were never missed)
  stats$Errors[is.na(stats$Errors)] <- 0
  
  # 5. Calculate Metrics
  stats$Accuracy <- 1 - (stats$Errors / stats$Total_Occurrences)
  stats$Error_Rate <- (stats$Errors / stats$Total_Occurrences)
  
  # 6. Sort by Error Rate (Desc) then Total Occurrences (Desc)
  # We usually care about alleles that appear at least a few times (e.g. > 5)
  stats <- stats %>%
    arrange(desc(Error_Rate), desc(Total_Occurrences))
  
  return(stats)
}

plot_confusion_heatmap <- function(incorrect_df, min_count = 2, gene, tool) {
  
  # Get specific pairs
  mismatches <- extract_mismatches(incorrect_df)
  
  # Count pairs
  pair_counts <- mismatches %>%
    group_by(True_Allele, Predicted_Allele) %>%
    summarise(Count = n(), .groups = "drop") %>%
    filter(Count >= min_count) # Filter out one-off errors to clean the plot
  
  ggplot(pair_counts, aes(x = Predicted_Allele, y = True_Allele, fill = Count)) +
    geom_tile(colour = "black") +
    scale_fill_gradient(low = "yellow", high = "red") +
    #theme_minimal() +
    theme_clean() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          legend.background = element_blank()) +
    labs(title = paste(tool, gene, "Allele Confusion Matrix", sep = " "),
         subtitle = paste("Mistakes occuring >=", min_count, "times"),
         y = "True Allele (Gold Standard)",
         x = "Tool Prediction")
}

get_allele_distance <- function(allele1, allele2, gene, hla_sequences) {
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

calculate_overall_stats <- function(summary_df) {
  
  # 1. Reconstruct raw counts from the percentages
  # We need raw numbers to sum them correctly across genes
  detailed_df <- summary_df %>%
    mutate(
      # Number of samples where the tool actually made a call
      Valid_Samples = Num_Samples - Num_Excluded_NA,
      
      # Total alleles evaluated (2 per sample)
      Total_Alleles_Evaluated = Valid_Samples * 2,
      
      # Number of alleles correct (Accuracy * Total Evaluated)
      # round() handles floating point noise
      Correct_Alleles = round(Accuracy * Total_Alleles_Evaluated)
    )
  
  # 2. Group by Tool and Sum across A, B, and C
  overall_stats <- detailed_df %>%
    group_by(Tool) %>%
    summarise(
      Gene = "Overall",
      
      # Aggregated Accuracy
      Total_Correct = sum(Correct_Alleles),
      Total_Evaluated = sum(Total_Alleles_Evaluated),
      Accuracy = Total_Correct / Total_Evaluated,
      
      # Aggregated Call Rate
      Total_Valid_Samples = sum(Valid_Samples),
      Total_Samples = sum(Num_Samples),
      Call_Rate = Total_Valid_Samples / Total_Samples,
      
      # Keep structure consistent
      Num_Samples = Total_Samples,
      Num_Excluded_NA = Total_Samples - Total_Valid_Samples
    ) %>%
    select(Gene, Tool, Accuracy, Call_Rate, Num_Samples, Num_Excluded_NA)
  
  # 3. Combine with original gene-specific rows
  final_df <- bind_rows(summary_df, overall_stats)
  
  return(final_df)
}

format_publication_table <- function(full_stats_df) {
  full_stats_df %>%
    mutate(
      # Format as percentage string
      Metric_Label = paste0(sprintf("%.1f", Accuracy * 100), "%")
    ) %>%
    select(Tool, Gene, Metric_Label) %>%
    pivot_wider(names_from = Gene, values_from = Metric_Label) %>%
    # Reorder columns to ensure A, B, C, Overall order
    select(Tool, A, B, C, Overall) %>%
    arrange(desc(Overall))
}

run_full_benchmark <- function(master_df, gold_standard, genes, tools) {
  
  # Initialize storage
  summary_stats <- data.frame() # For Accuracy/Call Rate table
  detailed_artifacts <- list()  # For specific error tables/heatmaps
  gold_standard_missing <- data.frame()   
  for (gene in genes) {
    print(paste("=== Processing HLA-", gene, " ===", sep=""))
    
    # 1. Prepare Gold Standard for this gene
    #    (Extract cols A1/A2, convert to list, Clean G-groups)
    cols <- c(paste0(gene, "1"), paste0(gene, "2"))
    gs_list_raw <- df_to_list_fix(gold_standard, cols = cols)
    print("dim(gold_standard)")
    print(dim(gold_standard))
    print("length(gs_list_raw)")
    print(length(gs_list_raw))
    gs_list_clean <- clean_list(gs_list_raw, gene)
    # Gold standard callable samples (non-NA alleles)
    gs_callable <- sapply(gs_list_clean, function(x) !any(is.na(x)))
    gs_callable_ids <- names(gs_callable)[gs_callable]
    n_gs_called <- length(gs_callable_ids)
    print(n_gs_called)
    gs_not_callable <- sapply(gs_list_clean, function(x) any(is.na(x)))
    sapply(gs_list_clean, function(x) if (any(is.na(x))){print(x)})
    gs_not_callable_ids <- names(gs_not_callable)[gs_not_callable]
    print(paste("non-callable gold standard:", gs_not_callable_ids, sep = " "))
    if (length(gs_not_callable_ids) > 0) {
      
      missing_df <- data.frame(
        Gene = gene,
        Sample = gs_not_callable_ids,
        GS_Allele1 = sapply(gs_list_clean[gs_not_callable_ids], function(x) x[1]),
        GS_Allele2 = sapply(gs_list_clean[gs_not_callable_ids], function(x) x[2]),
        stringsAsFactors = FALSE
      )
      
      gold_standard_missing <- rbind(gold_standard_missing, missing_df)
    }
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
      # n_total <- length(intersect(names(gs_list_clean), names(tool_list_clean)))
      # # Count how many in tool list are NOT NA
      # n_called <- sum(!sapply(tool_list_clean, function(x) any(is.na(x))))
      # call_rate <- n_called / n_total
      # Restrict to samples where GS made a call
      tool_sub <- tool_list_clean[gs_callable_ids]
      
      # Tool callable among GS-callable samples
      tool_callable <- sapply(tool_sub, function(x) !any(is.na(x)))
      
      n_total <- length(gs_callable_ids)        # gold-standard calls
      n_called <- sum(tool_callable)             # tool calls where GS is callable
      n_excluded_na <- n_total - n_called        # guaranteed ≥ 0
      
      call_rate <- n_called / n_total
      
      summary_stats <- rbind(summary_stats, data.frame(
        Gene = gene,
        Tool = tool,
        Accuracy = metrics$final_accuracy,
        Call_Rate = call_rate,
        Num_Samples = n_total,              # GS callable
        Num_GS_Called = n_gs_called,         # same as Num_Samples (explicit)
        Num_Tool_Called = n_called,
        Num_Excluded_NA = n_excluded_na
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
  
  return(list(summary = summary_stats, 
              details = detailed_artifacts,
              gold_standard_missing = gold_standard_missing
      )
  )
}

run_full_benchmark_concordance <- function(master_df, gold_standard, genes, tools) {
  
  # Initialize storage
  summary_stats <- data.frame() # For Accuracy/Call Rate table
  detailed_artifacts <- list()  # For specific error tables/heatmaps
  gold_standard_missing <- data.frame()   
  concordance_stats <- data.frame()
  all_discordant <- list()
  for (gene in genes) {
    print(paste("=== Processing HLA-", gene, " ===", sep=""))
    
    # 1. Prepare Gold Standard for this gene
    #    (Extract cols A1/A2, convert to list, Clean G-groups)
    cols <- c(paste0(gene, "1"), paste0(gene, "2"))
    gs_list_raw <- df_to_list_fix(gold_standard, cols = cols)
    print("dim(gold_standard)")
    print(dim(gold_standard))
    print("length(gs_list_raw)")
    print(length(gs_list_raw))
    gs_list_clean <- clean_list(gs_list_raw, gene)
    # Gold standard callable samples (non-NA alleles)
    gs_callable <- sapply(gs_list_clean, function(x) !any(is.na(x)))
    gs_callable_ids <- names(gs_callable)[gs_callable]
    n_gs_called <- length(gs_callable_ids)
    print(n_gs_called)
    gs_not_callable <- sapply(gs_list_clean, function(x) any(is.na(x)))
    sapply(gs_list_clean, function(x) if (any(is.na(x))){print(x)})
    gs_not_callable_ids <- names(gs_not_callable)[gs_not_callable]
    print(paste("non-callable gold standard:", gs_not_callable_ids, sep = " "))
    if (length(gs_not_callable_ids) > 0) {
      
      missing_df <- data.frame(
        Gene = gene,
        Sample = gs_not_callable_ids,
        GS_Allele1 = sapply(gs_list_clean[gs_not_callable_ids], function(x) x[1]),
        GS_Allele2 = sapply(gs_list_clean[gs_not_callable_ids], function(x) x[2]),
        stringsAsFactors = FALSE
      )
      
      gold_standard_missing <- rbind(gold_standard_missing, missing_df)
    }
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
      tool_sub <- tool_list_clean[gs_callable_ids]
      
      # Tool callable among GS-callable samples
      tool_callable <- sapply(tool_sub, function(x) !any(is.na(x)))
      
      n_total <- length(gs_callable_ids)        # gold-standard calls
      n_called <- sum(tool_callable)             # tool calls where GS is callable
      n_excluded_na <- n_total - n_called        # guaranteed ≥ 0
      
      call_rate <- n_called / n_total
      
      summary_stats <- rbind(summary_stats, data.frame(
        Gene = gene,
        Tool = tool,
        Accuracy = metrics$final_accuracy,
        Call_Rate = call_rate,
        Num_Samples = n_total,              # GS callable
        Num_GS_Called = n_gs_called,         # same as Num_Samples (explicit)
        Num_Tool_Called = n_called,
        Num_Excluded_NA = n_excluded_na
      ))
      
      # 7. Store Detailed Artifacts (nested list)
      # Structure: List$Gene$Tool$Object
      detailed_artifacts[[gene]][[tool]] <- list(
        metrics = metrics,
        error_types = error_types,
        difficult_alleles = diff_alleles
      )
    }
    tool_gene_data <- lapply(tools, function(tool) {
      df <- master_df %>% filter(tool == !!tool)
      tool_list_raw <- df_to_list_fix(df, cols = cols)
      clean_list(tool_list_raw, gene)
    })
    names(tool_gene_data) <- tools
    
    tool_pairs <- combn(tools, 2, simplify = FALSE)
    
    for (pair in tool_pairs) {
      t1 <- pair[1]
      t2 <- pair[2]
      
      x <- tool_gene_data[[t1]]
      y <- tool_gene_data[[t2]]
      
      common_ids <- intersect(names(x), names(y))
      
      if (length(common_ids) == 0) next
      
      x_sub <- x[common_ids]
      y_sub <- y[common_ids]
      
      both_called <- !sapply(x_sub, function(z) any(is.na(z))) &
        !sapply(y_sub, function(z) any(is.na(z)))
      
      usable_ids <- common_ids[both_called]
      
      if (length(usable_ids) == 0) next
      
      # agree <- sapply(usable_ids, function(id) {
      #   identical(sort(x[[id]]), sort(y[[id]]))
      # })
      discordant <- list()
      agree <- sapply(usable_ids, function(id) {
        
        a <- x[[id]]
        b <- y[[id]]
        
        a_sorted <- sort(a)
        b_sorted <- sort(b)
        
        is_match <- identical(a_sorted, b_sorted)
        
        # Debug print (only when mismatch or always if you prefer)
        # if (!is_match) {
        # 
        #   message("Mismatch in gene: ", gene,
        #           " | ID: ", id,
        #           " | Tool1: ", paste(a_sorted, collapse = ","), 
        #           " | Tool2: ", paste(b_sorted, collapse = ","))
        #   
        #   discordant[[length(discordant) + 1]] <<- data.frame(
        #     Gene = gene,
        #     Tool1 = t1,
        #     Tool2 = t2,
        #     Sample = id,
        #     Tool1_call = paste(a, collapse = ","),
        #     Tool2_call = paste(b, collapse = ","),
        #     stringsAsFactors = FALSE
        #   )
        # }
        if (!is_match) {
          
          gs_call <- gs_list_clean[[id]]
          
          discordant[[length(discordant) + 1]] <<- data.frame(
            Gene = gene,
            Tool1 = t1,
            Tool2 = t2,
            Sample = id,
            
            Tool1_call = paste(a_sorted, collapse = ","),
            Tool2_call = paste(b_sorted, collapse = ","),
            
            GS_call = paste(sort(gs_call), collapse = ","),
            
            Tool1_vs_GS_match = identical(a_sorted, sort(gs_call)),
            Tool2_vs_GS_match = identical(b_sorted, sort(gs_call)),
            
            Tool1_overlap = length(intersect(a_sorted, gs_call)) / 2,
            Tool2_overlap = length(intersect(b_sorted, gs_call)) / 2,
            
            stringsAsFactors = FALSE
          )
        }
        # } else {
        #   message("Match | Gene: ", gene,
        #           " | ID: ", id)
        # }
        
        is_match
      })
      
      # if (length(discordant) > 0) {
        discordant_df <- do.call(rbind, discordant)
        all_discordant[[length(all_discordant) + 1]] <- discordant_df

      #   
      #   message("Top discordant samples for ", gene, " | ", t1, " vs ", t2)
      #   
      #   print(
      #     head(
      #       discordant_df[order(discordant_df$Sample), ],
      #       10
      #     )
      #   )
      # }
      
      concordance_stats <- rbind(concordance_stats, data.frame(
        Gene = gene,
        Tool1 = t1,
        Tool2 = t2,
        Concordance = sum(agree) / length(agree),
        N = length(agree)
      ))
      
    }
  }
  all_discordant_df <- do.call(rbind, all_discordant)
  
  tool_long <- bind_rows(
    concordance_stats %>% select(Gene, Tool = Tool1, Tool2, Concordance),
    concordance_stats %>% select(Gene, Tool = Tool2, Tool2 = Tool1, Concordance)
  ) %>% dplyr::filter(Tool != "Optitype_ad")
  
  tool_gene_discordance <- tool_long %>%
    group_by(Gene, Tool) %>%
    summarise(
      mean_concordance = mean(Concordance),
      discordance = 1 - mean_concordance,
      n_partners = n(),
      .groups = "drop"
    )
  
  # all_discordant_df %>%
  #   dplyr::count(Sample, sort = TRUE)
  return(list(summary = summary_stats, 
              concordance = concordance_stats,
              details = detailed_artifacts,
              gold_standard_missing = gold_standard_missing,
              discordance = all_discordant_df,
              tool_gene_discordance = tool_gene_discordance
  )
  )
}
