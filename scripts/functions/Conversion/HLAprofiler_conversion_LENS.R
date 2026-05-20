# Read in the output and store as R dataframe
#' @export
toolOutputToR.HLAprofiler <- function(outputFolder, trim = FALSE) {
        box::use(
          dplyr[...],
          vroom[...],
          stringr[...],
          tidyr[...],
          tibble[...]
        )
        
        fileList <- list.files(
          path = outputFolder,
          pattern = "\\.HLATypes\\.txt$",
          full.names = TRUE
        )
        
        res <- vroom::vroom(
          fileList,
          id = "filename",
          show_col_types = FALSE,
          .name_repair = "minimal"
        )
        
        # Keep only A/B/C rows
        res <- res %>%
          dplyr::filter(stringr::str_detect(Allele1, "^[ABC]\\*"))
        
        # Extract sample ID
        res <- res %>%
          dplyr::mutate(
            sample_id = stringr::str_remove(basename(filename), "\\.HLATypes\\.txt$")
          )
        
        # Extract gene + clean allele
        res <- res %>%
          dplyr::mutate(
            gene = stringr::str_extract(Allele1, "^[A-Z]+"),
            Allele1 = stringr::str_remove(Allele1, "^[^*]+\\*"),
            Allele2 = stringr::str_remove(Allele2, "^[^*]+\\*")
          )
        
        clean_allele <- function(x) {
          x %>%
            stringr::str_remove("_[^:]+$")  # removes _novel, _updated, etc.
        }
        
        res <- res %>%
          dplyr::mutate(
            Allele1 = clean_allele(Allele1),
            Allele2 = clean_allele(Allele2)
          )
        
        # Trim to 2-field if requested
        if (trim) {
          res <- res %>%
            dplyr::mutate(
              dplyr::across(
                c(Allele1, Allele2),
                ~ stringr::str_extract(.x, "^[^:]+:[^:]+")
              )
            )
        }
        ambig <- res %>%
          dplyr::group_by(sample_id, gene) %>%
          dplyr::summarise(
            multi_call_ambiguous = dplyr::n() > 1,
            .groups = "drop"
          )
        res_best <- res %>%
          dplyr::group_by(sample_id, gene) %>%
          dplyr::slice_max(Pair_score, n = 1, with_ties = FALSE) %>%
          dplyr::ungroup()
        
        res_wide <- res_best %>%
          tidyr::pivot_longer(
            cols = c(Allele1, Allele2),
            names_to = "copy",
            values_to = "allele"
          ) %>%
          dplyr::mutate(copy = ifelse(copy == "Allele1", 1, 2)) %>%
          tidyr::pivot_wider(
            id_cols = sample_id,
            names_from = c(gene, copy),
            values_from = allele
          )
        
        ambig_wide <- ambig %>%
          tidyr::pivot_wider(
            names_from = gene,
            values_from = multi_call_ambiguous,
            names_glue = "{gene}_ambiguous"
          )
        
        res <- dplyr::left_join(res_wide, ambig_wide, by = "sample_id")
        print(res)
        # Ensure all columns exist
        expected_cols <- c("A_1","A_2","B_1","B_2","C_1","C_2", "A_ambiguous", "B_ambiguous", "C_ambiguous")
        
        for (col in expected_cols) {
          if (!col %in% colnames(res)) {
            res[[col]] <- NA
          }
        }
        
        res <- res[, c("sample_id", expected_cols)]
        
        # Rename 
        colnames(res)[1:7] <- c("sample_id", rep(c("A","B","C"), each = 2))
        
        res <- tibble::column_to_rownames(res, "sample_id")
        
        if (nrow(res) == 0) {
                sample_id <- stringr::str_remove(
                  basename(fileList[1]),
                  "\\.HLATypes\\.txt$"
                )
                
                empty <- data.frame(
                  sample_id = sample_id,
                  A_1 = NA, A_2 = NA,
                  B_1 = NA, B_2 = NA,
                  C_1 = NA, C_2 = NA,
                  A_ambiguous = NA,
                  B_ambiguous = NA,
                  C_ambiguous = NA,
                  stringsAsFactors = FALSE
                )
                rownames(empty) <- empty$sample_id
                empty$sample_id <- NULL
                colnames(empty)[1:6] <- rep(c("A","B","C"), each = 2)
                return(empty)
        }
        return(res)
}