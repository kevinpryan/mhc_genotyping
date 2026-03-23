library(tidyverse)

# HLA naming at the time of publication of Adams 2005
old_naming = read_delim(
  "downloads/HLA_nomenclature/Allelelist.2090.txt",
  delim = " ",
  col_names = c("AlleleID", "Allele"),
  col_types = cols(AlleleID = col_character(),
                   Allele = col_character())
)

new_naming = read_delim(
  "downloads/HLA_nomenclature/Allelelist.3440.txt",
  delim = ",",
  comment = "#",
  col_types = cols(
    AlleleID = col_character(),
    Allele = col_character()
  )
)

new_naming_3630 = read_delim(
  "downloads/HLA_nomenclature/Allelelist.3630.txt",
  delim = ",",
  comment = "#",
  col_types = cols(
    AlleleID = col_character(),
    Allele = col_character()
  )
)

make_allele_map <- function(old_naming, new_naming) {
  left_join(old_naming, new_naming, by="AlleleID",
            suffix=c(".old", ".new")) %>%
    select(Allele.old, Allele.new) %>%
    separate("Allele.old", c("gene.old", "allele.old"), sep="\\*") %>%
    separate("Allele.new", c("gene.new", "allele.new"), sep="\\*")
}

run_adams_pipeline <- function(allele_map_fullres,
                               allele_map_4d,
                               problematic_mapping_4d,
                               adams_html_glob = "downloads/pub/adams_2005/MHC*.html") {
  
  # ---- helper functions (identical logic to your script) ----
  
  read_html_table <- function(path) {
    xml2::read_html(path) %>%
      rvest::html_element(css = ".data") %>%
      rvest::html_table()
  }
  
  split_allele_col <- function(col) {
    colname <- colnames(col)
    gene <- stringr::str_remove(colname, regex(" Locus$", ignore_case = TRUE))
    
    if (gene == "DRβ1") gene <- "DRB1"
    
    tidyr::separate(col, colname, paste0("c", c(1, 2)), ", *", fill = "right") %>%
      dplyr::mutate(c2 = dplyr::if_else(is.na(c2), c1, c2)) %>%
      dplyr::rename_with(~ paste0(gene, ".", 1:2))
  }
  
  clean_allele_name <- function(col) {
    tidyr::gather(col, "key", "allele.old") %>%
      tidyr::separate("key", c("gene.old", "idx"), sep="\\.", remove = FALSE) %>%
      
      # enforce ≥4-digit resolution
      dplyr::mutate(
        allele.old = dplyr::if_else(
          stringr::str_detect(allele.old, "^([0-9]{4}).*"),
          allele.old,
          NA_character_
        )
      ) %>%
      
      # full-resolution mapping first
      dplyr::left_join(allele_map_fullres,
                       by = c("gene.old", "allele.old")) %>%
      dplyr::mutate(rownum = dplyr::row_number()) %>%
      dplyr::group_by(!is.na(allele.new)) %>%
      dplyr::group_modify(function(df, state) {
        
        if (!state) {
          
          df_trim <- df %>%
            dplyr::select(-dplyr::ends_with(".new")) %>%
            dplyr::mutate(
              allele.old = stringr::str_replace(
                allele.old, "^([0-9]{4}).*", "\\1"
              )
            )
          
          # hard fail on ambiguous 4-digit mappings
          stopifnot(
            nrow(
              dplyr::semi_join(
                df_trim,
                problematic_mapping_4d,
                by = c("gene.old", "allele.old")
              )
            ) == 0
          )
          
          dplyr::left_join(df_trim,
                           allele_map_4d,
                           by = c("gene.old", "allele.old"))
          
        } else {
          df
        }
      }) %>%
      dplyr::ungroup() %>%
      dplyr::select(-`!is.na(allele.new)`) %>%
      dplyr::arrange(rownum) %>%
      
      # collapse to 2-field resolution
      dplyr::mutate(
        allele.new = stringr::str_replace(
          allele.new, "^([0-9]+:[0-9]+).*", "\\1"
        )
      ) %>%
      dplyr::transmute(
        rownum,
        key = paste0(gene.new, ".", idx),
        value = allele.new
      ) %>%
      tidyr::spread(key, value) %>%
      dplyr::select(-rownum, -dplyr::starts_with("NA"))
  }
  
  # ---- parse Adams et al. tables ----
  
  tbl.orig <- Sys.glob(adams_html_glob) %>%
    tibble::tibble(path = .) %>%
    dplyr::mutate(html_obj = purrr::map(path, read_html_table)) %>%
    dplyr::pull(html_obj) %>%
    purrr::reduce(dplyr::full_join,
                  by = c("ID", "Cell Line", "Tissue")) %>%
    dplyr::transmute(
      dplyr::select(dplyr::cur_data(), !dplyr::ends_with("Locus")),
      alleles = purrr::lmap(
        dplyr::select(dplyr::cur_data(), dplyr::ends_with("Locus")),
        split_allele_col
      )
    )
  
  # ---- clean and remap allele names ----
  
  tbl <- tbl.orig %>%
    dplyr::mutate(purrr::lmap(alleles, clean_allele_name)) %>%
    dplyr::select(-alleles)
  
  return(tbl)
}
allele_map_fullres_3430 = make_allele_map(old_naming, new_naming)
allele_map_fullres_3630 = make_allele_map(old_naming, new_naming_3630)

list_not_one_to_one <- function(df) {
  df %>%
    distinct %>%
    group_by(gene.old, allele.old) %>%
    add_count() %>%
    filter(n != 1 | any(is.na(allele.new))) %>%
    distinct(gene.old, allele.old)
}

allele_map_4d_3430 = allele_map_fullres_3430 %>%
  mutate(allele.old=str_replace(allele.old, "([0-9][0-9])([0-9][0-9]).*", "\\1\\2"),
         allele.new=str_replace(allele.new, "([0-9]+):([0-9]+).*", "\\1:\\2")) %>%
  distinct %>%
  # Filter for the genes we need (MICA and MICB cause problems)
  filter(gene.old %in% c("A", "B", "Cw", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"))

# List for which alleles there is still no 1 -> 1 mapping at 4 digit level
problematic_mapping_4d_3430 = list_not_one_to_one(allele_map_4d_3430)

allele_map_4d_3630 = allele_map_fullres_3630 %>%
  mutate(allele.old=str_replace(allele.old, "([0-9][0-9])([0-9][0-9]).*", "\\1\\2"),
         allele.new=str_replace(allele.new, "([0-9]+):([0-9]+).*", "\\1:\\2")) %>%
  distinct %>%
  # Filter for the genes we need (MICA and MICB cause problems)
  filter(gene.old %in% c("A", "B", "Cw", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"))

# List for which alleles there is still no 1 -> 1 mapping at 4 digit level
problematic_mapping_4d_3630 = list_not_one_to_one(allele_map_4d_3630)


tbl_3430 <- run_adams_pipeline(
  allele_map_fullres = allele_map_fullres_3430,
  allele_map_4d      = allele_map_4d_3430,
  problematic_mapping_4d = problematic_mapping_4d_3430
)

tbl_3630 <- run_adams_pipeline(
  allele_map_fullres = allele_map_fullres_3630,
  allele_map_4d      = allele_map_4d_3630,
  problematic_mapping_4d = problematic_mapping_4d_3630
)

