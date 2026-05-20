library(vroom)
library(knitr)
file.exists("/hlamajority-paper/external/mhc_genotyping/scripts/functions/Conversion/arcasHLA_conversion_LENS.R")
module_path <- "/hlamajority-paper/external/mhc_genotyping/scripts/functions"
options(box.path = module_path)
box::use(Conversion/arcasHLA_conversion_LENS[
  toolOutputToR.arcasHLA_conversion_LENS
],
Conversion/seq2HLA_conversion_LENS[
  toolOutputToR.seq2HLA
],
Conversion/HLAprofiler_conversion_LENS[
  toolOutputToR.HLAprofiler
],
Conversion/Optitype_conversion_LENS[
  toolOutputToR.Optitype
],
Conversion/LENS_consensus_conversion[
  toolOutputToR.LENS_consensus
]
)
setwd("/hlamajority-paper/external/mhc_genotyping/")
source("scripts/functions/ggroup_mapper.R")
source("scripts/functions/evaluate_predictions_functions.R")
gold.standard.nci60 <- readRDS("data/gold_standard_nci60.rds")
colnames(gold.standard.nci60)[1] <- "sample"
colnames(gold.standard.nci60)[2:7] <- gsub(pattern = "\\.", replacement = "", colnames(gold.standard.nci60[2:7]))
gold.standard.nci60 <- as.data.frame(gold.standard.nci60)

sample_names_nci60 <- readRDS("data/sample_names_nci60_srx.rds") %>% rename(sample = "sample_id.srx")
base_dir <- "../../data/raw/cell-lines-after-polysolver-change/lens/v1.2-dev"
base_dir_1.8 <- "../../data/raw/cell-lines-after-polysolver-change/lens/v1.8"
# get all sample directories
samples <- list.dirs(base_dir, recursive = FALSE, full.names = TRUE)
samples_1.8 <- list.dirs(base_dir_1.8, recursive = FALSE, full.names = TRUE)


results_list <- lapply(samples, function(sample_dir) {
  
  # find the arcashla_genotype directory inside each sample
  genotype_dir <- list.dirs(sample_dir, recursive = TRUE, full.names = TRUE)
  genotype_dir <- genotype_dir[grepl("arcashla_genotype$", genotype_dir)]
  
  # skip if not found
  if (length(genotype_dir) == 0) return(NULL)
  #  print("dir...")
  #  print(genotype_dir[1])
  # run your function
  df <- toolOutputToR.arcasHLA_conversion_LENS(genotype_dir[1], trim = T)
  
  return(df)
})

# remove NULLs (samples that failed / missing)
results_list <- Filter(Negate(is.null), results_list)

# combine into one dataframe
arcasHLA_df <- do.call(rbind, results_list)
colnames(arcasHLA_df) <- c("A1", "A2", "B1", "B2", "C1", "C2")
srx_srr_ids <- stringr::str_extract(
  rownames(arcasHLA_df),
  "SRX[0-9]+_SRR[0-9]+"
)
cell_lines <- sub(".*_", "", rownames(arcasHLA_df))
cell_lines <- sub("^[^-]+-", "", cell_lines)
arcasHLA_df$sample <- cell_lines
arcasHLA_df$srx_srr <- srx_srr_ids

arcasHLA_df$tool <- "arcasHLA"
results_list_seq2HLA <- lapply(samples, function(sample_dir) {
  genotype_dir <- list.dirs(sample_dir, recursive = TRUE, full.names = TRUE)
  genotype_dir <- genotype_dir[grepl("seq2hla", genotype_dir)]
  
  # skip if not found
  if (length(genotype_dir) == 0) return(NULL)
  df <- toolOutputToR.seq2HLA(genotype_dir)
  
  return(df)
})
# remove NULLs (samples that failed / missing)
results_list_seq2HLA <- Filter(Negate(is.null), results_list_seq2HLA)

# combine into one dataframe
seq2HLA_df <- do.call(rbind, results_list_seq2HLA)
colnames(seq2HLA_df)[1:6] <- c("A1", "A2", "B1", "B2", "C1", "C2")
seq2HLA_df <- seq2HLA_df %>% dplyr::select("A1", "A2", "B1", "B2", "C1", "C2")
srx_srr_ids <- stringr::str_extract(
  rownames(seq2HLA_df),
  "SRX[0-9]+_SRR[0-9]+"
)
cell_lines <- sub(".*_", "", rownames(seq2HLA_df))
cell_lines <- sub("^[^-]+-", "", cell_lines)
seq2HLA_df$sample <- cell_lines
seq2HLA_df$srx_srr <- srx_srr_ids
seq2HLA_df$tool <- "seq2HLA"

# HLAprofiler
results_list_HLAprofiler <- lapply(samples, function(sample_dir) {
  genotype_dir <- list.dirs(sample_dir, recursive = TRUE, full.names = TRUE)
  genotype_dir <- genotype_dir[grepl("hlaprofiler_predict", genotype_dir)]
  
  # skip if not found
  if (length(genotype_dir) == 0) return(NULL)
  df <- toolOutputToR.HLAprofiler(genotype_dir, trim = T)
  
  return(df)
})
# remove NULLs (samples that failed / missing)
results_list_HLAprofiler <- Filter(Negate(is.null), results_list_HLAprofiler)

# combine into one dataframe
HLAprofiler_df <- do.call(rbind, results_list_HLAprofiler)
colnames(HLAprofiler_df)[1:6] <- c("A1", "A2", "B1", "B2", "C1", "C2")
HLAprofiler_df <- HLAprofiler_df %>% dplyr::select("A1", "A2", "B1", "B2", "C1", "C2", "A_ambiguous", "B_ambiguous", "C_ambiguous")
srx_srr_ids <- stringr::str_extract(
  rownames(HLAprofiler_df),
  "SRX[0-9]+_SRR[0-9]+"
)

cell_lines <- sub(".*_", "", rownames(HLAprofiler_df))
cell_lines <- sub("^[^-]+-", "", cell_lines)
HLAprofiler_df$sample <- cell_lines
HLAprofiler_df$srx_srr <- srx_srr_ids
HLAprofiler_df$tool <- "HLAprofiler"

# Optitype ar
results_list_Optitype_ar <- lapply(samples_1.8, function(sample_dir) {
  genotype_dir <- list.dirs(sample_dir, recursive = TRUE, full.names = TRUE)
  genotype_dir <- genotype_dir[grepl("optitype", genotype_dir)]
  genotype_dir <- genotype_dir[grepl("ar-", genotype_dir)]
  #print(genotype_dir)
  # skip if not found
  if (length(genotype_dir) == 0) return(NULL)
  df <- toolOutputToR.Optitype(genotype_dir)
  
  return(df)
})
# remove NULLs (samples that failed / missing)
results_list_Optitype_ar <- Filter(Negate(is.null), results_list_Optitype_ar)

# combine into one dataframe
Optitype_ar_df <- do.call(rbind, results_list_Optitype_ar)
Optitype_ar_df <- tibble::column_to_rownames(Optitype_ar_df, "sample_id")
colnames(Optitype_ar_df)[1:6] <- c("A1", "A2", "B1", "B2", "C1", "C2")
Optitype_ar_df <- Optitype_ar_df %>% dplyr::select("A1", "A2", "B1", "B2", "C1", "C2")
cell_lines <- stringr::str_extract(rownames(Optitype_ar_df), "(?<=-)[A-Za-z0-9-]+(?=\\.optitype_calls)")
srx_srr_ids <- stringr::str_extract(
  rownames(Optitype_ar_df),
  "SRX[0-9]+_SRR[0-9]+"
)
Optitype_ar_df$sample <- cell_lines
Optitype_ar_df$srx_srr <- srx_srr_ids

# manually read in 578T
Hs_578T_Optitype_ar <- toolOutputToR.Optitype("../../data/raw/cell-lines-after-polysolver-change/lens/v1.8/Hs_578T/ar-SRX3728595_SRR6755990-Hs_578T/optitype")
Hs_578T_Optitype_ar$sample_id <- "578T"
Hs_578T_Optitype_ar <- tibble::column_to_rownames(Hs_578T_Optitype_ar, "sample_id")
colnames(Hs_578T_Optitype_ar)[1:6] <- c("A1", "A2", "B1", "B2", "C1", "C2")
Hs_578T_Optitype_ar <- Hs_578T_Optitype_ar %>% dplyr::select("A1", "A2", "B1", "B2", "C1", "C2")
Hs_578T_Optitype_ar$sample <- "578T"
Hs_578T_Optitype_ar$srx_srr <- "SRX3728595_SRR6755990"
Optitype_ar_df <- rbind(Optitype_ar_df, Hs_578T_Optitype_ar)

# manually read in IMVI
IMVI_Optitype_ar <- toolOutputToR.Optitype("../../data/raw/cell-lines-after-polysolver-change/lens/v1.8/LOX_IMVI/ar-SRX3728612_SRR6755973-LOX_IMVI/optitype")
IMVI_Optitype_ar$sample_id <- "IMVI"
IMVI_Optitype_ar <- tibble::column_to_rownames(IMVI_Optitype_ar, "sample_id")
colnames(IMVI_Optitype_ar)[1:6] <- c("A1", "A2", "B1", "B2", "C1", "C2")
IMVI_Optitype_ar <- IMVI_Optitype_ar %>% dplyr::select("A1", "A2", "B1", "B2", "C1", "C2")
IMVI_Optitype_ar$sample <- "IMVI"
IMVI_Optitype_ar$srx_srr <- "SRX3728612_SRR6755973"
Optitype_ar_df <- rbind(Optitype_ar_df, IMVI_Optitype_ar)
Optitype_ar_df$tool <- "Optitype_ar"

# Optitype ad
results_list_Optitype_ad <- lapply(samples_1.8, function(sample_dir) {
  genotype_dir <- list.dirs(sample_dir, recursive = TRUE, full.names = TRUE)
  genotype_dir <- genotype_dir[grepl("optitype", genotype_dir)]
  genotype_dir <- genotype_dir[grepl("ad-", genotype_dir)]
  #print(genotype_dir)
  # skip if not found
  if (length(genotype_dir) == 0) return(NULL)
  df <- toolOutputToR.Optitype(genotype_dir)
  
  return(df)
})
# remove NULLs (samples that failed / missing)
results_list_Optitype_ad <- Filter(Negate(is.null), results_list_Optitype_ad)

# combine into one dataframe
Optitype_ad_df <- do.call(rbind, results_list_Optitype_ad)

Optitype_ad_df <- tibble::column_to_rownames(Optitype_ad_df, "sample_id")
colnames(Optitype_ad_df)[1:6] <- c("A1", "A2", "B1", "B2", "C1", "C2")
Optitype_ad_df <- Optitype_ad_df %>% dplyr::select("A1", "A2", "B1", "B2", "C1", "C2")
cell_lines <- stringr::str_extract(rownames(Optitype_ad_df), "(?<=-)[A-Za-z0-9-]+(?=\\.optitype_calls)")

srx_srr_ids <- stringr::str_extract(
  rownames(Optitype_ad_df),
  "SRX[0-9]+_SRR[0-9]+"
)
Optitype_ad_df$sample <- cell_lines
Optitype_ad_df$srx_srr <- srx_srr_ids

# manually read in 578T
Hs_578T_Optitype_ad <- toolOutputToR.Optitype("../../data/raw/cell-lines-after-polysolver-change/lens/v1.8/Hs_578T/ad-SRX4239543_SRR7366623-Hs_578T/optitype")
Hs_578T_Optitype_ad$sample_id <- "578T"
Hs_578T_Optitype_ad <- tibble::column_to_rownames(Hs_578T_Optitype_ad, "sample_id")
colnames(Hs_578T_Optitype_ad)[1:6] <- c("A1", "A2", "B1", "B2", "C1", "C2")
Hs_578T_Optitype_ad <- Hs_578T_Optitype_ad %>% dplyr::select("A1", "A2", "B1", "B2", "C1", "C2")
Hs_578T_Optitype_ad$sample <- "578T"
Hs_578T_Optitype_ad$srx_srr <- "SRX4239543_SRR7366623"
Optitype_ad_df <- rbind(Optitype_ad_df, Hs_578T_Optitype_ad)

# manually read in IMVI
IMVI_Optitype_ad <- toolOutputToR.Optitype("../../data/raw/cell-lines-after-polysolver-change/lens/v1.8/LOX_IMVI/ad-SRX4239573_SRR7366593-LOX_IMVI/optitype")
IMVI_Optitype_ad$sample_id <- "IMVI"
IMVI_Optitype_ad <- tibble::column_to_rownames(IMVI_Optitype_ad, "sample_id")
colnames(IMVI_Optitype_ad)[1:6] <- c("A1", "A2", "B1", "B2", "C1", "C2")
IMVI_Optitype_ad <- IMVI_Optitype_ad %>% dplyr::select("A1", "A2", "B1", "B2", "C1", "C2")
IMVI_Optitype_ad$sample <- "IMVI"
IMVI_Optitype_ad$srx_srr <- "SRX4239573_SRR7366593"
Optitype_ad_df <- rbind(Optitype_ad_df, IMVI_Optitype_ad)
Optitype_ad_df$tool <- "Optitype_ad"

Optitype_ad_df$tool <- "Optitype_ad"

# consensus HLA alleles
results_list_consensus_alleles <- lapply(samples_1.8, function(sample_dir) {
  genotype_dir <- list.dirs(sample_dir, recursive = TRUE, full.names = TRUE)
  genotype_dir <- genotype_dir[grepl("consensus_hla_alleles", genotype_dir)]
  #print(genotype_dir)
  # skip if not found
  if (length(genotype_dir) == 0) return(NULL)
  df <- toolOutputToR.LENS_consensus(genotype_dir)
  
  #return(df)
})

# combine into one dataframe
LENS_consensus_df <- do.call(rbind, results_list_consensus_alleles)
LENS_consensus_df <- tibble::column_to_rownames(LENS_consensus_df, "sample_id")
colnames(LENS_consensus_df)[1:6] <- c("A1", "A2", "B1", "B2", "C1", "C2")
LENS_consensus_df <- LENS_consensus_df %>% dplyr::select("A1", "A2", "B1", "B2", "C1", "C2")
cell_lines <- sub(".*v[0-9.]+-([A-Za-z0-9-]+)-.*", "\\1", rownames(LENS_consensus_df))
get_cell_line <- function(x) {
  x %>%
    stringr::str_remove("^cell-lines-lens-v[0-9.]+-") %>%
    stringr::str_split("-ad-", simplify = TRUE) %>%
    .[,1]
}
cell_lines <- get_cell_line(rownames(LENS_consensus_df))
srx_srr_ids <- stringr::str_extract(
  rownames(LENS_consensus_df),
  "SRX[0-9]+_SRR[0-9]+"
)
LENS_consensus_df$sample <- cell_lines
LENS_consensus_df$srx_srr <- srx_srr_ids
LENS_consensus_df$tool <- "LENS-v1.8-consensus"


standardise_schema <- function(df) {
  required_cols <- c(
    "A1","A2","B1","B2","C1","C2",
    "A_ambiguous","B_ambiguous","C_ambiguous",
    "sample", "srx_srr", "tool"
  )
  
  for (col in required_cols) {
    if (!col %in% names(df)) {
      df[[col]] <- NA
    }
  }
  
  df[, required_cols]
}

# combine all results

rbined_df <- rbind(
  standardise_schema(arcasHLA_df),
  standardise_schema(HLAprofiler_df),
  standardise_schema(seq2HLA_df),
  standardise_schema(Optitype_ar_df),
  standardise_schema(Optitype_ad_df),
  standardise_schema(LENS_consensus_df)
)

print("remove UO-31 - single-end")
rbined_df <- rbined_df %>% 
  dplyr::filter(sample != "UO-31")

# rename Hs_578T to 578T
rbined_df$sample <- gsub("Hs_578T", "578T", rbined_df$sample)
# rename LOX_IMVI to IMVI
rbined_df$sample <- gsub("LOX_IMVI", "IMVI", rbined_df$sample)
# ensure we have six rows for every sample
rows_per_sample <- rbined_df %>% group_by(sample) %>% summarise(n = n()) %>% arrange(n)
stopifnot(all(rows_per_sample$n == 6))
print("writing file...")
write.csv(rbined_df, file = "../../data/processed/cell-lines-after-polysolver-change/lens-majority-benchmark/combined-results/lens-v1.2-dev-v1.8-results-standardised.csv", quote = F, row.names = F)
