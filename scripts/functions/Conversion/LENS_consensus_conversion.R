# Read in the output and store as R dataframe
#' @export
toolOutputToR.LENS_consensus <- function(outputFolder) {
  box::use(
    dplyr[...],
    stringr[...],
    tibble[...],
    readr[read_lines]
  )
  
  fileList <- list.files(
    path = outputFolder,
    pattern = "supported_hla_alleles$",
    full.names = TRUE
  )
  
  if (length(fileList) == 0) {
    warning("No LENS consensus files found in: ", outputFolder)
    return(NULL)
  }
  
  parse_locus <- function(alleles, locus) {
    hits <- alleles[stringr::str_detect(alleles, paste0("^", locus))]
    
    if (length(hits) == 0) {
      return(c(NA, NA))
    }
    
    # remove HLA- prefix and locus prefix (A, B, C)
    clean <- stringr::str_remove(hits, "^HLA-[A-Z]")
    clean <- stringr::str_remove(clean, "^[A-Z]")
    
    if (length(clean) == 1) {
      return(c(clean, clean))  # homozygous duplication
    }
    
    if (length(clean) >= 2) {
      return(clean[1:2])
    }
    
    c(NA, NA)
  }
  
  results <- lapply(fileList, function(f) {
    
    alleles <- readr::read_lines(f)
    
    sample_id <- stringr::str_remove(
      basename(f),
      "\\.supported_hla_alleles$"
    )
    if (length(alleles) == 0 || alleles == "") {
      # return all NA for this sample
      results <- data.frame(
        sample_id = sample_id,
        A1 = NA, A2 = NA,
        B1 = NA, B2 = NA,
        C1 = NA, C2 = NA,
        stringsAsFactors = FALSE
      )
      print(results)
      #rownames(results) <- NULL
      #results <- tibble::column_to_rownames(results, "sample_id")
      return(results)
    }
    
    alleles <- strsplit(alleles, ",")[[1]]
    
    alleles <- stringr::str_remove(alleles, "^HLA-")
    
    A <- parse_locus(alleles, "A")
    B <- parse_locus(alleles, "B")
    C <- parse_locus(alleles, "C")
    
    data.frame(
      sample_id = sample_id,
      A1 = A[1], A2 = A[2],
      B1 = B[1], B2 = B[2],
      C1 = C[1], C2 = C[2],
      stringsAsFactors = FALSE
    )
  })
  
  results <- dplyr::bind_rows(results)
  #results <- tibble::column_to_rownames(results, "sample_id")
  
  return(results)
}