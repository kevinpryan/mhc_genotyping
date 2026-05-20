box::use(
  rjson[fromJSON],
  stringr[str_replace]
)

toolOutputToR.arcasHLA_conversion_LENS <- function(outputFolder, trim = T){
  # Load needed packages
  
  # Get a list of .json file names in the given outputFolder
  fileList <- list.files(outputFolder, pattern = "\\.genotype.json$", full.names = T)
  # Define the expected gene types
  geneTypes <- c("A","B","C")
  
  # Create a data frame to store the results in
  results <- data.frame(matrix(NA, nrow = length(fileList), ncol = length(geneTypes)*2))
  names(results) <- rep(geneTypes, each = 2)
  IDs <- c()
  
  # For every file in the given folder, extract and store sample ID and genotype output
  for(i in 1:length(fileList)){
    # Load in file as R object
    alleleList <- fromJSON(file = fileList[i])
    alleleList <- alleleList[geneTypes]
    alleleList[setdiff(geneTypes, names(alleleList))] <- list(NA)
    # Extract sample ID
    #IDs[i] <- basename(fileList[i]) %>% str_replace('\\.genotype.json$', '')
    #IDs[i] <- str_replace(
    #  basename(fileList[i]),
    #  "\\.genotype.json$",
    #  ""
    #)
    IDs[i] <- basename(dirname(dirname(fileList[i])))    
    # For every allele in the list, store the type of gene and allele result
    for (gene in geneTypes) {
      
      alleles <- alleleList[[gene]]
      
      colIndex <- which(names(results) == gene)
      
      # missing gene
      if (is.na(alleles)[1]) {
        results[i, colIndex[1]] <- NA
        results[i, colIndex[2]] <- NA
        next
      }
      
      # single allele
      if (length(alleles) == 1) {
        parts <- strsplit(alleles, "[*]")[[1]]
        if (trim == TRUE){
        trimmed_genotype <- paste(strsplit(parts[2], split = ":")[[1]][1:2], collapse = ":")
        trimmed_genotype <- gsub("[NPQ]$", "", trimmed_genotype)
        } else {
          trimmed_genotype <- parts[2]
        }
        results[i, colIndex[1]] <- trimmed_genotype
        results[i, colIndex[2]] <- trimmed_genotype
        next
      }
      
      # two alleles
      parts1 <- strsplit(alleles[1], "[*]")[[1]]
      if (trim == TRUE){
      trimmed_genotype_1 <-  paste(strsplit(parts1[2], split = ":")[[1]][1:2], collapse = ":")
      trimmed_genotype_1 <-  gsub("[NPQ]$", "", trimmed_genotype_1)
      parts2 <- strsplit(alleles[2], "[*]")[[1]]
      trimmed_genotype_2 <-  paste(strsplit(parts2[2], split = ":")[[1]][1:2], collapse = ":")
      trimmed_genotype_2 <-  gsub("[NPQ]$", "", trimmed_genotype_2)
      } else {
        trimmed_genotype_1 <- parts1[2]
        trimmed_genotype_2 <- parts2[2]
      }
      results[i, colIndex[1]] <- trimmed_genotype_1
      results[i, colIndex[2]] <- trimmed_genotype_2
    }
  }
  # Add the IDs as rownames to the data frame
  row.names(results) <- IDs
  
  return(results)
}
box::export(toolOutputToR.arcasHLA_conversion_LENS)