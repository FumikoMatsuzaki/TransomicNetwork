tnetMainfun <- function(){
  
  library(dplyr)
  library(tidyr)
  library(PerseusR)
  library(grid)
  library(openxlsx)
  
  source("./generateTnet.R")
  source("./extractTnet.R")
  source("./extractTnet_ts.R")
  source("./generateTnet_ts.R")
  source("./extractMdata_row.R")
  
  #Settings
  dirsave <- "."
  load(paste0(dirsave, "/settings"))
  
  #Generate the whole trans-omic network
  generateTnet(mdata_suffix = "_ir", arrowRelThick = c(100, 120), dirsave, prefix = "1_", suffix = "_Whl", settings)
  
  #Extract sub-networks from the whole trans-omic network
  extractTnet(mdata_suffix_extracted = "_InsSig", extract_categs = "KeggpathName", dirsave, 
              categ_targets_list = list("Insulin signaling pathway"), settings)
  extractTnet(mdata_suffix_extracted = "_GenInfo", extract_categs = "KeggClass", dirsave,
              categ_targets_list = list(c("Genetic Information Processing, Transcription",
                                          "Genetic Information Processing, Translation",
                                          "Genetic Information Processing, Folding, sorting and degradation",
                                          "Genetic Information Processing, Replication and repair")), settings)
  extractTnet(mdata_suffix_extracted = "_Metab", extract_categs = "KeggClass", dirsave,
              categ_targets_list = list(c("Metabolism, Carbohydrate metabolism",
                                          "Metabolism, Energy metabolism",
                                          "Metabolism, Lipid metabolism",
                                          "Metabolism, Nucleotide metabolism",
                                          "Metabolism, Amino acid metabolism",
                                          "Metabolism, Metabolism of other amino acids",
                                          "Metabolism, Glycan biosynthesis and metabolism",
                                          "Metabolism, Metabolism of cofactors and vitamins",
                                          "Metabolism, Metabolism of terpenoids and polyketides",
                                          "Metabolism, Biosynthesis of other secondary metabolites",
                                          "Metabolism, Xenobiotics biodegradation and metabolism")), settings)
  arrowRelThick <- c(400, 400)
  generateTnet(mdata_suffix = "_InsSig", arrowRelThick, dirsave, prefix = "2_", suffix = "_InsSig", settings)
  generateTnet(mdata_suffix = "_GenInfo", arrowRelThick, dirsave, prefix = "3_", suffix = "_GenInfo", settings)
  generateTnet(mdata_suffix = "_Metab", arrowRelThick, dirsave, prefix = "4_", suffix = "_Metab", settings)
  
  #Generate time-series of the trans-omic network
  extractTnet_ts(mdata_suffix = "_ir", dirsave)
  generateTnet_ts(arrowRelThick = c(120, 200), dirsave, prefix = "5_", settings)
  
  return()
}
