runTnetFun <- function(){
  
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
  source("./extractMdata_row_list.R")
  
  dirsave <- "."
  
  #generate the whole trans-omic network
  generateTnet(mdata_suffix = "_ir",
               arrowRelThick = c(100, 120),
               dirsave = dirsave, 
               prefix = "1_", suffix = "_Whl")
  
  #extract subnetworks from the whole trans-omic network
  extractTnet(mdata_suffix_extracted = "_InsSig",
              extract_categs = c("KeggpathName"), dirsave = dirsave,
              categ_targets_list = list(c("Insulin signaling pathway")))
  extractTnet(mdata_suffix_extracted = "_GenInfo",
              extract_categs = c("KeggClass"), dirsave = dirsave,
              categ_targets_list = list(c("Genetic Information Processing, Transcription",
                                          "Genetic Information Processing, Translation",
                                          "Genetic Information Processing, Folding, sorting and degradation",
                                          "Genetic Information Processing, Replication and repair")))
  extractTnet(mdata_suffix_extracted = "_Metab",
              extract_categs = c("KeggClass"), dirsave = dirsave,
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
                                          "Metabolism, Xenobiotics biodegradation and metabolism")))
  
  #generate the sub-networks
  arrowRelThick <- c(400, 400)
  generateTnet(mdata_suffix = "_InsSig", 
               arrowRelThick = arrowRelThick, dirsave = dirsave,
               prefix = "2_", suffix = "_InsSig")
  generateTnet(mdata_suffix = "_GenInfo",
               arrowRelThick = arrowRelThick, dirsave = dirsave,
               prefix = "3_", suffix = "_GenInfo")
  
  generateTnet(mdata_suffix = "_Metab", 
               arrowRelThick = arrowRelThick, dirsave = dirsave, 
               prefix = "4_", suffix = "_Metab")
  
  #generate time-series of the trans-omic network
  extractTnet_ts(mdata_suffix = "_ir", dirsave = dirsave)
  generateTnet_ts(arrowRelThick = c(120, 200), dirsave = dirsave, prefix = "5_")
  
  return()
}
