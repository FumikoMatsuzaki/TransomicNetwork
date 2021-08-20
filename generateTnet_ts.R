generateTnet_ts <- function(arrowRelThick = c(100, 100), dirsave, prefix = NULL, settings){
  #Generate a time series of trans-omics network using subsets of multi-omic datasets extracted at each time point.
  
  source("./generateTnet.R")
  
  exeGenerateTnet <- function(compTimeUnits, arrowRelThick, dirsave, prefix, settings){
    for(compTimeUnit in compTimeUnits){
      mdata_suffix_extracted <- paste0("_", compTimeUnit)
      generateTnet(mdata_suffix = mdata_suffix_extracted,
                   arrowRelThick = arrowRelThick,
                   dirsave = dirsave,
                   prefix = prefix,
                   suffix = mdata_suffix_extracted,
                   settings = settings)
    }
    return()
  }
  
  mdatalist <- base::get(load(paste0(dirsave, "/MdataFiles/mdata_ir_list")))
  compTimeUnits <- colnames(annotCols(mdatalist[[1]])[grepl("^Time_\\d+$", colnames(annotCols(mdatalist[[1]])))])
  exeGenerateTnet(compTimeUnits, arrowRelThick, dirsave, prefix, settings)
  
  return()
}
