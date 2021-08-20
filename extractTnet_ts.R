extractTnet_ts <- function(mdata_suffix, dirsave){
  #Extract subsets at each time point from multi-omic datasets.
  
  source("./extractMdata_row.R")
  
  extractMdatalist_peak <- function(mdata_suffix, mdatalist, compTimeUnits, dirsave){
    
    mdataNum <- length(mdatalist)
    mdatalist_extracted <- vector("list", mdataNum)
    for(compTimeUnit in compTimeUnits){
      for(md in 1:mdataNum){
        mdatalist_extracted[[md]] <- extractMdata_row(mdatalist[[md]], annotCol_name = compTimeUnit, categ_targets = TRUE)
      }
      save(mdatalist_extracted, file = paste0(dirsave, "/MdataFiles/mdata_", compTimeUnit, "_list"))
    }
    return()
  }
  #read data
  mdatalist <- base::get(load(paste0(dirsave, "/MdataFiles/mdata", mdata_suffix, "_list")))
  compTimeUnits <- colnames(annotCols(mdatalist[[1]])[grepl(paste0("^Time_\\d+$"), colnames(annotCols(mdatalist[[1]])))])
  extractMdatalist_peak(mdata_suffix, mdatalist, compTimeUnits, dirsave)
  
  return()
}
