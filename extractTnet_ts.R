extractTnet_ts <- function(mdata_suffix, dirsave){
  
  source("./extractMdata_row.R")
  
  extractMdatalist_peak <- function(mdata_suffix, mdatalist, compTimeUnits, dirsave){
    
    mdataNum <- length(mdatalist)
    mdatalist_extracted <- vector("list", mdataNum)
    for(compTimeUnit in compTimeUnits){
      for(md in 1:mdataNum){
        mdatalist_extracted[[md]] <- extractMdata_row(mdatalist[[md]], annotCol_name = compTimeUnit, categ_targets = TRUE)
      }
      save(mdatalist_extracted, file = paste0(dirsave, "/MdataFiles/mdata_", sub("Insulin_compTP_maxmin2_TimeUnit_", "time", compTimeUnit), "_list"))
    }
    return()
  }
  #read data
  mdatalist <- base::get(load(paste0(dirsave, "/MdataFiles/mdata_ir_list")))
  compTimeUnits <- colnames(annotCols(mdatalist[[1]])[grepl(paste0("Insulin_compTP_maxmin2_TimeUnit_"), colnames(annotCols(mdatalist[[1]])))])
  extractMdatalist_peak(mdata_suffix, mdatalist, compTimeUnits, dirsave)
  
  return()
}
