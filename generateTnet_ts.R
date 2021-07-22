generateTnet_ts <- function(mdata_suffix, arrowRelThick = c(100, 100), dirsave, prefix = NULL, suffix = NULL){

  source("./generateTnet.R")
  
  exeGenerateTnet <- function(mdatalist, compTimeUnits, arrowRelThick, dirsave, prefix, suffix){
    for(compTimeUnit in compTimeUnits){
      mdata_suffix_extracted <- paste0("_", sub("Insulin_compTP_maxmin2_TimeUnit_", "time", compTimeUnit))
      mdatalist2 <- base::get(load(paste0(dirsave, "/MdataFiles/mdata", mdata_suffix_extracted, "_list")))
      generateTnet(mdatalist2,
                mdata_suffix = mdata_suffix_extracted,
                arrowRelThick = arrowRelThick,
                dirsave = dirsave,
                prefix = prefix,
                suffix = mdata_suffix_extracted)
    }
    return()
  }
  
  #read data
  mdatalist <- base::get(load(paste0(dirsave, "/MdataFiles/mdata_ir_list")))
  compTimeUnits <- colnames(annotCols(mdatalist[[1]])[grepl(paste0("Insulin_compTP_maxmin2_TimeUnit_"), colnames(annotCols(mdatalist[[1]])))])
  exeGenerateTnet(mdatalist, compTimeUnits, arrowRelThick, dirsave, prefix, suffix)
  
  return()
}