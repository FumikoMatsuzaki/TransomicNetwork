extractMdata_row <- function(mdata, annotCol_name, categ_targets){

  library(PerseusR)
  library(dplyr)
  
  main <- main(mdata)
  annotCols <- annotCols(mdata)
  annotRows <- annotRows(mdata)
  
  if(!annotCol_name %in% colnames(annotCols)){
    main_extracted <- main[NULL,]
    annotCols_extracted <- annotCols[NULL,]
    mdata_extracted <- matrixData(main = main_extracted, annotCols = annotCols_extracted, annotRows = annotRows)
    return(mdata_extracted)
  }
  
  factor_flag <- 0
  if(class(annotCols[,annotCol_name]) == "factor"){
    factor_flag <- 1
    annotCols[,annotCol_name] <- as.character(annotCols[,annotCol_name])
  }
  row_log <- annotCols[,annotCol_name] %in% categ_targets
  
  if(factor_flag == 1){
    annotCols[,annotCol_name] <- as.factor(annotCols[,annotCol_name])
  }
  annotCols_extracted <- annotCols[row_log,]
  main_extracted <- main[row_log,]
  mdata_extracted <- matrixData(main = main_extracted, annotCols = annotCols_extracted, annotRows = annotRows)
  
  return(mdata_extracted)
}
