extractTnet <- function(mdata_suffix_extracted, extract_categs = NULL, categ_targets_list = NULL, dirsave = NULL){
  #Extract data with specific values (categ_targets_list) in the specified 'annotCols' columns (extract_categs)
  #as well as the data which is known to relate to those specified data (e.g. TF, NetWorKIN, etc.)
  
  library(dplyr)
  library(PerseusR)
  
  source("./extractMdata_row_list.R")
  source("./extractMdata_row.R")

  mdatalist <- base::get(load(paste0(dirsave, "/MdataFiles/mdata_ir_list")))
  dataTypes <- c("phosphoproteome", "transcriptome", "proteome", "metabolome")
  annotCol_names_list <- list()
  annotCol_names_list[[1]] <- c(phosphoproteome = "Protein_AApos",
                                transcriptome = "Symbol",
                                proteome = "Majority.protein.IDs",
                                metabolome = "Compound.Name")
  annotCol_names_list[[2]] <- c(phosphoproteome = "Symbol",
                                transcriptome = "EnsemblGeneId",
                                proteome = "Symbol",
                                metabolome = "KeggCpdId")
  annotCol_targets_list <- priorbased_targets_list <- vector("list", length(dataTypes))
  names(annotCol_targets_list) <- names(priorbased_targets_list) <- dataTypes
  for(i in seq_along(mdatalist)){
    mdata <- mdatalist[[i]]
    
    #Nodes
    mdata_tmp1 <- extractMdata_row_list(mdata, 
                                        extract_categs = extract_categs, 
                                        categ_targets_list = categ_targets_list, splitCol = TRUE)
    if(dataTypes[i] == "transcriptome"){
      extract_categs_parents <- paste0("TF_", extract_categs)
    }else{
      extract_categs_parents <- paste0("NKINparent_", extract_categs)
    }
    categ_targets_list_parents <- categ_targets_list
    mdata_tmp2 <- extractMdata_row_list(mdata, 
                                        extract_categs = extract_categs_parents, 
                                        categ_targets_list = categ_targets_list_parents, splitCol = TRUE)
    
    annotCol_targets_list[[i]] <- unique(c(annotCols(mdata_tmp1)[[annotCol_names_list[[1]][i]]], 
                                           annotCols(mdata_tmp2)[[annotCol_names_list[[1]][i]]]))
    
    #Extract data based on upstream extracted datasets ("TF_xxx" and "NKINparent_xxx" cols)
    annotCols <- annotCols(mdata_tmp2)
    if(dataTypes[i] == "phosphoproteome"){
      #proteomeへ
      tmp <- unlist(annotCols[,colnames(annotCols)[colnames(annotCols) == "Proteins"]])
      if(length(tmp) > 0){
        priorbased_targets_list[[which(dataTypes == "proteome")]] <- unique(c(priorbased_targets_list[[which(dataTypes == "proteome")]],　unlist(strsplit(as.character(tmp), ";"))))
      }
    }else if(dataTypes[i] == "transcriptome"){
      tmp <- unlist(annotCols[,colnames(annotCols)[colnames(annotCols) %in% c("Protein.IDs", "TF_Protein.IDs")]])
      if(length(tmp) > 0){
        priorbased_targets_list[[which(dataTypes == "proteome")]] <- unique(c(priorbased_targets_list[[which(dataTypes == "proteome")]],　unlist(strsplit(as.character(tmp), ";"))))
      }
    }else if(dataTypes[i] == "proteome"){
      if("EC4" %in% colnames(annotCols)){
        tmp <- unlist(as.character(annotCols$EC4))
        if(length(tmp) > 0){
          priorbased_targets_list[[which(dataTypes == "metabolome")]] <- unique(c(priorbased_targets_list[[which(dataTypes == "metabolome")]],　unlist(strsplit(as.character(tmp), ";"))))
        }
      }
    }
  }
  for(i in seq_along(mdatalist)){
    mdata <- mdatalist[[i]]
    annotCols <- annotCols(mdata)
    if(dataTypes[i] == "proteome"){
      if(annotCol_names_list[[1]][[i]] %in% c("Majority.protein.IDs","Protein.IDs")){
        toAdd <- annotCols[, annotCol_names_list[[1]][[i]]]
        toAdd <- toAdd[sapply(toAdd, function(x) any(unlist(strsplit(as.character(x), ";")) %in% priorbased_targets_list[[i]]))]
        toAdd <- as.character(toAdd)
        annotCol_targets_list[[i]] <- unique(c(annotCol_targets_list[[i]], toAdd))
      }else{
        toAdd <- annotCols[, c(annotCol_names_list[[1]][[i]], "Majority.protein.IDs")]
        toAdd <- toAdd[sapply(toAdd$Majority.protein.IDs, function(x) any(unlist(strsplit(x, ";")) %in% priorbased_targets_list[[i]])),]
        toAdd <- unique(toAdd[, annotCol_names_list[[1]][[i]]])
        annotCol_targets_list[[i]] <- unique(c(annotCol_targets_list[[i]], toAdd))
      }
    }else if(dataTypes[i] == "metabolome"){next}
    
    annotCol_targets <- annotCol_targets_list[[i]]
    mdatalist[[i]] <- extractMdata_row(mdata, 
                                       annotCol_name = annotCol_names_list[[1]][[i]], 
                                       categ_targets = annotCol_targets)
    if(dataTypes[i] == "proteome"){
      annotCols <- mdatalist[[i]]@annotCols
      if("EC4" %in% colnames(annotCols)){
        tmp <- unlist(as.character(annotCols$EC4))
        if(length(tmp) > 0){
          priorbased_targets_list[[which(dataTypes == "metabolome")]] <- unique(c(priorbased_targets_list[[which(dataTypes == "metabolome")]],　unlist(strsplit(as.character(tmp), ";"))))
        }
      }
    }
  }
  
  metIndex <- which(dataTypes == "metabolome")
  mdata <- mdatalist[[metIndex]]
  annotCols <- annotCols(mdata)
  if(annotCol_names_list[[1]][[metIndex]] == "extEC4"){
    toAdd <- annotCols$extEC4
    toAdd <- toAdd[sapply(toAdd, function(x) any(unlist(strsplit(as.character(x), ";")) %in% priorbased_targets_list[[metIndex]]))]
    toAdd <- as.character(toAdd)
    annotCol_targets_list[[metIndex]] <- unique(c(annotCol_targets_list[[metIndex]], toAdd))
  }else{
    toAdd <- annotCols[, c(annotCol_names_list[[1]][[metIndex]], "extEC4")]
    toAdd <- toAdd[sapply(toAdd$extEC4, function(x) any(unlist(strsplit(as.character(x), ";")) %in% priorbased_targets_list[[metIndex]])),]
    toAdd <- unique(toAdd[, annotCol_names_list[[1]][[metIndex]]])
    annotCol_targets_list[[metIndex]] <- unique(c(annotCol_targets_list[[metIndex]], toAdd))
  }
  annotCol_targets <- annotCol_targets_list[[metIndex]]
  mdatalist[[metIndex]] <- extractMdata_row(mdata, 
                                            annotCol_name = annotCol_names_list[[1]][[metIndex]], 
                                            categ_targets = annotCol_targets)
  save(mdatalist, file = paste0(dirsave, "/MdataFiles/mdata", mdata_suffix_extracted, "_list"))
  
  return()
}
