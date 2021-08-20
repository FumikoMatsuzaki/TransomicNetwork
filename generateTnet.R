generateTnet <- function(mdata_suffix, arrowRelThick = c(100, 100), dirsave, prefix = NULL, suffix = NULL, settings){
  #Generate trans-omics network ("tnet") using multi-omic datasets retained in matrixData objects within PerseusR.
  #arrowRelThick: To scale up or down the thicknesses of arrows from the default setting (100%) in the original and simplified diagrams, respectively.
  
  library(dplyr)
  library(tidyr)
  library(PerseusR)
  library(grid)
  library(openxlsx)
  
  mainFun <- function(mdatalist, mdata_suffix, arrowRelThick, dirsave, prefix, suffix, settings){
    #Settings
    dataTypes <- settings[["dataTypes"]]
    annotCol_names_list <- settings[["annotCol_names_list"]]
    sepGenes_vec_list <- settings[["sepGenes_vec_list"]]
    annotCol_names_tnet <- c(annotCol_names_list[[1]], 
                             phosphoproteome_protein = "Proteins")
    keycols_tnet<- c(ph1st = "NKINparent_desc", tr3st = "TF", tr4st = "miRNAname", pr5st = "miRNAname", 
                     tr6stLim = "TCreg_TrPrMatch_match", tr6stKey = "Protein.IDs", pr6enLim = "TCreg_TrPrMatch_match", 
                     ph7stLim = "timeVar_fcRatio_aucRatio_proteome_match", pr7enLim = "timeVar_fcRatio_aucRatio_phosphoproteome_match", 
                     ph8stLim = "extEC4_metabolome_match", ph8stKey = "EC4", me8enLim = "EC4_phosphoproteome_match", me8enKey = "extEC4", 
                     pr9stLim = "extEC4_metabolome_match", pr9stKey = "EC4", me9enLim = "EC4_proteome_match", me9enKey = "extEC4")
    #Create dataframe for trans-omic network
    df <- createTnetDF(mdatalist, annotCol_names_list, annotCol_names_tnet, keycols_tnet, dataTypes)
    #Make a list of the number of molecules at each point
    sumlists <- makeSumLists(df, mdatalist, annotCol_names_list, sepGenes_vec_list, dataTypes)
    #Save images
    saveImages(sumlists, arrowRelThick, dirsave)
    #Save data 
    dfout <- df[, c("InterlayerEdgeID","NodeLabel_start","NodeLabel_end","OriginalName_start","OriginalName_end","DataType_start","DataType_end","OriginalName_col_start","OriginalName_col_end")]
    colnames(dfout) <- c("Edge/Node Label","EdgeStart/Node Label","EdgeEnd Label","EdgeStart/Node ID","EdgeEnd ID","EdgeStart/Node SourceDataType","EdgeEnd SourceDataType","EdgeStart/Node IDtype","EdgeEnd IDtype")
    dfout <- unique(dfout)
    dfout <- dfout[order(dfout[,1],dfout[,2],dfout[,3],dfout[,4],dfout[,5]),]
    write.xlsx(dfout, paste0(dirsave,"/NetworkFiles/", prefix, "Tnet", suffix, ".xlsx"), col.names = TRUE, row.names = FALSE, overwrite = TRUE)
    
    return()
  }
  
  createTnetDF <- function(mdatalist, annotCol_names_list, annotCol_names_tnet, keycols_tnet, dataTypes){
    df <- as.data.frame(matrix(ncol = 23, nrow = 0))
    colnames(df) <- c("OriginalName_start","OriginalName_end","InterlayerEdgeID","LayerID_start","LayerID_end","DataType_start",
                      "DataType_end","OriginalName_col_start","OriginalName_col_end","LimitCol_start","LimitCol_end",
                      "annotColName_cols_start","annotColName_cols_end","annotColName_val1_start","annotColName_val2_start",
                      "annotColName_val1_end","annotColName_val2_end","rownames_start","rownames_end","VarName_start","VarName_end",
                      "NodeLabel_start","NodeLabel_end")
    #Add edges based on one omic dataset
    df <- addEdgeData_fromOne(df, mdatalist, dataType = "phosphoproteome", startName = keycols_tnet["ph1st"], endName = annotCol_names_tnet["phosphoproteome"], limitName = NULL, InterlayerID = 1, annotCol_names_list)
    df <- addEdgeData_fromOne(df, mdatalist, dataType = "phosphoproteome", startName = annotCol_names_tnet["phosphoproteome"], endName = annotCol_names_tnet["phosphoproteome_protein"], limitName = NULL, InterlayerID = 2, annotCol_names_list)
    df <- addEdgeData_fromOne(df, mdatalist, dataType = "transcriptome", startName = keycols_tnet["tr3st"], endName = annotCol_names_tnet["transcriptome"], limitName = NULL, InterlayerID = 3, annotCol_names_list)
    df <- addEdgeData_fromOne(df, mdatalist, dataType = "transcriptome", startName = keycols_tnet["tr4st"], endName = annotCol_names_tnet["transcriptome"], limitName = NULL, InterlayerID = 4, annotCol_names_list)
    df <- addEdgeData_fromOne(df, mdatalist, dataType = "proteome", startName = keycols_tnet["pr5st"], endName = annotCol_names_tnet["proteome"], limitName = NULL, InterlayerID = 5, annotCol_names_list)
    #Add edges using two omic datasets
    df <- addEdgeData_fromTwo(df, mdatalist, targetDataTypes = c("transcriptome","proteome"), 
                              startColNames = c(annotCol_names_tnet["transcriptome"], keycols_tnet["tr6stLim"], keycols_tnet["tr6stKey"]), 
                              endColNames = c(annotCol_names_tnet["proteome"], keycols_tnet["pr6enLim"], annotCol_names_tnet["proteome"]), 
                              InterlayerID = 6, annotCol_names_list)
    df <- addEdgeData_fromTwo(df, mdatalist, targetDataTypes = c("phosphoproteome","proteome"), 
                              startColNames = c(annotCol_names_tnet["phosphoproteome"], keycols_tnet["ph7stLim"], annotCol_names_tnet["phosphoproteome_protein"]), 
                              endColNames = c(annotCol_names_tnet["proteome"], keycols_tnet["pr7enLim"], annotCol_names_tnet["proteome"]), 
                              InterlayerID = 7, annotCol_names_list)
    df <- addEdgeData_fromTwo(df, mdatalist, targetDataTypes = c("phosphoproteome","metabolome"), 
                              startColNames = c(annotCol_names_tnet["phosphoproteome_protein"], keycols_tnet["ph8stLim"], keycols_tnet["ph8stKey"]), 
                              endColNames = c(annotCol_names_tnet["metabolome"], keycols_tnet["me8enLim"], keycols_tnet["me8enKey"]), 
                              InterlayerID = 8, annotCol_names_list)
    df <- addEdgeData_fromTwo(df, mdatalist, targetDataTypes = c("proteome","metabolome"), 
                              startColNames = c(annotCol_names_tnet["proteome"], keycols_tnet["pr9stLim"], keycols_tnet["pr9stKey"]), 
                              endColNames = c(annotCol_names_tnet["metabolome"], keycols_tnet["me9enLim"], keycols_tnet["me9enKey"]), 
                              InterlayerID = 9, annotCol_names_list)
    #Add isolated nodes
    df <- addIsolatedNodes(df, lay = 3, mdatalist, targetdatatype = "phosphoproteome", colname = annotCol_names_tnet["phosphoproteome"], annotCol_names_list)
    df <- addIsolatedNodes(df, lay = 4, mdatalist, targetdatatype = "transcriptome", colname = annotCol_names_tnet["transcriptome"], annotCol_names_list)
    df <- addIsolatedNodes(df, lay = 5, mdatalist, targetdatatype = "proteome", colname = annotCol_names_tnet["proteome"], annotCol_names_list)
    df <- addIsolatedNodes(df, lay = 6, mdatalist, targetdatatype = "metabolome", colname = annotCol_names_tnet["metabolome"], annotCol_names_list)
    #Add other data for graphic output
    df <- add2DdispData(df, mdatalist, annotCol_names_list, dataTypes)
    #Arrange
    df <- df[order(df$LayerID_start, df$LayerID_end),]
    roworder <- c(seq(1:8), "8_Pho", 9, "Isolated")
    roworder <- roworder[roworder %in% unique(df$InterlayerEdgeID)]
    df$InterlayerEdgeID <- as.factor(df$InterlayerEdgeID)
    df <- df[order(factor(df$InterlayerEdgeID, levels = roworder)),]
    if(nrow(df)>0){
      tmp <- data.frame(VarNames = c("PrioSIG","PrioEXP1","PrioEXP2","PrioEXP3","PhoSIG1","PhoSIG2","PhoSIG3",
                                     "TraEXP1","TraEXP2","TraEXP4","ExpSIG1","ExpSIG3",
                                     "ExpEXP1","ExpEXP2","ExpEXP3","ExpEXP5","MetMET1","MetMET2"),
                        NodeLabels = c("a","d","e","e","b","b","b","f","f","f","c","c","g","g","g","g","h","h"),
                        stringsAsFactors = FALSE)
      df$NodeLabel_start[is.na(df$NodeLabel_start)] <- tmp$NodeLabels[match(df$VarName_start[is.na(df$NodeLabel_start)], tmp$VarNames)]
      df$NodeLabel_end[is.na(df$NodeLabel_end)] <- tmp$NodeLabels[match(df$VarName_end[is.na(df$NodeLabel_end)], tmp$VarNames)]
    }
    return(df)
  }

  extAnnotCols <- function(annotCols, targetcol, unique = TRUE, sep = FALSE, othercols = NULL, limitcol = NULL){
    if(!all(c(targetcol, othercols, limitcol) %in% colnames(annotCols))){stop()}
    tmpdf <- as.data.frame(annotCols[,colnames(annotCols)[colnames(annotCols) %in% c(targetcol, othercols, limitcol)]])
    if(ncol(tmpdf) == 1){colnames(tmpdf) <- c(targetcol, othercols, limitcol)}
    tmpdf$rownames <- rownames(annotCols)
    tmpdf <- tmpdf[tmpdf[,targetcol] != "" & !is.na(tmpdf[,targetcol]),]
    if(!is.null(limitcol)){
      if(class(tmpdf[,limitcol]) == "logical"){
        tmpdf <- tmpdf[tmpdf[,limitcol],]
      }else{
        tmpdf <- tmpdf[tmpdf[,limitcol] != "" & !is.na(tmpdf[,limitcol]),]
      }
      if(!limitcol %in% c(targetcol, othercols)){
        tmpdf[,limitcol] <- NULL
      }
    }
    if(sep){
      colnames(tmpdf)[colnames(tmpdf) == targetcol] <- "targetcol"
      tmpdf <- tmpdf %>% dplyr::mutate(targetcol = strsplit(as.character(targetcol), ";")) %>% unnest(targetcol)
      colnames(tmpdf)[colnames(tmpdf) == "targetcol"] <- targetcol
    }
    if(unique){
      colnames(tmpdf)[colnames(tmpdf) == targetcol] <- "targetcol"
      tmpdf <- tmpdf %>%
        dplyr::group_by(targetcol) %>%
        summarise_all(list(~paste(unique(.), collapse=";")))
      colnames(tmpdf)[colnames(tmpdf) == "targetcol"] <- targetcol
    }
    return(tmpdf)
  }
  
  makeSumLists <- function(df, mdatalist, annotCol_names_list, sepGenes_vec_list, dataTypes){
    #Make sumCountList (summary of molecular counts at each point) and sumIDlist (summary of molecular IDs at each point) from tnetDf
    sumCountList <- sumIDlist <- list()
    tmpElem <- list(molcount_origCol = NA, molcount_1 = NA, molcount_2 = NA)
    tmpID <- list(annotCol_names1 = NA, annotCol_names2 = NA, annotCol_names = NA)
    nmElem <- c("PrioSIG","PrioEXP1","PrioEXP2","PrioEXP3","PhoSIG1","PhoSIG2","PhoSIG3","TraEXP1","TraEXP2","TraEXP3","TraEXP4","ExpSIG1","ExpSIG2","ExpSIG3","ExpEXP1","ExpEXP2","ExpEXP3","ExpEXP4","ExpEXP5","MetMET1","MetMET2","MetMET3")
    nmID <- c("PrioSIG","PrioEXP1","PrioEXP2","PrioEXP3","PhoSIG1","PhoSIG2","PhoSIG3","TraEXP1","TraEXP2","TraEXP3","TraEXP4","ExpSIG1","ExpSIG2","ExpSIG3","ExpEXP1","ExpEXP2","ExpEXP3","ExpEXP4","ExpEXP5","MetMET1","MetMET2","MetMET3")
    for(i in seq_along(nmElem)){sumCountList[[nmElem[i]]] <- tmpElem}
    for(i in seq_along(nmID)){sumIDlist[[nmID[i]]] <- tmpID}
    dfnode <- data.frame(OriginalName = c(df$OriginalName_start,df$OriginalName_end),
                         annotColName_cols = c(df$annotColName_cols_start,df$annotColName_cols_end),
                         annotColName_val1 = c(df$annotColName_val1_start,df$annotColName_val1_end), annotColName_val2 = c(df$annotColName_val2_start,df$annotColName_val2_end),
                         DataType = c(df$DataType_start,df$DataType_end),
                         rownames = c(df$rownames_start, df$rownames_end), VarName = c(df$VarName_start,df$VarName_end))
    nodenames <- c("PrioSIG","PrioEXP1","PrioEXP2","PrioEXP3","PhoSIG3","TraEXP4","ExpSIG3","ExpEXP5","ExpSIG2",
                   "PhoSIG1","TraEXP1","TraEXP2","ExpEXP3","ExpEXP1","ExpEXP2","MetMET1","MetMET2")
    for(nd in seq_along(nodenames)){
      nodename <- nodenames[nd]
      tmpdf <- dfnode[grep(paste0("^",nodename,"$"), dfnode$VarName),]
      tmpdf$VarName <- NULL
      datatype <- as.character(tmpdf$DataType[1])
      annotCol_names <- as.character(tmpdf$annotColName_cols[1])
      if(nrow(tmpdf) == 0){next}
      if(nodename %in% c("PrioSIG","PrioEXP1","PrioEXP2","PrioEXP3","TraEXP4","ExpSIG1","ExpSIG3","ExpEXP1","ExpEXP2")){
        tmpdf1 <- tmpdf %>% dplyr::group_by(OriginalName) %>% summarise_all(funs(paste(unique(.), collapse=";")))
      }else{
        tmpdf1 <- tmpdf %>% dplyr::group_by(rownames) %>% summarise_all(funs(paste(unique(.), collapse=";")))
      }
      tmpdf1$OriginalName <- tmpdf1$rownames <- NULL
      sumCountList[[nodename]][["molcount_origCol"]] <- nrow(tmpdf1)
      sumCountList[[nodename]][["molcount_1"]] <- nrow(tmpdf1[tmpdf1$annotColName_val1 != "" & !is.na(tmpdf1$annotColName_val1),])
      sumCountList[[nodename]][["molcount_2"]] <- nrow(tmpdf1[tmpdf1$annotColName_val2 != "" & !is.na(tmpdf1$annotColName_val2),])
      if(nrow(tmpdf1) != 0){
        sepGenes_vec <- c(sepGenes_vec_list[[1]][datatype], sepGenes_vec_list[[2]][datatype])
        if(sepGenes_vec[1]){
          annotCol_names1 <- paste(unique(as.character(unlist(sapply(tmpdf1$annotColName_val1, strsplit, ";")))), collapse = ";")
        }else{
          annotCol_names1 <- paste(unique(as.character(unlist(tmpdf1$annotColName_val1))), collapse = ";")
        }
        if(sepGenes_vec[2]){
          annotCol_names2 <- paste(unique(as.character(unlist(sapply(tmpdf1$annotColName_val2, strsplit, ";")))), collapse = ";")
        }else{
          annotCol_names2 <- paste(unique(as.character(unlist(tmpdf1$annotColName_val2))), collapse = ";")
        }
        sumIDlist[[nodename]] <- list(annotCol_names1 = annotCol_names1, annotCol_names2 = annotCol_names2, annotCol_names = annotCol_names)
      }
    }
    nodes <- c("PhoSIG2","TraEXP3","ExpEXP4","MetMET3","ExpSIG1")
    names(nodes) <- c("phosphoproteome","transcriptome","proteome","metabolome","phos_other")
    if(!any(c("PhoSIG1","PhoSIG2","PhoSIG3") %in% as.character(unique(dfnode$VarName)))){
      nodes <- nodes[nodes != "PhoSIG2"]
    }
    if(!any(c("ExpSIG1","ExpSIG2","ExpSIG3") %in% as.character(unique(dfnode$VarName)))){
      nodes <- nodes[nodes != "ExpSIG1"]
    }
    if(length(nodes) > 0){
      for(n in seq_along(nodes)){
        nodename <- nodes[n]
        if(nodename == "ExpSIG1"){
          dataType <- "phosphoproteome"
          ac <- annotCols(mdatalist[[which(dataTypes == "phosphoproteome")]])
          ac <- ac[,unique(c("Proteins", annotCol_names_list[[1]][dataType], annotCol_names_list[[2]][dataType]))]
          ac <- ac %>% dplyr::group_by(Proteins) %>% summarise_all(funs(paste(unique(.), collapse=";")))
        }else{
          dataType <- names(nodes)[n]
          ac <- annotCols(mdatalist[[which(dataTypes == dataType)]])
        }
        sumCountList[[nodename]][["molcount_origCol"]] <- nrow(ac)
        if(nrow(ac) == 0){
          sumCountList[[nodename]][["molcount_1"]] <- sumCountList[[nodename]][["molcount_2"]] <- 0
          annotCol_names1 <- annotCol_names2 <- ""
        }else{
          sumCountList[[nodename]][["molcount_1"]] <- nrow(ac[ac[,annotCol_names_list[[1]][dataType]] != "" & !is.na(ac[,annotCol_names_list[[1]][dataType]]),])
          sumCountList[[nodename]][["molcount_2"]] <- nrow(ac[ac[,annotCol_names_list[[2]][dataType]] != "" & !is.na(ac[,annotCol_names_list[[2]][dataType]]),])
          if(sepGenes_vec_list[[1]][dataType]){
            annotCol_names1 <- paste(unique(as.character(unlist(sapply(ac[,annotCol_names_list[[1]][dataType]], strsplit, ";")))), collapse = ";")
          }else{
            annotCol_names1 <- paste(unique(as.character(unlist(ac[,annotCol_names_list[[1]][dataType]]))), collapse = ";")
          }
          if(sepGenes_vec_list[[2]][dataType]){
            annotCol_names2 <- paste(unique(as.character(unlist(sapply(ac[,annotCol_names_list[[2]][dataType]], strsplit, ";")))), collapse = ";")
          }else{
            annotCol_names2 <- paste(unique(as.character(unlist(ac[,annotCol_names_list[[1]][dataType]]))), collapse = ";")
          }
        }
        sumIDlist[[nodename]] <- list(annotCol_names1 = annotCol_names1, annotCol_names2 = annotCol_names2, annotCol_names = paste0(annotCol_names_list[[1]][dataType],",",annotCol_names_list[[2]][dataType]))
      }
    }
    sumlists <- list(sumCountList = sumCountList, sumIDlist = sumIDlist)
    
    return(sumlists)
  }
  
  addEdgeData_fromOne <- function(df, mdatalist, dataType, startName, endName, limitName = NULL, InterlayerID, annotCol_names_list){
    dataTypes <- names(annotCol_names_list[[1]])
    mdata <- mdatalist[which(dataTypes == dataType)][[1]]
    annotCols <- annotCols(mdata)
    if(nrow(annotCols) == 0){return(df)}
    if(InterlayerID == 1){
      startlayerID <- 1; endlayerID <- 3
      varnameStart <- "PrioSIG"; varnameEnd <- "PhoSIG1"
      seps = c(TRUE,FALSE)
    }
    if(InterlayerID == 3){
      startlayerID <- 1; endlayerID <- 4
      varnameStart <- "PrioEXP1"; varnameEnd <- "TraEXP1"
      seps = c(TRUE,FALSE)
    }
    if(InterlayerID == 4){
      startlayerID <- 2; endlayerID <- 4
      varnameStart <- "PrioEXP2"; varnameEnd <- "TraEXP2"
      seps = c(TRUE,FALSE)
    }
    if(InterlayerID == 5){
      startlayerID <- 2; endlayerID <- 5
      varnameStart <- "PrioEXP3"; varnameEnd <- "ExpEXP3"
      seps = c(TRUE,FALSE)
    }
    if(InterlayerID %in% c(2,7)){
      startlayerID <- 3; endlayerID <- 5
      if(InterlayerID == 2){
        varnameStart <- "PhoSIG2"; varnameEnd <- "ExpSIG1"
      }else{
        varnameStart <- "PhoSIG3"; varnameEnd <- "ExpEXP1"
      }
      seps = c(FALSE,FALSE)
    }
    if((!is.null(limitName) & !all(c(startName, endName, limitName) %in% colnames(annotCols))) | (is.null(limitName) & !all(c(startName, endName) %in% colnames(annotCols)))){
      return(df)
    }
    tmpdf <- annotCols[,colnames(annotCols)[colnames(annotCols) %in% c(startName, endName, limitName)]]
    tmpdf$annotColName_val1_start <- annotCols[,annotCol_names_list[[1]][dataType]]
    tmpdf$annotColName_val2_start <- annotCols[,annotCol_names_list[[2]][dataType]]
    tmpdf$rownames_start <- tmpdf$rownames_end <- rownames(annotCols)
    tmpdf <- tmpdf[tmpdf[,startName] != "" & !is.na(tmpdf[,startName]),]
    tmpdf <- tmpdf[tmpdf[,endName] != "" & !is.na(tmpdf[,endName]),]
    if(nrow(tmpdf) == 0){return(df)}
    if(!is.null(limitName)){
      if(class(tmpdf[,limitName]) == "logical"){
        tmpdf <- tmpdf[tmpdf[,limitName],]
      }else{
        tmpdf <- tmpdf[tmpdf[,limitName] != "" & !is.na(tmpdf[,limitName]),]
      }
      tmpdf[,limitName] <- NULL
    }
    if(nrow(tmpdf) == 0){return(df)}
    colnames(tmpdf)[colnames(tmpdf) == startName] <- "StartNode"
    colnames(tmpdf)[colnames(tmpdf) == endName] <- "EndNode"
    if(seps[1]){tmpdf <- tmpdf %>% dplyr::mutate(StartNode = strsplit(as.character(StartNode), ";")) %>% unnest(StartNode)}
    if(seps[2]){tmpdf <- tmpdf %>% dplyr::mutate(EndNode = strsplit(as.character(EndNode), ";")) %>% unnest(EndNode)}
    tmpdf <- tmpdf[, c("StartNode", "EndNode", colnames(tmpdf)[!colnames(tmpdf) %in% c("StartNode", "EndNode")])]
    colnames(tmpdf)[1:2] <- c("OriginalName_start","OriginalName_end")
    tmpdf <- as.data.frame(tmpdf)
    # if(nrow(tmpdf) == 0){return(df)}
    if(is.null(limitName)){
      tmpdf$LimitCol_start <- tmpdf$LimitCol_end <- NA
    }else{
      tmpdf$LimitCol_start <- tmpdf$LimitCol_end <- rep(limitName, nrow(tmpdf))
    }
    tmpdf$annotColName_val1_end <- tmpdf$annotColName_val1_start
    tmpdf$annotColName_val2_end <- tmpdf$annotColName_val2_start
    tmpdfAdd <- data.frame(InterlayerEdgeID = InterlayerID, LayerID_start = startlayerID, LayerID_end = endlayerID,
                           DataType_start = dataType, DataType_end = dataType, OriginalName_col_start = startName, OriginalName_col_end = endName,
                           annotColName_cols_start = paste0(annotCol_names_list[[1]][dataType],";",annotCol_names_list[[2]][dataType]),
                           annotColName_cols_end = paste0(annotCol_names_list[[1]][dataType],";",annotCol_names_list[[2]][dataType]),
                           VarName_start = varnameStart, VarName_end = varnameEnd, stringsAsFactors=FALSE)
    rownames(tmpdfAdd) <- NULL
    tmpdf <- bind_cols(tmpdf, tmpdfAdd)
    #Order
    tmpdf <- tmpdf[,c("OriginalName_start","OriginalName_end","InterlayerEdgeID","LayerID_start","LayerID_end",
                      "DataType_start","DataType_end","OriginalName_col_start","OriginalName_col_end","LimitCol_start","LimitCol_end",
                      "annotColName_cols_start","annotColName_cols_end","annotColName_val1_start","annotColName_val2_start",
                      "annotColName_val1_end","annotColName_val2_end", "rownames_start", "rownames_end", "VarName_start", "VarName_end")]
    tmpdf <- unique(tmpdf)
    if(nrow(tmpdf) == 0){return(df)}
    df <- merge(df, tmpdf, all = TRUE)
    return(df)
  }
  
  addEdgeData_fromTwo <- function(df, mdatalist, targetDataTypes, startColNames, endColNames, InterlayerID, annotCol_names_list){
    dataTypes <- names(annotCol_names_list[[1]])
    if(InterlayerID %in% c(4,7)){
      startlayerID <- 3; endlayerID <- 5
      if(InterlayerID == 4){
        varnameStart <- "PhoSIG2"; varnameEnd <- "ExpSIG1"
      }else{
        varnameStart <- "PhoSIG3"; varnameEnd <- "ExpEXP1"
      }
      seps = c(FALSE,FALSE)
    }
    if(InterlayerID == 6){
      startlayerID <- 4; endlayerID <- 5 
      varnameStart <- "TraEXP4"; varnameEnd <- "ExpEXP2"
      seps = c(FALSE,FALSE)
    }
    if(InterlayerID %in% c(8,9)){
      startlayerID <- 5; endlayerID <- 6; 
      if(InterlayerID == 8){
        varnameStart <- "ExpSIG3"; varnameEnd <- "MetMET1"
        seps = c(FALSE,FALSE)
      }else{
        varnameStart <- "ExpEXP5"; varnameEnd <- "MetMET2"
        seps = c(FALSE,FALSE)
      }
    }
    
    #1(Start node)
    mdata1 <- mdatalist[which(dataTypes == targetDataTypes[1])][[1]]
    annotCols1 <- annotCols(mdata1)
    if(nrow(annotCols1) == 0){return(df)}
    if((startColNames[2] != FALSE & !all(startColNames %in% colnames(annotCols1))) | (startColNames[2] == FALSE & !all(startColNames[c(1,3)] %in% colnames(annotCols1)))){
      return(df)
    }
    tmpdf1 <- annotCols1[, startColNames[startColNames != "FALSE"]]
    tmpdf1$annotColName_val1_start <- annotCols1[,annotCol_names_list[[1]][targetDataTypes[1]]]
    tmpdf1$annotColName_val2_start <- annotCols1[,annotCol_names_list[[2]][targetDataTypes[1]]]
    tmpdf1$rownames_start <- rownames(annotCols1)
    tmpdf1 <- tmpdf1[tmpdf1[,startColNames[1]] != "" & !is.na(tmpdf1[,startColNames[1]]),]
    if(startColNames[2] != "FALSE"){
      if(class(tmpdf1[,startColNames[2]]) != "logical"){stop("The annotCols column assigned to 'startColNames[2]' is not 'logical.' Please check your data.")}
      tmpdf1 <- tmpdf1[tmpdf1[,startColNames[2]],]
      tmpdf1[,startColNames[2]] <- NULL
    }
    if(nrow(tmpdf1) == 0){return(df)}
    colnames(tmpdf1)[colnames(tmpdf1) == startColNames[1]] <- "StartNode"
    colnames(tmpdf1)[colnames(tmpdf1) == startColNames[3]] <- "key"
    tmpdf1 <- tmpdf1[, c("StartNode", "key", colnames(tmpdf1)[!colnames(tmpdf1) %in% c("StartNode", "key")])]
    colnames(tmpdf1)[1] <- "OriginalName_start"
    if(seps[1]){tmpdf1 <- tmpdf1 %>% dplyr::mutate(OriginalName_start = strsplit(as.character(OriginalName_start), ";")) %>% unnest(OriginalName_start)}
    tmpdf1 <- tmpdf1 %>% dplyr::mutate(key = strsplit(as.character(key), ";")) %>% unnest(key)
    if(nrow(tmpdf1) == 0){return(df)}
    
    #2(End node)
    mdata2 <- mdatalist[which(dataTypes == targetDataTypes[2])][[1]]
    annotCols2 <- annotCols(mdata2)
    if(nrow(annotCols2) == 0){return(df)}
    if((endColNames[2] != FALSE & !all(endColNames %in% colnames(annotCols2))) | (endColNames[2] == FALSE & !all(endColNames[c(1,3)] %in% colnames(annotCols2)))){
      return(df)
    }
    tmpdf2 <- annotCols2[,unique(endColNames[endColNames != "FALSE"])]
    tmpdf2$annotColName_val1_end <- annotCols2[,annotCol_names_list[[1]][targetDataTypes[2]]]
    tmpdf2$annotColName_val2_end <- annotCols2[,annotCol_names_list[[2]][targetDataTypes[2]]]
    tmpdf2$rownames_end <- rownames(annotCols2)
    tmpdf2 <- tmpdf2[tmpdf2[,endColNames[1]] != "" & !is.na(tmpdf2[,endColNames[1]]),]
    if(endColNames[2] != "FALSE"){
      if(class(tmpdf2[,endColNames[2]]) != "logical"){stop("The annotCols column assigned to 'endColNames[2]' is not 'logical.' Please check your data.")}
      tmpdf2 <- tmpdf2[tmpdf2[,endColNames[2]],]
      tmpdf2[,endColNames[2]] <- NULL
    }
    if(nrow(tmpdf2) == 0){return(df)}
    colnames(tmpdf2)[colnames(tmpdf2) == endColNames[1]] <- "EndNode"
    if(endColNames[1] == endColNames[3]){
      tmpdf2[,endColNames[3]] <- tmpdf2$EndNode
    }
    colnames(tmpdf2)[colnames(tmpdf2) == endColNames[3]] <- "key"
    tmpdf2 <- tmpdf2[, c("EndNode", "key", colnames(tmpdf2)[!colnames(tmpdf2) %in% c("EndNode", "key")])]
    colnames(tmpdf2)[1] <- "OriginalName_end"
    if(seps[2]){tmpdf2 <- tmpdf2 %>% dplyr::mutate(OriginalName_end = strsplit(as.character(OriginalName_end), ";")) %>% unnest(OriginalName_end)}
    tmpdf2 <- tmpdf2 %>% dplyr::mutate(key = strsplit(as.character(key), ";")) %>% unnest(key)
    if(nrow(tmpdf2) == 0){return(df)}

    #Combine two node dfs based on key values to explore edges.
    tmpdf <- merge(tmpdf1, tmpdf2, by = "key", all = TRUE)
    tmpdf$key <- NULL
    tmpdf <- tmpdf[tmpdf$OriginalName_start != "" & !is.na(tmpdf$OriginalName_start) & tmpdf$OriginalName_end != "" & !is.na(tmpdf$OriginalName_end),]
    tmpdfAdd <- data.frame(InterlayerEdgeID = InterlayerID, LayerID_start = startlayerID, LayerID_end=endlayerID, 
                           DataType_start = targetDataTypes[1], DataType_end = targetDataTypes[2], 
                           OriginalName_col_start = startColNames[1], OriginalName_col_end = endColNames[1], 
                           LimitCol_start = startColNames[2], LimitCol_end = endColNames[2], 
                           annotColName_cols_start = paste0(annotCol_names_list[[1]][targetDataTypes[1]],";",annotCol_names_list[[2]][targetDataTypes[1]]), 
                           annotColName_cols_end = paste0(annotCol_names_list[[1]][targetDataTypes[2]],";",annotCol_names_list[[2]][targetDataTypes[2]]), 
                           VarName_start = varnameStart, VarName_end = varnameEnd, stringsAsFactors=FALSE)
    rownames(tmpdfAdd) <- NULL
    tmpdf <- bind_cols(tmpdf, tmpdfAdd)
    tmpdf <- tmpdf[,c("OriginalName_start","OriginalName_end","InterlayerEdgeID","LayerID_start","LayerID_end","DataType_start","DataType_end",
                      "OriginalName_col_start","OriginalName_col_end","LimitCol_start","LimitCol_end",
                      "annotColName_cols_start","annotColName_cols_end","annotColName_val1_start","annotColName_val2_start",
                      "annotColName_val1_end","annotColName_val2_end", "rownames_start", "rownames_end", "VarName_start", "VarName_end")] #"DisplayName_start","DisplayName_end","DisplayName_col_start", "DisplayName_col_end", dispColNames,
    tmpdf <- unique(tmpdf)
    df <- bind_rows(df, tmpdf)
    return(df)
  }
  
  addIsolatedNodes <- function(df, lay, mdatalist, targetdatatype, colname, annotCol_names_list){
    dataTypes <- names(annotCol_names_list[[1]])
    mdata <- mdatalist[[which(dataTypes == targetdatatype)]]
    annotCols <- annotCols(mdata)
    if(nrow(annotCols) == 0){return(df)}
    annotCol_names <- c(annotCol_names_list[[1]][targetdatatype], annotCol_names_list[[2]][targetdatatype])

    removeRows <- unique(c(df$rownames_start[df$DataType_start == targetdatatype],df$rownames_end[df$DataType_end == targetdatatype]))
    annotCols_iso <- annotCols[rownames(annotCols)[!rownames(annotCols) %in% removeRows],]
    
    nr <- nrow(annotCols_iso)
    if(nr != 0){
      switch(targetdatatype,
             phosphoproteome = {nodelabel <- "b"},
             transcriptome = {nodelabel <- "f"},
             proteome = {nodelabel <- "g"},
             metabolome = {nodelabel <- "h"}
      )
      isonodedf <- data.frame(
        OriginalName_start = annotCols_iso[, colname], OriginalName_end = NA,
        InterlayerEdgeID = rep("Isolated", nr),
        LayerID_start = rep(lay, nr), LayerID_end = NA,
        DataType_start = rep(targetdatatype, nr), DataType_end = NA,
        OriginalName_col_start = rep(colname, nr), OriginalName_col_end = NA,
        LimitCol_start = NA, LimitCol_end = NA,
        annotColName_cols_start = rep(paste(annotCol_names, collapse = ";"), nr), annotColName_cols_end = NA,
        annotColName_val1_start = annotCols_iso[, annotCol_names[1]],
        annotColName_val2_start = annotCols_iso[, annotCol_names[2]],
        annotColName_val1_end = NA,
        annotColName_val2_end = NA,
        rownames_start = rownames(annotCols_iso), rownames_end = NA,
        VarName_start = NA, VarName_end = NA,
        NodeLabel_start = rep(nodelabel, nr), NodeLabel_end = NA,
        stringsAsFactors = FALSE)
      df$InterlayerEdgeID <- as.character(df$InterlayerEdgeID)
      df$NodeLabel_start <- as.character(df$NodeLabel_start)
      df <- bind_rows(df,isonodedf)
      df <- unique(df)
    }
    return(df)
  }

  add2DdispData <- function(df, mdatalist, annotCol_names_list, dataTypes){
    mdata <- mdatalist[which(dataTypes == "phosphoproteome")][[1]]
    mdata_met <- mdatalist[which(dataTypes == "metabolome")][[1]]
    annotCols <- annotCols(mdata)
    annotCol_names <- c(annotCol_names_list[[1]]["phosphoproteome"], annotCol_names_list[[2]]["phosphoproteome"])
    annotCols <- extAnnotCols(annotCols, targetcol = annotCol_names_list[[1]][["phosphoproteome"]],
                              unique = FALSE, sep = FALSE, othercols = c(annotCol_names, "EC4"), limitcol = "extEC4_metabolome_match")
    annotCols_met <- annotCols(mdata_met)
    annotCols_met <- extAnnotCols(annotCols_met, targetcol = "extEC4", unique = TRUE, sep = TRUE, 
                                  othercols = annotCol_names_list[[1]][["metabolome"]], 
                                  limitcol = "EC4_phosphoproteome_match")
    if(nrow(annotCols) == 0 | nrow(annotCols_met) == 0){return(df)}
    
    annotCols <- annotCols %>% dplyr::mutate(EC4 = strsplit(as.character(EC4), ";")) %>% unnest(EC4)
    colnames(annotCols_met)[colnames(annotCols_met) == "extEC4"] <- "EC4"
    colnames(annotCols_met)[colnames(annotCols_met) == "rownames"] <- "rownames_end"
    annotCols <- merge(annotCols, annotCols_met, by = "EC4")
    annotCols <- annotCols %>% dplyr::group_by(rownames) %>% dplyr::summarise_all(funs(paste(unique(.), collapse=";")))
    
    nr <- nrow(annotCols)
    if(nr == 0){return(df)}
    tmpdf <- data.frame(#DisplayName_start = NA, DisplayName_end = NA,
                        OriginalName_start = annotCols[,annotCol_names_list[[1]][["phosphoproteome"]]], OriginalName_end = NA,
                        InterlayerEdgeID = rep("8_Pho",nr), LayerID_start = NA, LayerID_end = NA,
                        DataType_start = rep("phosphoproteome", nr), DataType_end = rep("metabolome", nr),
                        OriginalName_col_start = rep(annotCol_names_list[[1]][["phosphoproteome"]], nr), OriginalName_col_end = NA,
                        LimitCol_start = NA, LimitCol_end = NA,
                        annotColName_cols_start = rep(paste(annotCol_names, collapse = ";"), nr), 
                        annotColName_cols_end = rep(annotCol_names_list[[1]][["metabolome"]], nr), 
                        annotColName_val1_start = annotCols[, annotCol_names[1]], 
                        annotColName_val2_start = annotCols[, annotCol_names[2]],
                        annotColName_val1_end = annotCols[,annotCol_names_list[[1]][["metabolome"]]], annotColName_val2_end = NA,
                        VarName_start = rep("ExpSIG2", nr), VarName_end = NA,
                        rownames_start = annotCols$rownames, rownames_end = annotCols$rownames_end)
    colnames(tmpdf)[c(1,14:16)] <- c("OriginalName_start","annotColName_val1_start", "annotColName_val2_start", "annotColName_val1_end")
    df <- merge(df, tmpdf, all = TRUE)
    
    return(df)
  }
  
  saveImages <- function(sumlists, arrowRelThick, dirsave){
    sumCountList <- sumlists[["sumCountList"]]
    thicklist <- sumCountList[names(sumCountList)]
    thickmax <- 13000 #Adjustment val (temp)
    th <- list()
    for(s in names(thicklist)){
      relthick <- as.numeric(thicklist[[s]][["molcount_origCol"]]) / thickmax
      th[[s]] <- relthick * 0.1
    }
    #make diagrams
    saveDiagram1(sumCountList, th, fac = arrowRelThick[1], dirsave)
    saveDiagram2(sumCountList, th, fac = arrowRelThick[2], dirsave)
  }
  
  makeGrobs <- function(th_start, th_end, rellen = 1, label = NULL){
    if(is.null(label)){nullflg <- TRUE}else{nullflg <- FALSE}
    label <- paste0(label, "\n")
    if(is.na(th_start)){
      th_start <- 0
      label <- paste0(label, " *NA Start")}
    if(is.na(th_end)){
      th_end <- 0
      label <- paste0(label, " *NA End")}
    pdf(file = NULL)
    txt <- grid.text(label, x = unit(2.5, "cm"), y = unit(0.5, "npc"), gp = gpar(fontsize = 9))
    arrow <- grid.polygon(x=c(2.5-th_start*5, 2.5-th_end*5, 2.5-(0.1+th_end*8), 2.5, 2.5+(0.1+th_end*8), 2.5+th_end*5, 2.5+th_start*5, 2.5-th_start*5), y=c(3.9*rellen, 0.2+th_end*3, 0.3+th_end*3, 0, 0.3+th_end*3, 0.2+th_end*3, 3.9*rellen, 3.9*rellen), default.units="cm", gp=gpar(fill="black"))
    
    grid.newpage()
    if(nullflg){
      g <- arrow
    }else{
      g <- gtable(widths = unit(5, 'cm'), heights = unit(c(1, 3.9*rellen+0.1), 'cm'))
      g <- gtable_add_grob(g, txt, t = 1, l=1, b=1, r=1)
      g <- gtable_add_grob(g, arrow, t = 2, l=1, b=2, r=1)
    }
    grid.draw(g)
    dev.off()
    return(g)
  }
  
  saveDiagram1 <- function(sumCountList, th, fac = 100, dirsave){
    fac <- fac/2
    colours <- c("#ffe367", "#ffc732", "#4ea24e", "#f1f132", "#aed75a", "#4a708b")
    cairo_pdf(file = paste0(dirsave, "/NetworkFiles/", prefix, "Tnet", suffix, ".pdf"), width = 11.69/2, height = 8.27/2, family = "Arial")
    grid.newpage()
    pushViewport(viewport(layout=grid.layout(8, 4, widths=unit(c(0, 0.4, 1.8, 1.8), c("inch", "inch", "inch", "inch")),
                                             heights=unit(c(0.2, 0.4, 0.65, 0.65, 0.65, 0.65, 0.65, 0.2), c("inch", "inch", "inch", "inch", "inch", "inch", "inch")))))
    #BG
    for(c in seq_along(colours)){
      grid.rect(gp=gpar(fill = colours[c], lwd = NA), height=unit(0.95, "npc"), width=unit(0.99, "npc"), vp=viewport(layout.pos.col=3:4, layout.pos.row=c+1))
    }
    grid.rect(gp=gpar(fill = "#656565", lwd = NA), y = unit(0.17, "npc"), height=unit(0.133, "npc"), width=unit(0.98, "npc"), vp=viewport(layout.pos.col=3, layout.pos.row=1))
    grid.rect(gp=gpar(fill = "#656565", lwd = NA), y = unit(0.17, "npc"), height=unit(0.133, "npc"), width=unit(0.98, "npc"), vp=viewport(layout.pos.col=4, layout.pos.row=1))#line above boxes (right; Expression regulation)
    grid.rect(gp=gpar(fill = "#656565", lwd = NA), y = unit(0.8, "npc"), height=unit(0.13, "npc"), width=unit(1, "npc"), vp=viewport(layout.pos.col=3:4, layout.pos.row=8))#line below boxes (Metabolic response)
    grid.rect(gp=gpar(fill = NA, lwd = 2, col="#656565"), height=unit(0.97, "npc"), width=unit(0.97, "npc"), vp=viewport(layout.pos.col=3, layout.pos.row=3:6))#outline of Signaling
    grid.rect(gp=gpar(fill = NA, lwd = 2, col="#656565"), height=unit(0.97, "npc"), width=unit(0.97, "npc"), vp=viewport(layout.pos.col=4, layout.pos.row=3:6))#outline of Expression regulation
    grid.rect(gp=gpar(fill = NA, lwd = 2, col="#656565"), height=unit(0.97, "npc"), width=unit(0.99, "npc"), vp=viewport(layout.pos.col=3:4, layout.pos.row=7))#outline of Metabolic response
    #BG Text
    grid.text("Phos-mediated signaling",
              y = unit(0.6, "npc"), 
              gp = gpar(fontsize = 9, fontface = "bold.italic", col = "black"), 
              vp = viewport(layout.pos.col = 3, layout.pos.row = 1))
    grid.text("Expression regulation", 
              y = unit(0.6, "npc"), 
              gp = gpar(fontsize = 9, fontface = "bold.italic", col = "black", lineheight = 0.9), 
              vp = viewport(layout.pos.col = 4, layout.pos.row = 1))
    grid.text("Metabolic response", 
              y = unit(0.4, "npc"), 
              gp = gpar(fontsize = 9, fontface = "bold.italic", col = "black", lineheight = 0.9), 
              vp = viewport(layout.pos.col = 3:4, layout.pos.row = 8))
    #Left column
    grid.text("Input", gp = gpar(fontsize = 8, col = "black"), 
              vp = viewport(layout.pos.col=1.95, layout.pos.row=2))
    grid.text("Effector",
              gp = gpar(fontsize = 8, col = "black", lineheight = 0.9), 
              vp = viewport(layout.pos.col = 1.95, layout.pos.row = 3))
    grid.text("Ph-ome",
              y = unit(0.65, "npc"), gp = gpar(fontsize = 8, col = "black"), 
              vp = viewport(layout.pos.col = 1.95, layout.pos.row = 4))
    grid.text("Tr-ome",
              y = unit(0.65, "npc"), gp = gpar(fontsize = 8, col = "black"), 
              vp=viewport(layout.pos.col=1.95, layout.pos.row=5))
    grid.text("Pr-ome",
              y = unit(0.6, "npc"), gp = gpar(fontsize = 8, col = "black"), 
              vp=viewport(layout.pos.col=1.95, layout.pos.row=6))
    grid.text("Met-ome",
              y = unit(0.6, "npc"), gp = gpar(fontsize = 8, col = "black"), 
              vp=viewport(layout.pos.col=1.95, layout.pos.row=7))
    #Flow chart
    grid.text("Blood insulin \nlevel",
              gp = gpar(fontsize = 7, col = "black", lineheight = 0.9), 
              vp = viewport(layout.pos.col = 3:4, layout.pos.row = 2))
    grid.text("Kin/Pase/Pdoms",
              gp = gpar(fontsize = 7, col = "black", lineheight = 0.9), 
              vp = viewport(layout.pos.col = 3, layout.pos.row = 3))
    if(!is.na(sumCountList[["PrioSIG"]][["molcount_origCol"]]) & sumCountList[["PrioSIG"]][["molcount_origCol"]] != 0){
      grid.text(sumCountList[["PrioSIG"]][["molcount_origCol"]], 
                x = unit(0.67, "npc"), y = unit(0.15, "npc"), 
                gp = gpar(fontsize = 5, col = "black"), 
                vp = viewport(layout.pos.col = 3, layout.pos.row = 3))
    }
    grid.text("TFs",
              x = unit(0.29, "npc"),
              gp = gpar(fontsize = 7, col = "black", lineheight = 0.9), 
              vp = viewport(layout.pos.col = 4, layout.pos.row = 3))
    if(!is.na(sumCountList[["PrioEXP1"]][["molcount_origCol"]]) & sumCountList[["PrioEXP1"]][["molcount_origCol"]] != 0){
      grid.text(sumCountList[["PrioEXP1"]][["molcount_origCol"]], 
                x = unit(0.2, "npc"), y = unit(0.35, "npc"), 
                gp = gpar(fontsize = 5, col = "black"), 
                vp = viewport(layout.pos.col = 4, layout.pos.row = 3))
    }
    grid.text("miRNAs", 
              x = unit(0.7, "npc"),
              gp = gpar(fontsize = 7, col = "black"), 
              vp = viewport(layout.pos.col = 4, layout.pos.row = 3))
    if(!is.na(sumCountList[["PrioEXP2"]][["molcount_origCol"]]) & sumCountList[["PrioEXP2"]][["molcount_origCol"]] != 0){
      grid.text(sumCountList[["PrioEXP2"]][["molcount_origCol"]], x = unit(0.56, "npc"), y = unit(0.35, "npc"), gp = gpar(fontsize = 5, col = "black"), vp=viewport(layout.pos.col=4, layout.pos.row=3))
    }
    if(!is.na(sumCountList[["PrioEXP3"]][["molcount_origCol"]]) & sumCountList[["PrioEXP3"]][["molcount_origCol"]] != 0){
      grid.text(sumCountList[["PrioEXP3"]][["molcount_origCol"]], x = unit(0.9, "npc"), y = unit(0.35, "npc"), gp = gpar(fontsize = 5, col = "black"), vp=viewport(layout.pos.col=4, layout.pos.row=3))
    }
    if(!is.na(sumCountList[["PhoSIG1"]][["molcount_origCol"]]) & sumCountList[["PhoSIG1"]][["molcount_origCol"]] != 0){
      grid.text(sumCountList[["PhoSIG1"]][["molcount_origCol"]], x = unit(0.67, "npc"), y = unit(0.9, "npc"), gp = gpar(fontsize = 5, col = "black"), vp=viewport(layout.pos.col=3, layout.pos.row=4))
    }
    if(!is.na(sumCountList[["PhoSIG2"]][["molcount_1"]]) & sumCountList[["PhoSIG2"]][["molcount_1"]] != 0){
      grid.text(paste0(sumCountList[["PhoSIG2"]][["molcount_1"]], " Phos-sites"),
                gp = gpar(fontsize = 7, col = "black", lineheight = 0.9), 
                vp = viewport(layout.pos.col = 3, layout.pos.row = 4))
    }
    if(!is.na(sumCountList[["PhoSIG3"]][["molcount_origCol"]]) & sumCountList[["PhoSIG3"]][["molcount_origCol"]] != 0){
      grid.text(sumCountList[["PhoSIG3"]][["molcount_origCol"]], x = unit(0.7, "npc"), y = unit(0.75, "npc"), gp = gpar(fontsize = 5, col = "black"), vp=viewport(layout.pos.col=3, layout.pos.row=5))
    }
    if(!is.na(sumCountList[["TraEXP1"]][["molcount_origCol"]]) & sumCountList[["TraEXP1"]][["molcount_origCol"]] != 0){
      grid.text(sumCountList[["TraEXP1"]][["molcount_origCol"]], x = unit(0.2, "npc"), y = unit(0.85, "npc"), gp = gpar(fontsize = 5, col = "black"), vp=viewport(layout.pos.col=4, layout.pos.row=5))
    }
    if(!is.na(sumCountList[["TraEXP2"]][["molcount_origCol"]]) & sumCountList[["TraEXP2"]][["molcount_origCol"]] != 0){
      if(th["TraEXP2"] > 0.05){#When the arrow is too wide, move the text slightly to the left. 
        grid.text(sumCountList[["TraEXP2"]][["molcount_origCol"]], x = unit(0.47, "npc"), y = unit(0.85, "npc"), gp = gpar(fontsize = 5, col = "black"), vp=viewport(layout.pos.col=4, layout.pos.row=5))
      }else{
        grid.text(sumCountList[["TraEXP2"]][["molcount_origCol"]], x = unit(0.54, "npc"), y = unit(0.85, "npc"), gp = gpar(fontsize = 5, col = "black"), vp=viewport(layout.pos.col=4, layout.pos.row=5))
      }
    }
    if(!is.na(sumCountList[["TraEXP3"]][["molcount_origCol"]]) & sumCountList[["TraEXP3"]][["molcount_origCol"]] != 0){
      grid.text(paste0(sumCountList[["TraEXP3"]][["molcount_origCol"]]," mRNAs"), y = unit(0.45, "npc"), gp = gpar(fontsize = 7, col = "black"), vp=viewport(layout.pos.col=4, layout.pos.row=5))
    }
    if(!is.na(sumCountList[["TraEXP4"]][["molcount_origCol"]]) & sumCountList[["TraEXP4"]][["molcount_origCol"]] != 0){
      grid.text(sumCountList[["TraEXP4"]][["molcount_origCol"]], x = unit(0.56, "npc"), y = unit(0.26, "npc"), gp = gpar(fontsize = 5, col = "black"), vp=viewport(layout.pos.col=4, layout.pos.row=5))
    }
    if(!is.na(sumCountList[["ExpSIG1"]][["molcount_origCol"]]) & sumCountList[["ExpSIG1"]][["molcount_origCol"]] != 0){
      grid.text(paste0(sumCountList[["ExpSIG1"]][["molcount_origCol"]], " Proteins"), y = unit(0.55, "npc"), gp = gpar(fontsize = 7, col = "black", lineheight = 0.9), vp=viewport(layout.pos.col=3, layout.pos.row=6))
    }
    if(!is.na(sumCountList[["ExpSIG2"]][["molcount_origCol"]]) & sumCountList[["ExpSIG2"]][["molcount_origCol"]] != 0){
      grid.text(paste0(sumCountList[["ExpSIG2"]][["molcount_origCol"]], " Phos-sites"),
                just = "right", x = unit(0.6, "npc"), y = unit(0.3, "npc"), 
                gp = gpar(fontsize = 5, col = "black"), 
                vp = viewport(layout.pos.col = 3, layout.pos.row = 6))
    }
    if(!is.na(sumCountList[["ExpSIG3"]][["molcount_origCol"]]) & sumCountList[["ExpSIG3"]][["molcount_origCol"]] != 0){
      grid.text(paste0(sumCountList[["ExpSIG3"]][["molcount_origCol"]]," Proteins"), 
                just = "right", x = unit(0.6, "npc"), y = unit(0.2, "npc"), 
                gp = gpar(fontsize = 5, col = "black"), 
                vp = viewport(layout.pos.col = 3, layout.pos.row = 6))
    }
    if(!is.na(sumCountList[["ExpEXP1"]][["molcount_origCol"]]) & sumCountList[["ExpEXP1"]][["molcount_origCol"]] != 0){
      grid.text(sumCountList[["ExpEXP1"]][["molcount_origCol"]], x = unit(0.295, "npc"), y = unit(0.62, "npc"), 
                gp = gpar(fontsize = 5, col = "black"), 
                vp = viewport(layout.pos.col = 4, layout.pos.row = 6))
    }
    if(!is.na(sumCountList[["ExpEXP2"]][["molcount_origCol"]]) & sumCountList[["ExpEXP2"]][["molcount_origCol"]] != 0){
      grid.text(sumCountList[["ExpEXP2"]][["molcount_origCol"]], x = unit(0.56, "npc"), y = unit(0.8, "npc"), gp = gpar(fontsize = 5, col = "black"), vp=viewport(layout.pos.col=4, layout.pos.row=6))
    }
    if(!is.na(sumCountList[["ExpEXP3"]][["molcount_origCol"]]) & sumCountList[["ExpEXP3"]][["molcount_origCol"]] != 0){
      grid.text(sumCountList[["ExpEXP3"]][["molcount_origCol"]], x = unit(0.9, "npc"), y = unit(0.8, "npc"), gp = gpar(fontsize = 5, col = "black"), vp=viewport(layout.pos.col=4, layout.pos.row=6))
    }
    if(!is.na(sumCountList[["ExpEXP4"]][["molcount_origCol"]]) & sumCountList[["ExpEXP4"]][["molcount_origCol"]] != 0){
      grid.text(paste0(sumCountList[["ExpEXP4"]][["molcount_origCol"]], " Proteins"), x = unit(0.65, "npc"), y = unit(0.5, "npc"), gp = gpar(fontsize = 7, col = "black"), vp=viewport(layout.pos.col=4, layout.pos.row=6))
    }
    if(!is.na(sumCountList[["ExpEXP5"]][["molcount_origCol"]]) & sumCountList[["ExpEXP5"]][["molcount_origCol"]] != 0){
      grid.text(sumCountList[["ExpEXP5"]][["molcount_origCol"]], x = unit(0.75, "npc"), y = unit(0.3, "npc"), gp = gpar(fontsize = 5, col = "black"), vp=viewport(layout.pos.col=4, layout.pos.row=6))
    }
    if(!is.na(sumCountList[["MetMET3"]][["molcount_origCol"]]) & sumCountList[["MetMET3"]][["molcount_origCol"]] != 0){
      grid.text(paste0(sumCountList[["MetMET3"]][["molcount_origCol"]], "\nMetabolites"), gp = gpar(fontsize = 7, col = "black", lineheight = 0.9), vp=viewport(layout.pos.col=3:4, layout.pos.row=7))
    }
    if(!is.na(sumCountList[["MetMET1"]][["molcount_origCol"]]) & sumCountList[["MetMET1"]][["molcount_origCol"]] != 0){
      grid.text(sumCountList[["MetMET1"]][["molcount_origCol"]], x = unit(0.35, "npc"), y = unit(0.65, "npc"), gp = gpar(fontsize = 5, col = "black"), vp=viewport(layout.pos.col=3:4, layout.pos.row=7))
    }
    if(!is.na(sumCountList[["MetMET2"]][["molcount_origCol"]]) & sumCountList[["MetMET2"]][["molcount_origCol"]] != 0){
      grid.text(sumCountList[["MetMET2"]][["molcount_origCol"]], x = unit(0.65, "npc"), y = unit(0.65, "npc"), gp = gpar(fontsize = 5, col = "black"), vp=viewport(layout.pos.col=3:4, layout.pos.row=7))
    }
    
    if(!is.na(sumCountList[["PrioSIG"]][["molcount_origCol"]]) & sumCountList[["PrioSIG"]][["molcount_origCol"]] != 0 && !is.na(sumCountList[["PhoSIG1"]][["molcount_origCol"]]) & sumCountList[["PhoSIG1"]][["molcount_origCol"]] != 0){
      drawArrowsStraight(poscols = 3, posrows = 3:4, 
                         rellen = 0.32, center = c(0.8, 0.64), 
                         th_start = th[["PrioSIG"]], th_end = th[["PhoSIG1"]], fac = fac)
    }
    if(!is.na(sumCountList[["PrioEXP1"]][["molcount_origCol"]]) & sumCountList[["PrioEXP1"]][["molcount_origCol"]] != 0 && !is.na(sumCountList[["TraEXP1"]][["molcount_origCol"]]) & sumCountList[["TraEXP1"]][["molcount_origCol"]] != 0){
      drawArrowsStraight(poscols = 4, posrows = 3:5, rellen = 0.6, center = c(0.5, 0.96), 
                         th_start = th[["PrioEXP1"]], th_end = th[["TraEXP1"]], fac = fac)
    }
    if(!is.na(sumCountList[["PrioEXP2"]][["molcount_origCol"]]) & sumCountList[["PrioEXP2"]][["molcount_origCol"]] != 0 && !is.na(sumCountList[["TraEXP2"]][["molcount_origCol"]]) & sumCountList[["TraEXP2"]][["molcount_origCol"]] != 0){
      drawArrowsStraight(poscols = 4, posrows = 3:5, rellen = 0.6, center = c(1.2, 0.96), th_start = th[["PrioEXP2"]], th_end = th[["TraEXP2"]], fac = fac)
    }
    if(!is.na(sumCountList[["PrioEXP3"]][["molcount_origCol"]]) & sumCountList[["PrioEXP3"]][["molcount_origCol"]] != 0 && !is.na(sumCountList[["ExpEXP3"]][["molcount_origCol"]]) & sumCountList[["ExpEXP3"]][["molcount_origCol"]] != 0){
      drawArrowsStraight(poscols = 4, posrows = 3:6, rellen = 0.645, center = c(1.45, 1.33), th_start = th[["PrioEXP3"]], th_end = th[["ExpEXP3"]], fac = fac)
    }
    if(!is.na(sumCountList[["PhoSIG2"]][["molcount_2"]]) & sumCountList[["PhoSIG2"]][["molcount_2"]] != 0 && !is.na(sumCountList[["ExpSIG1"]][["molcount_origCol"]]) & sumCountList[["ExpSIG1"]][["molcount_origCol"]] != 0){
      drawArrowsStraight(poscols = 3, posrows = 4:6, rellen = 0.52, center = c(0.8, 0.93), th_start = th[["PhoSIG2"]], th_end = th[["ExpSIG1"]], fac = fac)
    }
    if(!is.na(sumCountList[["TraEXP4"]][["molcount_origCol"]]) & sumCountList[["TraEXP4"]][["molcount_origCol"]] != 0 && !is.na(sumCountList[["ExpEXP2"]][["molcount_origCol"]]) & sumCountList[["ExpEXP2"]][["molcount_origCol"]] != 0){
      drawArrowsStraight(poscols = 4, posrows = 5:6, rellen = 0.3, center = c(0.9, 0.65), th_start = th[["TraEXP4"]], th_end = th[["ExpEXP2"]], fac = fac)
    }
    if(!is.na(sumCountList[["PhoSIG3"]][["molcount_origCol"]]) & sumCountList[["PhoSIG3"]][["molcount_origCol"]] != 0 && !is.na(sumCountList[["ExpEXP1"]][["molcount_origCol"]]) & sumCountList[["ExpEXP1"]][["molcount_origCol"]] != 0){
      drawArrowsSlanted(th_start = th[["PhoSIG3"]], th_end = th[["ExpEXP1"]], 
                        rellen = 0.45, angle = 55,
                        vp = viewport(x = unit(3.2, "inch"), y = unit(1.75, "inch")), fac = fac)
    }
    if(!is.na(sumCountList[["ExpSIG3"]][["molcount_origCol"]]) & sumCountList[["ExpSIG3"]][["molcount_origCol"]] != 0 && !is.na(sumCountList[["MetMET1"]][["molcount_origCol"]]) & sumCountList[["MetMET1"]][["molcount_origCol"]] != 0){
      drawArrowsSlanted(th_start = th[["ExpSIG3"]], th_end = th[["MetMET1"]], 
                        rellen = 0.16, angle = 43,
                        vp = viewport(x = unit(2.7, "inch"), y = unit(0.88, "inch")), fac = fac)
    }
    if(!is.na(sumCountList[["ExpEXP5"]][["molcount_origCol"]]) & sumCountList[["ExpEXP5"]][["molcount_origCol"]] != 0 && !is.na(sumCountList[["MetMET2"]][["molcount_origCol"]]) & sumCountList[["MetMET2"]][["molcount_origCol"]] != 0){
      drawArrowsSlanted(th_start = th[["ExpEXP5"]], th_end = th[["MetMET2"]], 
                        rellen = 0.25, angle = -62,
                        vp = viewport(x = unit(3.7, "inch"), y = unit(0.88, "inch")), fac = fac)
    }
    #Dotted arrows
    grid.lines(x = unit(c(3.1, 2.1), "inch"),
               y = unit(c(3.55, 3.3), "inch"), 
               gp = gpar(fill="black", lty=2),
               arrow = arrow(length = unit(0.07, "inch"), ends="last", type="closed"))
    grid.lines(x = unit(c(3.1, 3.55), "inch"),
               y = unit(c(3.55, 3.3), "inch"),
               gp = gpar(fill="black", lty=2),
               arrow = arrow(length = unit(0.07, "inch"), ends="last", type="closed"))
    grid.lines(x = unit(c(3.1, 4.2), "inch"),
               y = unit(c(3.55, 3.3), "inch"),
               gp = gpar(fill="black", lty=2),
               arrow = arrow(length = unit(0.07, "inch"), ends="last", type="closed"))
    dev.off()
    
  }
  
  saveDiagram2 <- function(sumCountList, th, fac = 100, dirsave){
    fac <- fac/20
    colours <- c("#ffc732", "#6ec26e", "#f2f252", "#cef77a", "#6a90ab")
    fontsizeNorm <- 4
    fontsizeLarge <- 7
    cairo_pdf(file = paste0(dirsave, "/NetworkFiles/", prefix, "Tnet_lite", suffix, ".pdf"), width = 8.27/6, height = 11.69/4, family = "Arial")
    drawDiagram3(sumCountList, th, fac, colours, fontsizeNorm, fontsizeLarge)
    dev.off()
    
    return()
  }
  
  drawDiagram2 <- function(sumCountList, th, fac, colours){
    grid.newpage()
    pushViewport(viewport(layout=grid.layout(16, 1, widths=unit(1, "npc"), 
                                             heights=unit(c(0.065, 0.065, 0.065, 0.045, 0.065, 0.065, 0.045, 0.065, 0.065, 0.045, 0.065, 0.065, 0.045, 0.065, 0.065, 0.065), rep("npc", 11)))))
    #BG
    grid.rect(gp=gpar(fill = "#656565", lwd = NA), x = unit(0.36, "npc"), y = unit(0.2, "npc"), height=unit(0.2, "npc"), width=unit(0.41, "npc"), vp=viewport(layout.pos.col=1, layout.pos.row=1))
    grid.rect(gp=gpar(fill = "#656565", lwd = NA), x = unit(0.785, "npc"), y = unit(0.2, "npc"), height=unit(0.2, "npc"), width=unit(0.41, "npc"), vp=viewport(layout.pos.col=1, layout.pos.row=1))
    grid.rect(gp=gpar(fill = "#656565", lwd = NA), x = unit(0.427, "npc"), y = unit(0.82, "npc"), height=unit(0.2, "npc"), width=unit(0.845, "npc"), vp=viewport(layout.pos.col=1, layout.pos.row=16))
    
    grid.text("Sig", x = unit(0.35, "npc"), y = unit(0.65, "npc"), gp = gpar(fontsize = 11, fontface = "bold.italic", col = "black"), vp=viewport(layout.pos.col=1, layout.pos.row=1))
    grid.text("Exp", x = unit(0.80, "npc"), y = unit(0.65, "npc"), gp = gpar(fontsize = 11, fontface = "bold.italic", col = "black"), vp=viewport(layout.pos.col=1, layout.pos.row=1))
    grid.text("MetRes", x = unit(0.42, "npc"), y = unit(0.35, "npc"), gp = gpar(fontsize = 11, fontface = "bold.italic", col = "black"), vp=viewport(layout.pos.col=1, layout.pos.row=16))
    
    for(c in seq_along(colours)){
      grid.polygon(x=c(0.15, 0, 0.85, 1), y=c(1, 0, 0, 1), default.units="npc", 
                   gp=gpar(fill=colours[c], lwd = NA), vp = viewport(layout.pos.col = 1, layout.pos.row = ((c-1)*3+2):((c-1)*3+3)))
    }
    #Arrows
    if(!is.na(sumCountList[["PrioSIG"]][["molcount_origCol"]]) & sumCountList[["PrioSIG"]][["molcount_origCol"]] != 0 && !is.na(sumCountList[["PhoSIG1"]][["molcount_origCol"]]) & sumCountList[["PhoSIG1"]][["molcount_origCol"]] != 0){
      drawArrowsStraightNPC(poscols = 1, posrows = 4:5, rellen = 1, center = c(0.3, 0.5), th_start = th[["PrioSIG"]]*0.7, th_end = th[["PhoSIG1"]]*0.7, fac = fac)
    }
    if(!is.na(sumCountList[["PrioEXP1"]][["molcount_origCol"]]) & sumCountList[["PrioEXP1"]][["molcount_origCol"]] != 0 && !is.na(sumCountList[["TraEXP1"]][["molcount_origCol"]]) & sumCountList[["TraEXP1"]][["molcount_origCol"]] != 0){
      drawArrowsStraightNPC(poscols = 1, posrows = 4:8, rellen = 1, center = c(0.62, 0.5), th_start = th[["PrioEXP1"]]*0.7, th_end = th[["TraEXP1"]]*0.7, fac = fac)
    }
    if(!is.na(sumCountList[["PrioEXP2"]][["molcount_origCol"]]) & sumCountList[["PrioEXP2"]][["molcount_origCol"]] != 0 && !is.na(sumCountList[["TraEXP2"]][["molcount_origCol"]]) & sumCountList[["TraEXP2"]][["molcount_origCol"]] != 0){
      drawArrowsStraightNPC(poscols = 1, posrows = 4:8, rellen = 1, center = c(0.73, 0.5), th_start = th[["PrioEXP2"]]*0.7, th_end = th[["TraEXP2"]]*0.7, fac = fac)
    }
    if(!is.na(sumCountList[["PrioEXP3"]][["molcount_origCol"]]) & sumCountList[["PrioEXP3"]][["molcount_origCol"]] != 0 && !is.na(sumCountList[["ExpEXP3"]][["molcount_origCol"]]) & sumCountList[["ExpEXP3"]][["molcount_origCol"]] != 0){
      drawArrowsStraightNPC(poscols = 1, posrows = 4:11, rellen = 1, center = c(0.84, 0.5), th_start = th[["PrioEXP3"]]*0.7, th_end = th[["ExpEXP3"]]*0.7, fac = fac)
    }
    if(!is.na(sumCountList[["PhoSIG2"]][["molcount_2"]]) & sumCountList[["PhoSIG2"]][["molcount_2"]] != 0 && !is.na(sumCountList[["ExpSIG1"]][["molcount_origCol"]]) & sumCountList[["ExpSIG1"]][["molcount_origCol"]] != 0){
      drawArrowsStraightNPC(poscols = 1, posrows = 7:11, rellen = 1, center = c(0.2, 0.5), th_start = th[["PhoSIG2"]]*0.7, th_end = th[["ExpSIG1"]]*0.7, fac = fac)
    }
    if(!is.na(sumCountList[["PhoSIG3"]][["molcount_origCol"]]) & sumCountList[["PhoSIG3"]][["molcount_origCol"]] != 0 && !is.na(sumCountList[["ExpEXP1"]][["molcount_origCol"]]) & sumCountList[["ExpEXP1"]][["molcount_origCol"]] != 0){
      drawArrowsSlanted(th_start = th[["PhoSIG3"]]*0.7, th_end = th[["ExpEXP1"]]*0.7, rellen = 0.21, angle = 20, vp = viewport(x = unit(-0.85, "npc"), y = unit(0.01, "npc")), fac = fac*4)
    }
    if(!is.na(sumCountList[["TraEXP4"]][["molcount_origCol"]]) & sumCountList[["TraEXP4"]][["molcount_origCol"]] != 0 && !is.na(sumCountList[["ExpEXP2"]][["molcount_origCol"]]) & sumCountList[["ExpEXP2"]][["molcount_origCol"]] != 0){
      drawArrowsStraightNPC(poscols = 1, posrows = 10:11, rellen = 1, center = c(0.675, 0.5), th_start = th[["TraEXP4"]]*0.7, th_end = th[["ExpEXP2"]]*0.7, fac = fac)
    }
    if(!is.na(sumCountList[["ExpSIG3"]][["molcount_origCol"]]) & sumCountList[["ExpSIG3"]][["molcount_origCol"]] != 0 && !is.na(sumCountList[["MetMET1"]][["molcount_origCol"]]) & sumCountList[["MetMET1"]][["molcount_origCol"]] != 0){
      drawArrowsSlanted(th_start = th[["ExpSIG3"]]*0.7, th_end = th[["MetMET1"]]*0.7, rellen = 0.1, angle = 33, vp = viewport(x = unit(-0.74, "npc"), y = unit(-0.39, "npc")), fac = fac*4)
    }
    if(!is.na(sumCountList[["ExpEXP5"]][["molcount_origCol"]]) & sumCountList[["ExpEXP5"]][["molcount_origCol"]] != 0 && !is.na(sumCountList[["MetMET2"]][["molcount_origCol"]]) & sumCountList[["MetMET2"]][["molcount_origCol"]] != 0){
      drawArrowsSlanted(th_start = th[["ExpEXP5"]]*0.7, th_end = th[["MetMET2"]]*0.7, rellen = 0.1, angle = -33, vp = viewport(x = unit(-0.97, "npc"), y = unit(0.445, "npc")), fac = fac*4)
    }
  }
  
  drawDiagram3 <- function(sumCountList, th, fac, colours, fontsizeNorm, fontsizeLarge){
    drawDiagram2(sumCountList, th, fac, colours)
    if(!is.na(sumCountList[["PrioSIG"]][["molcount_origCol"]]) & sumCountList[["PrioSIG"]][["molcount_origCol"]] != 0){
      grid.text(sumCountList[["PrioSIG"]][["molcount_origCol"]], x = unit(0.3, "npc"), y = unit(0.36, "npc"), gp = gpar(fontsize = fontsizeNorm, col = "black"), vp=viewport(layout.pos.col=1, layout.pos.row=3))
    }
    if(!is.na(sumCountList[["PrioEXP1"]][["molcount_origCol"]]) & sumCountList[["PrioEXP1"]][["molcount_origCol"]] != 0){
      grid.text(sumCountList[["PrioEXP1"]][["molcount_origCol"]], x = unit(0.6, "npc"), y = unit(0.36, "npc"), gp = gpar(fontsize = fontsizeNorm, col = "black"), vp=viewport(layout.pos.col=1, layout.pos.row=3))
    }
    if(!is.na(sumCountList[["PrioEXP2"]][["molcount_origCol"]]) & sumCountList[["PrioEXP2"]][["molcount_origCol"]] != 0){
      grid.text(sumCountList[["PrioEXP2"]][["molcount_origCol"]], x = unit(0.72, "npc"), y = unit(0.36, "npc"), gp = gpar(fontsize = fontsizeNorm, col = "black"), vp=viewport(layout.pos.col=1, layout.pos.row=3))
    }
    if(!is.na(sumCountList[["PrioEXP3"]][["molcount_origCol"]]) & sumCountList[["PrioEXP3"]][["molcount_origCol"]] != 0){
      grid.text(sumCountList[["PrioEXP3"]][["molcount_origCol"]], x = unit(0.85, "npc"), y = unit(0.36, "npc"), gp = gpar(fontsize = fontsizeNorm, col = "black"), vp=viewport(layout.pos.col=1, layout.pos.row=3))
    }
    if(!is.na(sumCountList[["PhoSIG1"]][["molcount_origCol"]]) & sumCountList[["PhoSIG1"]][["molcount_origCol"]] != 0){
      grid.text(sumCountList[["PhoSIG1"]][["molcount_origCol"]], x = unit(0.38, "npc"), y = unit(0.5, "npc"), gp = gpar(fontsize = fontsizeNorm, col = "black"), vp=viewport(layout.pos.col=1, layout.pos.row=5))
    }
    if(!is.na(sumCountList[["PhoSIG2"]][["molcount_1"]]) & sumCountList[["PhoSIG2"]][["molcount_1"]] != 0){
      grid.text(sumCountList[["PhoSIG2"]][["molcount_1"]], x = unit(0.3, "npc"),  y = unit(0.65, "npc"), gp = gpar(fontsize = fontsizeLarge, col = "black", lineheight = 0.9), vp=viewport(layout.pos.col=1, layout.pos.row=6))
    }
    if(!is.na(sumCountList[["PhoSIG3"]][["molcount_origCol"]]) & sumCountList[["PhoSIG3"]][["molcount_origCol"]] != 0){
      grid.text(sumCountList[["PhoSIG3"]][["molcount_origCol"]], x = unit(0.35, "npc"), y = unit(0.7, "npc"), gp = gpar(fontsize = fontsizeNorm, col = "black"), vp=viewport(layout.pos.col=1, layout.pos.row=7))
    }
    if(!is.na(sumCountList[["TraEXP1"]][["molcount_origCol"]]) & sumCountList[["TraEXP1"]][["molcount_origCol"]] != 0){
      grid.text(sumCountList[["TraEXP1"]][["molcount_origCol"]], x = unit(0.58, "npc"), y = unit(0.5, "npc"), gp = gpar(fontsize = fontsizeNorm, col = "black"), vp=viewport(layout.pos.col=1, layout.pos.row=8))
    }
    if(!is.na(sumCountList[["TraEXP2"]][["molcount_origCol"]]) & sumCountList[["TraEXP2"]][["molcount_origCol"]] != 0){
      if(th["TraEXP2"] > 0.05){
        grid.text(sumCountList[["TraEXP2"]][["molcount_origCol"]], x = unit(0.67, "npc"), y = unit(0.5, "npc"), gp = gpar(fontsize = fontsizeNorm, col = "black"), vp=viewport(layout.pos.col=1, layout.pos.row=8))
      }else{
        grid.text(sumCountList[["TraEXP2"]][["molcount_origCol"]], x = unit(0.69, "npc"), y = unit(0.5, "npc"), gp = gpar(fontsize = fontsizeNorm, col = "black"), vp=viewport(layout.pos.col=1, layout.pos.row=8))
      }
    }
    if(!is.na(sumCountList[["TraEXP3"]][["molcount_origCol"]]) & sumCountList[["TraEXP3"]][["molcount_origCol"]] != 0){
      grid.text(sumCountList[["TraEXP3"]][["molcount_origCol"]], x = unit(0.67, "npc"), y = unit(0.65, "npc"), gp = gpar(fontsize = fontsizeLarge, col = "black"), vp=viewport(layout.pos.col=1, layout.pos.row=9))
    }
    if(!is.na(sumCountList[["TraEXP4"]][["molcount_origCol"]]) & sumCountList[["TraEXP4"]][["molcount_origCol"]] != 0){
      grid.text(sumCountList[["TraEXP4"]][["molcount_origCol"]], just = "left", x = unit(0.69, "npc"), y = unit(0.2, "npc"), gp = gpar(fontsize = fontsizeNorm, col = "black"), vp=viewport(layout.pos.col=1, layout.pos.row=9))
    }
    if(!is.na(sumCountList[["ExpSIG1"]][["molcount_origCol"]]) & sumCountList[["ExpSIG1"]][["molcount_origCol"]] != 0){
      grid.text(sumCountList[["ExpSIG1"]][["molcount_origCol"]], x = unit(0.2, "npc"), y = unit(0.65, "npc"), gp = gpar(fontsize = fontsizeLarge, col = "black", lineheight = 0.9), vp=viewport(layout.pos.col=1, layout.pos.row=12))
    }
    if(!is.na(sumCountList[["ExpSIG2"]][["molcount_origCol"]]) & sumCountList[["ExpSIG2"]][["molcount_origCol"]] != 0){
      grid.text(paste0(sumCountList[["ExpSIG2"]][["molcount_origCol"]]," Pho"), just = "right", x = unit(0.3, "npc"), y = unit(0.75, "npc"), gp = gpar(fontsize = fontsizeNorm, col = "black"), vp=viewport(layout.pos.col=1, layout.pos.row=13))
    }
    if(!is.na(sumCountList[["ExpSIG3"]][["molcount_origCol"]]) & sumCountList[["ExpSIG3"]][["molcount_origCol"]] != 0){
      grid.text(paste0(sumCountList[["ExpSIG3"]][["molcount_origCol"]]," Pro"), just = "right", x = unit(0.3, "npc"), y = unit(0.3, "npc"), gp = gpar(fontsize = fontsizeNorm, col = "black"), vp=viewport(layout.pos.col=1, layout.pos.row=13))
    }
    if(!is.na(sumCountList[["ExpEXP1"]][["molcount_origCol"]]) & sumCountList[["ExpEXP1"]][["molcount_origCol"]] != 0){
      grid.text(sumCountList[["ExpEXP1"]][["molcount_origCol"]], x = unit(0.5, "npc"), y = unit(0.6, "npc"), gp = gpar(fontsize = fontsizeNorm, col = "black"), vp=viewport(layout.pos.col=1, layout.pos.row=11))
    }
    if(!is.na(sumCountList[["ExpEXP2"]][["molcount_origCol"]]) & sumCountList[["ExpEXP2"]][["molcount_origCol"]] != 0){
      grid.text(sumCountList[["ExpEXP2"]][["molcount_origCol"]], just = "left", x = unit(0.69, "npc"), y = unit(0.5, "npc"), gp = gpar(fontsize = fontsizeNorm, col = "black"), vp=viewport(layout.pos.col=1, layout.pos.row=11))
    }
    if(!is.na(sumCountList[["ExpEXP3"]][["molcount_origCol"]]) & sumCountList[["ExpEXP3"]][["molcount_origCol"]] != 0){
      grid.text(sumCountList[["ExpEXP3"]][["molcount_origCol"]], x = unit(0.9, "npc"), y = unit(0.5, "npc"), gp = gpar(fontsize = fontsizeNorm, col = "black"), vp=viewport(layout.pos.col=1, layout.pos.row=11))
    }
    if(!is.na(sumCountList[["ExpEXP4"]][["molcount_origCol"]]) & sumCountList[["ExpEXP4"]][["molcount_origCol"]] != 0){
      grid.text(sumCountList[["ExpEXP4"]][["molcount_origCol"]], x = unit(0.67, "npc"), y = unit(0.65, "npc"), gp = gpar(fontsize = fontsizeLarge, col = "black"), vp=viewport(layout.pos.col=1, layout.pos.row=12))
    }
    if(!is.na(sumCountList[["ExpEXP5"]][["molcount_origCol"]]) & sumCountList[["ExpEXP5"]][["molcount_origCol"]] != 0){
      grid.text(sumCountList[["ExpEXP5"]][["molcount_origCol"]], x = unit(0.77, "npc"), y = unit(0.2, "npc"), gp = gpar(fontsize = fontsizeNorm, col = "black"), vp=viewport(layout.pos.col=1, layout.pos.row=12))
    }
    if(!is.na(sumCountList[["MetMET1"]][["molcount_origCol"]]) & sumCountList[["MetMET1"]][["molcount_origCol"]] != 0){
      grid.text(sumCountList[["MetMET1"]][["molcount_origCol"]], x = unit(0.35, "npc"), y = unit(0.4, "npc"), gp = gpar(fontsize = fontsizeNorm, col = "black"), vp=viewport(layout.pos.col=1, layout.pos.row=14))
    }
    if(!is.na(sumCountList[["MetMET3"]][["molcount_origCol"]]) & sumCountList[["MetMET3"]][["molcount_origCol"]] != 0){
      grid.text(sumCountList[["MetMET3"]][["molcount_origCol"]], y = unit(0.65, "npc"), gp = gpar(fontsize = fontsizeLarge, col = "black"), vp=viewport(layout.pos.col=1, layout.pos.row=15))
    }
    if(!is.na(sumCountList[["MetMET2"]][["molcount_origCol"]]) & sumCountList[["MetMET2"]][["molcount_origCol"]] != 0){
      grid.text(sumCountList[["MetMET2"]][["molcount_origCol"]], x = unit(0.65, "npc"), y = unit(0.4, "npc"), gp = gpar(fontsize = fontsizeNorm, col = "black"), vp=viewport(layout.pos.col=1, layout.pos.row=14))
    }
  }
  
  drawArrowsStraight <- function(poscols, posrows, rellen = 1, center = NULL, th_start, th_end, fac = 100){
    fac <- fac/2
    zoom <- 0.5
    
    #Set in inch
    if(is.null(center)){
      rown <- length(posrows)
      centx <- 1.6*zoom; centy <- 0.7*zoom*rown
    }else{
      centx <- center[1]; centy <- center[2]
    }
    pushViewport(viewport(layout.pos.col = poscols, layout.pos.row = posrows))
    grid.polygon(x=c(centx-th_start*fac/2, centx-th_end*fac/2, centx-th_end*fac/2-(0.03+th_end/2*fac), centx, 
                     centx+th_end*fac/2+(0.03+th_end/2*fac), centx+th_end*fac/2, centx+th_start*fac/2),
                 y=c(centy*(1+rellen), centy*(1-rellen)+0.07+th_end*0.45, centy*(1-rellen)+0.09+th_end*0.5, centy*(1-rellen), 
                     centy*(1-rellen)+0.09+th_end*0.5, centy*(1-rellen)+0.07+th_end*0.45, centy*(1+rellen)),
                 default.units="inch", gp=gpar(fill="black"))
    
    upViewport()
  }
  
  drawArrowsStraightNPC <- function(poscols, posrows, rellen = 1, center = NULL, th_start, th_end, fac = 1){
    centx <- center[1]
    centy <- center[2]
    pushViewport(viewport(layout.pos.col = poscols, layout.pos.row = posrows))
    grid.polygon(x=c(centx-th_start*1*fac, centx-th_end*1*fac, centx-(0.01+th_end*2*fac), centx, centx+(0.01+th_end*2*fac), centx+th_end*1*fac, centx+th_start*1*fac, centx-th_start*1*fac),
                 y=c(centy*(1+rellen), centy*(1-rellen)+0.05+th_end*0.3, centy*(1-rellen)+0.07+th_end*0.3, centy*(1-rellen), centy*(1-rellen)+0.07+th_end*0.3, centy*(1-rellen)+0.05+th_end*0.3, centy*(1+rellen), centy*(1+rellen)),
                 default.units="npc", gp=gpar(fill="black"))
    upViewport()
  }
  
  drawArrowsSlanted <- function(th_start, th_end, rellen, angle = 0, vp = viewport(x = unit(1, "npc"), y = unit(1, "npc")), fac = 100){
    pushViewport(vp)
    fac <- fac/2
    
    #Set in inch (A4: width = 11.69, height = 8.27)
    centx <- 5.845/2; centy <- 4.135/2
    grid.polygon(x=c(centx-th_start*fac/2, centx-th_end*fac/2, centx-th_end*fac/2-(0.03+th_end/2*fac), centx, 
                     centx+th_end*fac/2+(0.03+th_end/2*fac), centx+th_end*fac/2, centx+th_start*fac/2),
                 y=c(centy*(1+rellen), centy*(1-rellen)+0.07+th_end*0.45, centy*(1-rellen)+0.09+th_end*0.5, centy*(1-rellen), 
                     centy*(1-rellen)+0.09+th_end*0.5, centy*(1-rellen)+0.07+th_end*0.45, centy*(1+rellen)),
                 default.units="inch",
                 gp=gpar(fill="black"), vp=viewport(angle = angle))
    upViewport()
  }
  
  dir.create(paste0(dirsave, "/NetworkFiles"), showWarnings = FALSE)
  mdatalist <- base::get(load(paste0(dirsave, "/MdataFiles/mdata", mdata_suffix, "_list")))
  mainFun(mdatalist, mdata_suffix, arrowRelThick, dirsave, prefix, suffix, settings)
  
  return()
}
