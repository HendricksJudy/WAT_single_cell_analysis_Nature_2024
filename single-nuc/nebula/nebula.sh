
nebula_diff_gex = function(data, idents, celltype, comparison, covariates, sample_id ="sample", offset=NULL, downsample = FALSE, downsamplenumber = NULL, method = "LN", reml = 0, resid = FALSE) {
  
  # @param data = seurat object with counts and meta-data
  # @param idents = cell subset identifiers (e.g. main_clusters, or sub_cluster_ids)
  # @param celltype = subcluster of interest
  # @param comparison = vector of comparisons eg. Obese v. Lean
  # @param covariates is a vector of covariates to include in the design.
  # @param sample_id, default "sample"
  # @param offset is a vector of the column name of the scaling factor, default is NULL
  # @param downsample is TRUE/FALSE to downsample to smallest condition, default is FALSE
  # @param method is nebula method - "LN" or "HL"
  # @param when ratio n subjects to n subject-level variables is small (<10) use restricted maximum likelihood (REML, reml = 1) estimate instead of maximum likelihood estimate (MLE, reml = 0).
  # @param return residuals or not
  
  require(nebula)  

  # set up data for nebula
  
  # subset cell types and conditions for comparisions

  Idents(data) = idents
  subsetdata = subset(data, idents = celltype)
  subsetdata = subset(subsetdata, subset = (condition == comparison[1]| condition == comparison[2]))
  
  # downsample 
  
  if(downsample) {
    
    if(is.null(downsamplenumber)){
      
      print("downssample to smallest condition")
      Idents(subsetdata) = "condition"
      min = min(table(Idents(subsetdata)))
      subsetdata = subset(subsetdata, downsample = min)
      
    }
    
    if(is.numeric(downsamplenumber)){
      
      print("downssample to selected number")
      #print(downsamplenumber)
      Idents(subsetdata) = "condition"
      subsetdata = subset(subsetdata, downsample = downsamplenumber)
      print(table(subsetdata$condition))
    }
  }
  
  pcts = FindMarkers(subsetdata, ident.1 = comparison[1], ident.2 = comparison[2], group.by="condition", logfc.threshold = 0,  min.pct = 0, only.pos = FALSE)[,c("pct.1", "pct.2")]
  pcts$pct_delta = pcts$pct.1 - pcts$pct.2
  colnames(pcts) = c(paste0("pct_",comparison[1]),paste0("pct_",comparison[2]),"pct_delta")
  pcts$gene = rownames(pcts)

  # set up data in nebula format

  count = round(subsetdata@assays$RNA@counts)
  id = unlist(subsetdata[[sample_id]])
  covariates = c(covariates,"condition")
  variables = subsetdata[[covariates]]

  if("eth" %in% colnames(variables)) {variables$eth = factor(variables$eth) } 
  if("sex" %in% colnames(variables)) {variables$sex = factor(variables$sex) } 

  variables[] <- lapply(variables, function(x) if(is.factor(x)) factor(x) else x) 	# remove missing factor levels after subsetting

  if(is.null(offset) == FALSE) { offs = unlist(subsetdata[[offset]]) }
  if(is.null(offset) == TRUE) { offs = NULL }  
  design = as.formula( paste(" ~ ",paste(covariates, collapse=" + "), sep="") )
  design = model.matrix(design, data = variables)
  
  # order data
  ordered_data = nebula::group_cell(count=count,id=id,pred=design,offset=offs)
  
  # run nebula ln and hl
  
  print("run nebula for celltype")
  results = nebula::nebula(count = ordered_data$count, id = ordered_data$id, pred = ordered_data$pred, offset = ordered_data$offs, method = method, reml = reml)
  results$pct = pcts
  
  if(resid) {
    
    residuals = nbresidual(results,count = ordered_data$count, id = ordered_data$id, pred = ordered_data$pred, offset = ordered_data$offs)
    results = list(results, residuals)
    names(results) = c("results", "residuals")
    return(results)
    
  } else {
    
    return(results)
    
  }

}
