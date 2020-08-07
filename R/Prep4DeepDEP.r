#' Main function of Prep4DeepDEP
#'
#'@description Prep4DeepDEP generates the genomic and gene fingerprint data tables from user's datasets. Prep4DeepDEP has two main modes:
#'\itemize{
#'\item The “Prediction” mode generates data for DeepDEP to predict gene dependency scores of unscreened CCLs or tumors. It extracts and orders the required genomic features from user’s genome-wide datasets, and generates the functional fingerprints of gene dependencies of interest (DepOIs) from a user-provided list or the default 1,298 genes we studied in the paper. For the copy number alteration (CNA) data, an embedded R function (PrepCNA) converts copy-number segments to bins (every 10k bases in the genome) and calculate per-bin CNA scores.
#'\item The “Training” mode generates data to train a new DeepDEP model using user's genome-wide genomic data and gene dependency scores from an in-house CRISPR screening experiment. It creates data tables of genomics and gene dependencies for all CCL-DepOI pairs (number of samples = number of CCLs x number of DepOIs). Functional fingerprints are generated based on the list of genes available in the gene dependency dataset.}
#'Please refer to the paper and DeepDEP package (https://codeocean.com/capsule/3348251/tree) about how to use the generated data tables for DeepDEP model training and prediction.
#'
#'@param exp.data Gene expression data (a data.frame object) of cell lines or tumors. Rows and columns of the data frame correspond to genes and samples, respectively. The data frame should contain sample names as column names and gene symbols (e.g., CCND1) as the first column. Row names are not used by this function. Expression levels are presented by log2(TPM+1) per gene.
#'@param mut.data Mutation data (a data.frame object) of cell lines or tumors. Rows and columns of the data frame correspond to genes and samples, respectively. The data frame should contain sample names as column names and gene symbols (e.g., TP53) as the first column. Row names are not used by this function. Mutations are represented by 0/1 binary values per gene, with 1s denoting missense and nonsense mutations, frameshift insertions and deletions, and splice-site mutations.
#'@param meth.data DNA methylation data (a data.frame object) of cell lines or tumors. Rows and columns of the data frame correspond to probes and samples, respectively. The data frame should contain sample names as column names and probe ID (e.g., cg00000292) as the first column. Row names are not used by this function. DNA methylation is measured by beta values per probe of Infinium® HumanMethylation27 or HumanMethylation450 BeadChips.
#'@param cna.data Copy number alteration (CNA) data (a data.frame object) of cell lines or tumors. CNA should be prepared as segmented copy-number profiles using the .seg file format against the reference genome hg19. Example of the CNA data can be downloaded from the CCLE portal (https://portals.broadinstitute.org/ccle/data). Rows and columns of the data frame correspond to CNA segments per sample and CNA information, respectively. The following columns are required: CCLE_name (sample name), Chromosome (numeric without ‘Chr’), Start (numeric), End (numeric), and Segment_Mean (in the log2(CN/2) scale).
#'@param dep.data Gene symbols of dependency genes of interest (DepOIs) with or without user’s in-house gene dependency scores.
#' For the “Training” mode, this argument is required and expects a data.frame object of which rows and columns correspond to DepOIs and samples, respectively. The data frame should contain sample names as column names and gene symbol (e.g., TP53) as the first column.
#' For the “Prediction” mode, this argument is optional and expects a data.frame object with a single column of gene symbols (e.g., TP53) of DepOIs that user would like to predict. If the argument is left NULL, the 1298 default genes as studied in the original paper will be used.

#'@param mode “Training” or “Prediction”.
#' The “Training” mode creates data tables of genomics and gene dependencies for all CCL-DepOI pairs (number of samples = number of CCLs x number of DepOIs). Functional fingerprints are generated based on the list of genes of “dep.data”.
#' The ‘Prediction’ mode generates data tables of genomics for all samples (number of samples = number of CCLs/tumors). Functional fingerprints are generated based on the genes of “dep.data”.

#'@param filename.out Path and prefix for the output files.
#'
#'@details
#'\itemize{
#'\item For each genomic data, Prep4DeepDEP extracts and orders the genomic features that are required to run the Python DeepDEP tool (4539 mutations, 6016 gene expressions, 7460 CNA bins, and 6617 methylation probes). At least one of the four genomics should be provided in order to run Prep4DeepDEP. If multiple genomic profiles are provided, only the samples contained in the first genomic profile will be analyzed across the provided genomic profiles. Please make sure the sample names (CCL name or tumor ID) are consistent across genomic profiles and dep.data.
#'\item For each DepOI, Prep4DeepDEP generates the binary status of 3115 functional fingerprints based on the chemical and genetic perturbation (CGP) signatures of the MSigDB database (https://www.gsea-msigdb.org/gsea/msigdb).
#'\item Output files of the function: A txt file for each genomic profile is written to the path/filename indicated by filename.out. Another txt file is outputted for the gene fingerprints. If running in the Training mode, a dependency score file is also outputted that contains the reformatted dependency scores from user-provided dep.data.
#'\item In the “training” mode, all output samples are generated by the CCL-DepOI combinations; i.e., samples are C1G1 (CCL1-DepOI1), C1G2, …, C2G1, C2G2, …
#'\item Missing values in exp.data are filled by the mean value of the corresponding genomic feature in the 278 CCLs used in our study. Missing values in mut.data are filled by the median of CCLs. Missing values in meth.data are filled by zero.}
#'
#'@note The “Prediction” mode can be slow and memory-heavy if huge numbers of samples and DepOIs are provided since the generated data and output files have the sample size of #samples x #DepOIs.
#'
#'@export

Prep4DeepDEP <- function(exp.data = NULL, mut.data = NULL,
                         meth.data = NULL,cna.data  = NULL,
                         dep.data = NULL, mode = c("Training","Prediction"),
                         filename.out = "data_out"){

check.cellNames <- NULL
cat("Mode",mode,"\n\n")
path <- system.file("extdata/",package = "Prep4DeepDEP")

#Check file
if(sum(is.null(exp.data),is.null(mut.data),is.null(meth.data),is.null(cna.data)) == 4){
  stop(c("All genomic profiles are missing. Please provide at least one of mut.data, exp.data, meth.data, and cna.data."))
}

if(is.null(dep.data) & tolower(mode) == "prediction") {
  cat("dep.data is not provided, running with the default 1298 DepOIs...","\n")
  load(paste0(path,"default_dep_genes_1298.RData"))
}

if(is.null(dep.data) & tolower(mode) == "training") {
  cat("dep.data is not provided. Please provide gene dependency scores for the training mode...", "\n")
}

if(ncol(dep.data) == 1 & tolower(mode) == "training"){
  stop(c("Only one column detected in dep.data. Please provide gene dependency symbols and scores for the training mode."),call. = FALSE)
}

#Check gene symbol
load(paste0(path,"gene_fingerprints_CGP.RData"))
list.genes <- .CheckGeneSymbol(dep.data = dep.data,filename.out= filename.out)
n <- nrow(list.genes)

#Expression
if(!is.null(exp.data)){
  if (!is.character(exp.data[1,1])) {
    stop("exp.data format error, please check!", call. = FALSE)
  }

  cat(c("Exp started..."),"\n")
  colnames(exp.data)[1] <- "Gene"
  ncell <- ncol(exp.data[,-1])
  if(is.null(check.cellNames)){
    check.cellNames <- colnames(exp.data[,-1])
  }else if(length(check.cellNames) != length( colnames(exp.data[,-1])) |
           sum(check.cellNames %in% colnames(exp.data[,-1])) != length(check.cellNames)){
    stop(c("Cell line names are inconsistent!"),call. = FALSE)
  }

  cat("Precessing",paste(length(check.cellNames),"cell lines..."),"\n")

  inputData <- exp.data[!duplicated(exp.data$Gene),]

  load(paste0(path,"ccle_exp_for_missing_value_6016.RData"))

  outputData <- merge(exp.index,inputData, by = "Gene", sort = FALSE, all.x = TRUE )
  Gene <- outputData$Gene
  rownames(outputData) <- outputData$Gene
  value_NA <-  rowSums(outputData[,-c(1,2)])

  cat(sum(is.na(value_NA)),"genes with NA values in exp.data. Substitute by mean values of CCLE.","\n")
  if(round((sum(is.na(value_NA))/nrow(outputData)),digits = 2)> 0.2){
    warning("NA found in >20% genes, please check if input format is correct!")
  }

  for (i in 1:nrow(outputData)) {
    if(is.na(value_NA[i])){
      outputData[i,is.na(outputData[i,])] <- outputData$Mean[i]
    }

  }
  outputData <- round(as.matrix(outputData[,-1]), digits = 4)
  outputData.final.exp  <- cbind(Gene,as.data.frame( outputData[,-1],stringsAsFactors = FALSE))
  outputData.final.exp  <- outputData.final.exp[,c("Gene", check.cellNames)]

    if(tolower(mode) == "prediction"){
    write.table(outputData.final.exp, file = paste(filename.out,"exp_prediction.txt",sep = "_"),
                sep = "\t",row.names = FALSE,col.names = TRUE, quote = FALSE)
    }

    if (tolower(mode) == "training") {
      k = 2
      rep_col <- do.call("cbind", replicate(n,outputData.final.exp[,2], simplify = FALSE))
      colnames(rep_col) <- paste0("C1G",seq(1,n,1))
      if(ncol(outputData.final.exp) >= 3){
      for (i in 3:ncol(outputData.final.exp)) {
        rep_col.1 <- do.call("cbind", replicate(n,outputData.final.exp[,i],
                                                simplify = FALSE))
        colnames(rep_col.1) <- paste0("C",k,"G",seq(1,n,1))
        k = k+1
        rep_col <- cbind(rep_col,rep_col.1)
      }

      rep_col <- cbind(rownames(outputData.final.exp),as.data.frame(rep_col,stringsAsFactors = FALSE))
      colnames(rep_col)[1] <- "Gene"
    }

      write.table(rep_col, file = paste(filename.out,"exp_training.txt",sep = "_"),
                  sep = "\t",row.names = FALSE,col.names = TRUE, quote = FALSE)

    }
  cat("Exp completed!","\n\n")
}


#Mutation
if(!is.null(mut.data)){
  if (!is.character(mut.data[1,1])) {
    stop("mut.data format error, please check!", call. = FALSE)
  }
  cat(c("Mut started..."),"\n")

  colnames(mut.data)[1] <- "Gene"
  ncell <- ncol(exp.data[,-1])
  if(is.null(check.cellNames)){
    check.cellNames <- colnames(mut.data[,-1])
  }else if(length(check.cellNames) != length( colnames(mut.data[,-1])) |
    sum(check.cellNames %in% colnames(mut.data[,-1])) != length(check.cellNames)){
    stop(c("Cell line names are inconsistent!"),call. = FALSE)
  }
  cat("Precessing",paste(length(check.cellNames),"cell lines..."),"\n")

  inputData  <- mut.data[!duplicated(mut.data$Gene),]

  load(paste0(path,"ccle_mut_for_missing_value_4539.RData"))

  outputData <- merge(mut.index,inputData, by = "Gene", sort = FALSE, all.x = TRUE )
  Gene <- outputData$Gene
  rownames(outputData) <- outputData$Gene
  value_NA <-  rowSums(outputData[,-c(1,2)])
  cat(sum(is.na(value_NA)),"genes with NA values in mut.data. Substitute by median values of CCLE.","\n")
  if(round((sum(is.na(value_NA))/nrow(outputData)),digits = 2)> 0.2){
    warning("NA found in >20% genes, please check if input format is correct!")
  }

  for (i in 1:nrow(outputData)) {
    if(is.na(value_NA[i])){
      outputData[i,is.na(outputData[i,])] <- outputData$Median[i]
    }

  }

  outputData.final.mut  <- cbind(Gene, as.data.frame(outputData[,-c(1,2)]))
  outputData.final.mut  <- outputData.final.mut[,c("Gene",check.cellNames)]
    if(tolower(mode) == "prediction"){
    write.table(outputData.final.mut, file = paste(filename.out,"mut_prediction.txt",sep = "_"),
                sep = "\t",row.names = FALSE,col.names = TRUE, quote = FALSE)
    }
  if (tolower(mode) == "training") {
    k = 2
    rep_col <- do.call("cbind", replicate(n,outputData.final.mut[,2]
                                          , simplify = FALSE))
    colnames(rep_col) <- paste0("C1G",seq(1,n,1))

    if(ncol(outputData.final.mut) >= 3){
    for (i in 3:ncol(outputData.final.mut)) {
      rep_col.1 <- do.call("cbind", replicate(n,outputData.final.mut[,i],
                                              simplify = FALSE))
      colnames(rep_col.1) <- paste0("C",k,"G",seq(1,n,1))
      k = k+1
      rep_col <- cbind(rep_col,rep_col.1)
    }

    rep_col <- cbind(rownames(outputData.final.mut),as.data.frame(rep_col,stringsAsFactors = FALSE))
    colnames(rep_col)[1] <- "Gene"
  }

    write.table(rep_col, file = paste(filename.out,"mut_training.txt",sep = "_"),
                sep = "\t",row.names = FALSE,col.names = TRUE, quote = FALSE)
  }
  cat("Mut completed!","\n\n")
  }


#Methylation
if(!is.null(meth.data)){
  if (!is.character(meth.data[1,1])) {
    stop("meth.data format error, please check!", call. = FALSE)
  }
  cat(c("Meth started..."),"\n")

  colnames(meth.data)[1] <- "Probe"
  ncell <- ncol(exp.data[,-1])
  if(is.null(check.cellNames)){
    check.cellNames <- colnames(meth.data[,-1])
  }else if(length(check.cellNames) != length( colnames(meth.data[,-1])) |
    sum(check.cellNames %in% colnames(meth.data[,-1])) != length(check.cellNames)){
    stop(c("Cell line names are inconsistent!"),call. = FALSE)
  }
  cat("Precessing",paste(length(check.cellNames),"cell lines..."),"\n")

  inputData  <- meth.data[!duplicated(meth.data$Probe),]

  load(paste0(path,"ccle_meth_for_missing_value_6617.RData"))

  outputData <- merge(meth.index,inputData, by = "Probe", sort = FALSE, all.x = TRUE )
  Probe <- outputData$Probe
  rownames(outputData) <- outputData$Probe
  value_NA <-  rowSums(outputData[,-c(1,2)])

  cat(sum(is.na(value_NA)),"genes with NA values in meth.data. Substitute by 0.","\n")
  if(round((sum(is.na(value_NA))/nrow(outputData)),digits = 2)> 0.2){
    warning("NA found in >20% genes, please check if input format is correct!")
  }

  for (i in 1:nrow(outputData)) {
    if(is.na(value_NA[i])){
      if(sum(is.na(outputData[i,])) == sum(ncol(outputData)-2)){
        # No values on all cell lines
        #stop(c(paste("Missing value",rownames(outputData)[i], sep = ":"),"\n"),call. = FALSE)
        outputData[i,is.na(outputData[i,])] <- 0
      }else{
        #No values on part of cell lines
      outputData[i,is.na(outputData[i,])] <- 0
      }
    }
  }

  outputData <- round(as.matrix(outputData[,-1]), digits = 4)
  outputData.final.meth  <- cbind(Probe, as.data.frame(outputData[,-1]))
  outputData.final.meth  <- outputData.final.meth[,c("Probe",check.cellNames)]

  if(tolower(mode) == "prediction"){
  write.table(outputData.final.meth, file = paste(filename.out,"meth_prediction.txt",sep = "_"),
              sep = "\t",row.names = FALSE,col.names = TRUE, quote = FALSE)
  }

  if (tolower(mode) == "training") {
    k = 2
    rep_col <- do.call("cbind", replicate(n,outputData.final.meth[,2]
                                          , simplify = FALSE))
    colnames(rep_col) <- paste0("C1G",seq(1,n,1))

    if(ncol(outputData.final.meth) >= 3){
    for (i in 3:ncol(outputData.final.meth)) {
      rep_col.1 <- do.call("cbind", replicate(n,outputData.final.meth[,i],
                                              simplify = FALSE))
      colnames(rep_col.1) <- paste0("C",k,"G",seq(1,n,1))
      k = k+1
      rep_col <- cbind(rep_col,rep_col.1)
    }

    rep_col <- cbind(rownames(outputData.final.meth),as.data.frame(rep_col,stringsAsFactors = FALSE))
    colnames(rep_col)[1] <- "Probe"
  }
    write.table(rep_col, file = paste(filename.out,"meth_training.txt",sep = "_"),
                sep = "\t",row.names = FALSE,col.names = TRUE, quote = FALSE)
  }
  cat("Meth completed!","\n\n")
}


#Fingerprints
if(!is.character(dep.data[1,1])){
  stop("dep.data format error, please check!", call. = FALSE)
}else{
  cat(c("Fingerprint started..."),"\n")
  #load(paste0(path,"gene_fingerprints_CGP.RData"))

  colnames(dep.data)[1] <- "Gene"
  idx <- which(fingerprint[1,] %in% c("GeneSet",list.genes$Gene))
  outputData <- fingerprint[,idx]


  if(tolower(mode) == "prediction"){
    write.table(outputData, file = paste(filename.out,"fingerprint_prediction.txt",sep = "_"),
              sep = "\t",row.names = FALSE,col.names = FALSE, quote = FALSE)
  }

  if(tolower(mode) == "training"){
    outputData.train <- cbind(outputData[,1],
                                  do.call("cbind", replicate(ncell,outputData[,-1], simplify = FALSE)))
    outputData.train[1,] <- c("GeneSet",paste0(paste0("C",rep(seq(1,ncell,1),each = n)),"G",seq(1,n,1)))

    write.table(outputData.train, file = paste(filename.out,"fingerprint_training.txt",sep = "_"),
                sep = "\t",row.names = FALSE,col.names = FALSE, quote = FALSE)
  }

cat("Fingerprint completed!","\n\n")
}

#Dependency scores (training only)
if (!is.null(dep.data) & tolower(mode) == "training") {
  if(is.null(check.cellNames)){
    check.cellNames <- colnames(dep.data[,-1])
  }else if(length(check.cellNames) != length( colnames(dep.data[,-1])) |
    sum(check.cellNames %in% colnames(dep.data[,-1])) != length(check.cellNames)){
    stop(c("Cell line names are inconsistent!"),call. = FALSE)
  }

   cat("Gene dependency scores (training mode) start...","\n")
   crispr.input <- dep.data[which(dep.data$Gene %in% list.genes$Gene),
                            which(colnames(dep.data) %in% c("Gene", check.cellNames))]
   crispr.output <- t(crispr.input[,2])
   colnames(crispr.output) <- paste0("C1G",seq(1,n,1))
   k=2
   for (i in 3:ncol(crispr.input)) {
     table <- t(crispr.input[,i])
     colnames(table) <- paste0("C",k,"G",seq(1,n,1))
     k = k+1
     crispr.output <- cbind(crispr.output,table)
   }

   crispr.output <- cbind("score",crispr.output)
   colnames(crispr.output)[1] <- "Dep_Score"
   write.table(crispr.output,file = paste(filename.out,"DepScore_training.txt",sep = "_"),
               sep = "\t",row.names = FALSE,col.names = TRUE, quote = FALSE)

   cat("Gene dependency scores (training) completed!","\n\n")
}

#Copy Number alteration
if(!is.null(cna.data)){
  if (sum(colnames(cna.data) %in% c("CCLE_name","Chromosome","Start","End","Num_Probe", "Segment_Mean"))!=5) {
    stop("cna.data format error, please check!", call. = FALSE)
  }
  cat(c("CNA started..."),"\n")
  ncell <- length(unique(cna.data$CCLE_name))

  if(is.null(check.cellNames)){
  outputData.cna <- .PrepCNA(cna.original = cna.data,filename.out ,exportTable = FALSE)
  }else{
    idx <- which(cna.data$CCLE_name %in% check.cellNames)
    if(length(check.cellNames) != length(unique(cna.data$CCLE_name[idx])) |
       sum(check.cellNames %in% unique(cna.data$CCLE_name[idx])) != length(check.cellNames)){
      stop(c("Cell line names are inconsistent!"),call. = FALSE)
    }
    outputData.cna <- .PrepCNA(cna.original = cna.data[idx,], filename.out ,exportTable = FALSE)
    outputData.cna <- outputData.cna[,c("CNA",check.cellNames)]
  }

  if(tolower(mode) == "prediction"){
    colnames(outputData.cna)[1] <- "Bin"
    write.table(outputData.cna, file = paste(filename.out,"cna_prediction.txt",sep = "_"),
                sep = "\t",row.names = FALSE,col.names = TRUE, quote = FALSE)
    }

  if (tolower(mode) == "training") {
    k = 2
    rep_col <- do.call("cbind", replicate(n,outputData.cna[,2]
                                          , simplify = FALSE))
    colnames(rep_col) <- paste0("C1G",seq(1,n,1))

    if(ncol(outputData.cna) >=3 ){
      for (i in 3:ncol(outputData.cna)) {
        rep_col.1 <- do.call("cbind", replicate(n,outputData.cna[,i],
                                                simplify = FALSE))
        colnames(rep_col.1) <- paste0("C",k,"G",seq(1,n,1))
        k = k+1
        rep_col <- cbind(rep_col,rep_col.1)
      }
      rep_col <- cbind(outputData.cna$CNA,as.data.frame(rep_col,stringsAsFactors = FALSE))
      colnames(rep_col)[1]<- "Bin"
    }

    write.table(rep_col, file = paste(filename.out,"cna_training.txt",sep = "_"),
                sep = "\t",row.names = FALSE,col.names = TRUE, quote = FALSE)
  }
    cat("CNA completed!","\n\n")
  }
}
