#'Check DepOIs in the fingerprint database
#'
#'@export

.CheckGeneSymbol <- function(dep.data,filename.out){
  path <- system.file("extdata/",package = "Prep4DeepDEP")

  load(paste0(path,"gene_fingerprints_CGP.RData"))
  colnames(dep.data)[1] <- "Gene"

  if(length(intersect(dep.data$Gene, fingerprint[1,])) != length(dep.data$Gene)){

    list.genes <- data.frame(Gene = dep.data[which(dep.data$Gene %in% intersect(dep.data$Gene, fingerprint[1,])),1],stringsAsFactors = FALSE)
    list.gene.exclude <- dep.data[which(!dep.data$Gene %in% intersect(dep.data$Gene, fingerprint[1,])),1]


    list.gene.exclude <- data.frame(GeneList = list.gene.exclude, stringsAsFactors = FALSE)
    write.table(list.gene.exclude, file = paste(filename.out,"fingerprint",nrow(list.gene.exclude),"excluded_genes.txt",sep = "_"),
                sep = "\t",row.names = FALSE,col.names = TRUE, quote = FALSE)
    #n.list.genes.exclude <- nrow(list.gene.exclude)
    #n.list.genes        <- nrow(list.genes)
  }else{

    list.genes <- data.frame( Gene = dep.data[,1],stringsAsFactors = FALSE)

  }

  return(list.genes)
}
