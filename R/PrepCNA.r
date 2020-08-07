#'Genomic binning and calculation of CNA scores
#'
#'@description PrepCNA extracts 7460 bins (each with 10k bases) required by Prep4DeepDEP and each CNA (in log scale) is mapped onto these bins weighted by the percentage covered by CNA. The function is called by the main function of this package, Prep4DeepDEP.

#'@export
.PrepCNA <- function(cna.original, filenames,exportTable = FALSE){
#Cell line
cellLine     <- unique(cna.original$CCLE_name)

#chr X and Y ---> chr 23 & 24
if(sum(tolower(colnames(cna.original)) %in% c("chr","chromosome")) == 1){

  col_idx <- which(tolower(colnames(cna.original)) %in% c("chr","chromosome"))
  cna.original[which(tolower(cna.original[,col_idx]) == "x"),col_idx] <- 23
  cna.original[which(tolower(cna.original[,col_idx]) == "y"),col_idx] <- 24
  colnames(cna.original)[col_idx] <- "Chromosome"
}
chr                                                     <- unique(cna.original[,col_idx])

#length of chr.
chrom_length <- c(249260000,243200000,198030000,191160000,180920000,171120000,
                 159140000,146370000,141220000,135540000,135010000,133860000,
                 115170000,107350000,102540000,90360000,81200000,78080000,59130000,
                 63030000,48130000,51310000,155280000,59380000)
chrom_bin    <- ceiling(chrom_length/10^4)


bigTable <- data.frame(matrix(data = 0, ncol = length(cellLine)+1,
                              nrow = sum(chrom_bin)), stringsAsFactors = FALSE)

colnames(bigTable) <- c("CNA", cellLine)
#saveTable <- bigTable
#bigTable <- saveTable
k = 1
for (i in 1:length(chr)) {
    bin_start <- seq(0, chrom_bin[i]-1,1)
    bin_end   <- seq(1, chrom_bin[i],1)
    bigTable$CNA[k:(k+chrom_bin[i]-1)] <- paste(paste0("chr",i),paste0(bin_start, "to",bin_end),
                                          "10k", sep = "_")
    k = chrom_bin[i]+k

}


#Filter
t1 <- proc.time()
#Default: 7460 cna
path <- system.file("extdata/",package = "Prep4DeepDEP")
load(paste0(path,"cna_table_7460.RData"))

#lst loop: chromosome
for (i in unique(filterTable$Chr)) {
  cna.filterTable <- filterTable[which(filterTable$Chr == i),]
  idx.chr <- which(cna.original$Chromosome == i)
  Table.chr <- cna.original[idx.chr,]
  cna.length.chr <- (Table.chr$End - Table.chr$Start)+1

  if(i == 1){
    l=0
  }else{
  l <- sum(chrom_bin[1:(i-1)])
  }

  #2nd loop: cell line
  for (j in 2:ncol(bigTable)) {
    idx.cell <- which(Table.chr$CCLE_name == colnames(bigTable)[j])
    cellTable.chr <- Table.chr[idx.cell,]
    cna.length <- cna.length.chr[idx.cell]


    #3rd: bin
    for (k in cna.filterTable$End) {

      end_matrix <- data.frame(matrix(data = 1,nrow = nrow(cellTable.chr))*10^4*k,
                               stringsAsFactors = FALSE)
      end_matrix$cellTable <- cellTable.chr$End
      start_matrix   <- data.frame(matrix(data = 1,nrow = nrow(cellTable.chr))*10^4*(k-1)+1,
                                   stringsAsFactors = FALSE)
      start_matrix$cellTable <- cellTable.chr$Start


      overlap.length <- (10^4) + cna.length -
        (apply(end_matrix,1,max) - apply(start_matrix,1,min)+1)

      overlap.length[overlap.length < 0] <- 0

      bigTable[l+k,j] <- round(sum(overlap.length*cellTable.chr$Segment_Mean)/10^4,
                               digits = 4)

    }

   # cat(paste("Cell line:", colnames(bigTable)[j],"Chr",i,"done",sep = " "),"\n")


  }

  #cat(paste("All Chr",i,"done",sep = " "),"\n")


}
#a <- bigTable[which(substr(bigTable$CNA,1,4) == "chr7"),1:4]
bigTable.filter <- merge(filterTable, bigTable, by = "CNA", all.x = TRUE, sort = FALSE)
bigTable        <- bigTable.filter[,-c(2:4)]

if(exportTable == TRUE){
write.table(bigTable, file = paste(filenames,nrow(bigTable.filter),"_CNA_filter.txt",sep = "_"),
            sep = "\t", col.names = TRUE,row.names = FALSE,quote = FALSE)
}

t2 <- proc.time()
t  <- round((t2-t1)/60,digits = 2)

print(c("Computation time (mins)",t[1:3]))

return(bigTable)
}
