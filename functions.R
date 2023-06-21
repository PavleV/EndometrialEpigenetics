# load libraries

library(tidyverse)
library(ggplot2)
library(GenomicRanges)

# load data

RNA_data_df <- readRDS("./Data/ALL_avg_celltype_RNA.rds")
Peak_data_df <- readRDS("./Data/ALL_avg_celltype_peaks.rds")
Key_GRCh38_p13_df <- readRDS("./Data/Key.GRCh38.p13.rds")

geneKey.ranges <- GRanges(seq=Key_GRCh38_p13_df$Chromosome,IRanges(start=Key_GRCh38_p13_df$Start, end=Key_GRCh38_p13_df$End), strand=Key_GRCh38_p13_df$Strand,mcols=Key_GRCh38_p13_df[,c("EnsemblID","GeneName","Description")])

#load ATAC count matrix and genomic coordinates

Peak_coords <- as.data.frame(str_split(row.names(Peak_data_df),pattern="-", simplify = T))
Peak_coords[,2] <- as.numeric(Peak_coords[,2])
Peak_coords[,3] <- as.numeric(Peak_coords[,3])

AllPeaks.granges <- GRanges(seqnames = Peak_coords[,1], ranges = IRanges(start=Peak_coords[,2], end=Peak_coords[,3]), mcols = data.frame(PeakID = row.names(Peak_data_df)))


tpm3 <- function(counts,len) {
  x <- counts/len
  return(t(t(x)*1e6/colSums(x)))
}

ATAC_TPMmatrix <- as.data.frame(tpm3(Peak_data_df,width(AllPeaks.granges)))




# Functions

# subset data according to gene name or ensembl ID
subset.RNA.FUN <- function(mydata=AH_EL_RNA_ALLREPS, geneName=NULL, ensemblID=NULL){

  if(!is.null(geneName) & !is.null(ensemblID)){
    return(subset(mydata, GeneName %in% geneName & EnsemblID == ensemblID))
  }
  if(!is.null(geneName) & is.null(ensemblID)){
    return(subset(mydata, GeneName %in% geneName))
  }
  if(is.null(geneName) & !is.null(ensemblID)){
    return(subset(mydata, EnsemblID == ensemblID))
  }

}

#arrange subsetted data for plotting

arrange.FUN <- function(mydata, geneName=NULL, ensemblID=NULL){

  mydata.subset <- subset.RNA.FUN(mydata=mydata, geneName = geneName, ensemblID=ensemblID)
  mydata.subset <-  pivot_longer(mydata.subset, cols = colnames(mydata)[2:40])

  mydata.subset <- cbind(mydata.subset,str_split(mydata.subset$name,pattern="_", simplify = T))
  colnames(mydata.subset)[3:6] <- c("Sample","TPM","Biopsy","Hours")

  return(mydata.subset)

}




plotRNA.FUN <- function(mydata, geneName=NULL, ensemblID=NULL, biopsies = c("S169","S170","S506","S508"), times = c(0,3,6,9,12,18,24,36,48,96)){

  if(is.null(geneName) & is.null(ensemblID)){
    return(
      print(
        ggplot()+theme_bw()
      )
    )
  }

  data.sub <- arrange.FUN(mydata=mydata, geneName=geneName, ensemblID=ensemblID)

  data.sub <- subset(data.sub, Biopsy %in% biopsies & Hours %in% times)

  if(length(geneName) == 1){

    p1 <- ggplot(subset(data.sub))+
      geom_line(aes(x=as.numeric(Hours),y=TPM,colour=Biopsy), size = 2)+
      stat_summary(aes(x=as.numeric(Hours),y=TPM), geom = "line", fun.y = mean, size = 2)+
      theme_bw()+theme(axis.title = element_text(size=20), axis.text = element_text(size=20), legend.text= element_text(size=15), title= element_text(size=20) )+
      scale_x_continuous(breaks=times, minor_breaks = NULL)+
      labs(y="Transcripts per Million", x="Hours")+ggtitle("RNA-seq")

  }

  if(length(geneName) >= 2){

    p1 <- ggplot(subset(data.sub))+
      stat_summary(aes(x=as.numeric(Hours),y=TPM,colour=GeneName), geom = "line", fun.y = mean, size = 2)+
      theme_bw()+theme(axis.title = element_text(size=20), axis.text = element_text(size=20), legend.text= element_text(size=15), title= element_text(size=20) )+
      scale_x_continuous(breaks=times, minor_breaks = NULL)+
      labs(y="Transcripts per Million", x="Hours")+ggtitle("RNA-seq")

  }


  print(p1)
}



#load ATAC count matrix and genomic coordinates

ATAC_countsmatrix_cleaned <- read.csv("./Data/ATAC_countsmatrix_cleaned_230616.csv", row.names=1)

Alex_ATAC_peaks_hg19_201110 <- read.delim("./Data/Alex_ATAC_peaks_hg19_201110.bed", header=FALSE)

AllPeaks.granges <- GRanges(seqnames = Alex_ATAC_peaks_hg19_201110[,1], ranges = IRanges(start=Alex_ATAC_peaks_hg19_201110[,2], end=Alex_ATAC_peaks_hg19_201110[,3]), mcols = data.frame(PeakID = Alex_ATAC_peaks_hg19_201110[,4]))


tpm3 <- function(counts,len) {
  x <- counts/len
  return(t(t(x)*1e6/colSums(x)))
}

ATAC_TPMmatrix <- as.data.frame(tpm3(ATAC_countsmatrix_cleaned,width(AllPeaks.granges)))

# function for extracting genomic coordinates from text input and generating genomic ranges object

extractCoord <- function(string.input){

  coordinates <- str_split_1(string.input, pattern=regex("[:,_-[:space:]]"))

  if (length(coordinates) != 3 || any(is.na(as.numeric(coordinates[2:3])))) {
    stop("Invalid genomic region provided. Please provide a valid genomic region in the format 'chr:start-end'.")
  }

  return(GRanges(seqnames = coordinates[1], ranges = IRanges(start=as.numeric(coordinates[2]), end=as.numeric(coordinates[3]))))


}

# extend coordinates

grexpand <- function(x, upstream=0, downstream=0){
  if (any(strand(x) == "*")){
    warning("'*' ranges were treated as '+'")}
  on_plus <- strand(x) == "+" | strand(x) == "*"
  new_start <- start(x) - ifelse(on_plus, upstream, downstream)
  new_end <- end(x) + ifelse(on_plus, downstream, upstream)
  ranges(x) <- IRanges(new_start, new_end)
  trim(x)
}

grcontract <- function(x, gene.region){
  if (any(strand(x) == "*")){
    warning("'*' ranges were treated as '+'")}
  on_plus <- strand(x) == "+" | strand(x) == "*"

  if(gene.region == "TSS"){
    new_start <- ifelse(on_plus, start(x), end(x)-1)
    new_end <- ifelse(on_plus,start(x)+1, end(x))
  }
  if(gene.region == "TTS"){

    new_start <- ifelse(on_plus, end(x)-1, start(x))
    new_end <- ifelse(on_plus, end(x), start(x)+1)

  }
  if(gene.region == "Whole"){
    new_start <- start(x)
    new_end <- end(x)
  }

  ranges(x) <- IRanges(new_start, new_end)
  trim(x)
}


# find coordinates based on gene names

extractCoordfromGene <- function(genes=NULL, geneCoordKey=geneKey.ranges, add.upstream = 0, add.downstream = 0, gene.region = "Whole"){

  coordinates <- subset(geneKey.ranges, mcols.GeneName %in% genes)
  coordinates <- grcontract(coordinates, gene.region = gene.region)

  coordinates <- grexpand(coordinates, upstream = add.upstream, downstream= add.downstream)

  return(coordinates)

}


# function for subsetting data matrix given a set of coordinates

subset.ATAC.FUN <- function(mydata = ATAC_TPMmatrix, coordinate.key = AllPeaks.granges, coordinates){

  if (length(coordinates) == 0 || any(is.na(coordinates$start)) || any(is.na(coordinates$end))) {
    stop("Invalid gene name or genomic location provided. Please provide a valid gene name or genomic location.")
  }

  mydata.subset <- mydata[subjectHits(findOverlaps(coordinates, coordinate.key)),]

  mydata.subset$PeakID <- row.names(mydata.subset)

  mydata.subset <-  pivot_longer(mydata.subset, cols = colnames(mydata))

  mydata.subset <- cbind(mydata.subset, str_split(mydata.subset$name,pattern="_", simplify = T))
  colnames(mydata.subset)[4:5] <- c("Biopsy","Hours")

  return(mydata.subset)

}


# function for plotting ATAC peaks based on set of valid genomic coordinates


plotATAC.FUN <- function(mydata = ATAC_TPMmatrix, coordinate.key = AllPeaks.granges, coordinates, times = c(0,3,6,9,12,18,24,36,48,96)){

  if (length(coordinates) == 0 || any(is.na(coordinates))) {
    stop("Invalid gene name or genomic location provided. Please provide a valid gene name or genomic location.")
  }

  plot.data <- subset.ATAC.FUN(mydata = mydata, coordinate.key = coordinate.key, coordinates = coordinates)

  p1 <- ggplot(plot.data)+
    stat_summary(aes(x=as.numeric(Hours),y=value,colour=PeakID), geom = "line", fun.y = mean, size = 2)+
    scale_x_continuous(breaks=times, minor_breaks = NULL)+
    theme_bw()+theme(axis.title = element_text(size=20), axis.text = element_text(size=20), legend.text= element_text(size=15), title= element_text(size=20) )+
    labs(y="RPM", x="Hours")+ggtitle("ATAC-seq")

  print(p1)

}


plot_error_message <- function(message) {
  plot(1, type = "n", xlim = c(0, 1), ylim = c(0, 1), xlab = "", ylab = "", main = "")
  text(x = 0.5, y = 0.5, label = message, cex = 1.5, col = "red", font = 2)
}
