#' @title  Create GenomicRanges object of given the bin size  for human genome.
#'
#' @description Create GenomicRanges object of given the bin size  for human genome.
#' @param assmblyName hg19 or hg38
#' @param binSize  size of the genomic block
#' @return GenomicRanges object of given the bin size
#' @export
#' @examples
#' hg_19_Bins.gr=getGbins(assmblyName='hg19', binSize=1000 )
#' hg_19_Bins.gr
#'
#'
getGbins=function(assmblyName, binSize=1000){
  if(assmblyName=='hg19'){
    library(BSgenome.Hsapiens.UCSC.hg19)
    BS_genome=BSgenome.Hsapiens.UCSC.hg19
  }else if(assmblyName=='hg38'){
    library(BSgenome.Hsapiens.UCSC.hg38)
    BS_genome=BSgenome.Hsapiens.UCSC.hg38
  }
  fgbins=genomeBlocks(BS_genome, chrs = seqnames(BS_genome), width = binSize) # library("Repitools")
  fgbins=addGrInformation(fgbins,'hg19')
  fgBinSeqs=getSeq(x=BS_genome, fgbins)
  fgbins$CpG_counts=2*oligonucleotideFrequency(DNAStringSet(fgBinSeqs), width=2, step=1)[,'CG']
  return(fgbins)
}
