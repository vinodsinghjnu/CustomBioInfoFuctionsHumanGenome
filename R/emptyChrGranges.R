#' @title  Create empty charomosome GenomicRange object for a given human genome assembly
#'
#' @description Create empty charomosome GenomicRange object for a given human genome assembly (for standard chromosomes )
#' @param assmblyName hg19 or hg38
#' @return empty charomosome GenomicRange
#' @export
#' @examples
#' hg_19_Chr.gr=emptyChrGranges('hg19')
#' hg_19_Chr.gr
#'
#'
emptyChrGranges=function(assmblyName){
  if(assmblyName=='hg19'){
    library(BSgenome.Hsapiens.UCSC.hg19)
    BS_genome=BSgenome.Hsapiens.UCSC.hg19
  }else if(assmblyName=='hg38'){
    library(BSgenome.Hsapiens.UCSC.hg38)
    BS_genome=BSgenome.Hsapiens.UCSC.hg38
  }
  emptyChrGr=GRanges(names(BS_genome), IRanges(start=1, end=seqlengths(BS_genome)), strand = '*' )

  emptyChrGr=addGrInformation(emptyChrGr,assmblyName)

  chrSeqs=getSeq(BS_genome, emptyChrGr)
  names(chrSeqs)=seqnames(emptyChrGr)
  emptyChrGr$Seqs=chrSeqs
  return(emptyChrGr)
}

