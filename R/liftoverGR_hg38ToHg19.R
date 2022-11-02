#' @title  Liftover genomic range from hg38 to hg19 (Human genome)
#'
#' @description  Liftover genomic range from hg38 to hg19 (Human genome) and will also add hg19 assembly information to the output genomic ranges. (filter out non-standard chromosomes)
#' @param gr.f A genomic range.
#' @return Liftover of hg38 genomic loci to hg19 genomic loci.
#' @export
#' @examples
#' data(hg_38_gr)
#' outGr_h19=liftoverGR_hg38ToHg19(gr.f=hg_38_gr)
#' outGr_h19
#'
#'
liftoverGR_hg38ToHg19=function(gr.f){
  seqlevelsStyle(gr.f) = "UCSC"  # necessary
  path = system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain") #library(rtracklayer); library(liftOver)
  ch = import.chain(path)
  gr.f = unlist(liftOver(gr.f, ch))
  gr.f=addGrInformation(gr.f,'hg19')
}
