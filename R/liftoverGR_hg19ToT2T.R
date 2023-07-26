#' @title  Liftover genomic range from hg19 to t2t (Human genome)
#'
#' @description  Liftover genomic range from hg19 to t2t (Human genome) and will also add t2t assembly information to the output genomic ranges. (filter out non-standard chromosomes)
#' @param gr.f A genomic range.
#' @return Liftover of hg19 genomic loci to t2t genomic loci.
#' @export
#' @examples
#' data(hg_38_gr)
#' outGr_h19=liftoverGR_hg38ToHg19(gr.f=hg_38_gr)
#' outGr_t2t=liftoverGR_hg19ToT2T(gr.f=outGr_h19)
#' outGr_t2t
#'
#'
liftoverGR_hg19ToT2T=function(gr.f){
  seqlevelsStyle(gr.f) = "UCSC"  # necessary
  ch.url='https://hgdownload.gi.ucsc.edu/hubs/GCA/009/914/755/GCA_009914755.4/liftOver/hg19-chm13v2.over.chain.gz'
  td = tempdir() # # create a temporary directory
  tf = tempfile(tmpdir=td, fileext=".gz") # # create the placeholder file
  download.file(ch.url, tf) # # download into the placeholder file
  fname = gsub("[.]gz$", "", tf) 
  gunzip(filename=tf, destname = fname , overwrite = FALSE, remove = TRUE) # # unzip the file to the temporary directory
  ch = import.chain(fname)
  gr.f = unlist(liftOver(gr.f, ch))
  gr.f=addGrInformation(gr.f,'t2t')
}