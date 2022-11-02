#' @title  Generate all possible DNA sequences from a \code{\link[https://genomevolution.org/wiki/index.php/Ambiguous_nucleotide]{Ambiguous nucleotide sequence}}
#'
#' @description Generate all possible DNA sequences from a \code{\link[https://genomevolution.org/wiki/index.php/Ambiguous_nucleotide]{Ambiguous nucleotide sequence}}
#' @param pat Ambiguous nucleotide sequence
#' @return all possible DNA sequences for a given Ambiguous nucleotide sequence
#' @export
#' @examples
#' DNA_seqs=DNASeqsForPattern(pat='NYYN')
#' DNA_seqs
#'
#'
DNASeqsForPattern=function(pat){
  possibleCombination=as.vector(unlist(Disambiguate(DNAStringSet(pat))))
  return(possibleCombination)
}
