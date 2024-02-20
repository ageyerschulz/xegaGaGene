#
# (c) 2023 Andreas Geyer-Schulz
#     Simple Genetic Algorithm in R. V0.1
#     Layer: Gene-Level Functions 
#            For a binary gene representation.
#     Package: xegaGaGene
#

#' Generate a random binary gene.
#'
#' @description \code{xegaGaInitGene()} generates a random binary gene 
#'              with a given length.
#'
#' @param lF   The local configuration of the genetic algorithm.
#'
#' @return A binary gene (a named list):
#'         \itemize{
#'         \item \code{$evaluated}: FALSE. See package \code{xegaSelectGene}.
#'         \item \code{$evalFail}:  FALSE. Set by the error handler(s)
#'                                  of the evaluation functions
#'                                  in package \code{xegaSelectGene} 
#'                                  in the case of failure.
#'         \item \code{$fit}:       Fitness.
#'         \item \code{$gene1}:     Binary gene.
#'         }
#'
#' @family Gene Generation.
#'
#' @examples
#' xegaGaInitGene(lFxegaGaGene)
#'
#' @export
xegaGaInitGene<-function(lF)
{gene1<-sample(0:1, lF$penv$genelength(), replace=TRUE)
return(list(evaluated=FALSE, evalFail=FALSE, fit=0, gene1=gene1))
}


