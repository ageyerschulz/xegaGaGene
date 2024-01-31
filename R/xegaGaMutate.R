#
# (c) 2023 Andreas Geyer-Schulz
#     Simple Genetic Algorithm in R. V0.1
#     Layer: Gene-Level Functions 
#            For a binary gene representation.
#     Package: xegaGaGene
#

#' Mutate a gene.
#'
#' @description \code{xegaGaMutateGene} mutates a binary gene.
#'               The per-bit mutation rate is given by MutationRate().
#'
#' @param gene   A binary gene.
#' @param lF     The local configuration of the genetic algorithm
#'
#' @return A binary gene.
#'
#' @family Mutation
#'
#' @examples
#' parm<-function(x) {function() {return(x)}}
#' lFxegaGaGene$BitMutationRate1<-parm(1.0)
#' gene1<-xegaGaInitGene(lFxegaGaGene)
#' xegaGaDecodeGene(gene1, lFxegaGaGene)
#' lFxegaGaGene$BitMutationRate1()
#' gene<-xegaGaMutateGene(gene1, lFxegaGaGene)
#' xegaGaDecodeGene(gene, lFxegaGaGene)
#' @importFrom stats runif 
#' @export
xegaGaMutateGene<-function(gene, lF)
{
ng<-gene
gene1<-gene$gene1
ng$gene1<-1*xor(gene1,stats::runif(length(gene1), 0, 1)<lF$BitMutationRate1())
if ((ng$evaluated==TRUE) && all(gene1==ng$gene1)) 
	{ng$evaluated<-TRUE} else {ng$evaluated<-FALSE}
return(ng) 
}

#' Individually variable adaptive mutation of a gene.
#'
#' @description \code{xegaGaIVAdaptiveMutateGene} mutates a binary gene.
#'              Two mutation rates (\code{lF$MutationRate()} 
#'              and \code{lF$MutationRate2()} which is higher than the first)
#'              are used depending on the relative fitness of the gene.
#'              \code{lF$CutoffFit} and \code{lF$CBestFitness} are used
#'              to determine the relative fitness of the gene.
#'              The rationale is that mutating genes having a low fitness
#'              with a higher probability rate improves the performance
#'              of a genetic algorithm, because the gene gets a higher 
#'              chance to improve.
#'
#' @details This principle is a candidate for a more abstract implementation,
#'          because it applies to all variants of evolutionary algorithms.
#'
#'          The goal is to separate the threshold code and the 
#'          representation-dependent part and 
#'          to combine them in the factory properly.
#'
#' @references
#'   Stanhope, Stephen A. and Daida, Jason M. (1996)
#'   An Individually Variable Mutation-rate Strategy for Genetic Algorithms.
#'   In: Koza, John (Ed.)
#'   Late Breaking Papers at the Genetic Programming 1996 Conference.
#'   Stanford University Bookstore, Stanford, pp. 177-185.
#'   (ISBN:0-18-201-031-7)
#'
#' @param gene     A binary gene.
#' @param lF       The local configuration of the genetic algorithm.
#'
#' @return A binary gene
#' 
#' @family Mutation
#'
#' @examples
#' parm<-function(x) {function() {return(x)}}
#'   lFxegaGaGene$BitMutationRate1<-parm(1.0)
#'   lFxegaGaGene$BitMutationRate2<-parm(0.5)
#' gene1<-xegaGaInitGene(lFxegaGaGene)
#' xegaGaDecodeGene(gene1, lFxegaGaGene)
#' gene<-xegaGaIVAdaptiveMutateGene(gene1, lFxegaGaGene)
#' xegaGaDecodeGene(gene, lFxegaGaGene)
#' @importFrom stats runif 
#' @export
xegaGaIVAdaptiveMutateGene<-function(gene, lF)
{
	ng<-gene
        gene1<-gene$gene1
	f<-gene$fit
	### A threshold based rate adaptation function.
	if (f>(lF$CutoffFit()*lF$CBestFitness())) 
		{MutRate<-lF$BitMutationRate1()} 
	else
		{MutRate<-lF$BitMutationRate2()} 
	###
	ng$gene1<-1*xor(gene1, stats::runif(length(gene1), 0, 1)<MutRate)
	if ((ng$evaluated==TRUE) && all(gene1==ng$gene1)) 
		{ng$evaluated<-TRUE} else {ng$evaluated<-FALSE}
	return(ng) 
}

#' Configure the mutation function of a genetic algorithm.
#'
#' @description \code{xegaGaMutationFactory} implements the selection
#'              of one of the mutation functions in this
#'              package by specifying a text string.
#'              The selection fails ungracefully (produces
#'              a runtime error) if the label does not match.
#'              The functions are specified locally.
#'
#'              Current support:
#'
#'              \enumerate{
#'              \item "MutateGene" returns \code{xegaGaMutateGene}.
#'              \item "IVMGene" returns \code{xegaGaIVAdaptiveMutateGene}.
#'              }
#'
#' @param method    A string specifying the mutation function.
#'
#' @return A mutation function for genes.
#'
#' @family Configuration
#'
#' @examples
#' parm<-function(x) {function() {return(x)}}
#' lFxegaGaGene$BitMutationRate1<-parm(1.0)
#' Mutate<-xegaGaMutationFactory("MutateGene")
#' gene1<-xegaGaInitGene(lFxegaGaGene)
#' gene1
#' Mutate(gene1, lFxegaGaGene)
#' @export
xegaGaMutationFactory<-function(method="MutateGene") {
if (method=="MutateGene") {f<- xegaGaMutateGene}
if (method=="IVM") {f<- xegaGaIVAdaptiveMutateGene}
if (!exists("f", inherits=FALSE))
        {stop("sga Mutation label ", method, " does not exist")}
return(f)
}

