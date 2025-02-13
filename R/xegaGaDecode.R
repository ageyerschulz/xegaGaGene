
# (c) 2023 Andreas Geyer-Schulz
#     Simple Genetic Algorithm in R. V0.1
#     Layer: Gene-Level Functions 
#            For a binary gene representation.
#     Package: xegaGaGene
#

#
# TODO: Bin2IEEE754
# See Goldberg, 1991, CSUR, 	
# What Every Computer Scientist Should Know About Floating-point Arithmetic 
# IEEE 754 Standard.
# 1 bit:  Sign
# m bit:  Mantissa
# r bit:  Exponent
#
# b: Base 2 or 10
# Bias: 2^r-1-1  (for unsigned exponents.
#

#' Map the bit strings of a binary gene to an identical bit vector.
#'
#' @description \code{xegaGaGeneMapIdentity()} maps the bit strings 
#'              of a binary vector
#'              to an identical binary vector.
#'              Faster for all problems with single-bit coding.
#'              Examples: Knapsack, Number Partitioning into 2 partitions.
#'
#' @param gene    A binary gene (the genotype).
#' @param penv    A problem environment.
#'
#' @return A binary gene (the phenotype).
#'
#' @family Decoder
#'
#' @examples
#' gene<-xegaGaInitGene(lFxegaGaGene)
#' xegaGaGeneMapIdentity(gene$gene1, lFxegaGaGene$penv)
#'
#' @export
xegaGaGeneMapIdentity<-function(gene, penv)
{ 
return(gene) 
}

#' Map the bit strings of a binary gene to parameters in an interval.
#'
#' @description \code{xegaGaGeneMap()} maps the bit strings of a binary string 
#'              to parameters in an interval.
#'              Bit vectors are mapped into equispaced numbers in the interval.
#'              Examples: Optimization of problems with real-valued 
#'                        parameter vectors.
#'
#' @param gene    A binary gene (the genotype).
#' @param penv    A problem environment.
#'
#' @return The decoded gene (the phenotype).
#'
#' @family Decoder
#'
#' @examples
#' gene<-xegaGaInitGene(lFxegaGaGene)
#' xegaGaGeneMap(gene$gene1, lFxegaGaGene$penv)
#'
#' @export
xegaGaGeneMap<-function(gene, penv)
{ nparm<-length(penv$bitlength())
  parm<-rep(0, nparm)
  s<-1
  for (i in 1:nparm)
  { bv<-gene[s:(s-1+penv$bitlength()[i])]
	s<-s+penv$bitlength()[i]
	parm[i]<-penv$lb()[i]+(penv$ub()[i]-penv$lb()[i])*
		(sum(bv*2^((length(bv)-1):0))/(2^length(bv)-1)) }
return(parm) 
}

#' Map Gray code to binary.
#' 
#' @details Start with the highest order bit, and 
#'          \code{r[k-i]<- xor(n[k], n[k-1])}.
#'
#' @param x   Gray code (boolean vector).
#' 
#' @return   Binary code (boolean vector).
#'
#' @references 
#' Gray, Frank (1953):
#' Pulse Code Communication. US Patent 2 632 058.
#'
#'@examples
#'  Gray2Bin(c(1, 0, 0, 0))
#'  Gray2Bin(c(1, 1, 1, 1))
#'@export
Gray2Bin<-function(x)
{
        r<-rep(0, length(x))
        r[1]<-x[1]
        for (i in 2:length(x))
        {r[i]<-xor(r[i-1], x[i])}
        return(r)
}


#' Map the bit strings of a gray-coded gene to parameters in an interval.
#'
#' @description \code{xegaGaGeneMapGray()} maps the bit strings of 
#'              a binary string 
#'              interpreted as Gray codes to parameters in an interval.
#'              Bit vectors are mapped into equispaced numbers in the interval.
#'              Examples: Optimization of problems with real-valued 
#'                        parameter vectors.
#'
#'
#' @param gene    A binary gene (the genotype).
#' @param penv    A problem environment.
#'
#' @return The decoded gene (the phenotype).
#'
#' @family Decoder
#'
#' @examples
#' gene<-xegaGaInitGene(lFxegaGaGene)
#' xegaGaGeneMapGray(gene$gene1, lFxegaGaGene$penv)
#'
#' @export
xegaGaGeneMapGray<-function(gene, penv)
{ nparm<-length(penv$bitlength())
  parm<-rep(0, nparm)
  s<-1
  for (i in 1:nparm)
  { bv<-Gray2Bin(gene[s:(s-1+penv$bitlength()[i])])
	s<-s+penv$bitlength()[i]
	parm[i]<-penv$lb()[i]+(penv$ub()[i]-penv$lb()[i])*
		(sum(bv*2^((length(bv)-1):0))/(2^length(bv)-1)) }
return(parm) 
}

#' Returns elements of
#' vector \code{x} without elements in \code{y}.
#'
#' @param x  A vector.
#' @param y  A vector.
#'
#' @return A vector.
#'
#' @family Utility
#'
#' @examples
#' a<-sample(1:15,15, replace=FALSE)
#' b<-c(1, 3, 5)
#' without(a, b)
#' @export
without<-function(x, y)
{
        x[!unlist(lapply(x, function(z) is.element(z, y)))]
}

#' Map the bit strings of a binary gene to a permutation.
#'
#' @description \code{xegaGaGeneMapPerm()} maps 
#'              the bit strings of a binary string 
#'              to a permutation of integers.
#'              Example: Traveling Salesman Problem (TSP).
#'
#' @param gene    A binary gene (the genotype).
#' @param penv    A problem environment.
#'
#' @return A permutation (the decoded gene (the phenotype))
#'
#' @family Decoder
#'
#' @examples
#' gene<-xegaGaInitGene(lFxegaGaGene)
#' xegaGaGeneMapPerm(gene$gene1, lFxegaGaGene$penv)
#'
#' @export
xegaGaGeneMapPerm<-function(gene, penv)
{ nparm<-length(penv$bitlength())
  p<-1:nparm
  parm<-rep(0, nparm)
  s<-1
  for (i in 1:nparm)
  { bv<-gene[s:(s-1+penv$bitlength()[i])]
	s<-s+penv$bitlength()[i]
	parm[i]<-p[(sum(bv*2^((length(bv)-1):0))%%(nparm+1-i))+1]
        p<-without(p, parm[i])        
  }
return(parm) 
}

#' Configure the gene map function of a genetic algorithm.
#'
#' @description \code{xegaGaGeneMapFactory()} implements the selection
#'              of one of the GeneMap functions in this
#'              package by specifying a text string.
#'              The selection fails ungracefully (produces
#'              a runtime error) if the label does not match.
#'              The functions are specified locally.
#'
#'              Current support:
#'
#'              \enumerate{
#'              \item "Bin2Dec" returns \code{xegaGaGeneMap()}. (Default).
#'              \item "Gray2Dec" returns \code{xegaGaGeneMapGray()}.
#'              \item "Identity" returns \code{xegaGaGeneMapIdentity()}.
#'              \item "Permutation" returns \code{xegaGaGeneMapPerm()}.
#'              }
#'
#' @param method A string specifying the GeneMap function.
#'
#' @return A gene map function for genes.
#'
#' @family Configuration
#'
#' @examples
#' XGene<-xegaGaGeneMapFactory("Identity")
#' gene<-xegaGaInitGene(lFxegaGaGene)
#' XGene(gene$gene1, lFxegaGaGene$penv)
#' @export
xegaGaGeneMapFactory<-function(method="Bin2Dec") {
if (method=="Bin2Dec") {f<-xegaGaGeneMap}
if (method=="Gray2Dec") {f<-xegaGaGeneMapGray}
if (method=="Identity") {f<-xegaGaGeneMapIdentity}
if (method=="Permutation") {f<-xegaGaGeneMapPerm}
if (!exists("f", inherits=FALSE))
        {stop("sga GeneMap label ", method, " does not exist")}
return(f)
}

#' Decode a gene.
#'
#' @description \code{xegaGaDecodeGene()} decodes a binary gene.
#'
#' @param gene   A binary gene (the genotype).
#' @param lF     The local configuration of the genetic algorithm.
#'
#' @return The decoded gene (the phenotype).
#'
#' @family Decoder
#'
#' @examples
#' gene<-xegaGaInitGene(lFxegaGaGene)
#' xegaGaDecodeGene(gene, lFxegaGaGene)
#'
#' @export
xegaGaDecodeGene<-function(gene, lF)
{
  lF$GeneMap(gene$gene1, lF$penv)
}

