# (c) 2025 Andreas Geyer-Schulz
#     xegaOperatorPipelines.R
#

#
# Pipeline 1: No mutation, no crossover
# 

#' Converts a gene into a genetic operator pipeline (a function closure).
#'
#' @param  g   A gene.
#' @param  lF  The local function configuration.
#'
#' @return Closure of genetic operator pipeline  
#'         without mutation and crossover.
#'         The argument of the closure \code{lF} 
#'         configures the behavior of the pipeline.
#' 
#' @family Genetic Operator Pipelines
#' 
#' @examples
#' lFxegaGaGene$CrossGene<-xegaGaCrossGene
#' lFxegaGaGene$MutationRate<-function(fit, lF) {0.5}
#' lFxegaGaGene$CrossRate<-function(fit, lF) {0.5}
#' lFxegaGaGene$Accept<-function(OperatorPipeline, gene, lF) {gene}
#' g<-xegaGaInitGene(lFxegaGaGene)
#' a<-newPipeline(g, lFxegaGaGene)
#' print(a)
#' a(lFxegaGaGene)
#' @export
newPipeline<-function(g, lF)
{ gene<-g
  Pipeline<-function(lF) 
  {  lF$EvalGene(gene, lF) }
  # force
  a<-gene
  return(Pipeline)
}

#
# Pipeline 2: Mutation.
#

#' Converts a gene into a genetic operator pipeline (a function closure).
#'
#' @param  g   A gene.
#' @param  lF  The local function configuration.
#'
#' @return Closure of genetic operator pipeline 
#'         with mutation only.
#'         The argument of the closure \code{lF} 
#'         configures the behavior of the pipeline.
#' 
#' @family Genetic Operator Pipelines
#' 
#' @examples
#' lFxegaGaGene$CrossGene<-xegaGaCrossGene
#' lFxegaGaGene$MutationRate<-function(fit, lF) {0.5}
#' lFxegaGaGene$CrossRate<-function(fit, lF) {0.5}
#' lFxegaGaGene$Accept<-function(OperatorPipeline, gene, lF) {gene}
#' g<-xegaGaInitGene(lFxegaGaGene)
#' a<-newMutPipeline(g, lFxegaGaGene)
#' print(a)
#' a(lFxegaGaGene)
#' @export
newMutPipeline<-function(g, lF)
{
gene<-g
Pipeline<-function(lF) 
{  lF$EvalGene(lF$Accept(lF$MutateGene, gene, lF), lF) }
# force
a<-gene
return(Pipeline)
}

#
# Pipeline 3: Crossover
#

#' Converts a gene into a genetic operator pipeline (a function closure).
#'
#' @param  g   A gene.
#' @param  g1   A gene.
#' @param  lF  The local function configuration.
#'
#' @return Closure of genetic operator pipeline 
#'         with crossover only.
#'         The argument of the closure \code{lF} 
#'         configures the behavior of the pipeline.
#' 
#' @family Genetic Operator Pipelines
#' 
#' @examples
#' lFxegaGaGene$CrossGene<-xegaGaCrossGene
#' lFxegaGaGene$MutationRate<-function(fit, lF) {0.5}
#' lFxegaGaGene$CrossRate<-function(fit, lF) {0.5}
#' lFxegaGaGene$Accept<-function(OperatorPipeline, gene, lF) {gene}
#' g<-xegaGaInitGene(lFxegaGaGene)
#' g1<-xegaGaInitGene(lFxegaGaGene)
#' a<-newCrossPipeline(g, g1, lFxegaGaGene)
#' print(a)
#' a(lFxegaGaGene)
#' @export
newCrossPipeline<-function(g, g1, lF)
{
gene<-g
gene1<-g1
Pipeline<-function(lF) 
{  OPpip<-function(g, lF)
    { lF$CrossGene(gene, gene1, lF)[[1]]}
   lF$EvalGene(lF$Accept(OPpip, gene, lF), lF) }
# force
a<-gene
a<-gene1
return(Pipeline)
}

#
# Pipeline 4: Crossover and Mutation
#


#' Converts a gene into a genetic operator pipeline (a function closure).
#'
#' @param  g   A gene.
#' @param  g1   A gene.
#' @param  lF  The local function configuration.
#'
#' @return Closure of genetic operator pipeline 
#'         with mutation and crossover.
#'         The argument of the closure \code{lF} 
#'         configures the behavior of the pipeline.
#' 
#' @family Genetic Operator Pipelines
#' 
#' @examples
#' lFxegaGaGene$CrossGene<-xegaGaCrossGene
#' lFxegaGaGene$MutationRate<-function(fit, lF) {0.5}
#' lFxegaGaGene$CrossRate<-function(fit, lF) {0.5}
#' lFxegaGaGene$Accept<-function(OperatorPipeline, gene, lF) {gene}
#' g<-xegaGaInitGene(lFxegaGaGene)
#' g1<-xegaGaInitGene(lFxegaGaGene)
#' a<-newCrossMutPipeline(g, g1, lFxegaGaGene)
#' print(a)
#' a(lFxegaGaGene)
#' @export
newCrossMutPipeline<-function(g, g1, lF)
{
gene<-g
gene1<-g1
Pipeline<-function(lF) 
{  OPpip<-function(g, lF)
    { g1<-lF$CrossGene(gene, gene1, lF)[[1]]
      lF$MutateGene(g1, lF)} 
   lF$EvalGene(lF$Accept(OPpip, gene, lF), lF) }
# force
a<-gene
a<-gene1
return(Pipeline)
}

