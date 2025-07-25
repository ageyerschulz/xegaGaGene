# (c) 2025 Andreas Geyer-Schulz
#     xegaOperatorPipelines.R
#

#
# Pipeline 1: No mutation, no crossover
# 

#' Converts a gene into a genetic operator pipeline (a function closure).
#'
#' @description The pipeline is \code{evaluate(gene)}.
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

#' Converts a gene into a genetic operator pipeline with mutation (a function closure).
#'
#' @description The pipeline is \code{evaluate(accept(mutate, gene))}.
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

#' Converts two genes into a genetic operator pipeline with crossover (1 kid).
#'
#' @description The pipeline is \code{evaluate(accept(crossover, gene, gene1))}.
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

#' Converts two genes into a genetic operator pipeline with crossover (2 kids).
#'
#' @description The pipeline is \code{evaluate(accept(crossover, gene, gene1))}.
#'              The execution of this pipeline produces two genes.
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
#' lFxegaGaGene$CrossGene<-xegaGaCross2Gene
#' lFxegaGaGene$MutationRate<-function(fit, lF) {0.5}
#' lFxegaGaGene$CrossRate<-function(fit, lF) {0.5}
#' lFxegaGaGene$Accept<-function(OpPipeline, gene, lF) {OpPipeline(gene, lF)}
#' g<-xegaGaInitGene(lFxegaGaGene)
#' g1<-xegaGaInitGene(lFxegaGaGene)
#' a<-newCross2Pipeline(g, g1, lFxegaGaGene)
#' print(a)
#' a(lFxegaGaGene)
#' @export
newCross2Pipeline<-function(g, g1, lF)
{
gene<-g
gene1<-g1
Pipeline<-function(lF) 
{   g1g2<-lF$CrossGene(gene, gene1, lF)
    OPpip1<-function(g, lF)
    {g1g2[[1]]}
    OPpip2<-function(g, lF)
    {g1g2[[2]]}
   list(lF$EvalGene(lF$Accept(OPpip1, gene, lF), lF), 
   lF$EvalGene(lF$Accept(OPpip2, gene, lF), lF)) 
     }
# force
a<-gene
a<-gene1
return(Pipeline)
}

#
# Pipeline 4: Crossover and Mutation
#

#' Converts a gene into a genetic operator pipeline with crossover and mutation (a function closure).
#'
#' @description The pipeline is \code{evaluate(accept((crossover o mutation), gene, gene1))}.
#'              The symbol \code{o} is short for \code{mutation(crossover(gene, gene1))}
#'              in the accept function.
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

#
# Pipeline 5: Crossover and Mutation
#

#' Converts two genes into a pipeline with crossover and mutation for both kids.
#'
#' @description The pipeline is \code{evaluate(accept((crossover o mutation), gene, gene1))}.
#'              Mutation is applied to both kids.
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
#' lFxegaGaGene$CrossGene<-xegaGaCross2Gene
#' lFxegaGaGene$MutationRate<-function(fit, lF) {0.5}
#' lFxegaGaGene$CrossRate<-function(fit, lF) {0.5}
#' lFxegaGaGene$Accept<-function(OpPipeline, gene, lF) {OpPipeline(gene, lF)}
#' g<-xegaGaInitGene(lFxegaGaGene)
#' g1<-xegaGaInitGene(lFxegaGaGene)
#' a<-newCrossMut2Pipeline(g, g1, lFxegaGaGene)
#' print(a)
#' a(lFxegaGaGene)
#' @export
newCrossMut2Pipeline<-function(g, g1, lF)
{
gene<-g
gene1<-g1
Pipeline<-function(lF) 
{   g1g2<-lF$CrossGene(gene, gene1, lF)
    OPpip1<-function(g, lF)
      {lF$MutateGene(g1g2[[1]], lF)} 
    OPpip2<-function(g, lF)
      {lF$MutateGene(g1g2[[2]], lF)} 
   list(lF$EvalGene(lF$Accept(OPpip1, gene, lF), lF), 
   lF$EvalGene(lF$Accept(OPpip2, gene, lF), lF)) 
     }
# force
a<-gene
a<-gene1
return(Pipeline)
}

#' Converts two genes into a pipeline with crossover (2 kids) and mutation for first kid.
#'
#' @description The pipeline is \code{evaluate(accept((crossover o mutation), gene, gene1))}.
#'              Mutation is applied to the first kid.
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
#' lFxegaGaGene$CrossGene<-xegaGaCross2Gene
#' lFxegaGaGene$MutationRate<-function(fit, lF) {0.5}
#' lFxegaGaGene$CrossRate<-function(fit, lF) {0.5}
#' lFxegaGaGene$Accept<-function(OpPipeline, gene, lF) {OpPipeline(gene, lF)}
#' g<-xegaGaInitGene(lFxegaGaGene)
#' g1<-xegaGaInitGene(lFxegaGaGene)
#' a<-newCross2Mut1Pipeline(g, g1, lFxegaGaGene)
#' print(a)
#' a(lFxegaGaGene)
#' @export
newCross2Mut1Pipeline<-function(g, g1, lF)
{
gene<-g
gene1<-g1
Pipeline<-function(lF) 
{   g1g2<-lF$CrossGene(gene, gene1, lF)
    OPpip1<-function(g, lF)
      {lF$MutateGene(g1g2[[1]], lF)} 
    OPpip2<-function(g, lF)
      {g1g2[[2]]} 
   list(lF$EvalGene(lF$Accept(OPpip1, gene, lF), lF), 
   lF$EvalGene(lF$Accept(OPpip2, gene, lF), lF)) 
     }
# force
a<-gene
a<-gene1
return(Pipeline)
}

#' Converts two genes into a pipeline with crossover (2 kids) and mutation for second kid.
#'
#' @description The pipeline is \code{evaluate(accept((crossover o mutation), gene, gene1))}.
#'              Mutation is applied to the second kid.
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
#' lFxegaGaGene$CrossGene<-xegaGaCross2Gene
#' lFxegaGaGene$MutationRate<-function(fit, lF) {0.5}
#' lFxegaGaGene$CrossRate<-function(fit, lF) {0.5}
#' lFxegaGaGene$Accept<-function(OpPipeline, gene, lF) {OpPipeline(gene, lF)}
#' g<-xegaGaInitGene(lFxegaGaGene)
#' g1<-xegaGaInitGene(lFxegaGaGene)
#' a<-newCross2Mut1Pipeline(g, g1, lFxegaGaGene)
#' print(a)
#' a(lFxegaGaGene)
#' @export
newCross2Mut2Pipeline<-function(g, g1, lF)
{
gene<-g
gene1<-g1
Pipeline<-function(lF) 
{   g1g2<-lF$CrossGene(gene, gene1, lF)
    OPpip1<-function(g, lF)
      {g1g2[[1]]} 
    OPpip2<-function(g, lF)
      {lF$MutateGene(g1g2[[2]], lF)} 
   list(lF$EvalGene(lF$Accept(OPpip1, gene, lF), lF), 
   lF$EvalGene(lF$Accept(OPpip2, gene, lF), lF)) 
     }
# force
a<-gene
a<-gene1
return(Pipeline)
}
