# (c) 2025 Andreas Geyer-Schulz
#     xegaOperatorPipelinesInGene.R
#

#
# Pipeline 1: No mutation, no crossover
# 

#' Converts a gene into a gene with genetic operator pipeline. 
#'
#' @description The pipeline is \code{evaluate(gene)}.
#'
#' @details \code{newPipelineG} is a constructor 
#'          which integrates a genetic operator pipeline which 
#'          contains the evaluation of a gene into a gene.
#'
#' @param  g   A gene.
#'
#' @return A gene with embedded genetic operator pipeline  
#'         without mutation and crossover.
#'         The argument \code{lF} of the function \code{$Pipeline}
#'         configures the behavior of the pipeline.
#' 
#' @family Genetic Operator Pipelines in Gene
#' 
#' @examples
#' g<-xegaGaInitGene(lFxegaGaGene)
#' a<-newPipelineG(g)
#' print(a)
#' b<-a$Pipeline(a, lFxegaGaGene)
#' print(b)
#' @export
newPipelineG<-function(g)
{ g$Pipeline<-function(g, lF) 
  {  ng<-lF$EvalGene(g, lF) 
     ng$Pipeline<-NULL
     return(ng) }
  environment(g$Pipeline)<-globalenv()
  return(g)
}

#
# Pipeline 2: Mutation.
#

#' Converts a gene into a genetic operator pipeline embedded in a gene with mutation.
#'
#' @description The embedded pipeline is \code{evaluate(accept(mutate, gene))}.
#'
#' @param  g   A gene.
#'
#' @return A gene with embedded genetic operator pipeline 
#'         with mutation only.
#'         The argument \code{lF} of the function \code{$Pipeline}
#'         configures the behavior of the pipeline.
#' 
#' @family Genetic Operator Pipelines in Gene
#' 
#' @examples
#' lFxegaGaGene$MutationRate<-function(fit, lF) {0.5}
#' lFxegaGaGene$BitMutationRate1<-function(fit, lF) {0.5}
#' lFxegaGaGene$CrossRate<-function(fit, lF) {0.5}
#' lFxegaGaGene$Accept<-function(OperatorPipeline, gene, lF) {OperatorPipeline(gene, lF)}
#' g<-xegaGaInitGene(lFxegaGaGene)
#' print((lFxegaGaGene$EvalGene(g, lFxegaGaGene))$gene1)
#' a<-newMutPipelineG(g)
#' print(a)
#' b<-a$Pipeline(a, lFxegaGaGene)
#' @export
newMutPipelineG<-function(g)
{
g$Pipeline<-function(g, lF) 
{  ng<-lF$EvalGene(lF$Accept(lF$MutateGene, g, lF), lF) 
   ng$Pipeline<-NULL
   return(ng)
   }
  environment(g$Pipeline)<-globalenv()
return(g)
}

#
# Pipeline 3: Crossover
#

#' Converts two genes into a genetic operator pipeline in a gene with crossover (1 kid).
#'
#' @description The embedded pipeline is \code{evaluate(accept(crossover, gene, gene1))}.
#'
#' @param  g   A gene.
#' @param  g1   A gene.
#'
#' @return A gene with embedded genetic operator pipeline 
#'         with crossover only.
#'         The argument \code{lF} of the function \code{$Pipeline}
#'         configures the behavior of the pipeline.
#' 
#' @family Genetic Operator Pipelines in Gene
#' 
#' @examples
#' lFxegaGaGene$CrossGene<-xegaGaCrossGene
#' lFxegaGaGene$MutationRate<-function(fit, lF) {0.5}
#' lFxegaGaGene$BitMutationRate1<-function(fit, lF) {0.9}
#' lFxegaGaGene$CrossRate<-function(fit, lF) {0.5}
#' lFxegaGaGene$Accept<-function(OperatorPipeline, gene, lF) {OperatorPipeline(gene, lF)}
#' g<-xegaGaInitGene(lFxegaGaGene)
#' g1<-xegaGaInitGene(lFxegaGaGene)
#' a<-newCrossPipelineG(g, g1)
#' print(a)
#' a$Pipeline(a, lFxegaGaGene)
#' @export
newCrossPipelineG<-function(g, g1)
{
g$gmate<-g1
g$Pipeline<-function(g, lF) 
{  OPpip<-function(g, lF)
    {lF$CrossGene(g, g$gmate, lF)[[1]]}
   ng<-lF$EvalGene(lF$Accept(OPpip, g, lF), lF) 
   ng$gmate<-NULL; ng$Pipeline<-NULL
   return(ng) }
  environment(g$Pipeline)<-globalenv()
return(g)
}

#' Converts two genes into a genetic operator pipeline embedded in a gene with crossover (2 kids).
#'
#' @description The embedded pipeline is \code{evaluate(accept(crossover, gene, gene1))}.
#'              The execution of this pipeline produces two genes.
#'
#' @param  g   A gene.
#' @param  g1   A gene.
#'
#' @return A gene with embedded genetic operator pipeline 
#'         with crossover only.
#'         The argument \code{lF} of function \code{$Pipeline()} 
#'         configures the behavior of the pipeline.
#' 
#' @family Genetic Operator Pipelines in Gene
#' 
#' @examples
#' lFxegaGaGene$CrossGene<-xegaGaCross2Gene
#' lFxegaGaGene$MutationRate<-function(fit, lF) {0.5}
#' lFxegaGaGene$CrossRate<-function(fit, lF) {0.5}
#' lFxegaGaGene$Accept<-function(OpPipeline, gene, lF) {OpPipeline(gene, lF)}
#' g<-xegaGaInitGene(lFxegaGaGene)
#' g1<-xegaGaInitGene(lFxegaGaGene)
#' a<-newCross2PipelineG(g, g1)
#' print(a)
#' a$Pipeline(a, lFxegaGaGene)
#' @export
newCross2PipelineG<-function(g, g1)
{
g$gmate<-g1
g$Pipeline<-function(g, lF) 
{   g1g2<-lF$CrossGene(g, g$gmate, lF)
    OPpip1<-function(g1g2, lF) {g1g2[[1]]}
    OPpip2<-function(g1g2, lF) {g1g2[[2]]}
   ng1<-lF$EvalGene(lF$Accept(OPpip1, g1g2, lF), lF)
   ng1$gmate<-NULL; ng1$Pipeline<-NULL
   ng2<-lF$EvalGene(lF$Accept(OPpip2, g1g2, lF), lF)
   ng2$gmate<-NULL; ng2$Pipeline<-NULL
   return(list(ng1, ng2)) }
  environment(g$Pipeline)<-globalenv()
return(g)
}

#
# Pipeline 4: Crossover and Mutation
#

#' Converts a gene into a gene with embedded genetic operator pipeline  with crossover and mutation (a function closure).
#'
#' @description The embedded pipeline is \code{evaluate(accept((crossover o mutation), gene, gene1))}.
#'              The symbol \code{o} is short for \code{mutation(crossover(gene, gene1))}
#'              in the accept function.
#'
#' @param  g   A gene.
#' @param  g1   A gene.
#'
#' @return A gene with embedded genetic operator pipeline 
#'         with mutation and crossover.
#'         The argument \code{lF} of the function \code{$Pipeline()} 
#'         configures the behavior of the pipeline.
#' 
#' @family Genetic Operator Pipelines in Gene
#' 
#' @examples
#' lFxegaGaGene$CrossGene<-xegaGaCrossGene
#' lFxegaGaGene$MutationRate<-function(fit, lF) {0.5}
#' lFxegaGaGene$BitMutationRate1<-function(fit, lF) {0.2}
#' lFxegaGaGene$CrossRate<-function(fit, lF) {0.5}
#' lFxegaGaGene$Accept<-function(OperatorPipeline, gene, lF) {OperatorPipeline(gene, lF)}
#' g<-xegaGaInitGene(lFxegaGaGene)
#' g1<-xegaGaInitGene(lFxegaGaGene)
#' a<-newCrossMutPipelineG(g, g1)
#' print(a)
#' a$Pipeline(a, lFxegaGaGene)
#' @export
newCrossMutPipelineG<-function(g, g1)
{  g$gmate<-g1
   g$Pipeline<-function(g, lF) 
{  OPpip<-function(g, lF)
    { g2<-lF$CrossGene(g, g$gmate, lF)[[1]]
      lF$MutateGene(g2, lF)} 
   ng<-lF$EvalGene(lF$Accept(OPpip, g, lF), lF) 
   ng$gmate<-NULL
   ng$Pipeline<-NULL
   return(ng)
   }
 environment(g$Pipeline)<-globalenv()
return(g) }

#
# Pipeline 5: Crossover and Mutation
#

#' Converts two genes into a pipeline embedded in a gene with crossover and mutation for both kids.
#'
#' @description The embedded pipeline is \code{evaluate(accept((crossover o mutation), gene, gene1))}.
#'              Mutation is applied to both kids.
#'
#' @param  g   A gene.
#' @param  g1   A gene.
#'
#' @return A gene with embedded genetic operator pipeline 
#'         with crossover with two kids and mutation on both kids.
#'         The argument \code{lF} of the function \code{$Pipeline()} 
#'         configures the behavior of the pipeline.
#' 
#' @family Genetic Operator Pipelines in Gene
#' 
#' @examples
#' lFxegaGaGene$CrossGene<-xegaGaCross2Gene
#' lFxegaGaGene$MutationRate<-function(fit, lF) {0.5}
#' lFxegaGaGene$BitMutationRate1<-function(fit, lF) {0.2}
#' lFxegaGaGene$CrossRate<-function(fit, lF) {0.5}
#' lFxegaGaGene$Accept<-function(OpPipeline, gene, lF) {OpPipeline(gene, lF)}
#' g<-xegaGaInitGene(lFxegaGaGene)
#' g1<-xegaGaInitGene(lFxegaGaGene)
#' a<-newCrossMut2PipelineG(g, g1)
#' print(a)
#' a$Pipeline(a, lFxegaGaGene)
#' @export
newCrossMut2PipelineG<-function(g, g1)
{
g$gmate<-g1
g$Pipeline<-function(g, lF) 
{   g1g2<-lF$CrossGene(g, g$gmate, lF)
    OPpip1<-function(g1g2, lF) {lF$MutateGene(g1g2[[1]], lF)} 
    OPpip2<-function(g1g2, lF) {lF$MutateGene(g1g2[[2]], lF)} 
   ng1<-lF$EvalGene(lF$Accept(OPpip1, g1g2, lF), lF) 
   ng1$gmate<-NULL; ng1$Pipeline<-NULL
   ng2<-lF$EvalGene(lF$Accept(OPpip2, g1g2, lF), lF) 
   ng2$gmate<-NULL; ng2$Pipeline<-NULL
   return(list(ng1,ng2))  }
   environment(g$Pipeline)<-globalenv()
return(g)
}

#' Converts two genes into a pipeline with crossover (2 kids) and mutation for first kid.
#'
#' @description The pipeline is \code{evaluate(accept((crossover o mutation), gene, gene1))}.
#'              Mutation is applied to the first kid.
#'
#' @param  g   A gene.
#' @param  g1   A gene.
#'
#' @return Closure of genetic operator pipeline 
#'         with crossover with 2 kids and mutation on the first kid only.
#'         The argument of the closure \code{lF} 
#'         configures the behavior of the pipeline.
#' 
#' @family Genetic Operator Pipelines in Gene
#' 
#' @examples
#' lFxegaGaGene$CrossGene<-xegaGaCross2Gene
#' lFxegaGaGene$MutationRate<-function(fit, lF) {0.5}
#' lFxegaGaGene$CrossRate<-function(fit, lF) {0.5}
#' lFxegaGaGene$Accept<-function(OpPipeline, gene, lF) {OpPipeline(gene, lF)}
#' g<-xegaGaInitGene(lFxegaGaGene)
#' g1<-xegaGaInitGene(lFxegaGaGene)
#' a<-newCross2Mut1PipelineG(g, g1)
#' print(a)
#' a$Pipeline(a, lFxegaGaGene)
#' @export
newCross2Mut1PipelineG<-function(g, g1)
{
g$gmate<-g1
g$Pipeline<-function(g, lF) 
{   g1g2<-lF$CrossGene(g, g$gmate, lF)
    OPpip1<-function(g1g2, lF) {lF$MutateGene(g1g2[[1]], lF)} 
    OPpip2<-function(g1g2, lF) {g1g2[[2]]} 
   ng1<-lF$EvalGene(lF$Accept(OPpip1, g1g2, lF), lF)
   ng1$gmate<-NULL; ng1$Pipeline<-NULL
   ng2<-lF$EvalGene(lF$Accept(OPpip2, g1g2, lF), lF) 
   ng2$gmate<-NULL; ng2$Pipeline<-NULL
   return(list(ng1, ng2))}
environment(g$Pipeline)<-globalenv()
return(g)
}

#' Converts two genes into a pipeline embedded in a gene with crossover (2 kids) and mutation for second kid.
#'
#' @description The embdded pipeline is \code{evaluate(accept((crossover o mutation), gene, gene1))}.
#'              Mutation is applied to the second kid.
#'
#' @param  g   A gene.
#' @param  g1   A gene.
#'
#' @return A gene with embedded genetic operator pipeline 
#'         with crossover with 2 kids and mutation on the second kid only.
#'         The argument \code{lF} of the function \code{Pipeline}
#'         configures the behavior of the pipeline.
#' 
#' @family Genetic Operator Pipelines in Gene
#' 
#' @examples
#' lFxegaGaGene$CrossGene<-xegaGaCross2Gene
#' lFxegaGaGene$MutationRate<-function(fit, lF) {0.5}
#' lFxegaGaGene$BitMutationRate1<-function(fit, lF) {1.0}
#' lFxegaGaGene$CrossRate<-function(fit, lF) {0.5}
#' lFxegaGaGene$Accept<-function(OpPipeline, gene, lF) {OpPipeline(gene, lF)}
#' g<-xegaGaInitGene(lFxegaGaGene)
#' g1<-xegaGaInitGene(lFxegaGaGene)
#' a<-newCross2Mut1PipelineG(g, g1)
#' print(a)
#' a$Pipeline(a, lFxegaGaGene)
#' @export
newCross2Mut2PipelineG<-function(g, g1)
{
g$gmate<-g1
g$Pipeline<-function(g, lF) 
{   g1g2<-lF$CrossGene(g, g$gmate, lF)
    OPpip1<-function(g1g2, lF) {g1g2[[1]]} 
    OPpip2<-function(g1g2, lF) {lF$MutateGene(g1g2[[2]], lF)} 
   ng1<-lF$EvalGene(lF$Accept(OPpip1, g, lF), lF)
   ng1$gmate<-NULL; ng1$Pipeline<-NULL
   ng2<-lF$EvalGene(lF$Accept(OPpip2, g, lF), lF) 
   ng2$gmate<-NULL; ng2$Pipeline<-NULL
   return(list(ng1, ng2))}
environment(g$Pipeline)<-globalenv()
return(g)
}
