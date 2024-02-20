#
# (c) 2021 Andreas Geyer-Schulz
#          Simple Genetic Algorithm in R. V 0.1
#          Layer: Gene/Population-level functions.
#                 Independent of gene representation.
#                 The replication mechanism and its variants
#          Package: xegaGaGene
#

#' Replicates a gene.
#'
#' @description \code{xegaGaReplicate2Gene()} replicates a gene
#'              by 2 random experiments which determine if a mutation 
#'              operator (boolean variable \code{mut})  and/or 
#'              a crossover operator (boolean variable \code{cross} 
#'              should be applied. For each of the 4 cases, the 
#'              appropriate code is executed.
#'
#' @details \code{xegaGaReplicate2Gene()} implements the control flow 
#'          by case distinction which  depends
#'          on the random choices for mutation and crossover:
#' \enumerate{
#'   \item A gene \code{g} is selected and the boolean variables \code{mut}
#'         and \code{cross} are set to \code{runif(1)<rate}. 
#'         \code{rate} is given by 
#'         \code{lF$MutationRate()} or \code{lF$CrossRate()}. 
#'   \item The truth values of \code{cross} and \code{mut} determine 
#'         the code that is executed:
#'   \enumerate{      
#'   \item \code{(cross==TRUE) & (mut==TRUE)}: 
#'           Mate selection,  crossover, mutation. 
#'   \item \code{(cross==TRUE) & (mut==FALSE)}: 
#'           Mate selection, crossover. 
#'   \item \code{(cross==FALSE) & (mut==TRUE)}: 
#'           Mutation. 
#'   \item \code{(cross==FALSE) & (mut==FALSE)} is implicit: 
#'          Returns a gene list. 
#'   }
#'   }
#'
#' @param pop    A population of binary genes.
#' @param fit    Fitness vector.
#' @param lF     The local configuration of the genetic algorithm.
#'
#' @return A list of either 1 or 2 binary genes.
#'
#' @family Replication
#'
#' @examples
#' lFxegaGaGene$CrossGene<-xegaGaCross2Gene
#' lFxegaGaGene$MutationRate<-function(fit, lF) {0.001}
#' names(lFxegaGaGene)
#' pop10<-lapply(rep(0,10), function(x) xegaGaInitGene(lFxegaGaGene))
#' epop10<-lapply(pop10, lFxegaGaGene$EvalGene, lF=lFxegaGaGene)
#' fit10<-unlist(lapply(epop10, function(x) {x$fit}))
#' newgenes<-xegaGaReplicate2Gene(pop10, fit10, lFxegaGaGene)
#'
#' @importFrom stats runif
#' @export
xegaGaReplicate2Gene<- function(pop, fit, lF)
{
    g<-pop[[lF$SelectGene(fit, lF)]]
    mut<-runif(1)<lF$MutationRate(g$fit, lF)
    cross<-runif(1)<lF$CrossRate(g$fit, lF)
     
   if (cross && mut) 
    { return(lapply(lF$CrossGene(g, pop[[lF$SelectMate(fit, lF)]], lF), 
	       lF$MutateGene, lF=lF))}
  if ((cross) && (!mut)) { 
    return(lF$CrossGene(g, pop[[lF$SelectMate(fit, lF)]], lF)) }
  if ((!cross) && (mut)) { return(list(lF$MutateGene(g, lF))) }
  return(list(g))
}

#' Replicates a gene.
#'
#' @description \code{xegaGaReplicateGene()} replicates a gene
#'              by applying a gene reproduction pipeline 
#'              which uses crossover and
#'              mutation and finishes with an acceptance rule.
#'              The control flow starts
#'              by selecting a gene from the population
#'              followed by the case distinction:
#'              \itemize{
#'               \item
#'              Check if the mutation operation should be applied.
#'              (\code{mut} is \code{TRUE} with a probability of \code{lF$MutationRate()}).
#'              \item
#'              Check if the crossover operation should be applied.
#'              (\code{cross} is \code{TRUE} with a probability of \code{lF$CrossRate()}).
#'              }
#'              The state distinction determines which genetic operations are performed.
#'
#' @details \code{xegaGaReplicateGene()} implements the control flow 
#'          by a dynamic definition of the operator pipeline depending
#'          on the random choices for mutation and crossover:
#' \enumerate{
#'   \item A gene \code{g} is selected and the boolean variables \code{mut}
#'         and \code{cross} are set to \code{runif(1)<rate}. 
#'   \item The local function for the operator pipeline \code{OPpip(g, lF)}
#'         is defined by the truth values of \code{cross} and \code{mut}:
#'   \enumerate{      
#'   \item \code{(cross==FALSE) & (mut==FALSE)}: 
#'           Identity function. 
#'   \item \code{(cross==TRUE) & (mut==TRUE)}: 
#'           Mate selection,  crossover, mutation. 
#'   \item \code{(cross==TRUE) & (mut==FALSE)}: 
#'           Mate selection, crossover. 
#'   \item \code{(cross==FALSE) & (mut==TRUE)}: 
#'           Mutation. 
#'   }
#'   \item  Perform the operator pipeline and accept the result.
#'          The acceptance step allows the combination of a genetic algorithm
#'          with other heuristic algorithms like simulated 
#'          annealing by executing an acceptance rule. 
#'          For the genetic algorithm, the identity function is used.            
#'   }
#'
#' @param pop    Population of binary genes.
#' @param fit    Fitness vector.
#' @param lF     Local configuration of the genetic algorithm.
#'
#' @return A list of one gene.
#'
#' @family Replication
#'
#' @examples
#' lFxegaGaGene$CrossGene<-xegaGaCrossGene
#' lFxegaGaGene$MutationRate<-function(fit, lF) {0.001}
#' lFxegaGaGene$Accept<-function(OperatorPipeline, gene, lF) {gene}
#' pop10<-lapply(rep(0,10), function(x) xegaGaInitGene(lFxegaGaGene))
#' epop10<-lapply(pop10, lFxegaGaGene$EvalGene, lF=lFxegaGaGene)
#' fit10<-unlist(lapply(epop10, function(x) {x$fit}))
#' newgenes<-xegaGaReplicateGene(pop10, fit10, lFxegaGaGene)
#' @importFrom stats runif
#' @export
xegaGaReplicateGene<- function(pop, fit, lF)
{
    g<-pop[[lF$SelectGene(fit, lF)]]
    mut<-runif(1)<lF$MutationRate(g$fit, lF)
    cross<-runif(1)<lF$CrossRate(g$fit, lF)
    OPpip<-function(g, lF) {g}

    if (cross && mut)
    {OPpip<-function(g, lF)
    { g1<-lF$CrossGene(g, pop[[lF$SelectMate(fit, lF)]], lF)
    lF$MutateGene(g1[[1]], lF)
    }}

  if ((cross) && (!mut))
    {OPpip<-function(g, lF)
    { lF$CrossGene(g, pop[[lF$SelectMate(fit, lF)]], lF)[[1]]}}

  if ((!cross) && (mut)) 
    {OPpip<-function(g, lF) { lF$MutateGene(g, lF)}}
    
  return(list(lF$Accept(OPpip, g, lF)))
}

#' Configure the replication function of a genetic algorithm.
#'
#' @description \code{xegaGaReplicationFactory()} implements the selection
#'              of a replication method. 
#'
#'              Current support:
#'
#'              \enumerate{
#'              \item "Kid1" returns \code{xegaGaReplicateGene()}.
#'              \item "Kid2" returns \code{xegaGaReplicate2Gene()}.
#'              }
#'
#' @param method     A string specifying the replication function.
#'
#' @return A replication function for genes.
#'
#' @family Configuration
#'
#' @examples
#' lFxegaGaGene$CrossGene<-xegaGaCrossGene
#' lFxegaGaGene$MutationRate<-function(fit, lF) {0.001}
#' lFxegaGaGene$Accept<-function(OperatorPipeline, gene, lF) {gene}
#' Replicate<-xegaGaReplicationFactory("Kid1")
#' pop10<-lapply(rep(0,10), function(x) xegaGaInitGene(lFxegaGaGene))
#' epop10<-lapply(pop10, lFxegaGaGene$EvalGene, lF=lFxegaGaGene)
#' fit10<-unlist(lapply(epop10, function(x) {x$fit}))
#' newgenes1<-Replicate(pop10, fit10, lFxegaGaGene)
#' lFxegaGaGene$CrossGene<-xegaGaCross2Gene
#' Replicate<-xegaGaReplicationFactory("Kid2")
#' newgenes2<-Replicate(pop10, fit10, lFxegaGaGene)
#' @export
xegaGaReplicationFactory<-function(method="Kid1") {
if (method=="Kid1") {f<- xegaGaReplicateGene}
if (method=="Kid2") {f<- xegaGaReplicate2Gene}
if (!exists("f", inherits=FALSE))
        {stop("xegaGaGene Replication label ", method, " does not exist")}
return(f)
}

