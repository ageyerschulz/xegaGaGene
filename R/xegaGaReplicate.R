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
#' @description \code{xegaGaReplicate2Gene} replicates a gene
#'              by applying a gene reproduction pipeline 
#'              which uses crossover and
#'              mutation. The control flow is as follows:
#'              \itemize{
#'              \item A gene is selected from the population.
#'              Check if the crossover operation should be applied.
#'              (The check is \code{TRUE} with a probability of \code{crossrate}).
#'              If the check is \code{TRUE}:
#'              \itemize{
#'                \item Select a mating gene from the population.
#'                \item Perform the crossover operation.
#'                \item Apply mutation with a probability of \code{mutrate}.
#'                \item Return a list with both genes.
#'                            }
#'                \item Apply mutation with a probability of \code{mutrate}.
#'                \item Return a list with a single gene.
#'              }
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
#' @description \code{xegaGaReplicateGene} replicates a gene
#'              by applying a gene reproduction pipeline 
#'              which uses crossover and
#'              mutation. 
#'              The control flow may have the following steps:
#'              \itemize{
#'              \item A gene is selected from the population.
#'              Check if the crossover operation should be applied.
#'              (The check is \code{TRUE} with a probability of \code{crossrate}).
#'              If the check is \code{TRUE}:
#'              \itemize{
#'                \item Select a mating gene from the population.
#'                \item Perform the crossover operation.
#'                \item Apply mutation with a probability of \code{mutrate}.
#'                \item Return a list one gene.
#'                            }
#'                \item Apply mutation with a probability of \code{mutrate}.
#'                \item Accept gene. For genetic algorithms: Identity.
#'                \item Return a list with a single gene.
#'              }
#'
#' @details \code{xegaGaReplicateGene} implements the control flow 
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
#' @description \code{ReplicationFactory} implements the selection
#'              of a replication method. 
#'
#'              Current support:
#'
#'              \enumerate{
#'              \item "Kid1" returns \code{ReplicateGene}.
#'              \item "Kid2" returns \code{Replicate2Gene}.
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

