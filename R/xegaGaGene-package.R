
#' Genetic operations for binary coded genetic algorithms.
#' 
#' For an introduction to this class of algorithms, see Goldberg, D. (1989).
#'
#' For binary-coded genes, the \code{xegaGaGene} package provides
#' \itemize{
#' \item Gene initiatilization.
#' \item Decoding of parameters as well as a function factory for configuration.
#' \item Mutation functions as well as a function factory for configuration.
#' \item Crossover functions as well as a function factory for configuration.
#'       We provide two families of crossover functions:
#'  \enumerate{
#' \item Crossover functions with two kids:
#'       Crossover preserves the genetic information in the gene pool.
#' \item Crossover functions with one kid:
#'       These functions allow the construction of gene evaluation pipelines.
#'       One advantage of this is a simple control structure 
#'       at the population level.
#' \item Gene replication functions as well as a function factory for 
#'       configuration. The replication functions implement control flows
#'       for sequences of gene operations. For \code{xegaReplicateGene}, 
#'       an acceptance step has been added. Simulated annealing algorithms
#'       can be configured e.g. by configuring uniform random selection combined
#'       with a Metropolis Acceptance Rule and a suitable cooling schedule.
#' }
#' }
#' 
#' @section Binary Gene Representation:
#'            
#' A binary gene is a named list:
#'   \itemize{
#'    \item $gene1      the gene must be a binary vector.
#'    \item $fit        the fitness value of the gene
#'                      (for EvalGeneDet and EvalGeneU) or
#'                      the mean fitness (for stochastic functions
#'                      evaluated with EvalGeneStoch).
#'    \item $evaluated  has the gene been evaluated?
#'    \item $evalFail   has the evaluation of the gene failed?
#'    \item $var        the cumulative variance of the fitness 
#'                      of all evaluations of a gene.
#'                      (For stochastic functions)
#'    \item $sigma      the standard deviation of the fitness of 
#'                      all evaluations of a gene.
#'                      (For stochastic functions)
#'    \item $obs        the number of evaluations of a gene.
#'                      (For stochastic functions)
#'   }
#'
#' @section Abstract Interface of Problem Environment:
#'
#' A problem environment \code{penv} must provide:
#'   \itemize{
#'     \item \code{$f(parameters, gene, lF)}: 
#'   Function with a real parameter vector as first argument 
#'   which returns a gene 
#'   with evaluated fitness.
#'   
#'   \item $genelength(): The number of bits of the binary-coded
#'                        real parameter vector. Used in \code{InitGene}.
#'     \item $bitlength(): A vector specifying the number of bits 
#'                        used for coding each real parameter.
#'                        If \code{penv$bitlength()[1]} is \code{20}, 
#'                        then \code{parameters[1]} is coded by 20 bits.
#'           Used in \code{GeneMap}.
#'     \item $lb(): The lower bound vector of each parameter.
#'           Used in \code{GeneMap}.
#'     \item $ub(): The upper bound vector of each parameter.
#'           Used in \code{GeneMap}.
#'   } 
#'
#' @section Abstract Interface of Mutation Functions:
#'
#' Each mutation function has the following function signature:
#'
#'     newGene<-Mutate(gene, lF) 
#'
#' All local parameters of the mutation function configured are 
#' expected in the local function list lF.
#' 
#' @section Local Constants of Mutation Functions:
#'
#' The local constants of a mutation function determine 
#' the behavior of the function. The default values in the 
#' table below are set in \code{lFxegaGaGene}.
#'
#' \tabular{rcl}{ 
#' \strong{Constant} \tab \strong{Default} \tab \strong{Used in} \cr 
#' lF$BitMutationRate1() \tab 0.01       \tab xegaGaMutateGene() \cr 
#'                     \tab              \tab xegaGaIVAdaptiveMutateGene() \cr 
#' lF$BitMutationRate2() \tab 0.20       \tab xegaGaIVAdaptiveMutateGene() \cr 
#' lF$CutoffFit()        \tab 0.5        \tab xegaGaIVADaptiveMutateGene() \cr
#' }
#'
#' @section Abstract Interface of Crossover Functions:
#'
#' The signatures of the abstract interface to the 2 families 
#' of crossover functions are:
#'
#'     ListOfTwoGenes<-Crossover2(gene1, gene2, lF) 
#'
#'     ListOfOneGene<-Crossover(gene1, gene2, lF) 
#'
#' All local parameters of the crossover function configured are 
#' expected in the local function list lF.
#'
#' @section Local Constants of Crossover Functions:
#'
#' The local constants of a crossover function determine the 
#' the behavior of the function. 
#'
#' \tabular{rcl}{ 
#' \strong{Constant} \tab \strong{Default} \tab \strong{Used in} \cr 
#' lF$UCrossSwap()     \tab 0.2              \tab UPCross2Gene() \cr 
#'                     \tab                \tab UPCrossGene() \cr 
#' }
#'
#' @section Abstract Interface of Gene Replication Functions:
#'
#' The signatures of the abstract interface to the 2 
#' gene replication functions are:
#'
#'     ListOfTwoGenes<-Replicate2Gene(gene1, gene2, lF) 
#'
#'     ListOfOneGene<-ReplicateGene(gene1, gene2, lF) 
#'
#' @section Configuration and Constants of Replication Functions:
#'
#' \strong{Configuration for ReplicateGene (1 Kid, Default).}
#'
#' \tabular{rcl}{ 
#' \strong{Function} \tab \strong{Default} \tab Configured By \cr 
#' lF$SelectGene()   \tab SelectSUS()      \tab SelectGeneFactory() \cr 
#' lF$SelectMate()   \tab SelectSUS()      \tab SelectGeneFactory() \cr 
#' lF$CrossGene()    \tab CrossGene()      \tab xegaGaCrossoverFactory() \cr
#' lF$MutateGene()   \tab MutateGene()     \tab xegaGaMutationFactory() \cr
#' lF$Accept()       \tab AcceptNewGene()  \tab AcceptFactory() \cr
#' }
#'
#' \strong{Configuration for Replicate2Gene (2 Kids).}
#'
#' \tabular{rcl}{ 
#' \strong{Function} \tab \strong{Default} \tab Configured By \cr 
#' lF$SelectGene()   \tab SelectSUS()        \tab SelectGeneFactory() \cr 
#' lF$SelectMate()   \tab SelectSUS()        \tab SelectGeneFactory() \cr 
#' lF$CrossGene()    \tab CrossGene()        \tab xegaGaCrossoverFactory() \cr
#' lF$MutateGene()   \tab MutateGene()       \tab xegaGaMutationFactory() \cr
#' }
#'
#' \strong{Global Constants.}
#'
#' Global constants specify the probability that a mutation or
#' crossover operator is applied to a gene.
#' In the xega-architecture, these rates can be configured to be 
#' adaptive.
#'
#' \tabular{rcl}{ 
#' \strong{Constant} \tab \strong{Default} \tab \strong{Used in} \cr 
#' lF$MutationRate() \tab 1.0 (static)     \tab xegaGaReplicateGene() \cr 
#'                   \tab                  \tab xegaGaReplicate2Gene() \cr 
#' lF$CrossRate()    \tab 0.2 (static)     \tab xegaGaReplicateGene() \cr 
#'                   \tab                  \tab xegaGaReplicate2Gene() \cr 
#' }
#'
#' \strong{Local Constants.}
#'
#' \tabular{rcl}{ 
#' \strong{Constant} \tab \strong{Default} \tab \strong{Used in} \cr 
#' lF$BitMutationRate1() \tab 0.01    \tab xegaGaMutateGene() \cr 
#'                       \tab         \tab xegaGaIVAdaptiveMutateGene() \cr 
#' lF$BitMutationRate2() \tab 0.20    \tab xegaGaIVAdaptiveMutateGene() \cr 
#' lF$CutoffFit()        \tab 0.5     \tab xegaGaIVADaptiveMutateGene() \cr
#' lF$UCrossSwap()       \tab 0.2     \tab xegaGaUPCross2Gene() \cr 
#'                       \tab         \tab xegaGaUPCrossGene() \cr 
#' }
#'
#' In the xega-architecture, these rates can be configured to be 
#' adaptive.
#'
#' @section The Architecture of the xegaX-Packages:
#' 
#' The xegaX-packages are a family of R-packages which implement 
#' eXtended Evolutionary and Genetic Algorithms (xega).  
#' The architecture has 3 layers, 
#' namely the user interface layer,
#' the population layer, and the gene layer: 
#' 
#' \itemize{
#' \item
#' The user interface layer (package \code{xega}) 
#' provides a function call interface and configuration support
#' for several algorithms: genetic algorithms (sga), 
#' permutation-based genetic algorithms (sgPerm), 
#' derivation-free algorithms as e.g. differential evolution (sgde), 
#' grammar-based genetic programming (sgp) and grammatical evolution
#' (sge). 
#'
#' \item
#' The population layer (package \code{xegaPopulation}) contains
#' population-related functionality as well as support for 
#' population statistics dependent adaptive mechanisms and parallelization.
#'
#' \item 
#' The gene layer is split into a representation-independent and 
#' a representation-dependent part:
#' \enumerate{
#' \item 
#'  The representation indendent part (package \code{xegaSelectGene})
#'  is responsible for variants of selection operators, evaluation 
#'  strategies for genes, as well as profiling and timing capabilities.        
#' \item 
#'  The representation dependent part consists of the following packages: 
#' \itemize{
#' \item \code{xegaGaGene} for binary coded genetic algorithms.
#' \item \code{xegaPermGene} for permutation-based genetic algorithms.
#' \item \code{xegaDfGene} for derivation-free algorithms as e.g. 
#'                         differential evolution.
#' \item \code{xegaGpGene} for grammar-based genetic algorithms.
#' \item \code{xegaGeGene} for grammatical evolution algorithms.
#' }
#' The packages \code{xegaDerivationTrees} and \code{xegaBNF} support
#' the last two packages:
#' \code{xegaBNF} essentially provides a grammar compiler, and 
#' \code{xegaDerivationTrees} is an abstract data type for derivation trees.
#' }} 
#'
#' @references
#' Goldberg, David E. (1989)
#' Genetic Algorithms in Search, Optimization and Machine Learning.
#' Addison-Wesley, Reading. 
#' (ISBN:0-201-15767-5)
#'
#' @family Package Description
#'
#' @name xegaGaGene
#' @aliases xegaGaGene
#' @docType package
#' @title Package xegaGaGene.
#' @author Andreas Geyer-Schulz
#' @section Copyright: (c) 2023 Andreas Geyer-Schulz
#' @section License: MIT
#' @section URL: <https://github.com/ageyerschulz/xegaGaGene>
#' @section Installation: From CRAN by \code{install.packages('xegaGaGene')}
"_PACKAGE"

