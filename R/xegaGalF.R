#
# (c) 2023 Andreas Geyer-Schulz
#     Simple Genetic Algorithm in R. V0.1
#     Layer: Gene-Level Functions 
#            For a binary gene representation.
#     Package: xegaGaGene
#

#' The local function list lFxegaGaGene.
#'
#' @description
#' We enhance the configurability of our code by introducing 
#'  a function factory. The  function factory contains
#'  all the functions that are needed for defining
#'  local functions in genetic operators. The local function  
#'  list keeps the signatures of functions (e.g. mutation functions)
#'  uniform and small. At the same time, variants of functions
#'  can use different local functions. 
#'
#' @details
#'    We use the local function list for 
#'    \enumerate{
#'    \item
#'       replacing all constants by constant functions.
#'       
#'       Rationale: We need one formal argument (the local function list lF)
#'       and we can dispatch multiple functions. E.g.  \code{lF$verbose()}
#'   \item    
#'       dynamically binding a local function with a definition from a
#'       proper function factory. E.g., the selection methods 
#'       \code{lF$SelectGene()} and \code{lF$SelectMate()}.
#'       
#'  \item gene representations which require special functions to handle them:
#'        \code{lF$InitGene()}, \code{lF$DecodeGene()}, \code{lF$EvalGene()},
#'        \code{lF$ReplicateGene()}, ...
#'       
#'  } 
#'
#' @family Configuration
#'
#' @importFrom xegaSelectGene Parabola2DFactory
#' @importFrom xegaSelectGene SelectGeneFactory
#' @importFrom xegaSelectGene parm
#' @importFrom xegaSelectGene EvalGeneFactory
#' @export 
lFxegaGaGene<-list(
penv=xegaSelectGene::Parabola2DFactory(),
replay=xegaSelectGene::parm(0),
verbose=xegaSelectGene::parm(4),
CutoffFit=xegaSelectGene::parm(0.5),
CBestFitness=xegaSelectGene::parm(100),
CWorstFitness=xegaSelectGene::parm(-100),
MutationRate1=xegaSelectGene::parm(0.01),
MutationRate2=xegaSelectGene::parm(0.20),
BitMutationRate1=xegaSelectGene::parm(0.01),
BitMutationRate2=xegaSelectGene::parm(0.20),
MutationRate=function(fit, lF) {1.0},
MutateGene=xegaGaMutationFactory(),
CrossRate=function(fit, lF) {0.2},
UCrossSwap=xegaSelectGene::parm(0.2),
CrossGene=xegaGaCrossoverFactory(),
Max=xegaSelectGene::parm(1),
Offset=xegaSelectGene::parm(1),
Eps=xegaSelectGene::parm(0.01),
Elitist=xegaSelectGene::parm(TRUE),
TournamentSize=xegaSelectGene::parm(2),
GeneMap=xegaGaGeneMapFactory(method="Bin2Dec"),
SelectGene=xegaSelectGene::SelectGeneFactory(method="Proportional"),
SelectMate=xegaSelectGene::SelectGeneFactory(method="Uniform"),
Accept=function(op, gene, lF) {op(gene, lF)},
ReportEvalErrors=function() {TRUE},
Pipeline=xegaSelectGene::parm(FALSE),
InitGene=xegaGaInitGene,
DecodeGene=xegaGaDecodeGene,
EvalGene=xegaSelectGene::EvalGeneFactory(method="EvalGeneU"),
SelectionContinuation=xegaSelectGene::parm(TRUE),
Verbose=xegaSelectGene::parm(4),
lapply=base::lapply
)

