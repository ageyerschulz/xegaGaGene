#
# (c) 2023 Andreas Geyer-Schulz
#     Simple Genetic Algorithm in R. V0.1
#     Layer: Gene-Level Functions 
#            For a binary gene representation.
#     Package: xegaGaGene
#

#' One point crossover of 2 genes.
#'
#' @description \code{xegaGaCross2Gene()} randomly determines a cut point.
#'              It combines the bits before the cut point of the first gene
#'              with the bits after the cut point from the second gene (kid 1).
#'              It combines the bits before the cut point of the second gene
#'              with the bits after the cut point from the first gene (kid 2).
#'              It returns 2  genes.
#'
#' @param gg1     A binary gene.
#' @param gg2     A binary gene.
#' @param lF      The local configuration of the genetic algorithm.
#'
#' @return A list of 2 binary genes.
#'
#' @family Crossover (Returns 2 Kids)
#'
#' @examples
#' gene1<-xegaGaInitGene(lFxegaGaGene)
#' gene2<-xegaGaInitGene(lFxegaGaGene)
#' xegaGaDecodeGene(gene1, lFxegaGaGene)
#' xegaGaDecodeGene(gene2, lFxegaGaGene)
#' newgenes<-xegaGaCross2Gene(gene1, gene2, lFxegaGaGene)
#' xegaGaDecodeGene(newgenes[[1]], lFxegaGaGene)
#' xegaGaDecodeGene(newgenes[[2]], lFxegaGaGene)
#' @importFrom utils head
#' @importFrom utils tail
#' @export
xegaGaCross2Gene<-function(gg1, gg2, lF)
{
    l<-length(gg1$gene1)
    cut<-sample(1:max(1, l-1), 1)
    ng1<-gg1; ng2<-gg2
    ng1$gene1<-c(utils::head(gg1$gene1,cut), utils::tail(gg2$gene1, l-cut))
    # Tests on equality?
    ng1$evaluated<-FALSE
    ng2$gene1<-c(utils::head(gg2$gene1,cut), utils::tail(gg1$gene1, l-cut))
    # Tests on equality?
    ng2$evaluated<-FALSE
    return(list(ng1, ng2))
}

#' Uniform crossover of 2 genes.
#'
#' @description \code{xegaGaUCross2Gene()} swaps alleles of both genes
#'              with a probability of 0.5. It generates a random 
#'              mask which is used to build the new genes.
#'              It returns 2 genes.
#'
#' @param gg1     A binary gene.
#' @param gg2     A binary gene.
#' @param lF      The local configuration of the genetic algorithm.
#'
#' @return A list of 2 binary genes.
#'
#' @references 
#'   Syswerda, Gilbert (1989):
#'   Uniform Crossover in Genetic Algorithms. 
#'   In: Schaffer, J. David (Ed.)
#'   Proceedings of the Third International Conference on Genetic Algorithms,
#'   Morgan Kaufmann Publishers, Los Altos, California, pp. 2-9.
#'   (ISBN:1-55860-066-3)
#'
#' @family Crossover (Returns 2 Kids)
#'
#' @examples
#' gene1<-xegaGaInitGene(lFxegaGaGene)
#' gene2<-xegaGaInitGene(lFxegaGaGene)
#' xegaGaDecodeGene(gene1, lFxegaGaGene)
#' xegaGaDecodeGene(gene2, lFxegaGaGene)
#' newgenes<-xegaGaUCross2Gene(gene1, gene2, lFxegaGaGene)
#' xegaGaDecodeGene(newgenes[[1]], lFxegaGaGene)
#' xegaGaDecodeGene(newgenes[[2]], lFxegaGaGene)
#' @importFrom stats runif 
#' @export
xegaGaUCross2Gene<-function(gg1, gg2, lF)
{
    ng1<-gg1; ng2<-gg2
    n1<-g1<-gg1$gene1; n2<-g2<-gg2$gene1
    mask<-0.5>runif(rep(1, length(g1)))  
    n1[mask]<-g2[mask]; n2[mask]<-g1[mask]
    # Tests on equality?
    ng1$evaluated<-FALSE; ng2$evaluated<-FALSE
    ng1$gene1<-n1; ng2$gene1<-n2
    return(list(ng1, ng2))
}

#' Parameterized uniform crossover of 2 genes.
#'
#' @description \code{xegaGaUP2CrossGene()} swaps alleles of both genes
#'              with a probability of \code{lF$UCrossSwap()}. 
#'              It generates a random 
#'              mask which is used to build the new gene.
#'              It returns 2  genes.
#' @references 
#'   Spears William and De Jong, Kenneth (1991):
#'   On the Virtues of Parametrized Uniform Crossover. 
#'   In: Belew, Richar K. and Booker, Lashon B. (Ed.)
#'   Proceedings of the Fourth International Conference on Genetic Algorithms,
#'   Morgan Kaufmann Publishers, Los Altos, California, pp. 230-236.
#'   (ISBN:1-55860-208-9)
#'
#' @param gg1     A binary gene.
#' @param gg2     A binary gene.
#' @param lF      The local configuration of the genetic algorithm.
#'
#' @return A list of 2 binary genes.
#'
#' @family Crossover (Returns 2 Kids)
#'
#' @examples
#' gene1<-xegaGaInitGene(lFxegaGaGene)
#' gene2<-xegaGaInitGene(lFxegaGaGene)
#' xegaGaDecodeGene(gene1, lFxegaGaGene)
#' xegaGaDecodeGene(gene2, lFxegaGaGene)
#' newgenes<-xegaGaUPCross2Gene(gene1, gene2, lFxegaGaGene)
#' xegaGaDecodeGene(newgenes[[1]], lFxegaGaGene)
#' xegaGaDecodeGene(newgenes[[2]], lFxegaGaGene)
#' @importFrom stats runif 
#' @export
xegaGaUPCross2Gene<-function(gg1, gg2, lF)
{
    ng1<-gg1; ng2<-gg2
    n1<-g1<-gg1$gene1; n2<-g2<-gg2$gene1
    mask<-lF$UCrossSwap()>runif(rep(1, length(n1)))  
    n1[mask]<-g2[mask]; n2[mask]<-g1[mask]
    # Tests on equality?
    ng1$evaluated<-FALSE; ng2$evaluated<-FALSE
    ng1$gene1<-n1; ng2$gene1<-n2
    return(list(ng1, ng2))
}

#' One point crossover of 2 genes.
#'
#' @description \code{xegaGaCrossGene()} randomly determines a cut point.
#'              It combines the bits before the cut point of the first gene
#'              with the bits after the cut point from the second gene (kid 1).
#'
#' @param gg1     A binary gene.
#' @param gg2     A binary gene.
#' @param lF      The local configuration of the genetic algorithm.
#'
#' @return A list of one binary gene.
#'
#' @family Crossover (Returns 1 Kid)
#'
#' @examples
#' gene1<-xegaGaInitGene(lFxegaGaGene)
#' gene2<-xegaGaInitGene(lFxegaGaGene)
#' xegaGaDecodeGene(gene1, lFxegaGaGene)
#' xegaGaDecodeGene(gene2, lFxegaGaGene)
#' gene3<-xegaGaCrossGene(gene1, gene2, lFxegaGaGene)
#' xegaGaDecodeGene(gene3[[1]], lFxegaGaGene)
#' @importFrom utils head
#' @importFrom utils tail
#' @export
xegaGaCrossGene<-function(gg1, gg2, lF)
{
    g1<-gg1$gene1; g2<-gg2$gene1
    cut<-sample(1:max(1,(length(g1)-1)), 1)
    ng<-gg1
    ng$gene1<-c(utils::head(g1,cut), utils::tail(g2, length(g2)-cut))
    # Tests on equality?
    ng$evaluated<-FALSE
    return(list(ng))
}

#' Uniform crossover of 2 genes.
#'
#' @description \code{xegaGaUCrossGene()} swaps alleles of both genes
#'              with a probability of 0.5. It generates a random 
#'              mask which is used to build the new gene.
#'
#' @references 
#'   Syswerda, Gilbert (1989):
#'   Uniform Crossover in Genetic Algorithms. 
#'   In: Schaffer, J. David (Ed.)
#'   Proceedings of the Third International Conference on Genetic Algorithms,
#'   Morgan Kaufmann Publishers, Los Altos, California, pp. 2-9.
#'   (ISBN:1-55860-066-3)
#'
#' @param gg1     A binary gene.
#' @param gg2     A binary gene.
#' @param lF      The local configuration of the genetic algorithm.
#'
#' @return A list of one binary gene.
#'
#' @family Crossover (Returns 1 Kid)
#'
#' @examples
#' gene1<-xegaGaInitGene(lFxegaGaGene)
#' gene2<-xegaGaInitGene(lFxegaGaGene)
#' xegaGaDecodeGene(gene1, lFxegaGaGene)
#' xegaGaDecodeGene(gene2, lFxegaGaGene)
#' gene3<-xegaGaUCrossGene(gene1, gene2, lFxegaGaGene)
#' xegaGaDecodeGene(gene3[[1]], lFxegaGaGene)
#' @importFrom stats runif 
#' @export
xegaGaUCrossGene<-function(gg1, gg2, lF)
{
    ng1<-gg1; ng2<-gg2
    n1<-gg1$gene1; g2<-gg2$gene1
    mask<-0.5>runif(rep(1, length(n1)))  
    n1[mask]<-g2[mask]
    # Tests on equality?
    ng1$evaluated<-FALSE
    ng1$gene1<-n1
    return(list(ng1))
}

#' Parameterized uniform crossover of 2 genes.
#'
#' @description \code{xegaGaUPCrossGene()} swaps alleles of both genes
#'              with a probability of \code{lF$UCrossSwap()}. 
#'              It generates a random 
#'              mask which is used to build the new gene.
#'
#' @references 
#'   Spears William and De Jong, Kenneth (1991):
#'   On the Virtues of Parametrized Uniform Crossover. 
#'   In: Belew, Richar K. and Booker, Lashon B. (Ed.)
#'   Proceedings of the Fourth International Conference on Genetic Algorithms,
#'   Morgan Kaufmann Publishers, Los Altos, California, pp. 230-236.
#'   (ISBN:1-55860-208-9)
#'
#' @param gg1     A binary gene.
#' @param gg2     A binary gene.
#' @param lF      The local configuration of the genetic algorithm.
#'
#' @return A list of one binary gene.
#'
#' @family Crossover (Returns 1 Kid)
#'
#' @examples
#' gene1<-xegaGaInitGene(lFxegaGaGene)
#' gene2<-xegaGaInitGene(lFxegaGaGene)
#' xegaGaDecodeGene(gene1, lFxegaGaGene)
#' xegaGaDecodeGene(gene2, lFxegaGaGene)
#' gene3<-xegaGaUPCrossGene(gene1, gene2, lFxegaGaGene)
#' xegaGaDecodeGene(gene3[[1]], lFxegaGaGene)
#' @importFrom stats runif 
#' @export
xegaGaUPCrossGene<-function(gg1, gg2, lF)
{
    ng1<-gg1; ng2<-gg2
    n1<-gg1$gene1; g2<-gg2$gene1
    mask<-lF$UCrossSwap()>runif(rep(1, length(n1)))  
    n1[mask]<-g2[mask]
    # Tests on equality?
    ng1$evaluated<-FALSE
    ng1$gene1<-n1
    return(list(ng1))
}

#' Configure the crossover function of a genetic algorithm.
#'
#' @description \code{xegaGaCrossoverFactory()} implements the selection
#'              of one of the crossover functions in this
#'              package by specifying a text string.
#'              The selection fails ungracefully (produces
#'              a runtime error) if the label does not match.
#'              The functions are specified locally.
#'
#'              Current support:
#'
#'              \enumerate{
#'              \item Crossover functions with two kids:
#'              \enumerate{
#'              \item "Cross2Gene" returns \code{xegaGaCross2Gene()}.
#'              \item "UCross2Gene" returns \code{xegaGaUCross2Gene()}.
#'              \item "UPCross2Gene" returns \code{xegaGaUPCross2Gene()}.
#'              }
#'              \item Crossover functions with one kid:
#'              \enumerate{
#'              \item "CrossGene" returns \code{xegaGaCrossGene()}.
#'              \item "UCrossGene" returns \code{xegaGaUCrossGene()}.
#'              \item "UPCrossGene" returns \code{xegaGaUPCrossGene()}.
#'              }
#'              }
#'
#' @details Crossover operations which return 2 kids preserve the genetic
#'          material in the population. However, because we work with fixed 
#'          size populations, genes with 2 offsprings fill two slots in the
#'          new population with their genetic material. 
#'
#' @param method     A string specifying the crossover function.
#'
#' @return A crossover function for genes.
#'
#' @family Configuration
#'
#' @examples
#' XGene<-xegaGaCrossoverFactory("Cross2Gene")
#' gene1<-xegaGaInitGene(lFxegaGaGene)
#' gene2<-xegaGaInitGene(lFxegaGaGene)
#' XGene(gene1, gene2, lFxegaGaGene)
#' @export
xegaGaCrossoverFactory<-function(method="Cross2Gene") {
if (method=="Cross2Gene") {f<- xegaGaCross2Gene}
if (method=="UCross2Gene") {f<- xegaGaUCross2Gene}
if (method=="UPCross2Gene") {f<- xegaGaUPCross2Gene}
if (method=="CrossGene") {f<- xegaGaCrossGene}
if (method=="UCrossGene") {f<- xegaGaUCrossGene}
if (method=="UPCrossGene") {f<- xegaGaUPCrossGene}
if (!exists("f", inherits=FALSE))
        {stop("sga Crossover label ", method, " does not exist")}
return(f)
}

