

library(testthat)
library(xegaSelectGene)
library(xegaGaGene)


test_that("lFxegaGaGene OK",
          {
	   expect_identical(lFxegaGaGene$penv$name(), "Parabola2D")
           expect_equal(lFxegaGaGene$replay(), 0)
           expect_equal(lFxegaGaGene$verbose(), 4)
           expect_equal(lFxegaGaGene$CutoffFit(), 0.5)
           expect_equal(lFxegaGaGene$CBestFitness(), 100)
           expect_equal(lFxegaGaGene$MutationRate1(), 0.01)
           expect_equal(lFxegaGaGene$MutationRate2(), 0.20)
           expect_equal(lFxegaGaGene$CrossRate(), 0.5)
           expect_equal(lFxegaGaGene$UCrossSwap(), 0.2)
           expect_equal(lFxegaGaGene$Max(), 1)
           expect_equal(lFxegaGaGene$Offset(), 1)
           expect_equal(lFxegaGaGene$Eps(), 0.01)
           expect_identical(lFxegaGaGene$Elitist(), TRUE)
           expect_equal(lFxegaGaGene$TournamentSize(), 2)
          }
)

#
# Cross2Gene
#

test_that("Cross2Gene OK",
{
   # Random number seed fixed.
   # It might happen that both genes are equal (with very small probability).
   set.seed(5)
   g1<-xegaGaInitGene(lFxegaGaGene)
   g1<-EvalGene(g1, lFxegaGaGene)
   g2<-xegaGaInitGene(lFxegaGaGene)
   g2<-EvalGene(g2, lFxegaGaGene)
   ng<-xegaGaCross2Gene(g1, g2, lFxegaGaGene)
   ng1<-ng[[1]]
   ng2<-ng[[2]]
   expect_identical(all(ng1$gene1==g1$gene1), all(ng2$gene1==g2$gene1))
   expect_identical(all(ng1$gene1==g1$gene1), FALSE)
   expect_identical(all(ng2$gene1==g2$gene1), FALSE)
   expect_identical(ng1$evaluated, FALSE)
   expect_identical(ng2$evaluated, FALSE)
   }
)

test_that("Cross2Gene 0/1 genes OK",
{  set.seed(21)
   g1<-xegaGaInitGene(lFxegaGaGene)
   g1<-EvalGene(g1, lFxegaGaGene)
   g1$gene1<-rep(0, length(g1$gene1))
   g2<-xegaGaInitGene(lFxegaGaGene)
   g2<-EvalGene(g2, lFxegaGaGene)
   g2$gene1<-rep(1, length(g2$gene1))
   ng<-xegaGaCross2Gene(g1, g2, lFxegaGaGene)
   ng1<-ng[[1]]
   ng2<-ng[[2]]
   expect_identical(all(ng1$gene1==g1$gene1), FALSE)
   expect_identical(all(ng2$gene1==g2$gene1), FALSE)
   expect_identical(ng1$evaluated, FALSE)
   expect_identical(ng2$evaluated, FALSE)
   }
)

#
# UCross2Gene
#

test_that("UCross2Gene OK",
{  set.seed(13)
   g1<-xegaGaInitGene(lFxegaGaGene)
   g1<-EvalGene(g1, lFxegaGaGene)
   g2<-xegaGaInitGene(lFxegaGaGene)
   g2<-EvalGene(g2, lFxegaGaGene)
   ng<-xegaGaUCross2Gene(g1, g2, lFxegaGaGene)
   ng1<-ng[[1]]
   ng2<-ng[[2]]
   expect_identical(all(ng1$gene1==g1$gene1), all(ng2$gene1==g2$gene1))
   expect_identical(ng1$evaluated, FALSE)
   expect_identical(ng2$evaluated, FALSE)
   }
)

test_that("UCross2Gene 0/1 genes OK",
{  set.seed(299)
   g1<-xegaGaInitGene(lFxegaGaGene)
   g1<-EvalGene(g1, lFxegaGaGene)
   g1$gene1<-rep(0, length(g1$gene1))
   g2<-xegaGaInitGene(lFxegaGaGene)
   g2<-EvalGene(g2, lFxegaGaGene)
   g2$gene1<-rep(1, length(g2$gene1))
   ng<-xegaGaUCross2Gene(g1, g2, lFxegaGaGene)
   ng1<-ng[[1]]
   ng2<-ng[[2]]
   expect_identical(all(ng1$gene1==g1$gene1), FALSE)
   expect_identical(all(ng2$gene1==g2$gene1), FALSE)
   expect_identical(ng1$evaluated, FALSE)
   expect_identical(ng2$evaluated, FALSE)
   }
)

#
# UPCross2Gene
#

test_that("UPCross2Gene OK",
{  set.seed(139)
   g1<-xegaGaInitGene(lFxegaGaGene)
   g1<-EvalGene(g1, lFxegaGaGene)
   g2<-xegaGaInitGene(lFxegaGaGene)
   g2<-EvalGene(g2, lFxegaGaGene)
   ng<-xegaGaUPCross2Gene(g1, g2, lFxegaGaGene)
   ng1<-ng[[1]]
   ng2<-ng[[2]]
   expect_identical(all(ng1$gene1==g1$gene1), all(ng2$gene1==g2$gene1))
   expect_identical(ng1$evaluated, FALSE)
   expect_identical(ng2$evaluated, FALSE)
   }
)

test_that("UPCross2Gene 0/1 genes OK",
{  set.seed(1023)
   g1<-xegaGaInitGene(lFxegaGaGene)
   g1<-EvalGene(g1, lFxegaGaGene)
   g1$gene1<-rep(0, length(g1$gene1))
   g2<-xegaGaInitGene(lFxegaGaGene)
   g2<-EvalGene(g2, lFxegaGaGene)
   g2$gene1<-rep(1, length(g2$gene1))
   ng<-xegaGaUPCross2Gene(g1, g2, lFxegaGaGene)
   ng1<-ng[[1]]
   ng2<-ng[[2]]
   expect_identical(all(ng1$gene1==g1$gene1), FALSE)
   expect_identical(all(ng2$gene1==g2$gene1), FALSE)
   expect_identical(ng1$evaluated, FALSE)
   expect_identical(ng2$evaluated, FALSE)
   }
)

#
# CrossGene
#

test_that("CrossGene OK",
{
   g1<-xegaGaInitGene(lFxegaGaGene)
   g1<-EvalGene(g1, lFxegaGaGene)
   g2<-xegaGaInitGene(lFxegaGaGene)
   g2<-EvalGene(g2, lFxegaGaGene)
   ng<-xegaGaCrossGene(g1, g2, lFxegaGaGene)
   expect_identical(ng[[1]]$evaluated, FALSE)
   }
)

test_that("CrossGene 0/1 genes OK",
{
   g1<-xegaGaInitGene(lFxegaGaGene)
   g1<-EvalGene(g1, lFxegaGaGene)
   g1$gene1<-rep(0, length(g1$gene1))
   g2<-xegaGaInitGene(lFxegaGaGene)
   g2<-EvalGene(g2, lFxegaGaGene)
   g2$gene1<-rep(1, length(g2$gene1))
   ng<-xegaGaCrossGene(g1, g2, lFxegaGaGene)
   expect_identical(all(ng[[1]]$gene1==g1$gene1), FALSE)
   expect_identical(all(ng[[1]]$gene1==g2$gene1), FALSE)
   expect_identical(ng[[1]]$evaluated, FALSE)
   }
)

#
# UCrossGene
#

test_that("UCrossGene OK",
{
   g1<-xegaGaInitGene(lFxegaGaGene)
   g1<-EvalGene(g1, lFxegaGaGene)
   g2<-xegaGaInitGene(lFxegaGaGene)
   g2<-EvalGene(g2, lFxegaGaGene)
   ng<-xegaGaUCrossGene(g1, g2, lFxegaGaGene)
   expect_identical(ng[[1]]$evaluated, FALSE)
   }
)

test_that("UCrossGene 0/1 genes OK",
{  set.seed(331)
   g1<-xegaGaInitGene(lFxegaGaGene)
   g1<-EvalGene(g1, lFxegaGaGene)
   g1$gene1<-rep(0, length(g1$gene1))
   g2<-xegaGaInitGene(lFxegaGaGene)
   g2<-EvalGene(g2, lFxegaGaGene)
   g2$gene1<-rep(1, length(g2$gene1))
   ng<-xegaGaUCrossGene(g1, g2, lFxegaGaGene)
   expect_identical(all(ng[[1]]$gene1==g1$gene1), FALSE)
   expect_identical(all(ng[[1]]$gene1==g2$gene1), FALSE)
   expect_identical(ng[[1]]$evaluated, FALSE)
   }
)

#
# UPCrossGene
#

test_that("UPCrossGene OK",
{
   g1<-xegaGaInitGene(lFxegaGaGene)
   g1<-EvalGene(g1, lFxegaGaGene)
   g2<-xegaGaInitGene(lFxegaGaGene)
   g2<-EvalGene(g2, lFxegaGaGene)
   ng<-xegaGaUPCrossGene(g1, g2, lFxegaGaGene)
   expect_identical(ng[[1]]$evaluated, FALSE)
   }
)

test_that("UPCrossGene 0/1 genes OK",
{  set.seed(997)
   g1<-xegaGaInitGene(lFxegaGaGene)
   g1<-EvalGene(g1, lFxegaGaGene)
   g1$gene1<-rep(0, length(g1$gene1))
   g2<-xegaGaInitGene(lFxegaGaGene)
   g2<-EvalGene(g2, lFxegaGaGene)
   g2$gene1<-rep(1, length(g2$gene1))
   ng<-xegaGaUPCrossGene(g1, g2, lFxegaGaGene)
   expect_identical(all(ng[[1]]$gene1==g1$gene1), FALSE)
   expect_identical(all(ng[[1]]$gene1==g2$gene1), FALSE)
   expect_identical(ng[[1]]$evaluated, FALSE)
   }
)

#
# xegaGaCrossoverFactory
#

test_that("xegaGaCrossoverFactory HUGO OK",
          {
           expect_error(xegaGaCrossoverFactory("HUGO"))
          }
)

test_that("xegaGaCrossoverFactory() OK",
{ set.seed(771)
   g1<-xegaGaInitGene(lFxegaGaGene)
   g1<-EvalGene(g1, lFxegaGaGene)
   g2<-xegaGaInitGene(lFxegaGaGene)
   g2<-EvalGene(g2, lFxegaGaGene)
   cCrossover<-xegaGaCrossoverFactory()
   ng<-cCrossover(g1, g2, lFxegaGaGene)
   ng1<-ng[[1]]
   ng2<-ng[[2]]
   expect_identical(all(ng1$gene1==g1$gene1), all(ng2$gene1==g2$gene1))
   expect_identical(ng1$evaluated, FALSE)
   expect_identical(ng2$evaluated, FALSE)
   }
)

test_that("xegaGaCrossoverFactory Cross2Gene OK",
{  set.seed(1719)
   g1<-xegaGaInitGene(lFxegaGaGene)
   g1<-EvalGene(g1, lFxegaGaGene)
   g2<-xegaGaInitGene(lFxegaGaGene)
   g2<-EvalGene(g2, lFxegaGaGene)
   cCrossover<-xegaGaCrossoverFactory(method="Cross2Gene")
   ng<-cCrossover(g1, g2, lFxegaGaGene)
   ng1<-ng[[1]]
   ng2<-ng[[2]]
   expect_identical(all(ng1$gene1==g1$gene1), all(ng2$gene1==g2$gene1))
   expect_identical(ng1$evaluated, FALSE)
   expect_identical(ng2$evaluated, FALSE)
   }
)

test_that("xegaGaCrossoverFactory UCross2Gene OK",
{  set.seed(2521)
   g1<-xegaGaInitGene(lFxegaGaGene)
   g1<-EvalGene(g1, lFxegaGaGene)
   g2<-xegaGaInitGene(lFxegaGaGene)
   g2<-EvalGene(g2, lFxegaGaGene)
   cCrossover<-xegaGaCrossoverFactory("UCross2Gene")
   ng<-cCrossover(g1, g2, lFxegaGaGene)
   ng1<-ng[[1]]
   ng2<-ng[[2]]
   expect_identical(all(ng1$gene1==g1$gene1), all(ng2$gene1==g2$gene1))
   expect_identical(ng1$evaluated, FALSE)
   expect_identical(ng2$evaluated, FALSE)
   }
)

test_that("xegaGaCrossoverFactory UPCross2Gene OK",
{  set.seed(3119)
   g1<-xegaGaInitGene(lFxegaGaGene)
   g1<-EvalGene(g1, lFxegaGaGene)
   g2<-xegaGaInitGene(lFxegaGaGene)
   g2<-EvalGene(g2, lFxegaGaGene)
   cCrossover<-xegaGaCrossoverFactory("UPCross2Gene")
   ng<-cCrossover(g1, g2, lFxegaGaGene)
   ng1<-ng[[1]]
   ng2<-ng[[2]]
   expect_identical(all(ng1$gene1==g1$gene1), all(ng2$gene1==g2$gene1))
   expect_identical(ng1$evaluated, FALSE)
   expect_identical(ng2$evaluated, FALSE)
   }
)

test_that("xegaGaCrossoverFactory CrossGene OK",
{  set.seed(2591)
   g1<-xegaGaInitGene(lFxegaGaGene)
   g1<-EvalGene(g1, lFxegaGaGene)
   g2<-xegaGaInitGene(lFxegaGaGene)
   g2<-EvalGene(g2, lFxegaGaGene)
   cCrossover<-xegaGaCrossoverFactory("CrossGene")
   ng<-cCrossover(g1, g2, lFxegaGaGene)
   expect_identical(ng[[1]]$evaluated, FALSE)
   }
)

test_that("xegaGaCrossoverFactory UCrossGene OK",
{  set.seed(19)
   g1<-xegaGaInitGene(lFxegaGaGene)
   g1<-EvalGene(g1, lFxegaGaGene)
   g2<-xegaGaInitGene(lFxegaGaGene)
   g2<-EvalGene(g2, lFxegaGaGene)
   cCrossover<-xegaGaCrossoverFactory("UCrossGene")
   ng<-cCrossover(g1, g2, lFxegaGaGene)
   expect_identical(ng[[1]]$evaluated, FALSE)
   }
)

test_that("xegaGaCrossoverFactory UPCrossGene OK",
{  set.seed(17)
   g1<-xegaGaInitGene(lFxegaGaGene)
   g1<-EvalGene(g1, lFxegaGaGene)
   g2<-xegaGaInitGene(lFxegaGaGene)
   g2<-EvalGene(g2, lFxegaGaGene)
   cCrossover<-xegaGaCrossoverFactory("UPCrossGene")
   ng<-cCrossover(g1, g2, lFxegaGaGene)
   expect_identical(ng[[1]]$evaluated, FALSE)
   }
)

test_that("xegaGaCrossoverFactory Cross2Gene OK",
 {
 f<-xegaGaCrossoverFactory(method="Cross2Gene")
 expect_identical(body(f), body(xegaGaGene::xegaGaCross2Gene))
}
)

test_that("xegaGaCrossoverFactory UCross2Gene OK",
 {
 f<-xegaGaCrossoverFactory(method="UCross2Gene")
 expect_identical(body(f), body(xegaGaGene::xegaGaUCross2Gene))
}
)

test_that("xegaGaCrossoverFactory UPCross2Gene OK",
 {
 f<-xegaGaCrossoverFactory(method="UPCross2Gene")
 expect_identical(body(f), body(xegaGaGene::xegaGaUPCross2Gene))
}
)

test_that("xegaGaCrossoverFactory CrossGene OK",
 {
 f<-xegaGaCrossoverFactory(method="CrossGene")
 expect_identical(body(f), body(xegaGaGene::xegaGaCrossGene))
}
)

test_that("xegaGaCrossoverFactory UCrossGene OK",
 {
 f<-xegaGaCrossoverFactory(method="UCrossGene")
 expect_identical(body(f), body(xegaGaGene::xegaGaUCrossGene))
}
)

test_that("xegaGaCrossoverFactory UPCrossGene OK",
 {
 f<-xegaGaCrossoverFactory(method="UPCrossGene")
 expect_identical(body(f), body(xegaGaGene::xegaGaUPCrossGene))
}
)

