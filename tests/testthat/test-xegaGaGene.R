

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

test_that("InitGene OK",
          { set.seed(325709)
	   g<-xegaGaInitGene(lFxegaGaGene)
	   expect_identical(g$evaluated, FALSE)
	   expect_identical(g$evalFail, FALSE)
           expect_equal(g$fit, 0)
           expect_equal(length(g$gene1), 40)
           expect_identical(1 %in% g$gene1, TRUE)
           expect_identical(0 %in% g$gene1, TRUE)
          }
)

test_that("GeneMap OK",
          {
	   g<-xegaGaInitGene(lFxegaGaGene)
           p1<- xegaGaGeneMap(g$gene1, lFxegaGaGene$penv)
           expect_equal(length(p1), 2)
           expect_identical(p1>=lFxegaGaGene$penv$lb(), rep(TRUE,2))
           expect_identical(p1<=lFxegaGaGene$penv$ub(), rep(TRUE,2))
          }
)

test_that("DecodeGene OK",
          {
	   g<-xegaGaInitGene(lFxegaGaGene)
           p1<-xegaGaDecodeGene(g, lFxegaGaGene)
           expect_equal(length(p1), 2)
           expect_identical(p1>=lFxegaGaGene$penv$lb(), rep(TRUE,2))
           expect_identical(p1<=lFxegaGaGene$penv$ub(), rep(TRUE,2))
          }
)

