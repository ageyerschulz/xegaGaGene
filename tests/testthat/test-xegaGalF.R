

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
           expect_equal(lFxegaGaGene$CrossRate(), 0.2)
           expect_equal(lFxegaGaGene$UCrossSwap(), 0.2)
           expect_equal(lFxegaGaGene$Max(), 1)
           expect_equal(lFxegaGaGene$Offset(), 1)
           expect_equal(lFxegaGaGene$Eps(), 0.01)
           expect_identical(lFxegaGaGene$Elitist(), TRUE)
           expect_equal(lFxegaGaGene$TournamentSize(), 2)
          }
)

