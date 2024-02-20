

library(testthat)
library(xegaSelectGene)
library(xegaGaGene)

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

