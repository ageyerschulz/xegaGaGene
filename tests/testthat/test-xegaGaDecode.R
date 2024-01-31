

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
          {
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

test_that("GeneMapIdentity OK",
          {
	   g<-xegaGaInitGene(lFxegaGaGene)
           p1<- xegaGaGeneMapIdentity(g$gene1, lFxegaGaGene$penv)
           expect_identical(p1, g$gene1)
          }
)

test_that("Gray2Bin OK",
          {
           expect_identical(Gray2Bin(c(1, 0, 0, 0)), c(1, 1, 1, 1))
           expect_identical(Gray2Bin(c(0, 1, 0)), c(0, 1, 1))
           expect_identical(Gray2Bin(c(1, 0)), c(1, 1))
           expect_identical(Gray2Bin(c(0, 0, 0, 1, 0)), c(0, 0, 0, 1, 1))
          }
)

test_that("GeneMapGray OK",
          {
	   g<-xegaGaInitGene(lFxegaGaGene)
           p1<- xegaGaGeneMapGray(g$gene1, lFxegaGaGene$penv)
           expect_equal(length(p1), 2)
           expect_identical(p1>=lFxegaGaGene$penv$lb(), rep(TRUE,2))
           expect_identical(p1<=lFxegaGaGene$penv$ub(), rep(TRUE,2))
          }
)

test_that("without OK",
          {
	   a<-sample(1:15,15, replace=FALSE)
	   b<-sample(1:15, 3, replace=FALSE)
	   c<-without(a, b)
	   d<-without(c, b)
	   e<-without(c, a)
           expect_identical(sum(c), (sum(a)-sum(b)))
           expect_identical(d, c)
           expect_identical(e, integer(0))
          }
)

test_that("GeneMapPerm OK",
          {
           newlF<-lFxegaGaGene
	   parm<-function(x){function() {return(x)}}
           newlF$penv$bitlength<-parm(rep(40, 10))	  
           newlF$penv$genelength<-parm(sum(rep(40, 10)))	  
	   g<-xegaGaInitGene(newlF)
           p1<- xegaGaGeneMapPerm(g$gene1, newlF$penv)
           l<-length(newlF$penv$bitlength())
           expect_equal(length(p1), l)
           expect_equal(sort(p1, decreasing =FALSE), (1:l))
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

test_that("xegaGaGeneMapFactory OK",
          {lF<-lFxegaGaGene
          GM<-xegaGaGeneMapFactory()
          expect_identical(body(GM), body(xegaGaGeneMap))
          GM<-xegaGaGeneMapFactory("Bin2Dec")
          expect_identical(body(GM), body(xegaGaGeneMap))
          GM<-xegaGaGeneMapFactory("Gray2Dec")
          expect_identical(body(GM), body(xegaGaGeneMapGray))
          GM<-xegaGaGeneMapFactory("Identity")
          expect_identical(body(GM), body(xegaGaGeneMapIdentity))
          GM<-xegaGaGeneMapFactory("Permutation")
          expect_identical(body(GM), body(xegaGaGeneMapPerm))
          expect_error(xegaGaGeneMapFactory("Stchastic"))
          }
)

