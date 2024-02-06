library(testthat)
library(xegaSelectGene)
library(xegaGaGene)

test_that("Replicate2Gene OK",
          {
          lFxegaGaGene$CrossGene<-xegaGaCross2Gene
          # Case 0: P(Mutation)=P(Crossover)=0.5
          lFxegaGaGene$MutationRate<-function(fit, lF) {0.5}
          lFxegaGaGene$CrossRate<-function(fit, lF) {0.5}
          pop10<-lapply(rep(0,10), function(x) xegaGaInitGene(lFxegaGaGene))
          epop10<-lapply(pop10, lFxegaGaGene$EvalGene, lF=lFxegaGaGene)
          fit10<-unlist(lapply(epop10, function(x) {x$fit}))
          newgenes<-xegaGaReplicate2Gene(pop10, fit10, lFxegaGaGene)
	  expect_equal(((length(newgenes)==1)|(length(newgenes)== 2)), TRUE)
          newgenes<-xegaGaReplicate2Gene(pop10, fit10, lFxegaGaGene)
	  expect_equal(((length(newgenes)==1)|(length(newgenes)== 2)), TRUE)
          newgenes<-xegaGaReplicate2Gene(pop10, fit10, lFxegaGaGene)
	  expect_equal(((length(newgenes)==1)|(length(newgenes)== 2)), TRUE)
          # Case 1: P(Mutation)=0.0 P(Crossover)=1.0
          lFxegaGaGene$MutationRate<-function(fit, lF) {0.0}
          lFxegaGaGene$CrossRate<-function(fit, lF) {1.0}
          newgenes<-xegaGaReplicate2Gene(pop10, fit10, lFxegaGaGene)
	  expect_equal(length(newgenes), 2)
          # Case 2: P(Mutation)=1.0 P(Crossover)=1.0
          lFxegaGaGene$MutationRate<-function(fit, lF) {1.0}
          lFxegaGaGene$CrossRate<-function(fit, lF) {1.0}
          newgenes<-xegaGaReplicate2Gene(pop10, fit10, lFxegaGaGene)
	  expect_equal(length(newgenes), 2)
          # Case 3: P(Mutation)=1.0 P(Crossover)=0.0
          lFxegaGaGene$MutationRate<-function(fit, lF) {1.0}
          lFxegaGaGene$CrossRate<-function(fit, lF) {0.0}
          newgenes<-xegaGaReplicate2Gene(pop10, fit10, lFxegaGaGene)
	  expect_equal(length(newgenes), 1)
          # Case 4: P(Mutation)=0.0 P(Crossover)=0.0
          lFxegaGaGene$MutationRate<-function(fit, lF) {0.0}
          lFxegaGaGene$CrossRate<-function(fit, lF) {0.0}
          newgenes<-xegaGaReplicate2Gene(pop10, fit10, lFxegaGaGene)
	  expect_equal(length(newgenes), 1)
          }
)

test_that("ReplicateGene OK",
{
AcceptNewGene<-function(OperatorPipeline, gene, lF)
{
newGene<-OperatorPipeline(gene, lF)
  if (newGene$evalFail) {
          cat("AcceptNewGene: Evaluation failed. Gene: \n")
          print(newGene)
  }
return(newGene)}
          lFxegaGaGene$CrossGene<-xegaGaCrossGene
          pop10<-lapply(rep(0,10), function(x) xegaGaInitGene(lFxegaGaGene))
          epop10<-lapply(pop10, lFxegaGaGene$EvalGene, lF=lFxegaGaGene)
          fit10<-unlist(lapply(epop10, function(x) {x$fit}))
          lFxegaGaGene$Accept<-AcceptNewGene
          # Case 0: P(Mutation)=0.5 P(Crossover)=0.5
          lFxegaGaGene$MutationRate<-function(fit, lF) {0.5}
          lFxegaGaGene$CrossRate<-function(fit, lF) {0.5}
          newgenes<-xegaGaReplicateGene(pop10, fit10, lFxegaGaGene)
	  expect_equal(((length(newgenes)==1)|(length(newgenes)== 2)), TRUE)
          newgenes<-xegaGaReplicateGene(pop10, fit10, lFxegaGaGene)
	  expect_equal(((length(newgenes)==1)|(length(newgenes)== 2)), TRUE)
          newgenes<-xegaGaReplicateGene(pop10, fit10, lFxegaGaGene)
	  expect_equal(((length(newgenes)==1)|(length(newgenes)== 2)), TRUE)
          # Case 1: P(Mutation)=0.0 P(Crossover)=1.0
          lFxegaGaGene$MutationRate<-function(fit, lF) {0.0}
          lFxegaGaGene$CrossRate<-function(fit, lF) {1.0}
          newgenes<-xegaGaReplicateGene(pop10, fit10, lFxegaGaGene)
	  expect_equal(((length(newgenes)==1)|(length(newgenes)== 2)), TRUE)
          # Case 2: P(Mutation)=1.0 P(Crossover)=1.0
          lFxegaGaGene$MutationRate<-function(fit, lF) {1.0}
          lFxegaGaGene$CrossRate<-function(fit, lF) {1.0}
          newgenes<-xegaGaReplicateGene(pop10, fit10, lFxegaGaGene)
	  expect_equal(length(newgenes), 1)
          # Case 3: P(Mutation)=1.0 P(Crossover)=0.0
          lFxegaGaGene$MutationRate<-function(fit, lF) {1.0}
          lFxegaGaGene$CrossRate<-function(fit, lF) {0.0}
          newgenes<-xegaGaReplicateGene(pop10, fit10, lFxegaGaGene)
	  expect_equal(length(newgenes), 1)
          # Case 4: P(Mutation)=0.0 P(Crossover)=0.0
          lFxegaGaGene$MutationRate<-function(fit, lF) {0.0}
          lFxegaGaGene$CrossRate<-function(fit, lF) {0.0}
          newgenes<-xegaGaReplicateGene(pop10, fit10, lFxegaGaGene)
	  expect_equal(length(newgenes), 1)
          }
)

test_that("ReplicationFactory OK",
          {
          Fun<-xegaGaReplicationFactory()
          expect_identical(body(Fun), body(xegaGaReplicateGene))
          Fun<-xegaGaReplicationFactory("Kid1")
          expect_identical(body(Fun), body(xegaGaReplicateGene))
          Fun<-xegaGaReplicationFactory("Kid2")
          expect_identical(body(Fun), body(xegaGaReplicate2Gene))
          expect_error(xegaGaReplicationFactory("Stchastic"))
          }
)

