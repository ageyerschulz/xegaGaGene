

library(testthat)
library(xegaSelectGene)
library(xegaGaGene)

#
# MutateGene
#

test_that("MutateGene: g:  not evaluated, not mutated. OK",
{ 
   g<-xegaGaInitGene(lFxegaGaGene)
   lFxegaGaGene$BitMutationRate1<-parm(0.0)
   ng<-xegaGaMutateGene(g, lFxegaGaGene)
   expect_identical(all(g$gene1==ng$gene1), TRUE)
   expect_equal(ng$fit, 0)
   expect_equal(ng$evaluated, FALSE)
  }
)

test_that("MutateGene: g:  not evaluated, mutated. OK",
  { set.seed(3257)
   g<-xegaGaInitGene(lFxegaGaGene)
   lFxegaGaGene$BitMutationRate1<-parm(1.0)
   ng<-xegaGaMutateGene(g, lFxegaGaGene)
   expect_identical(all(g$gene1==ng$gene1), FALSE)
   expect_equal(ng$fit, 0)
   expect_identical(ng$evaluated, FALSE)
  }
)

test_that("MutateGene: g:  evaluated, mutated. OK",
  {  set.seed(32057)
   g<-xegaGaInitGene(lFxegaGaGene)
   g<-EvalGene(g, lFxegaGaGene)
   lFxegaGaGene$BitMutationRate1<-parm(1.0)
   ng<-xegaGaMutateGene(g, lFxegaGaGene)
   expect_identical(all(g$gene1==ng$gene1), FALSE)
   expect_identical(!ng$fit==0, TRUE)
   expect_identical(ng$evaluated, FALSE)
  }
)

test_that("MutateGene: g:  evaluated, not mutated. OK",
  { 
   g<-xegaGaInitGene(lFxegaGaGene)
   g<-EvalGene(g, lFxegaGaGene)
   lFxegaGaGene$BitMutationRate1<-parm(0.0)
   ng<-xegaGaMutateGene(g, lFxegaGaGene)
   expect_identical(all(g$gene1==ng$gene1), TRUE)
   expect_identical(!ng$fit==0, TRUE)
   expect_identical(ng$evaluated, TRUE)
  }
)

#
# Individually Variable Adaptive Mutate Gene
#

test_that("IVAdaptiveMutateGene: g: not eval, not mutated, low fit. OK",
{ 
   g<-xegaGaInitGene(lFxegaGaGene)
   lFxegaGaGene$BitMutationRate1<-parm(0.0)
   lFxegaGaGene$BitMutationRate2<-parm(0.0)
   ng<-xegaGaIVAdaptiveMutateGene(g, lFxegaGaGene)
   expect_identical(all(g$gene1==ng$gene1), TRUE)
   expect_equal(ng$fit, 0)
   expect_equal(ng$evaluated, FALSE)
  }
)

test_that("IVAdaptiveMutateGene: g: not eval, mutated, low fit. OK",
{ set.seed(2833) 
   g<-xegaGaInitGene(lFxegaGaGene)
   lFxegaGaGene$BitMutationRate1<-parm(0.0)
   lFxegaGaGene$BitMutationRate2<-parm(1.0)
   ng<-xegaGaIVAdaptiveMutateGene(g, lFxegaGaGene)
   expect_identical(all(g$gene1==ng$gene1), FALSE)
   expect_equal(ng$fit, 0)
   expect_identical(ng$evaluated, FALSE)
  }
)

test_that("IVAdaptiveMutateGene: g: eval, mutated, low fit. OK",
{  set.seed(4013)
   g<-xegaGaInitGene(lFxegaGaGene)
   g<-EvalGene(g, lFxegaGaGene)
   lFxegaGaGene$BitMutationRate<-parm(0.0)
   lFxegaGaGene$BitMutationRate2<-parm(1.0)
   ng<-xegaGaIVAdaptiveMutateGene(g, lFxegaGaGene)
   expect_identical(all(g$gene1==ng$gene1), FALSE)
   expect_identical(!(ng$fit==0), TRUE)
   expect_identical(ng$evaluated, FALSE)
  }
)

test_that("IVAdaptiveMutateGene: g: eval, not mutated, low fit. OK",
{  
   g<-xegaGaInitGene(lFxegaGaGene)
   g<-EvalGene(g, lFxegaGaGene)
   lFxegaGaGene$BitMutationRate1<-parm(0.0)
   lFxegaGaGene$BitMutationRate2<-parm(0.0)
   ng<-xegaGaIVAdaptiveMutateGene(g, lFxegaGaGene)
   expect_identical(all(g$gene1==ng$gene1), TRUE)
   expect_identical(!(ng$fit==0), TRUE)
   expect_identical(ng$evaluated, TRUE)
  }
)

#### Problems! Not stable
test_that("IVAdaptiveMutateGene: g: eval, not mutated, high fit. OK",
{ set.seed(5651)
   g<-xegaGaInitGene(lFxegaGaGene)
   g<-EvalGene(g, lFxegaGaGene)
   lFxegaGaGene$CutoffFit<-parm(0.000001)
   lFxegaGaGene$BitMutationRate1<-parm(0.0)
   lFxegaGaGene$BitMutationRate2<-parm(0.5)
   ng<-xegaGaIVAdaptiveMutateGene(g, lFxegaGaGene)
   expect_identical(all(g$gene1==ng$gene1), TRUE)
   expect_identical(!(ng$fit==0), TRUE)
   expect_identical(ng$evaluated, TRUE)
  }
)

#### Problems! Not stable
test_that("IVAdaptiveMutateGene: g: eval, mutated, high fit. OK",
{  set.seed(499)
   g<-xegaGaInitGene(lFxegaGaGene)
   g<-EvalGene(g, lFxegaGaGene)
   lFxegaGaGene$CutoffFit<-parm(0.000001)
   lFxegaGaGene$BitMutationRate1<-parm(1.0)
   lFxegaGaGene$BitMutationRate2<-parm(0.5)
   ng<-xegaGaIVAdaptiveMutateGene(g, lFxegaGaGene)
   expect_identical(all(g$gene1==ng$gene1), FALSE)
   expect_identical(!(ng$fit==0), TRUE)
   expect_identical(ng$evaluated, FALSE)
  }
)

#
# xegaGaMutationFactory
#

test_that("xegaGaMutationFactory() OK",
{  set.seed(29)
   g<-xegaGaInitGene(lFxegaGaGene)
   g<-EvalGene(g, lFxegaGaGene)
   lFxegaGaGene$BitMutationRate1<-parm(1.0)
   cMutate<-xegaGaMutationFactory()
   ng<-cMutate(g, lFxegaGaGene)
   expect_identical(all(g$gene1==ng$gene1), FALSE)
   expect_identical(!ng$fit==0, TRUE)
   expect_identical(ng$evaluated, FALSE)
   }
)

test_that("xegaGaMutationFactory  MutateGene OK",
{  set.seed(84871)
   g<-xegaGaInitGene(lFxegaGaGene)
   g<-EvalGene(g, lFxegaGaGene)
   lFxegaGaGene$BitMutationRate1<-parm(1.0)
   cMutate<-xegaGaMutationFactory(method="MutateGene")
   ng<-cMutate(g, lFxegaGaGene)
   expect_identical(all(g$gene1==ng$gene1), FALSE)
   expect_identical(!ng$fit==0, TRUE)
   expect_identical(ng$evaluated, FALSE)
   }
)

test_that("xegaGaMutationFactory  IVM OK",
{  set.seed(2053)
   g<-xegaGaInitGene(lFxegaGaGene)
   g<-EvalGene(g, lFxegaGaGene)
   lFxegaGaGene$BitMutationRate1<-parm(1.0)
   lFxegaGaGene$BitMutationRate2<-parm(1.0)
   cMutate<-xegaGaMutationFactory(method="IVM")
   ng<-cMutate(g, lFxegaGaGene)
   expect_identical(all(g$gene1==ng$gene1), FALSE)
   expect_identical(!ng$fit==0, TRUE)
   expect_identical(ng$evaluated, FALSE)
   }
)

test_that("xegaGaMutationFactory HUGO OK",
          {
           expect_error(xegaGaMutationFactory("HUGO"))
          }
)

