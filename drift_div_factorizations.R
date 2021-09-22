library(dplyr)
library(flashier)
library(tidyverse)
library(ebnm)
library(fields)
library(dyno)
library(dyntoy)
library(dequer)

#just a helper for debugging
get_sums_by_label <- function(tree,loading){
  summary <- as.data.frame(table(tree$csv$Labels))
  summary$NumActivated <- 0
  idx <- 1
  for(label in rownames(summary)){
    start_idx <- idx
    idx <- idx + summary[label,"Freq"]
    summary[label,"NumActivated"] <- sum(loading[start_idx:(idx-1)])
  }
  return(summary)
}

add_factor <- function(dat,loading,fl,prior,Fprior){
  K <- fl$n.factors
  ls.soln  <- t(crossprod(loading,  dat - fitted(fl))/sum(loading))
  EF <- list(loading, ls.soln)
  next_fl <- fl %>%
    flash.init.factors(
      EF,
      prior.family = c(prior,Fprior)
    ) %>%
    flash.fix.loadings(kset = K + 1, mode = 1L, is.fixed = (EF[[1]] == 0)) %>%
    flash.backfit(kset = K + 1)
  return(next_fl)
}

get_divergence_factor <- function(dat,loading,fl,divprior,Fprior){
  K <- fl$n.factors
  ls.soln  <- t(crossprod(loading,  dat - fitted(fl))/sum(loading))
  EF <- list(loading, ls.soln)
  next_fl <- fl %>%
    flash.init.factors(
      EF,
      prior.family = c(divprior,Fprior)
    ) %>%
    flash.fix.loadings(kset = K + 1, mode = 1L, is.fixed = (EF[[1]] == 0)) %>%
    flash.backfit(kset = K + 1)
  return(next_fl$loadings.pm[[1]][,K+1])
}

drift_fit <- function(tree,dat,filename,
                    divprior = prior.point.laplace(),
                    driftprior = as.prior(ebnm.fn = ebnm_point_exponential,sign=1),
                    Fprior = prior.normal(),
                    Kmax = Inf,
                    min_pve = 0,
                    verbose.lvl = 0,
                    eps=1e-2) {
  #the first loading will be the all-ones vector
  ones <- matrix(1, nrow = nrow(dat), ncol = 1)
  #first factor will be least sq soln: argmin_f ||Y - ones t(f)||_F^2
  ls.soln <- t(crossprod(ones, dat)/nrow(dat))

  #create flash object with initial drift loading and initial divergence loading
  fl <- flash.init(dat) %>%
    flash.set.verbose(verbose.lvl) %>%
    #initialize L to be the ones vector, and F to be the least squares solution
    flash.init.factors(list(ones, ls.soln),
                       prior.family = c(driftprior,Fprior)) %>%
    #only fixing the first factor, and the we want to fix row loadings, so mode=1
    flash.fix.loadings(kset = 1, mode = 1) %>%
    #backfit to match the priors
    flash.backfit() %>%
    #add initial divergence loading
    flash.add.greedy(
      Kmax = 1,
      prior.family = c(divprior, Fprior)
    )

  #add divergence factor to a queue (Breadth-first)
  divergence_queue <- queue()
  pushback(divergence_queue,fl$loadings.pm[[1]][,2])
  #remove divergence factor
  fl <- fl %>%
    flash.remove.factors(kset = 2)
  K <- 1

  while(length(divergence_queue) > 0 && K < Kmax) {
    #pop the first divergence off the queue
    current_divergence <- pop(divergence_queue)
    #split into loadings for positive and negative parts (1-0 indicator vectors)
    splus <- matrix(1L * (current_divergence > eps), ncol = 1)
    if (sum(splus) > 0) {
      if (verbose.lvl > 0) {
        cat(get_sums_by_label(tree,splus))
      }
      #add drift loading
      next_fl <- add_factor(dat,splus,fl,driftprior,Fprior)
      if (next_fl$pve[K + 1] > min_pve){
        fl <- next_fl
        K <- fl$n.factors
        #enqueue new divergence
        new_div <- get_divergence_factor(dat,splus,fl,divprior,Fprior)
        pushback(divergence_queue,new_div)
      }
    }
    sminus <- matrix(1L * (current_divergence < -eps), ncol = 1)
    if (sum(sminus) > 0 && K < Kmax) {
      if (verbose.lvl > 0) {
        cat(get_sums_by_label(tree,splus))
      }
      #add drift loading
      next_fl <- add_factor(dat,sminus,fl,driftprior,Fprior)
      if (next_fl$pve[K + 1] > min_pve){
        fl <- next_fl
        K <- fl$n.factors
        #enqueue new divergence
        new_div <- get_divergence_factor(dat,sminus,fl,divprior,Fprior)
        pushback(divergence_queue,new_div)
      }
    }

    if (verbose.lvl > 0) {
      cat("K:", K, "\n")
    }

    image.plot(fl$loadings.pm[[1]],xlab = fl$n.factors)
  }

  L <- fl$loadings.pm[[1]]
  F <- fl$loadings.pm[[2]]
  scale <- fl$loadings.scale

  write.table(L,file=paste(filename,"L.csv",sep=''),sep=',')
  write.table(F,file=paste(filename,"F.csv",sep=''),sep=',')
  write.table(scale,file=paste(filename,"scale.csv",sep=''),sep=',')
  write.table(fl$pve,file=paste(filename,"pve.csv",sep=''),sep=',')

  return(fl)
}

div_fit <- function(tree,dat,filename,
                      divprior = prior.point.laplace(),
                      driftprior = as.prior(ebnm.fn = ebnm_point_exponential,sign=1),
                      Fprior = prior.normal(),
                      Kmax = Inf,
                      min_pve = 0,
                      verbose.lvl = 0,
                      eps=1e-2) {
  #the first loading will be the all-ones vector
  ones <- matrix(1, nrow = nrow(dat), ncol = 1)
  #first factor will be least sq soln: argmin_f ||Y - ones t(f)||_F^2
  ls.soln <- t(crossprod(ones, dat)/nrow(dat))

  #create flash object with initial drift loading and initial divergence loading
  fl <- flash.init(dat) %>%
    flash.set.verbose(verbose.lvl) %>%
    #initialize L to be the ones vector, and F to be the least squares solution
    flash.init.factors(list(ones, ls.soln),
                       prior.family = c(driftprior,Fprior)) %>%
    #only fixing the first factor, and the we want to fix row loadings, so mode=1
    flash.fix.loadings(kset = 1, mode = 1) %>%
    #backfit to match the priors
    flash.backfit() %>%
    #add initial divergence loading
    flash.add.greedy(
      Kmax = 1,
      prior.family = c(divprior, Fprior)
    )

  #add divergence factor to a queue (Breadth-first)
  divergence_queue <- queue()
  pushback(divergence_queue,fl$loadings.pm[[1]][,2])

  K <- fl$n.factors

  while(length(divergence_queue) > 0 && K < Kmax) {
    #pop the first divergence off the queue
    current_divergence <- pop(divergence_queue)
    
    #add drift loading
    snonzero <- matrix(1L * (abs(current_divergence) > eps), ncol = 1)
    if (sum(snonzero) > 0 && (K != 2)) {
      if (verbose.lvl > 0) {
        cat(get_sums_by_label(tree,snonzero))
      }
      #add drift loading
      next_fl <- add_factor(dat,snonzero,fl,driftprior,Fprior)
      if (next_fl$pve[K + 1] > min_pve){
        fl <- next_fl
        K <- fl$n.factors
      }
    }
    #try to add a divergence within the positive loadings of the current divergence
    splus <- matrix(1L * (current_divergence > eps), ncol = 1)
    if (sum(splus) > 0 && K < Kmax) {
      #add divergence loading
      next_fl <- add_factor(dat,splus,fl,divprior,Fprior)
      if (next_fl$pve[K + 1] > min_pve){
        fl <- next_fl
        K <- fl$n.factors
        new_div <- fl$loadings.pm[[1]][,K]
        #enqueue new divergence
        pushback(divergence_queue,new_div)
      }
    }
    #try to add a divergence within the negative loadings of the current divergence
    sminus <- matrix(1L * (current_divergence < -eps), ncol = 1)
    if (sum(sminus) > 0 && K < Kmax) {
      #add divergence loading
      next_fl <- add_factor(dat,sminus,fl,divprior,Fprior)
      if (next_fl$pve[K + 1] > min_pve){
        fl <- next_fl
        K <- fl$n.factors
        new_div <- fl$loadings.pm[[1]][,K]
        #enqueue new divergence
        pushback(divergence_queue,new_div)
      }
    }

    if (verbose.lvl > 0) {
      cat("K:", K, "\n")
    }
    image.plot(fl$loadings.pm[[1]],xlab = fl$n.factors)
  }

  L <- fl$loadings.pm[[1]]
  F <- fl$loadings.pm[[2]]
  scale <- fl$loadings.scale

  write.table(L,file=paste(filename,"L.csv",sep=''),sep=',')
  write.table(F,file=paste(filename,"F.csv",sep=''),sep=',')
  write.table(scale,file=paste(filename,"scale.csv",sep=''),sep=',')
  write.table(fl$pve,file=paste(filename,"pve.csv",sep=''),sep=',')

  return(fl)
}

form_tree_from_file <- function(filename){
  tree <- vector(mode="list")
  tree$csv <- read.csv(filename,row.names=1)
  tree$raw <- tree$csv %>%
    select(Raw0:Raw499)
  tree$matrix <- as.matrix(tree$raw)
  tree$dimred <- tree$csv %>%
    select(tsne0:tsne1)
  #As input, dynwrap requires raw counts and normalized (log2) expression data.
  tree$counts <- round(2**(tree$raw)-1)
  tree$counts[tree$counts<0] <- 0
  tree$dataset <- wrap_expression(
    expression = tree$matrix,
    counts = as.matrix(tree$counts)
  )
  tree$dataset <- add_prior_information(
    tree$dataset,
    start_id = "Row0",
    start_n = 1,
    end_n = 4,
  )
  tree$trajectory <- vector(mode="list")
  return(tree)
}

form_tree_from_dataset <- function(dataset){
  tree <- vector(mode="list")
  tree$dataset <- dataset
  tree$matrix <- as.matrix(dataset$expression)
  tree$trajectory <- vector(mode="list")
  return(tree)
}

run_methods <- function(tree,outfile,Kmax,eps){
  #drift factorization method
  tree$trajectory$drift <- drift_fit(tree,
                                     tree$matrix,
                                     paste(outfile,"drift",sep=""),
                                     Kmax = Kmax,
                                     eps = eps)
  tree$trajectory$div <- div_fit(tree,
                                     tree$matrix,
                                     paste(outfile,"div",sep=""),
                                     Kmax = Kmax,
                                     eps = eps)
  return(tree)
}

#### FUNCTION DEFINITIONS STOP HERE
#trajectory trees
tree1 <- form_tree_from_file('TITesting.nosync/tree1/tree1.csv')
tree1 <- run_methods(tree1,'TITesting.nosync/tree1/EBMFfactors/',Kmax = 10,eps = 1e-2)
# image.plot(cov(t(fitted(tree1$trajectory$drift))))
# image.plot(cov(t(tree1$matrix)))

tree2 <- form_tree_from_file('TITesting.nosync/tree2/tree2.csv')
tree2 <- run_methods(tree2,'TITesting.nosync/tree2/EBMFfactors/',Kmax = 10,eps = 1e-2)

tree3 <- form_tree_from_file('TITesting.nosync/tree3/tree3.csv')
tree3 <- run_methods(tree3,'TITesting.nosync/tree3/EBMFfactors/',Kmax=20,eps = 1e-2)

tree4 <- form_tree_from_file('TITesting.nosync/tree4/tree4.csv')
tree4 <- run_methods(tree4,'TITesting.nosync/tree4/EBMFfactors/',Kmax=10,eps = 1e-2)

#node trees
nodetree1 <- form_tree_from_file('TITesting.nosync/NodeTree1/NodeTree1.csv')
nodetree1 <- run_methods(nodetree1,'TITesting.nosync/NodeTree1/EBMFfactors/',Kmax=10,eps = 2e-2)

nodetree2 <- form_tree_from_file('TITesting.nosync/NodeTree2/NodeTree2.csv')
nodetree2 <- run_methods(nodetree2,'TITesting.nosync/NodeTree2/EBMFfactors/',Kmax=10,eps = 2e-2)

nodetree3 <- form_tree_from_file('TITesting.nosync/NodeTree3/NodeTree3.csv')
nodetree3 <- run_methods(nodetree3,'TITesting.nosync/NodeTree3/EBMFfactors/',Kmax=10,eps = 2e-2)
#stops early

nodetree4 <- form_tree_from_file('TITesting.nosync/NodeTree4/NodeTree4.csv')
nodetree4 <- run_methods(nodetree4,'TITesting.nosync/NodeTree4/EBMFfactors/',Kmax=20,eps = 2e-2)

nodetree5 <- form_tree_from_file('TITesting.nosync/NodeTree5/NodeTree5.csv')
nodetree5 <- run_methods(nodetree5,'TITesting.nosync/NodeTree5/EBMFfactors/',Kmax=20,eps = 2e-2)

nodetree6 <- form_tree_from_file('TITesting.nosync/NodeTree6/NodeTree6.csv')
nodetree6 <- run_methods(nodetree6,'TITesting.nosync/NodeTree6/EBMFfactors/',Kmax=20,eps = 2e-2)
#stops early with eps 2e-2, gets one more factor with 1e-2 but still stops early

nodetree7 <- form_tree_from_file('TITesting.nosync/NodeTree7/NodeTree7.csv')
nodetree7 <- run_methods(nodetree7,'TITesting.nosync/NodeTree7/EBMFfactors/',Kmax=40,eps = 1e-2)

nodetree8 <- form_tree_from_file('TITesting.nosync/NodeTree8/NodeTree8.csv')
nodetree8 <- run_methods(nodetree8,'TITesting.nosync/NodeTree8/EBMFfactors/',Kmax=40,eps = 1e-2)

nodetree9 <- form_tree_from_file('TITesting.nosync/NodeTree9/NodeTree9.csv')
nodetree9 <- run_methods(nodetree9,'TITesting.nosync/NodeTree9/EBMFfactors/',Kmax=40,eps = 1e-2)

#dynverse dataset
dynverse_tree <- generate_dataset(
  model = model_bifurcating(),
  num_cells = 700,
  num_features = 500
)
dynverse_tree <- form_tree_from_dataset(dynverse_tree)
write.table(dynverse_tree$matrix,file='TITesting.nosync/dynversetree/dynversetree.csv',row.names=FALSE,sep=',')
dynverse_tree <- run_methods(dynverse_tree,'TITesting.nosync/dynversetree/EBMFfactors/',Kmax=40,eps = 1e-2)

#issues/questions with drift_fit
  #1 what to choose eps to be? choose eps adaptively? Better way to force more sparsity than brute forcing this?
  # probably makes it error on the side of assuming a cell is less specialized than it actually is-- good or bad
  # sometimes eps is still too small, not zeroing out enough
  #2 stopping point-- it doesn't stop adding factors when it shouldn't, unfortunately. 
  # It does with node data, but not trajectory data
  # maybe some sort of spatial way of determining stopping? Like is this divergence you drawing a line in the sand,
  # or is the split point close enough to the ancestral point that it's a reasonable split
  #3 what happens with non-tree-shaped data?
  #4 seems to be larger numbers for more specialized samples for each drift factor. Is that what we want?
  # It doesn't quite match up with our proofs, but maybe another proof can give us that
  #5 works well with node datasets up to about 3 levels of splits. Breaks down a bit for 4 levels
  #6 when and how to backfit?
#issues with div_fit
  # 1 doesn't find the grandchildren in NodeTree3, i.e. when standard deviation is high, can only find one split
  # divergence fit seems less reliable in general, not exactly sure why yet
