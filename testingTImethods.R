setwd("~/Documents/CAM/Year1/GeneticTrees")
library(dplyr)
library(flashier)
library(tidyverse)
library(ebnm)
library(dyno)
library(dyntoy)

#Jason's EBMF methods
div_cov_fit <- function(covmat, filename, prior = prior.point.laplace(), Kmax = 1000) {
  fl <- div_fit(covmat, filename, prior, Kmax)
  s2 <- max(0, mean(diag(covmat) - diag(fitted(fl))))
  s2_diff <- Inf
  while(s2 > 0 && abs(s2_diff - 1) > 1e-4) {
    covmat_minuss2 <- covmat - diag(rep(s2, ncol(covmat)))
    fl <- div_fit(covmat_minuss2,filename, prior, Kmax)
    old_s2 <- s2
    s2 <- max(0, mean(diag(covmat) - diag(fitted(fl))))
    s2_diff <- s2 / old_s2
  }

  fl$ebcovmf_s2 <- s2

  return(fl)
}

div_fit <- function(dat,filename, prior = prior.point.laplace(), Kmax = Inf, min_pve = 0, verbose.lvl = 0) {
  #the first loading will be the all-ones vector
  ones <- matrix(1, nrow = nrow(dat), ncol = 1)
  #first factor will be least sq soln: argmin_f ||Y - ones t(f)||_F^2
  ls.soln <- t(crossprod(ones, dat)/nrow(dat))

  fl <- flash.init(dat) %>%
    flash.set.verbose(verbose.lvl) %>%
    #initialize L to be the ones vector, and F to be the least squares solution
    flash.init.factors(list(ones, ls.soln)) %>%
    #only fixing the first factor, and the we want to fix row loadings, so mode=1
    flash.fix.loadings(kset = 1, mode = 1) %>%
    #backfit to match the priors
    flash.backfit() %>%
    #add anoter factor
    flash.add.greedy(
      Kmax = 1,
      #specified prior on L, and a normal distribution on F
      prior.family = c(prior, prior.normal())
    )

  current_k <- 2
  K <- 2

  while(current_k <= K && K < Kmax) {
    print(current_k)
    #split into loadings for positive and negative parts (1-0 indicator vectors)
    splus <- matrix(1L * (fl$loadings.pm[[1]][, current_k] > 0), ncol = 1)
    sminus <- matrix(1L * (fl$loadings.pm[[1]][, current_k] < 0), ncol = 1)

    if (sum(splus) > 0 && sum(sminus) > 0) {
      #lst sq soln for positive and negative factors:
      # argmin_f ||(Y-sum lk fk) - splus t(f)||_F^2
      # argmin_f ||(Y-sum lk fk) - sminus t(f)||_F^2
      ls.soln.plus  <- t(crossprod(splus,  dat - fitted(fl))/sum(splus))
      ls.soln.minus <- t(crossprod(sminus, dat - fitted(fl))/sum(sminus))

      #initializations of new loadings
      EF <- list(cbind(splus, sminus), cbind(ls.soln.plus, ls.soln.minus))

      next_fl <- fl %>%
        #initialize new loadings
        flash.init.factors(EF) %>%
        flash.fix.loadings(kset = K + 1:2, mode = 1L, is.fixed = (EF[[1]] == 0)) %>%
        flash.backfit(kset = K + 1:2)
      if (any(next_fl$pve[K + 1:2] > min_pve)) {
        fl <- next_fl
      }
    }

    current_k <- current_k + 1
    K <- fl$n.factors
    if (verbose.lvl > 0) {
      cat("K:", K, "\n")
    }
  }

  fl$loadings.lfsr[[1]][, 1] <- 0
  fl$loadings.lfsr[[1]][is.na(fl$loadings.lfsr[[1]])] <- 1

  L <- fl$loadings.pm[[1]]
  F <- fl$loadings.pm[[2]]
  scale <- fl$loadings.scale

  write.table(L,file=paste(filename,"L.csv",sep=''),sep=',')
  write.table(F,file=paste(filename,"F.csv",sep=''),sep=',')
  write.table(scale,file=paste(filename,"scale.csv",sep=''),sep=',')
  write.table(fl$pve,file=paste(filename,"pve.csv",sep=''),sep=',')

  return(fl)
}

# Reproduces the guidelines as created in the shiny app
answers <- dynguidelines::answer_questions(
  multiple_disconnected = FALSE,
  expect_topology = TRUE,
  expected_topology = "tree",
  n_cells = 696,
  n_features = 500,
  time = "10m",
  memory = "5GB",
  prior_information = c("start_id", "end_n", "start_n"),
  method_selection = "fixed_n_methods",
  fixed_n_methods = 10,
  docker = TRUE
)
guidelines <- dynguidelines::guidelines(answers = answers)

run_methods <- function(tree,outfile){
  tree$trajectory <- vector(mode="list")
  #jason's tree method
  tree$trajectory$ebmf_cov <- div_cov_fit(cov(t(tree$matrix)),paste(outfile,'cov',sep=''),Kmax = 30)
  tree$trajectory$ebmf <- div_fit(tree$matrix,outfile,Kmax = 30)
  #recommended methods from dynverse that worked
  tree$trajectory$mst <- infer_trajectory(tree$dataset,"mst",verbose=TRUE)
  tree$trajectory$slingshot <- infer_trajectory(tree$dataset,"slingshot",verbose=TRUE)
  tree$trajectory$elpigraph <- infer_trajectory(tree$dataset,"elpigraph",verbose=TRUE)
  tree$trajectory$sincell <- infer_trajectory(tree$dataset,"sincell",verbose=TRUE)
  tree$trajectory$slice <- infer_trajectory(tree$dataset,"slice",verbose=TRUE)
  #other methods from dynverse that worked
  tree$trajectory$slicer <- infer_trajectory(tree$dataset,"slicer",verbose=TRUE)
  tree$trajectory$tscan <- infer_trajectory(tree$dataset,"tscan",verbose=TRUE)
  tree$trajectory$waterfall <- infer_trajectory(tree$dataset,"waterfall",verbose=TRUE)

  # #recommended methods from dynverse that didn't work
  # tree$trajectory$paga_tree <- infer_trajectory(tree$dataset,"paga_tree",verbose=TRUE) #TODO index error
  # tree$trajectory$paga <- infer_trajectory(tree$dataset,"paga",verbose=TRUE) #TODO index error
  # tree$trajectory$pcreode <- infer_trajectory(tree$dataset,"pcreode",verbose=TRUE) #TODO runs very slowly
  # tree$trajectory$monocle_ica <- infer_trajectory(tree$dataset,"monocle_ica",verbose=TRUE) #TODO error undefined column #monocle 1
  # tree$trajectory$raceid_stemid <- infer_trajectory(tree$dataset,"raceid_stemid",verbose=TRUE) #TODO error execution halted
  # #other methods from dynverse that didn't work
  # tree$trajectory$urd <- infer_trajectory(tree$dataset,"urd",verbose=TRUE)#TODO error in length of dimnames
  # tree$trajectory$scuba <- infer_trajectory(tree$dataset,"scuba",verbose=TRUE)#TODO error in size of an array
  # tree$trajectory$mpath <- infer_trajectory(tree$dataset,"mpath",verbose=TRUE)#TODO needs group cluster labels as a prior
  # tree$trajectory$wishbone <- infer_trajectory(tree$dataset,"wishbone",verbose=TRUE)#TODO python value error
  # tree$trajectory$wanderlust <- infer_trajectory(tree$dataset,"wanderlust",verbose=TRUE)#TODO python value error
  # tree$trajectory$scoup <- infer_trajectory(tree$dataset,"scoup",verbose=TRUE)#TODO needs group cluster labels as a prior
  # tree$trajectory$monocle_ddrtree <- infer_trajectory(tree$dataset,"monocle_ica",verbose=TRUE)#TODO error execution halter #monocle 2

  return(tree)
}

form_tree_from_file <- function(filename){
  tree <- vector(mode="list")
  tree$csv <- read.csv(filename,row.names=1)
  tree$raw <- tree$csv %>%
    select(Raw0:Raw499)
  tree$matrix <- as.matrix(tree$raw)
  tree$dimred <- tree$csv %>%
    select(tsne0:tsne1)
  #As input, dynwrap requires raw counts and normalised (log2) expression data.
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

  return(tree)
}

form_tree_from_dataset <- function(dataset){
  tree <- vector(mode="list")
  tree$dataset <- dataset
  tree$matrix <- as.matrix(dataset$expression)
    return(tree)
}

dynverse_tree <- generate_dataset(
  model = model_bifurcating(),
  num_cells = 700,
  num_features = 500
)
dynverse_tree <- form_tree_from_dataset(dynverse_tree)
write.table(dynverse_tree$matrix,file='TITesting.nosync/dynversetree/dynversetree.csv',row.names=FALSE,sep=',')
dynverse_tree <- run_methods(dynverse_tree,'TITesting.nosync/dynversetree/EBMFfactors/')
tree1 <- form_tree_from_file('TITesting.nosync/tree1/tree1.csv')
tree1 <- run_methods(tree1,'TITesting.nosync/tree1/EBMFfactors/')
tree2 <- form_tree_from_file('TITesting.nosync/tree2/tree2.csv')
tree2 <- run_methods(tree2,'TITesting.nosync/tree2/EBMFfactors/')
tree3 <- form_tree_from_file('TITesting.nosync/tree3/tree3.csv')
tree3 <- run_methods(tree3,'TITesting.nosync/tree3/EBMFfactors/')
tree4 <- form_tree_from_file('TITesting.nosync/tree4/tree4.csv')
tree4 <- run_methods(tree4,'TITesting.nosync/tree4/EBMFfactors/')

plot_dimred(tree1$trajectory$mst,
            dimred = tree1$dimred,
            grouping = group_onto_nearest_milestones(tree1$trajectory$mst))
plot_dimred(tree1$trajectory$slingshot,
            dimred = tree1$dimred,
            grouping = group_onto_nearest_milestones(tree1$trajectory$slingshot))
plot_dimred(tree1$trajectory$elpigraph,
            dimred = tree1$dimred,
            grouping = group_onto_nearest_milestones(tree1$trajectory$elpigraph))
plot_dimred(tree1$trajectory$sincell,
            dimred = tree1$dimred,
            grouping = group_onto_nearest_milestones(tree1$trajectory$sincell))
plot_dimred(tree1$trajectory$slice,
            dimred = tree1$dimred,
            grouping = group_onto_nearest_milestones(tree1$trajectory$slice))
plot_dimred(tree1$trajectory$slicer,
            dimred = tree1$dimred,
            grouping = group_onto_nearest_milestones(tree1$trajectory$slicer))
plot_dimred(tree1$trajectory$tscan,
            dimred = tree1$dimred,
            grouping = group_onto_nearest_milestones(tree1$trajectory$tscan))
plot_dimred(tree1$trajectory$waterfall,
            dimred = tree1$dimred,
            grouping = group_onto_nearest_milestones(tree1$trajectory$waterfall))
