## Description: HIV transmission between MSM, transgender women (TW), and partners of TW (PTW) in Lima, Peru.
## Author: Diana M. Tordoff
## Date: 24 August 2021

library(data.table)
library(tibble)
library(dplyr)
library(Hmisc)
library(memisc)
library(ggplot2)
library(grid)
library(gridExtra)
library(ape)
library(seqinr)
library(phytools)
library(ggtree)
library(treeio)
library(treedater)
library(phydynR)
library(BayesianTools)

#load alignment, tree, and metadata
alignment <- read.fasta( )
tree <- read.newick( )
df <- read.csv( ) 


## STEP ONE: Defined Mathematical Model

demes <- c("msm", "tw", "ptw", "src")
m <- 4

# we used a fixed population size 
SRCSIZE = 3.2e7 #population of Peru
msm.pop = 5e5 #5% of the population
tw.pop = 5e4 #0.5% of the population
ptw.pop = tw.pop #probably an underestimate

#define parameters
params0 <- list( 
  lambda = 1.6/100, #incidence among msm (1.6 per 100 person-years)
  mu = 3/100,       #incidence among tw (made up)
  nu = 1/100,       #incidence among ptw (made up)
  p.msm2tw = 0.1,    #proportion of transmissions from msm -> tw
  p.msm2ptw = 0.1,    #proportion of transmissions from msm -> ptw
  p.tw2msm = 0.2,    #proportion of transmissions from tw -> msm
  p.tw2ptw = 0.5,    #proportion of transmissions from tw -> ptw
  p.ptw2tw = 0.8,    #proportion of transmissions from ptw -> tw
  p.ptw2msm = 0.1,    #proportion of transmissions from ptw -> tw
  delta = 1/30,     #mortality rate for infected
  psi = 1.2/100,    #importation rate
  Ne = 1.1/100      #effective source population size
)

#define initial conditions
X0 <- c(msm  = 0.206*msm.pop, #20.6% prevalence in Sabes
        tw   = 0.197*tw.pop,  #19.7% prevalence in Sabes
        ptw  = 0.081*ptw.pop, #8.1% prevalence in Sabes
        src  = 0.001*SRCSIZE) #0.1% prevalence in Peru

#create birth matrix:
births <- matrix( 0, nrow=m, ncol=m)
rownames(births) = colnames(births) <- demes
births['msm', 'msm'] <- '(1 - params0$p.msm2tw - params0$p.msm2ptw) * params0$lambda * msm'
births['msm', 'tw']  <- 'params0$p.msm2tw * params0$lambda * msm'
births['msm', 'ptw'] <- 'params0$p.msm2ptw * params0$lambda * msm'
births['tw', 'msm']  <- 'params0$p.tw2msm * params0$mu * tw'
births['tw', 'tw']   <- '(1 - params0$p.tw2msm - params0$p.tw2ptw) * params0$mu * tw'
births['tw', 'ptw']  <- 'params0$p.tw2ptw * params0$mu * tw'
births['ptw', 'msm'] <- 'params0$p.ptw2msm * params0$nu * ptw'
births['ptw', 'tw']  <- 'params0$p.ptw2tw * params0$nu * ptw'
births['ptw', 'ptw']  <- '(1 - params0$p.ptw2msm - params0$p.ptw2tw) * params0$nu * ptw'
births['src', 'src'] <- '0.5 * SRCSIZE^2 / params0$Ne'

#create migration matrix:
migs <- matrix( 0, nrow=m, ncol=m) 
rownames(migs) = colnames(migs) <- demes 
migs['msm', 'src'] <- 'params0$psi * msm'
migs['tw', 'src']  <- 'params0$psi * tw'
migs['ptw', 'src'] <- 'params0$psi * ptw'
migs['src', 'msm'] <- 'params0$psi * msm'
migs['src', 'tw']  <- 'params0$psi * tw'
migs['src', 'ptw'] <- 'params0$psi * ptw'

#create death matrix:
deaths <- setNames( rep(0, m), demes)
deaths['msm'] <- 'params0$delta * msm'
deaths['tw']  <- 'params0$delta * tw'
deaths['ptw'] <- 'params0$delta * ptw'
deaths['src'] <- 'params0$delta * src'

#build demographic model:
dm <- build.demographic.process(births, migrations=migs, death=deaths, parameter=names(params0), rcpp=FALSE)


## STEP TWO: Set up metadata for ``phydynR`` input
  
#create vector of sample times: 
sampleTimes <- sts
head(sampleTimes) 

#create vector of sample states:
msm <- tw <- ptw <- src <- rep(0, length(tree$tip.label))
msm[metadata$group == "MSM"] <- 1 
tw[metadata$group == "TW"] <- 1 
ptw[metadata$group == "PTW"] <- 1 
src[metadata$group == "LANL"] <- 1 
sampleStates <- cbind(msm, tw, ptw, src) 
rownames(sampleStates) <- c(metadata$label)
head(sampleStates)

#date tree
seqlen <- 700
dated <- dater(tree, sts, seqlen)
rootToTipRegressionPlot( dated )
dated.tree <- phydynR::DatedTree(phylo = dated, sampleTimes = sampleTimes, sampleStates = sampleStates, minEdgeLength = 0, tol = 0.1)


## STEP THREE: Estimate Parameters using MCMC

#define parameters to be estimates:
obj_fun <- function(parameters){
  parameters <- unname(parameters)
  THETA <- params0
  THETA$lambda <- parameters[1] 
  THETA$mu <- parameters[2] 
  THETA$nu <- parameters[3] 
  THETA$p.msm2tw <- parameters[4] 
  THETA$p.msm2ptw <- parameters[5] 
  THETA$p.tw2msm <- parameters[6] 
  THETA$p.tw2ptw <- parameters[7] 
  THETA$p.ptw2tw <- parameters[8] 
  THETA$p.ptw2msm <- parameters[9] 
  THETA$delta <- parameters[10] 
  THETA$psi <- parameters[11] 
  THETA$Ne <- parameters[12] 
  mll <- phydynR::colik(tree = dated.tree, 
                        theta = THETA,
                        demographic.process.model = dm,
                        x0 = X0,
                        t0 = 1975,
                        res = 1e3, 
                        timeOfOriginBoundaryCondition = FALSE, 
                        AgtY_penalty = 1,
                        maxHeight = 20)
  return(mll) 
}

#informative priors (density + sampler)
densities <- function(par){
  d1 = dgamma(par[1], shape = 2, rate = 100, log=TRUE) 
  d2 = dgamma(par[2], shape = 2, rate = 100, log=TRUE) 
  d3 = dgamma(par[3], shape = 1.5, rate = 15, log=TRUE)
  d4 = dbeta(par[4], shape1 = 2, shape2 = 25, log=TRUE) 
  d5 = dbeta(par[5], shape1 = 2, shape2 = 15, log=TRUE) 
  d6 = dbeta(par[6], shape1 = 5, shape2 = 45, log=TRUE) 
  d7 = dbeta(par[7], shape1 = 30, shape2 = 8, log=TRUE) 
  d8 = dbeta(par[8], shape1 = 22, shape2 = 10, log=TRUE) 
  d9 = dbeta(par[9], shape1 = 8, shape2 = 20, log=TRUE) 
  d10 = dunif(par[10], min = 1/31, max = 1/29, log=TRUE) 
  d11 = dexp(par[11], rate = 30, log = TRUE) 
  d12 = dexp(par[12], rate = 20, log = TRUE) 
  return(d1 + d2 + d3 + d4 + d5 + d6 + d7 + d8 + d9 + d10 + d11 + d12 )
}
sampler <- function(n=1){
  d1 = rgamma(n, shape = 2, rate = 100)
  d2 = rgamma(n, shape = 2, rate = 100) 
  d3 = rgamma(n, shape = 1.5, rate = 15)
  d4 = rbeta(n, shape1 = 2, shape2 = 25) 
  d5 = rbeta(n, shape1 = 2, shape2 = 15)
  d6 = rbeta(n, shape1 = 5, shape2 = 45) 
  d7 = rbeta(n, shape1 = 30, shape2 = 8)
  d8 = rbeta(n, shape1 = 22, shape2 = 10) 
  d9 = rbeta(n, shape1 = 8, shape2 = 20) 
  d10 = runif(n, min = 1/31, max = 1/29) 
  d11 = rexp(n, rate = 30) 
  d12 = rexp(n, rate = 20) 
  return(cbind(d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12))
}
prior <- createPrior(density = densities, sampler = sampler,
                     lower = c(0, 0, 0, 0, 0, 0, 0, 0, 0,  1/31, 0, 0.0001), 
                     upper = c(.1, .1, .2, 1, 1, 1, 1, 1, 1, 1/29, 0.15, 0.15))


#estimate parameters using MCMC
niter <- 100000
bayesianSetup = createBayesianSetup(likelihood = obj_fun, prior=prior)
settings = list(iterations = niter, message=FALSE) 
start_time <- Sys.time()
out <- runMCMC(bayesianSetup = bayesianSetup, sampler = "DEzs", settings = settings)
end_time <- Sys.time()
end_time-start_time

#output 
summary(out)
tracePlot(out)
correlationPlot(out)
marginalPlot(out, prior = T)
gelmanDiagnostics(out, plot = T)
