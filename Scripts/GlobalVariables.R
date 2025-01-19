# Title: GLOBAL VARIABLES
# Author: Ted Monyak
# Description: This file contains default parameters for use in all breeding simulations

# POPULATION
n.popSize = 1000 # Num. individuals in the founder population
n.subPopSize = 100 # Num individuals in the independent sub-populations
n.ne = n.popSize # Effective population size
n.segSites = 1000 # Initial # segregating alleles per chromosome in the population
n.markers = 1000 # Num. loci per chromosome to include on the SNP chip 
n.chr <- 10 # Num. chromosomes (sorghum has 10)
addSnpChip <- TRUE # If true, will add a SNP chip with n.markers markers
basicPop <- TRUE # If true, uses runMacs. If false, uses runMacs2 with more custom parameters

# TRAITS
n.qtlPerChr = 5 # Number of qtl per chromosome, per trait
n.h2 <- 0.6 # Narrow-sense heritability for each trait
n.initTraitVal <- 10
n.var <- 1 # Initial variance for each trait. This is used as the 'rate' parameter in a gamma distribution
n.shape <- 1 # Initial shape for the gamma distribution for each trait
n.allele <- 2 # This is the allele for which to track frequency over time 

# ADAPTIVE WALKS
n.margin <- 0 # Populations will terminate their adaptive walks once they reach this fitness value
n.burnInSelProp <- 0.95 # % of the population to advance during burn-in
n.gens <- 100 # maximum number of generations for an adaptive walk
n.burnInGens <- 5 # number of burn-in generations for founder population
n.nPops <- 2 # number of independent subpopulations to create
n.selProp <- 0.3 # % of the population to advance during main adaptive walk
n.r <- 1 # 0.9 # the r value to use in the geometric series for a decaying selection intensity. Set to 1 for a non-decaying intensity, and decrease this to increase the rate of decay

# SIMULATIONS
n.popResets <- 10 # number of times to reset the founder population in a simulation
n.sims <- 100 # number of monte carlo simulations to run per set of parameters

# Mapping populations
n.RILFams <- 200 # number of RIL families to create
n.indPerRILFam <- 4 # number of replicates in each RIL family

# QTL Mapping Parameters
n.mappingMethod <- "hk" # One of 'hk' (Haley-Knott), 'em', 'imp' (imputation), or others. See ?scanone
n.errorProb <- 0.001 # assumed genotyping error rate. See ?calc.genoprob
n.step <- 1 # max distance b/t calculated genotype probabilities. See ?calc.genoprob
n.cores <- 8 # Number of cores to use for QTL mapping and population generation
saveQtlPlots <- FALSE # If true, will write each qtl plot to disk. Should be off for monte carlo simulations

# Function for returning all parameters (including ones that have been updated)
# for writing to the params.txt file along with graphs
getParams <- function() {
  n.df <- data.frame(
    popSize=n.popSize,
    subPopSize=n.subPopSize,
    ne=n.ne,
    segSites=n.segSites,
    markers=n.markers,
    chr=n.chr,
    basicPop=basicPop,
    qtlPerChr=n.qtlPerChr,
    h2=n.h2,
    initTraitVal=n.initTraitVal,
    var=n.var,
    shape=n.shape,
    allel=n.allele,
    margin=n.margin,
    burnInSelProp=n.burnInSelProp,
    selProp=n.selProp,
    gens=n.gens,
    burnInGens=n.burnInGens,
    sims=n.sims,
    popResets=n.popResets,
    RILFams=n.RILFams,
    indPerRILFam=n.indPerRILFam,
    nPops=n.nPops,
    mappingMethod=n.mappingMethod,
    errorProb=n.errorProb,
    step=n.step,
    cores=n.cores
  )
  return (t(n.df))
}
