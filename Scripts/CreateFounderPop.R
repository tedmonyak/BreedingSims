# Title: CREATE FOUNDER POPULATION
# Author: Ted Monyak
# Description: This script creates a founder population to use in subsequent simulations

# For reproducibility
set.seed(123)

# Create initial population
if (basicPop) {
  # default founder parameters
  founders = runMacs(
    nInd=n.popSize,
    nChr=n.chr,
    segSites=n.segSites,
    nThreads=n.cores
  )
} else {
  # This is from a Brian Rice script, and is not currently being used
  founders = runMacs(nInd=n.ne,
                     nChr=10,
                     segSites=n.segSites,
                     inbred=TRUE,
                     manualCommand = paste(
                       "1000000000 -t", #Physical length 1e8 base pairs
                       2.5/1E8*(4*n.ne), #Mutation rate adjusted for Ne
                       "-r",1/1E8*(4*n.ne), #Recombination rate adjusted for Ne
                       "-eN",10/(4*n.ne),100/n.ne), #Modeling Ne=100 at 10 generations ago
                     manualGenLen=1, nThreads=n.cores) #Genetic length 1 Morgan
}


SP <- SimParam$new(founders)
# Allow markers and segmentation sites to overlap
SP$restrSegSites(overlap = T)

# Create two generic traits with normal distributions of effect sizes (or gamma, if !normalDist)
if (normalDist) {
  if (randParams) {
    SP$addTraitA(mean=n.initTraitVal, var=runif(n=1,min=0,max=0.1), nQtlPerChr=n.qtlPerChr)
    SP$addTraitA(mean=n.initTraitVal, var=runif(n=1,min=0,max=0.1), nQtlPerChr=n.qtlPerChr)
  } else {
    SP$addTraitA(mean=n.initTraitVal, var=n.var, nQtlPerChr=n.qtlPerChr)
    SP$addTraitA(mean=n.initTraitVal, var=n.var, nQtlPerChr=n.qtlPerChr)
  }
} else {
  SP$addTraitA(mean=n.initTraitVal, var=n.var, nQtlPerChr=n.qtlPerChr, gamma=TRUE, shape=n.shape)
  SP$addTraitA(mean=n.initTraitVal, var=n.var, nQtlPerChr=n.qtlPerChr, gamma=TRUE, shape=n.shape)
}

# If we want to use more realistic trait values, use the following for height, flowering time, yield
#SP$addTraitA(mean=55, var=20, nQtlPerChr=c(1,1,1,1,0,0,0,0,0,0), gamma = TRUE, shape = n.shape) # height
#SP$addTraitA(mean=150, var=50, nQtlPerChr=c(1,1,1,1,0,0,0,0,0,0), gamma = TRUE, shape = n.shape) # flowering time
#SP$addTraitA(mean=950, var=200, nQtlPerChr = 100, gamma = TRUE, shape = 1) # yield

# If using random parameters, set heritability to be a random value
if (randParams) {
  n.h2 <- runif(n=1,min=0,max=0.5)
}
# Set heritability for each of the traits
SP$setVarE(h2=c(n.h2, n.h2))


# Add a SNP chip with n.markers*n.chr markers
if (addSnpChip) SP$addSnpChip(nSnpPerChr=n.markers)

# Create base population
founderPop <- newPop(founders, simParam = SP)
