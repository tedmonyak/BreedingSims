# Title: INDEPENDENT POPULATIONS WITH GENE FLOW
# Author: Ted Monyak
# Description: This script creates n.nPops independent sub-populations from an initial founder
# population and has them follow independent adaptive walks to a fitness optimum.
# It allows for gene flow under a stepping stone model, under which n.m individuals can move between 
# subpopulations each generation

# Assumes that CreateFounderPop.R has been run already
# Each simulation should create a new save_dir, where this data is stored

# founderPop is created in CreateFounderPop.R
pop <- founderPop
# Check if # independent pops * # individuals per pop does not exceed the
# number of individuals in the founder pop
if (n.nPops*n.subPopSize > nInd(pop)) {
  stop(paste0("Population of size ", nInd(pop), " not large enough to create ",
              n.nPops, " subpopulations of size ", n.subPopSize, "."))
}

# Burn-in
for (gen in 1:n.burnInGens) {
  pop <- selectCross(pop, trait=twoTraitFitFunc, nInd=nInd(pop)*n.burnInSelProp, nCrosses=nInd(pop))
}

# Create a random vector of size n.pops, with a random order of sub-population ids
randVec <- sample(rep(c(1:n.nPops), times=n.popSize/n.nPops))

# Create n.nPops sub populations
pops <- vector(mode="list", length=n.nPops)
for (p in 1:n.nPops) {
  pops[[p]] <- selectInd(pop, trait=selectSubPop, selectTop=TRUE, nInd=n.subPopSize, idx=p, randVec=randVec)
}

# Iterate through the number of generations
for (gen in 1:n.gens) {
  for (p in 1:length(pops)) {
    pop <- pops[[p]]
    # Select n.m*2 individuals to migrate
    migratingInds <- selectInd(pop, use="rand", nInd=n.m*2)
    migrateLeft <- migratingInds[1:n.m] # Select half of the individuals to move right
    migrateRight <- migratingInds[n.m+1:n.m*2] # Select half of the individuals to move left
    if (p > 1) {
      pops[[p-1]] <- c(pops[[p-1]], migrateLeft)
    }
    if (p < length(pops)) {
      pops[[p+1]] <- c(pops[[p+1]], migrateRight)
    }
    pop <- selectCross(pop, trait=threeTraitFitFunc, nInd=nInd(pop)*n.selProp, nCrosses=nInd(pop))
    pops[[p]] <- pop
  }
}

write.table(fit.df, file.path(subpop_dir, "fitness.csv"), col.names=TRUE, quote=FALSE, sep=",")

