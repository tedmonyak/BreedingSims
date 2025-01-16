# Title: CREATE INDEPENDENT POPULATIONS
# Author: Ted Monyak
# Description: This script creates n.nPops independent sub-populations from an initial founder
# population and has them follow independent adaptive walks to a fitness optimum

# Assumes that CreateFounderPop.R has been run already

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

# Create n.nPops sub populations
pops <- vector(mode="list", length=n.nPops)
for (p in 1:n.nPops) {
  pops[[p]] <- selectInd(pop, nInd=n.subPopSize, use="rand")
}

# Each population follows an adaptive walk for a maximum of n.gens generations
# Each will terminate once it is within n.margin of the fitness optimum
for (gen in 1:n.gens) {
  for (p in 1:length(pops)) {
    pop <- pops[[p]]
    if (mean(twoTraitFitFunc(pheno(pop))) < n.margin) {
      # At each stage, select the top individuals according to how close each 
      # is from the fitness optimum
      pops[[p]] <- selectCross(pop, trait=twoTraitFitFunc, nInd=nInd(pop)*n.selProp, nCrosses=nInd(pop))
    }
  }
}

