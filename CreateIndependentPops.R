pop <- founderPop
if (n.nPops*n.subPopSize > nInd(pop)) {
  stop(paste0("Population of size ", nInd(pop), " not large enough to create ",
              n.nPops, " subpopulations of size ", n.subPopSize, "."))
}

# BURN-IN
for (gen in 1:n.burnInGens) {
  pop <- selectCross(pop, trait=twoTraitFitFunc, nInd=nInd(pop)*n.burnInSelProp, nCrosses=nInd(pop))
}

# Create n.nPops sub populations
pops <- vector(mode="list", length=n.nPops)
for (p in 1:n.nPops) {
  pops[[p]] <- selectInd(pop, nInd=n.subPopSize, use="rand")
}

for (gen in 1:n.gens) {
  for (p in 1:length(pops)) {
    pop <- pops[[p]]
    if (mean(twoTraitFitFunc(pheno(pop))) < n.margin) {
      pops[[p]] <- selectCross(pop, trait=twoTraitFitFunc, nInd=nInd(pop)*n.selProp, nCrosses=nInd(pop))
    }
  }
}

