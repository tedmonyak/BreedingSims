# Title: CREATE INDEPENDENT POPULATIONS
# Author: Ted Monyak
# Description: This script creates n.nPops independent sub-populations from an initial founder
# population and has them follow independent adaptive walks to a fitness optimum

# Assumes that CreateFounderPop.R has been run already

# founderPop is created in CreateFounderPop.R
pop <- founderPop
n.nPops <- 10
n.subPopSize <- 100
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
#fig <- plot_ly()
# Each population follows an adaptive walk for a maximum of n.gens generations
# Each will terminate once it is within n.margin of the fitness optimum
for (p in 1:length(pops)) {
  fit_df <- data.frame(gen=c(),
                       fitness=c(),
                       traitValA=c(),
                       traitValB=c())
  pop <- pops[[p]]
  for (gen in 1:n.gens) {
    if (mean(twoTraitFitFunc(pheno(pop))) < n.margin) {
      # At each stage, select the top individuals according to how close each 
      # is from the fitness optimum
      meanFitness <- mean(twoTraitFitFunc(pheno(pop)))
      selRat <- selectionRatio(meanFitness)
      fit_df <- rbind(fit_df, data.frame(gen=gen,
                                         fitness=meanFitness,
                                         traitValA=meanP(pop)[1],
                                         traitValB=meanP(pop)[2]))
      pop <- selectCross(pop, trait=twoTraitFitFunc, nInd=nInd(pop)*selRat, nCrosses=nInd(pop))
    }
  }
  
  # Plot the adaptive walks
  #fig <- add_trace(
  #  fig,
  #  fit_df,
  #  name = p,
  #  x = fit_df$traitValA,
  #  y = fit_df$traitValB,
  #  z = fit_df$fitness,
  #  type = 'scatter3d',
  #  mode = 'lines',
  #  opacity = 1,
  #  color = p,
  #  line = list(width = 5)
  #)
}

#fig
