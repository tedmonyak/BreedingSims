# Title: MAPPING POPULATIONS
# Author: Ted Monyak
# Description: Create Mapping populations
# Assumes that CreateFounderPop and CreateIndepdendentPops have been called

# Creates a biparental recombinant inbred line (RIL) population with n.RILFams families
# interPop: if TRUE, will select a random individual from popA and popB. if FALSE,
# will select two random individuals from popA
# Returns: a population where the first two individuals are parentA and parentB,
# and the rest is the RIL
# Assumes that popA and popB exist already
createRIL <- function(interPop=TRUE) {
  aIdx <- sample.int(nInd(popA),2)
  parentA <- popA[aIdx[1]]
  
  if (interPop) {
    parentB <- popB[sample.int(nInd(popA),1)]
  } else {
    parentB <- popA[aIdx[2]]
  }
  # First, make the parents inbreds by selfing them for 10 generations
  for (f in 1:10) {
    parentA <- self(parentA)
    parentB <- self(parentB)
  }
  
  # Cross parent A with parent B, and create n.RILFams progeny
  F1 <- randCross2(parentA,
                   parentB,
                   nCrosses=1,
                   nProgeny=n.RILFams)
  
  # Create F10s of each RIL family with SSD
  F2 <- self(F1)
  F3 <- self(F2)
  F4 <- self(F3)
  F5 <- self(F4)
  F6 <- self(F5)
  F7 <- self(F6)
  F8 <- self(F7)
  F9 <- self(F8)
  
  # Each RIL family has n.indPerRILFam replicates
  RIL <- self(F9, nProgeny=n.indPerRILFam)
  #BC1 <- randCross2(parentA, F1, nCrosses=200)
  return (c(parentA, parentB, RIL))
}

# Creates a nested association mapping (NAM) population, to be used for GWAS
# Assumes that pops - a list of populations, exists
createNAM <- function() {
  # The "reference" population is the first population
  refPop <- pops[[1]]
  # Randomly select an individual from the reference population
  refLine <- refPop[sample.int(nInd(refPop),1)]
  # Make the reference line inbred by selfing it for 10 generations
  for (f in 1:10) {
    refLine <- self(refLine)
  }
  # The rest of the founder lines will come from the other populations in pops
  founderLines <- c()
  # Iterate through the rest of the populations and select a random individual
  # to become a founder line
  for (p in 2:n.nPops) {
    subPop <- pops[[p]]
    ind <- subPop[sample.int(nInd(subPop),1)]
    for (f in 1:10) {
      ind <- self(ind)
    }
    founderLines <- append(founderLines, ind)
  }
  
  # Cross each of the founder lines by the reference line
  crossPlan = matrix(c(rep(1:10), rep(1,10)), nrow=10, ncol=2)
  
  F1 <- makeCross2(founderLines, refLine, crossPlan, simParam=SP)
  # Create F10s with SSD
  F2 <- self(F1)
  F3 <- self(F2)
  F4 <- self(F3)
  F5 <- self(F4)
  F6 <- self(F5)
  F7 <- self(F6)
  F8 <- self(F7)
  F9 <- self(F8)
  # Create n.indPerRILFam replicates per NAM family
  NAM <- self(F9, nProgeny=n.indPerRILFam)
  return (NAM)
}

# Creates a diversity panel population by merging all of the subpopulations
# Returns: a population
# Assumes that pops - a list of populations, exists
createDP <- function() {
  DP <- c()
  # Iterate through all of the populations and merge them together
  for (p in 1:n.nPops) {
    DP <- append(DP, pops[[p]])
  }
  return (DP)
}

