# Title: MAPPING POPULATIONS
# Author: Ted Monyak
# Description: Create Mapping populations
# Assumes that CreateFounderPop and CreateIndepdendentPops have been called

# Creates a biparental recombinant inbred line (RIL) population with n.RILFams families
# popA: the first population to sample an individual from
# popB: the second population to sample an individual from. If left as NA, select two
# individuals from popA
# Returns: a population where the first two individuals are parentA and parentB,
# and the rest is the RIL
createRIL <- function(popA, popB, save_dir) {
  # Select 2 random individuals from popA
  aIdx <- sample.int(nInd(popA),2)
  # parentA always comes from popA
  parentA <- popA[aIdx[1]]
  
  # 
  if (!is.na(popB)) {
    # if it is an inter-population cross, select an individual from popB
    parentB <- popB[sample.int(nInd(popB),1)]
  } else {
    # if it is an intra-population cross, take the other random individual from popA
    parentB <- popA[aIdx[2]]
  }
  # Make the parents inbreds by selfing them for 10 generations
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
  
  if (!is.na(popB)) {
    trait1.df <- as.data.frame(cbind(pheno(popA)[,1], pheno(popB)[,1], pheno(RIL)[,1]))
    colnames(trait1.df) <- c("popA", "popB", "RIL")
    trait1.df <- trait1.df %>%
      pivot_longer(c("popA", "popB", "RIL"), names_to="pop", values_to="pheno")
    t1 <- ggplot(trait1.df, aes(pheno, fill=pop, color=pop)) +
      geom_density(alpha=0.1) +
      labs(title="Trait 1")
    
    trait2.df <- as.data.frame(cbind(pheno(popA)[,2], pheno(popB)[,2], pheno(RIL)[,2]))
    colnames(trait2.df) <- c("popA", "popB", "RIL")
    trait2.df <- trait2.df %>%
      pivot_longer(c("popA", "popB", "RIL"), names_to="pop", values_to="pheno")
    t2 <- ggplot(trait2.df, aes(pheno, fill=pop, color=pop)) +
      geom_density(alpha=0.1) +
      labs(title="Trait 2")
    
    
    (t1|t2)
    ggplot2::ggsave(filename = "trait_distributions.pdf",
                    path=save_dir,
                    device = "pdf",
                    width=20,
                    height=7)
    (plotTraitArchitecture(popA, "Additive", "popA") | plotTraitArchitecture(popB, "Additive", "popB"))
    ggplot2::ggsave(filename = "popA_popB_traitarchitecture.pdf",
                    path=save_dir,
                    device = "pdf",
                    width=20,
                    height=7)
  }
  

  plotTraitArchitecture(RIL, "Additive", "RIL")
  ggplot2::ggsave(filename = "RIL_traitarchitecture.pdf",
                  path=save_dir,
                  device = "pdf",
                  width=10,
                  height=7)
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

