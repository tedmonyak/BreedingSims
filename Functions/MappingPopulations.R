# Title: MAPPING POPULATIONS
# Author: Ted Monyak
# Description: Create Mapping populations
# Assumes that CreateFounderPop and CreateIndepdendentPops have been called

# Creates a biparental recombinant inbred line (RIL) population with n.RILFams families
# popA: the first population to sample an individual from
# popB: the second population to sample an individual from.
# save_dir: directory to write plots t
# inter: if TRUE, will cross an individual from popA and popB.
# if FALSE, will cross two individuals from popA
# Returns: a population where the first two individuals are parentA and parentB,
# and the rest is the RIL
createRIL <- function(popA, popB, save_dir, inter=TRUE) {
  # Develop elite lines
  popA <- makeElite(popA)
  popB <- makeElite(popB)
  # Select 2 random individuals from popA
  aIdx <- sample.int(nInd(popA),2)
  # parentA always comes from popA
  parentA <- popA[aIdx[1]]
  if (inter) {
    # if it is an inter-population cross, select an individual from popB
    parentB <- popB[sample.int(nInd(popB),1)]
  } else {
    # if it is an intra-population cross, take the other random individual from popA
    parentB <- popA[aIdx[2]]
  }

  # Cross parent A with parent B, and create n.RILFams progeny
  F1 <- randCross2(parentA,
                   parentB,
                   nCrosses=1,
                   nProgeny=n.RILFams)

  # Create F10s of each RIL family with SSD
  F2 <- self(F1, nProgeny=n.indPerRILFam)
  F3 <- self(F2)
  F4 <- self(F3)
  F5 <- self(F4)
  F6 <- self(F5)
  F7 <- self(F6)
  F8 <- self(F7)
  F9 <- self(F8)
  
  # Each RIL family has n.indPerRILFam replicates
  RIL <- self(F9)
  
  if (inter) {
    if (saveTraitPlots) {
      # Get the max number of individuals
      n <- max(nInd(popA), nInd(popB), nInd(RIL))
      
      # Normalize all of the phenotype vectors to have the same length (filling with NA)
      # to allow for cbind()
      phenoAT1 <- pheno(popA)[,1]
      phenoAT2 <- pheno(popA)[,2]
      length(phenoAT1) <- n
      length(phenoAT2) <- n
      
      phenoBT1 <- pheno(popB)[,1]
      phenoBT2 <- pheno(popB)[,2]
      length(phenoBT1) <- n
      length(phenoBT2) <- n
      
      phenoRilT1 <- pheno(RIL)[,1]
      phenoRilT2 <- pheno(RIL)[,2]
      length(phenoRilT1) <- n
      length(phenoRilT2) <- n

      theme <- theme(
        plot.title = element_text(family="Helvetica", size=22, hjust = 0.5),
        axis.title.x = element_text(family="Helvetica", size=20, vjust=-0.5),
        axis.title.y = element_text(family="Helvetica", size=20),
        axis.text.x = element_text(angle = 0, hjust=1, size=16),
        axis.text.y = element_text(angle = 0, hjust=1, size=16),
        legend.text = element_text(family="Helvetica", size=18),
        legend.title = element_text(family="Helvetica", size=18),
        legend.key = element_rect(linewidth=0.05),
        legend.spacing.y = unit(20, 'pt'),
        plot.caption = element_text(family="Helvetica", size=10, hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", color = "black"),
        aspect.ratio = 1)
      
      trait1.df <- as.data.frame(cbind(phenoAT1,phenoBT1,phenoRilT1))
      colnames(trait1.df) <- c("Subpopulation 1", "Subpopulation 2", "RIL Family")
      trait1.df <- trait1.df %>%
        pivot_longer(c("Subpopulation 1", "Subpopulation 2", "RIL Family"), names_to="Population", values_to="pheno") %>%
        drop_na()
      t1 <- ggplot(trait1.df, aes(pheno, fill=Population, color=Population)) +
        scale_color_manual(values=c("#674ED7", "#CC0000", "#3C78D8")) +
        geom_density(size=1, alpha=0.1) +
        xlim(-1,1) +
        ylim(0,5) +
        labs(title="Trait A", x="Phenotype", y="Density") +
        theme
      
      trait2.df <- as.data.frame(cbind(phenoAT2,phenoBT2,phenoRilT2))
      colnames(trait2.df) <- c("Subpopulation 1", "Subpopulation 2", "RIL Family")
      trait2.df <- trait2.df %>%
        pivot_longer(c("Subpopulation 1", "Subpopulation 2", "RIL Family"), names_to="Population", values_to="pheno") %>%
        drop_na()
      t2 <- ggplot(trait2.df, aes(pheno, fill=Population, color=Population)) +
        scale_color_manual(values=c("#674ED7", "#CC0000", "#3C78D8")) +
        geom_density(size=1, alpha=0.1) +
        xlim(-1,1) +
        ylim(0,5) +
        labs(title="Trait B", x="Phenotype", y="Density") +
        theme
      
      (t1|t2) + plot_layout(guides='collect', axes='collect')
      ggplot2::ggsave(filename = "trait_distributions.pdf",
                      path=save_dir,
                      device = "pdf",
                      width=15,
                      height=7)
      (plotTraitArchitecture(popA, "Additive", "popA") | plotTraitArchitecture(popB, "Additive", "popB"))
      ggplot2::ggsave(filename = "popA_popB_traitarchitecture.pdf",
                      path=save_dir,
                      device = "pdf",
                      width=20,
                      height=7)
      
      
    } # saveTraitPlots
    if (saveFitnessPlots) {
      fname <- file.path(save_dir, "3DFitness.html")
      fig <- plot3dPopulationFitnessTwoPops(popA, popB)
      htmlwidgets::saveWidget(as_widget(fig), fname)
    }
  } # inter

  if (saveTraitPlots) {
    plotTraitArchitecture(RIL, "Additive", "RIL")
    ggplot2::ggsave(filename = "RIL_traitarchitecture.pdf",
                    path=save_dir,
                    device = "pdf",
                    width=10,
                    height=7)
  }
  return (c(parentA, parentB, RIL))
}

# Simulate a population going through a breeding program, under purifying selection
# First, the top landrace individuals are purified
# pop: the landrace
# Returns: an F8 population that has undergone selection and inbreeding
makeElite <- function(pop) {
  purifiedLandraces <- selectInd(pop, trait=twoTraitFitFunc, nInd=n.landraces)
  # Purify each landrace
  for (f in 1:10) {
    purifiedLandraces <- self(purifiedLandraces)
  }
  # Create n.landrace x n.landrace F1s
  F1 <- randCross(purifiedLandraces, nCrosses=(n.landraces^2), simParam=SP)
  # Create enough F2s to select n.F2 to advance
  F2 <- self(F1, nProgeny=ceiling(n.F2/nInd(F1)))
  F2 <- selectInd(F2, trait=twoTraitFitFunc, nInd=n.F2)
  F3 <- self(F2)
  F3 <- selectInd(F3, trait=twoTraitFitFunc, nInd=n.F3)
  F4 <- self(F3)
  F4 <- selectInd(F4, trait=twoTraitFitFunc, nInd=n.F4)
  F5 <- self(F4)
  F5 <- selectInd(F5, trait=twoTraitFitFunc, nInd=n.F5)
  F6 <- self(F5)
  F6 <- selectInd(F6, trait=twoTraitFitFunc, nInd=n.F6)
  F7 <- self(F6)
  F7 <- selectInd(F7, trait=twoTraitFitFunc, nInd=n.F7)
  F8 <- self(F7)
  F8 <- selectInd(F8, trait=twoTraitFitFunc, nInd=n.F8)
  return (F8)
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

