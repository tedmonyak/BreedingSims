# Title: QTL Monte Carlo
# Author Ted Monyak
# Description:
# This will run n.popResets * n.sims simulations, and for each simulation,
# create two subpopulations, then two biparental RIL populations:
# 1 "inter" (parents come from different subpopulations),
# 1 "intra" (parents come from the same sub population)
# For each RIL, calculate and store the number of significant LOD peaks in the linkage map

library(AlphaSimR)
library(ade4)
library(devtools)
library(dplyr)
library(ggpubr)
library(ggplot2)
library(factoextra)
library(patchwork)
library(plotly)
library(purrr)
library(qtl)
library(qtlTools)
library(reshape2)
library(rrBLUP)
library(snow)
library(tibble)
library(tidyr)
library(viridis)
rm(list = ls())

set.seed(123)

# Set to false if running on Alpine
runLocal = TRUE

if (runLocal) {
  setwd("~/Documents/CSU/R/BreedingSims")
  output_dir <- file.path(getwd(), "Output")
  n.cores <- 16
} else {
  setwd("/pl/active/Morris_CSU/Ted_Monyak/BreedingSims")
  output_dir <- ("/scratch/alpine/c837220672@colostate.edu/Output")
  n.cores <- 4
}

if (!dir.exists(output_dir)) dir.create(output_dir)

source("Functions/Fitness.R")
source("Functions/GenoConversions.R")
source("Functions/MappingPopulations.R")
source("Functions/PopGen.R")
source("Functions/QtlMapping.R")
source("Functions/TraitArchitecture.R")
source("Scripts/GlobalParameters.R")

# Number of founder populations to simulate
n.popResets <- 4
# Number of adaptive walk simulations per pair of subpopulations
n.sims <- 25

saveQtlPlots <- FALSE
saveTraitPlots <- FALSE
saveAllelePlots <- FALSE
saveFitnessPlots <- FALSE
randParams <- FALSE

n.h2 <- 0.2
n.selProp <- 0.1
#n.gens <- 200
n.var <- 0.01

qtl_vec <- c(2,20)
pop_vec = c(50,500)
gen_vec <- c(50,100,150,200)

for (gx in 1:length(gen_vec)) {
  n.gens <- gen_vec[gx]
  print(paste0("Gens: ", n.gens))
  for (px in 1:length(pop_vec)){
    n.subPopSize <- pop_vec[px]
    print(paste0("Pop Size: ", n.subPopSize))
    for (qx in 1:length(qtl_vec)) {
      n.qtlPerChr <- qtl_vec[qx]
      print(paste0("QTL: ", n.qtlPerChr))
      eff_size.df <- data.frame(orderFixed=c(),
                                effectSize=c())
      
      # Result dataframe
      res.df <- data.frame(pop=c(),
                           sim=c(),
                           type=c(),
                           nSigQtl=c(),
                           fst=c(),
                           gens=c())
      
      # Tidy dataframe to store the order in which each allele is fixed, and the effect size
      eff_size.df <- data.frame(orderFixed=c(),
                                effectSize=c())
      
      base_dir <- file.path(output_dir, "QtlMonteCarlo")
      if (!dir.exists(base_dir)) dir.create(base_dir)
      base_fname <- paste0(paste0("Ne_", n.subPopSize, "_qtl_", n.qtlPerChr,
                                  "_selProp_", n.selProp, "_h2_", n.h2, "_gens_", n.gens, "_var_", n.var, "_"))
      base_dir <- file.path(base_dir, paste0(base_fname, format(Sys.time(), "%F_%H_%M")))
      if (!dir.exists(base_dir)) dir.create(base_dir)
      
      # Reset the founder population n.popResets times
      for (r in 1:n.popResets) {
        # Update selProp to be a random value if using random parameters
        if (randParams) {
          n.selProp <- runif(n=1, min=0.1, max=0.5)
        }
        pop_dir <- file.path(base_dir, paste0("FounderPopulation", r))
        if (!dir.exists(pop_dir)) dir.create(pop_dir)
        print(paste0("Pop Reset ", r))
        source("Scripts/CreateFounderPop.R")
        # For each founder population, create independent populations n.sims times
        for (s in 1:n.sims) {
          print(paste0("Sim ", s))
          save_dir <- file.path(pop_dir, paste0("Sim", s))
          dir.create(save_dir)
          # Landrace heritability
          SP$setVarE(h2=c(n.h2, n.h2))
          source("Scripts/CreateIndependentPops.R")
          fst <- FST(pops)
          # Breeding heritability
          SP$setVarE(h2=c(n.h2Breeding, n.h2Breeding))
          qtl_dir <- file.path(save_dir, "Inter")
          dir.create(qtl_dir)
          # Create a biparental RIL population by sampling one individual from each subpopulation
          res <- createRIL(popA=pops[[1]], popB=pops[[2]], save_dir=qtl_dir, inter=TRUE)
          # Select the parents from the result, and the RIL population
          RIL <- res[-(1:2)]
          parentA <- res[1]
          parentB <- res[2]
          # Determine the number of significant QTL
          nQtl <- getSigQtl(RIL, parentA, parentB, qtl_dir)
          res.df <- rbind(res.df, data.frame(pop=r,
                                             sim=s,
                                             type="Inter",
                                             nSigQtl=nQtl,
                                             fst=fst,
                                             gens=n.gens))
          
          qtl_dir <- file.path(save_dir, "Intra")
          dir.create(qtl_dir)
          # Create a biparental RIL population by sampling two individuals from the same subpopulation
          res <- createRIL(popA=pops[[1]], popB=pops[[2]], save_dir=qtl_dir, inter=FALSE)
          RIL <- res[-(1:2)]
          parentA <- res[1]
          parentB <- res[2]
          nQtl <- getSigQtl(RIL, parentA, parentB, qtl_dir)
          # Update the result dataframe
          res.df <- rbind(res.df, data.frame(pop=r,
                                             sim=s,
                                             type="Intra",
                                             nSigQtl=nQtl,
                                             fst=NA,
                                             gens=n.gens))
        } # end n.sims
      } # end n.popResets
      
      g <- ggplot(data=res.df, aes(x=type, y=nSigQtl, fill=type)) +
        geom_boxplot() +
        stat_compare_means(method="anova", label.x = 1.3, label.y= 13) +
        scale_fill_manual(name = "Parental Cross Type",
                          values = c("Inter" = "orange",
                                     "Intra" = "darkgreen"),
                          labels = c("Inter-population", "Intra-population")) +
        labs(title="# of Significant QTL from RIL Populations",
             x="# Significant QTL",
             y="Biparental Cross Type")
      
      fname <- file.path(base_dir, paste0(base_fname, "significant_qtl.pdf"))
      ggplot2::ggsave(filename = fname,
                      device = "pdf")
  
      gf <- ggplot(data=res.df, aes(fst)) +
        geom_density() +
        labs(title="Average FST", x="FST", y="Density") +
        theme
      fname <- file.path(base_dir, paste0(base_fname, "meanfst.pdf"))
      ggplot2::ggsave(filename = fname,
                      device = "pdf")
  
      if (saveAllelePlots) {
        # Determine the average additive effect size at each 'step'
        eff_size.df <- eff_size.df %>%
          group_by(orderFixed) %>%
          summarize(meanEffectSize = mean(effectSize))
        
        g2 <- ggplot(eff_size.df, aes(x=orderFixed, y=meanEffectSize)) +
          geom_bar(stat="identity")
        
        fname <- file.path(base_dir, paste0(base_fname, "average_effect_size_fixed.pdf"))
        ggplot2::ggsave(filename = fname,
                        device = "pdf")
      }
      
      write.table(res.df, file.path(base_dir, "result_dataframe.csv"), col.names=TRUE, quote=FALSE, sep=",")
      write.table(getParams(), file.path(base_dir, "params.txt"), col.names=FALSE, quote=FALSE, sep=":\t")
    } # end qtl_vec
  } # end pop_vec
} # end gen_vec
