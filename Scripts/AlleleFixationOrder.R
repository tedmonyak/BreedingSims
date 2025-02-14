# Title: ALLELE FIXATION ORDER
# Author: Ted Monyak
# Description: This block runs n.nPops * n.sims Monte Carlo simulations of different
# populations to determine the change in the average effect size of alleles that are
# fixed along the adaptive walk, saving the order in which the alleles are fixed
# n.nPops populations are created, and each one undergoes n.sims unique adaptive walks

setwd("~/Documents/CSU/R/BreedingSims")

source("Scripts/GlobalParameters.R")
saveAllelePlots <- TRUE

if (saveFitnessPlots) {
  fig <- plot_ly()
}

# Tidy dataframe to store the order in which each allele is fixed, and the effect size
effSize.df <- data.frame(orderFixed=c(),
                         gen=c(),
                         fitness=c(),
                         effectSizeA=c(),
                         effectSizeF=c())
# Create a new population n.nPops times
for (p in 1:n.popResets) {
  print(paste0("Pop Reset: ", p))
  source("Scripts/CreateFounderPop.R")
  
  # Get the effect sizes of each qtl
  qtlEff.df <- getQtlEffectSizes(founderPop)
  
  # Get the names of all the QTLs
  qtl <- colnames(getUniqueQtl(founderPop))
  for (s in 1:n.sims) {
    pop <- founderPop
    if (saveFitnessPlots) {
      pop_df <- data.frame(gen=1:n.gens,
                           fitness=numeric(n.gens),
                           traitValA=numeric(n.gens),
                           traitValB=numeric(n.gens))
    }
    
    
    print(paste0("Sim: ", s))
    # idx is the order in which an allele is fixed along an adaptive walk
    idx <- 1
    # whether or not to increment the idx counter. Multiple alleles may be fixed
    # in the same generation, so this cannot be incremented until each locus
    # has been examined
    inc <- FALSE
    
    # Burn-in
    for (gen in 1:n.burnInGens) {
      pop <- selectCross(pop, trait=twoTraitFitFunc, nInd=nInd(pop)*n.burnInSelProp, nCrosses=nInd(pop))
    }
    
    # Select a subpopulation
    pop <- selectInd(pop, nInd=n.subPopSize, use="rand")
    
    # Main simulation
    for (gen in 1:n.gens) {
      meanFitness <- mean(twoTraitFitFunc(pheno(pop)))
      selRat <- selectionRatio(meanFitness)
      if (saveFitnessPlots) {
        pop_df$fitness[gen] <- meanFitness
        pop_df$traitValA[gen] <- meanP(pop)[1]
        pop_df$traitValB[gen] <- meanP(pop)[2]
      }
      # Keep track of the current population
      prevPop <- pop
      prevGeno <- getUniqueQtl(prevPop)
      # Advance the population based on fitness
      pop <- selectCross(pop, trait=twoTraitFitFunc, nInd=nInd(pop)*selRat, nCrosses=nInd(pop))
      # Get the qtl genotypes from the current and new populations, so we can compare them
      # Determine whether each locus is segregating
      prevHet <- as.data.frame(apply(getUniqueQtl(prevPop), MARGIN=2, FUN=hetLocus))
      curHet <- as.data.frame(apply(getUniqueQtl(pop), MARGIN=2, FUN=hetLocus))
      # Join the two sets of genetic data
      loci <- cbind(prevHet, curHet)
      colnames(loci) <- c("prev", "cur")
      # Add a column called "fixed" which is true if the previous genotype was segregating
      # and the current genotype is not
      loci <- loci %>%
        mutate(fixed=mapply(function(p,c) (p && !c), prev, cur))
      
      # Iterate through all loci, and if the locus was fixed this generation,
      # add it to the result dataframe
      for (l in 1:length(qtl)) {
        id <- qtl[l]
        if (loci[id, "fixed"] == TRUE) {
          # Increment the order counter after this generation
          inc <- TRUE
          effectSizeF <- getEffectSize()
          
          effSize.df <- rbind(effSize.df,
                              data.frame(orderFixed=c(idx),
                                         gen=c(gen),
                                         fitness=c(meanFitness),
                                         effectSizeA=c(qtlEff.df[id,1])),
                                         effectSizeF=c(getEffectSize(locus=prevGeno[l],
                                                                     id=id,
                                                                     pop=prevPop,
                                                                     methodType="Fitness")))
        }
      }
      # Check whether to increment the order counter and reset 'inc'
      if (inc) {
        idx <- idx + 1
        inc <- FALSE
      }
      # If all alleles are fixed, terminate the simulation
      if(!any(curHet)) {
        break
      }
      # If the population is within n.margin, terminate the simulation
      if (mean(twoTraitFitFunc(pheno(pop))) >= n.margin) {
        break
      }
      
    }
    # Add the adaptive walk of the sub-population
    if (saveFitnessPlots) {
      fig <- add_trace(
        fig,
        pop_df,
        name = s,
        x = pop_df$traitValA,
        y = pop_df$traitValB,
        z = pop_df$fitness,
        type = 'scatter3d',
        mode = 'lines',
        opacity = 1,
        color = s,
        line = list(width = 2)
      )
    }
  }
}

save_dir <- file.path(output_dir, "AverageEffectSize")
if (!dir.exists(save_dir)) dir.create(save_dir)
save_dir <- file.path(save_dir, format(Sys.time(), "%F_%H_%M_%S"))
if (!dir.exists(save_dir)) dir.create(save_dir)
write.table(getParams(), file.path(save_dir, "params.txt"), col.names=FALSE, quote=FALSE, sep=":\t")
write.table(effSize.df, file.path(save_dir, "effSize.csv"), col.names=TRUE, quote=FALSE, sep=",")

order.df <- effSize.df %>%
  group_by(orderFixed) %>%
  #filter(!(abs(effectSize - median(effectSize)) > 2*sd(effectSize))) %>%
  summarize(meanEffectSize = mean(effectSizeA), n=n()) %>%
  filter(n > 100)

gO <- ggplot(order.df, aes(x=orderFixed, y=meanEffectSize)) +
  geom_point() +
  theme +
  labs(title="", x="Order Fixed", y="Mean Effect Size") +
  geom_smooth(formula = y ~ x, method="loess", color="black")

ggplot2::ggsave(filename = paste0("average_effect_size_order.jpg"),
                path=save_dir,
                device = "jpg",
                width=3, height=3, units="in")

# Create a plot with the adaptive walks
if (saveFitnessPlots) {
  p <- fig %>%
    layout(legend=list(title=list(text='Population')),
           showlegend=FALSE,
           scene = list(xaxis = list(title = "Trait A"),
                        yaxis = list(title = "Trait B"),
                        zaxis = list(title = "Fitness"),
                        aspectmode='cube')) %>% hide_colorbar()
  fname <- file.path(save_dir, "adaptivewalk.html")
  htmlwidgets::saveWidget(as_widget(p), fname)
}