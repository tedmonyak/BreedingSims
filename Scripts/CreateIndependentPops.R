# Title: CREATE INDEPENDENT POPULATIONS
# Author: Ted Monyak
# Description: This script creates n.nPops independent sub-populations from an initial founder
# population and has them follow independent adaptive walks to a fitness optimum

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

# Each population follows an adaptive walk for a maximum of n.gens generations
# Each will terminate once it is within n.margin of the fitness optimum
for (p in 1:length(pops)) {
  subpop_dir <- file.path(save_dir, paste0("Subpop_", p))
  if (!dir.exists(subpop_dir)) dir.create(subpop_dir)
  fig <- plot_ly()
  pop <- pops[[p]]

  # Get the effect sizes of each qtl
  qtlEff.df <- getQtlEffectSizes(pop)
  
  # Get the names of all the QTLs
  qtl <- rownames(qtlEff.df)

  # Create a dataframe of all zeros where the columns are the QTL ids, and the # rows is the # of generations
  fit.df <- data.frame(matrix(ncol=length(qtl)+4, nrow=0))
  colnames(fit.df) <- c("gen", "fitness", "traitValA", "traitValB", qtl)
  
  for (gen in 1:n.gens) {
    # Get the qtl genotype data
    qtlGeno <- getUniqueQtl(pop)
    alleleFreq <- data.frame(matrix(0, nrow=1, ncol=length(qtl)))
    colnames(alleleFreq) <- qtl
    # Get the frequency of the '2' allele at each locus
    for (l in 1:length(qtl)) {
      # id is the name of the qtl (chr_site)
      id <- qtl[l]
      # A list of genotype data for each individual in the population at that locus
      locus <- qtlGeno[,l]
      # Calculate the allele frequency as the frequency of homozygous individuals (for 'allele') +
      # 1/2 * frequency of heterozygous individuals (assumes the locus is biallelic)
      alleleFreq[1,id] <- (sum(locus==n.allele)/n.popSize) + ((sum(locus==1)/n.popSize)/2)
    }
    if (mean(twoTraitFitFunc(pheno(pop))) < n.margin) {
      # At each stage, select the top individuals according to how close each 
      # is from the fitness optimum
      meanFitness <- mean(twoTraitFitFunc(pheno(pop)))
      # Calculate selection ratio
      selRat <- selectionRatio(meanFitness)
      newRow <- data.frame(gen=gen,
                           fitness=meanFitness,
                           traitValA=meanP(pop)[1],
                           traitValB=meanP(pop)[2])
      # Join the new row with the newly calculated allele frequencies
      newRow <- cbind(newRow, alleleFreq)
      fit.df <- rbind(fit.df, newRow)
      pop <- selectCross(pop, trait=twoTraitFitFunc, nInd=nInd(pop)*selRat, nCrosses=nInd(pop))
    }
  }
  # Update the population in the list
  pops[[p]] <- pop
  
  # All the following graphing is for each subpopulation
  # Plot the adaptive walks
  fig <- add_trace(
    fig,
    fit.df,
    name = p,
    x = fit.df$traitValA,
    y = fit.df$traitValB,
    z = fit.df$fitness,
    type = 'scatter3d',
    mode = 'lines',
    opacity = 1,
    color = p,
    line = list(width = 5)
  )
  fig <- fig %>% layout(legend=list(title=list(text='Population Size')),
         scene = list(xaxis = list(title = "Trait A"),
                      yaxis = list(title = "Trait B"),
                      zaxis = list(title = "Fitness"),
                      aspectmode='cube')) %>% hide_colorbar()
  
  fname <- file.path(subpop_dir, "adaptivewalk.html")
  htmlwidgets::saveWidget(as_widget(fig), fname)
  
  plotTraitArchitecture(pop, "Additive", paste0("Pop ", p))
  
  ggplot2::ggsave(filename = "traitarchitecture.pdf",
                  path=subpop_dir,
                  device = "pdf",
                  width=10,
                  height=7)

  # Make the dataframe tidy
  freq.df <- fit.df[-c(2:4)]
  freq.df<- melt(freq.df, id="gen", variable.name="id", value.name="freq")
  
  # Add the qtl effect size data to the dataframe
  freq.df <- merge(freq.df, qtlEff.df, by.x="id", by.y="row.names", all.x=TRUE)
  
  # Create a line plot for the change in frequency of alleles over time
  # Each line's opacity is a function of its effect size (higher=darker)
  ggplot(freq.df, aes(x=gen, y=freq, color=id, alpha=eff_size)) +
    geom_line(size=0.7, show.legend=FALSE)
    
  
  ggplot2::ggsave(filename = "allelefrequencies.pdf",
                  path=subpop_dir,
                  device = "pdf",
                  width=10,
                  height=7)
  write.table(fit.df, file.path(subpop_dir, "fitness.csv"), col.names=TRUE, quote=FALSE, sep=",")
  
}
