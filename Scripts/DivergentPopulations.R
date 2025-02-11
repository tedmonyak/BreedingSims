# Title: PLOT DIVERGENT POPULATIONS
# Author: Ted Monyak
# Description: This script will create a base population, from which two subpopulations
# are randomly selected. Those subpopulations then follow an adaptive walk.
# The walks are then plotted on a 2d contour graph.

source("Scripts/GlobalParameters.R")
source("Scripts/CreateFounderPop.R")

fit_df <- data.frame(gen=1:n.burnInGens,
                     fitness=numeric(n.burnInGens),
                     traitValA=numeric(n.burnInGens),
                     traitValB=numeric(n.burnInGens))
pop <- founderPop

# Burn-in generations
for (gen in 1:n.burnInGens) {
  fit_df$fitness[gen] <- mean(twoTraitFitFunc(pheno(pop)))
  fit_df$traitValA[gen] <- meanP(pop)[1]
  fit_df$traitValB[gen] <- meanP(pop)[2]
  pop <- selectCross(pop, trait=twoTraitFitFunc, nInd=nInd(pop)*n.burnInSelProp, nCrosses=nInd(pop))
}

# Create a random vector of size n.pops, with a random order of sub-population ids
randVec <- sample(rep(c(1:n.nPops), times=n.popSize/n.nPops))

# Select all of the "1" indexed individuals
popA <- selectInd(pop, trait=selectSubPop, selectTop=TRUE, nInd=n.subPopSize, idx=1, randVec=randVec)
# Select all of the "2" indexed individuals
popB <- selectInd(pop, trait=selectSubPop, selectTop=TRUE, nInd=n.subPopSize, idx=2, randVec=randVec)

# Create dataframes for each subpopulation, initializing with current values
popA_df <- data.frame(gen=c(1),
                      fitness=c(mean(twoTraitFitFunc(pheno(popA)))),
                      traitValA=c(meanP(popA)[1]),
                      traitValB=c(meanP(popA)[2]))
popB_df <- data.frame(gen=c(1),
                      fitness=c(mean(twoTraitFitFunc(pheno(popB)))),
                      traitValA=c(meanP(popB)[1]),
                      traitValB=c(meanP(popB)[2]))

# Iterate through the generations
for (gen in 1:n.gens) {
  # If popA is within the margin of the fitness optimum, don't progress it any further
  if (mean(twoTraitFitFunc(pheno(popA))) < n.margin) {
    # Advance the population
    meanFitness <- calculateFitnessTwoTrait(meanP(popA)[1], meanP(popA[2]))
    # Get a selection ratio based on fitness
    selRat <- n.selProp
    popA <- selectCross(popA, trait=twoTraitFitFunc, nInd=nInd(popA)*selRat, nCrosses=nInd(popA))
    # Update the dataframe with new values
    popA_df <- rbind(popA_df, data.frame(gen=gen,
                                         fitness=meanFitness,
                                         traitValA=meanP(popA)[1],
                                         traitValB=meanP(popA)[2]))
    
  }
  # If popB is within the margin of the fitness optimum, don't progress it any further
  if (mean(twoTraitFitFunc(pheno(popB))) < n.margin) {
    # If popA is within the margin of the fitness optimum, don't progress it any further
    meanFitness <- calculateFitnessTwoTrait(meanP(popB)[1], meanP(popB[2]))
    # Get a selection ratio based on fitnes
    selRat <- n.selProp
    popB <- selectCross(popB, trait=twoTraitFitFunc, nInd=nInd(popB)*selRat, nCrosses=nInd(popB))
    # Update the dataframe with new values
    popB_df <- rbind(popB_df, data.frame(gen=gen,
                                         fitness=meanFitness,
                                         traitValA=meanP(popB)[1],
                                         traitValB=meanP(popB)[2]))
  }
}
# Update rownames
rownames(popA_df) <- 1:nrow(popA_df)
rownames(popB_df) <- 1:nrow(popB_df)

# Plot the adaptive walks
# Create a matrix of fitness values, with a small increment along the x and y axes.
fitness_x = seq(n.initTraitVal*-1,n.initTraitVal*1, by=(n.initTraitVal/100))
fitness_y = seq(n.initTraitVal*-1,n.initTraitVal*1, by=(n.initTraitVal/100))
fitness_z = outer(fitness_x,fitness_y,calculateFitnessTwoTrait)

f <- list(family="Arial", size=30)

fig <- plot_ly() %>%
  layout(legend=list(text= "Subpopulation"),
         xaxis = list(title = "Trait A", theme=theme, constrain = "domain"),
         yaxis = list(title = "Trait B", theme=theme, scaleanchor="x"),
         font=f) %>%
  add_trace(
    fig,
    x=fitness_x,
    y=fitness_y,
    z=fitness_z,
    type='contour',
    colors = viridis(n=10),
    colorbar=list(title = "Fitness"),
    line = list(color = 'black', width = 1),
    opacity=1) %>%
  add_trace(
    fig,
    fit_df,
    name = "Founder Population",
    x = fit_df$traitValA,
    y = fit_df$traitValB,
    type='scatter',
    mode = 'lines',
    line = list(color = '#E69138', width = 6, dash = 'solid'),
    opacity = 1) %>%
  add_trace(
    fig,
    popA_df,
    name = "Subpopulation 1",
    x = popA_df$traitValA,
    y = popA_df$traitValB,
    type='scatter',
    mode = 'lines',
    line = list(color = '#CC0000', width = 6, dash = 'solid'),
    opacity = 1) %>%
  add_trace(
    fig,
    popB_df,
    name = "Subpopulation 2",
    x = popB_df$traitValA,
    y = popB_df$traitValB,
    type='scatter',
    mode = 'lines',
    line = list(color = '#3C78D8', width = 6, dash = 'solid'),
    opacity = 1)

save_dir <- file.path(output_dir, "DivergingPopulations")
if (!dir.exists(save_dir)) dir.create(save_dir)
save_dir <- file.path(save_dir, format(Sys.time(), "%F_%H_%M_%S"))
if (!dir.exists(save_dir)) dir.create(save_dir)
fname <- file.path(save_dir, "adaptivewalks.html")
htmlwidgets::saveWidget(as_widget(fig), fname)

fname <- file.path(save_dir, "2PopulationFitness.html")
p <- plot3dPopulationFitnessTwoPops(popA, popB)
htmlwidgets::saveWidget(as_widget(p), fname)
write.table(getParams(), file.path(save_dir, "params.txt"), col.names=FALSE, quote=FALSE, sep=":\t")

fname <- file.path(save_dir, "traitarchitecture.pdf")
pdf(fname)

p1 <- plotTraitArchitecture(popA, "Fitness", "popA")
p2 <- plotTraitArchitecture(popB, "Fitness", "popB")
p3 <- plotTraitArchitecture(popA, "Additive", "popA")
p4 <- plotTraitArchitecture(popB, "Additive", "popB")

(p1|p2)/(p3|p4)
dev.off()

# Create a density plot of trait 1
trait1.df <- as.data.frame(cbind(pheno(popA)[,1], pheno(popB)[,1]))
colnames(trait1.df) <- c("popA", "popB")
trait1.df <- trait1.df %>%
  pivot_longer(c("popA", "popB"), names_to="pop", values_to="pheno")
t1 <- ggplot(trait1.df, aes(pheno, fill=pop, color=pop)) +
  geom_density(alpha=0.1) +
  labs(title="Trait 1")

# Create a density plot of trait 2
trait2.df <- as.data.frame(cbind(pheno(popA)[,2], pheno(popB)[,2]))
colnames(trait2.df) <- c("popA", "popB")
trait2.df <- trait2.df %>%
  pivot_longer(c("popA", "popB"), names_to="pop", values_to="pheno")
t2 <- ggplot(trait2.df, aes(pheno, fill=pop, color=pop)) +
  geom_density(alpha=0.1) +
  labs(title="Trait 2")


(t1|t2)
ggplot2::ggsave(filename = "trait_distributions.pdf",
                path=save_dir,
                device = "pdf")
