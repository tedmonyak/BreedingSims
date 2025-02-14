#Title: DIFFERENT POPULATION SIZES
#Author: Ted Monyak
#Description: Simulate several adaptive walks with different population sizes

setwd("~/Documents/CSU/R/BreedingSims")
output_dir <- file.path(getwd(), "Output")

# Different starting population sizes
subPopSizes <- c(50,500)
n.selProp <- 0.1
n.h2 <- 0.2
n.var <- 0.05
n.qtlPerChr <- 2
n.gens <- 50

# Don't use a SNP chip for these simulations
addSnpChip <- FALSE
source("Scripts/GlobalParameters.R")
source("Scripts/CreateFounderPop.R")

f <- list(family="Arial", size=16)

fig <- plot_ly() %>%
  layout(legend=list(title=list(text='Subpopulation Size', font=f)),
         scene = list(xaxis = list(title = "Trait A", font=f),
                      yaxis = list(title = "Trait B", font=f),
                      zaxis = list(title = "Fitness", font=f),
                      aspectmode='cube'))

pop <- founderPop
for (gen in 1:n.burnInGens) {
  pop <- selectCross(pop, trait=twoTraitFitFunc, nInd=nInd(pop)*n.burnInSelProp, nCrosses=nInd(pop))
}

# Iterate through each population size
for (p in 1:length(subPopSizes)) {
  n.subPopSize <- subPopSizes[p]
  fit_df <- data.frame(gen=1:n.gens,
                       fitness=numeric(n.gens),
                       traitValA=numeric(n.gens),
                       traitValB=numeric(n.gens))
  subPop <- selectInd(pop, nInd=n.subPopSize, use="rand")
  # Iterate through the generations, update the result dataframe, and advance
  # progeny based on the two trait fitness funcion
  for(gen in 1:n.gens) {
    meanFitness <- mean(twoTraitFitFunc(pheno(subPop)))
    fit_df$fitness[gen] <- calculateFitnessTwoTraitModified(meanP(subPop)[1], meanP(subPop)[2])
    fit_df$traitValA[gen] <- meanP(subPop)[1]
    fit_df$traitValB[gen] <- meanP(subPop)[2]
    selRat <- selectionRatio(meanFitness)
    subPop <- selectCross(subPop, trait=twoTraitFitFunc, nInd=n.subPopSize*selRat, nCrosses=n.subPopSize)
  }
  # Add a trace to represent this population's adaptive walk
  fig <- fig %>% add_trace(
    fig,
    fit_df,
    name = n.subPopSize,
    x = fit_df$traitValA,
    y = fit_df$traitValB,
    z = fit_df$fitness,
    type = 'scatter3d',
    mode = 'lines',
    opacity = 1,
    line = list(autocolorscale=FALSE, color=p, which=2, width = 10)
  )
}

# Create a matrix of fitness values, with a small increment along the x and y axes.
fitness_x = seq(n.initTraitVal*-1,n.initTraitVal*1, by=(n.initTraitVal/100))
fitness_y = seq(n.initTraitVal*-1,n.initTraitVal*1, by=(n.initTraitVal/100))
fitness_z = outer(fitness_x,fitness_y,calculateFitnessTwoTrait)

fig <- fig %>%
  add_trace(
    fig,
    x=fitness_x,
    y=fitness_y,
    z=fitness_z,
    type='surface',
    colorbar=list(title = "Fitness"),
    colors = viridis(n=10),
    opacity=1.0)

save_dir <- file.path(output_dir, "DifferentPopulationSizes")
if (!dir.exists(save_dir)) dir.create(save_dir)
save_dir <- file.path(save_dir, format(Sys.time(), "%F_%H_%M_%S"))
if (!dir.exists(save_dir)) dir.create(save_dir)
write.table(getParams(), file.path(save_dir, "params.txt"), col.names=FALSE, quote=FALSE, sep=":\t")

fname <- file.path(save_dir, paste0("populationSizes_", subPopSizes[1], ",", subPopSizes[2], ".html"))
htmlwidgets::saveWidget(as_widget(fig), fname)
fig