# Title: FITNESS
# Author: Ted Monyak
# Description: This file contains functions for calculating fitness and
# plotting fitness landscapes

# Calculates fitness based on an optimum value of zero for one trait
# w = -(x^2)
oneTraitFitFunc <- function(x) {
  res <- -(x)^2
  return (res)
}

# Calculates fitness based on an optimum value of zero for each trait
# w = -(x^2) + -(y^2)
twoTraitFitFunc <- function(x) {
  res <- -((x[,1])^2) - ((x[,2])^2)
  return (res)
}

# Calculates fitness based on an optimum value of zero for each trait
# w = -(x^2) + -(y^2)
calculateFitnessTwoTrait <- function(x,y) {
  res <- -((x)^2) - (y^2)
  return (res)
}

# Calculate a decaying selection ratio based on the distance from the fitness optimum
# Uses a geometric series to determine the result, where a=(1-n.selProp),
# r is set as an initial parameter (n.r), and n is a function of the distance from the initial fitness
# w: the current fitness
# Returns: a ratio between 0 and 1 which determines what percentage of individuals to advance
selectionRatio <- function(w) {
  # Based on the simulation parameters, this is the starting fitness value
  initFit <- calculateFitnessTwoTrait(n.initTraitVal,n.initTraitVal)
  # The initial selection ratio
  a <- 1-n.selProp
  # The "n" term in the geometric series increases as the the distance from the initial fitness increases
  n <- initFit/w
  # Return a geometrically increasing value (which increases with n)
  return (1-(a * n.r^(n-1)))
}

# This section will return 1s for all of the indices in randVec that match 'idx',
# specifying which individuals to select. This ensures there are no overlaps.
# idx: The index of the subpopulation being selected
# randVec: should be created previously with this line: sample(rep(c(1:n.nPops), times=n.popSize/n.nPops))
# Returns: a vector of size randVec where all the values matching idx are 1, and all others are 0
selectSubPop <- function(x, idx, randVec) {
  as.numeric(lapply(randVec, idx, FUN=setequal))
}

# Set theme elements
theme <- theme(plot.background = ggplot2::element_blank(),
               panel.background = ggplot2::element_blank(),
               axis.line = ggplot2::element_line(linewidth = 0.2),
               plot.title = ggplot2::element_text(hjust = 0.5,
                                                  face = 'bold',
                                                  size = 12),
               axis.text = ggplot2::element_text(size  = 12,
                                                 color = 'black'),
               axis.title = ggplot2::element_text(size  = 12,
                                                  face = 'bold',
                                                  color = 'black'))



# Create a plot_ly figure of a fitness landscape, with an adaptive walk
# overlaid on it.
# df: a dataframe with 2 columns: traitValA, traitValB (for type "CONTOUR"),
# with a 3rd column (fitness) (for type "SURFACE")
# type: one of 'CONTOUR' (for a 2D landscape) or 'SURFACE' (for a 3D landscape)
# fit calc: the function for determining fitness based on two trait values
# Also, supply the min and max trait values for the fitness landscape
# TODO: add a fixed value to df, because phenotypic data is 'under' the fitness curve
overlayWalkOnLandscape <- function(df,
                                   type="CONTOUR",
                                   fitCalc,
                                   trait1Min=-1,
                                   trait1Max=1,
                                   trait2Min=-1,
                                   trait2Max=1) {
  
  # Create a matrix of fitness values, with a small increment along the x and y axes.
  fitness_x = seq(trait1Min,trait1Max, by=(trait1Max-trait1Min)/40)
  fitness_y = seq(trait2Min,trait2Max, by=(trait2Max-trait2Min)/40)
  fitness_z = outer(fitness_x,fitness_y,fitCalc)
  
  if (type == "CONTOUR") {
    fig <- plot_ly() %>%
      layout(xaxis = list(title = "Trait 1", constrain = "domain"),
             yaxis = list(title = "Trait 2", scaleanchor="x")) %>%
      add_trace(
        fig,
        x=fitness_x,
        y=fitness_y,
        z=fitness_z,
        type='contour',
        colorscale = "RdBu",
        colorbar=list(title = "Fitness"),
        line = list(color = 'black', width = 1),
        opacity=1) %>%
      add_trace(
        fig,
        df,
        name = n.popSize,
        x = df$traitValA,
        y = df$traitValB,
        type='scatter',
        mode = 'lines',
        line = list(color = 'black', width = 4, dash = 'solid'),
        opacity = 1)
    return (fig)
  } else if (type == "SURFACE"){
    fig <- plot_ly() %>%
      layout(scene = list(xaxis = list(title = "Trait 1"),
                          yaxis = list(title = "Trait 2"),
                          zaxis = list(title = "Fitness"),
                          aspectmode='cube')) %>%
      add_trace(
        fig,
        df,
        name = n.popSize,
        x = df$traitValA,
        y = df$traitValB,
        z = df$fitness,
        type = 'scatter3d',
        mode = 'lines',
        opacity = 1,
        line = list(width = 10)
      ) %>%
      add_trace(
        fig,
        x=fitness_x,
        y=fitness_y,
        z=fitness_z,
        type='surface',
        colorbar=list(title = "Fitness"),
        colors = "PuBuGn",
        opacity=0.9)
    return (fig)
  } else {
    print("Type Not Supported")
  }
}


# Plot a population on a a 3D fitness surface (only works for populations with 2 traits)
# pop: The population to plot
# fitCalc: The function for calculating fitness based on two traits
# TODO: figure out how to add colorbar label
plot3dPopulationFitness <- function(pop, fitCalc) {
  popSize <- nInd(pop)
  df <- data.frame(traitA=numeric(popSize),
                   traitB=numeric(popSize),
                   fitness=numeric(popSize))
  # Calculate the fitness for each individual in the population
  for (i in 1:popSize) {
    df$traitA[i] <- gv(pop)[i,1]
    df$traitB[i] <- gv(pop)[i,2]
    df$fitness[i] <- fitCalc(df$traitA[i], df$traitB[i])
  }
  fig <- plot_ly()
  fig <- plot_ly() %>%
    layout(scene = list(xaxis = list(title = "Trait 1"),
                        yaxis = list(title = "Trait 2"),
                        zaxis = list(title = "Fitness"),
                        aspectmode='cube')) %>%
    add_trace(
      fig,
      df,
      name = popSize,
      x = df$traitA,
      y = df$traitB,
      z = df$fitness,
      type = 'scatter3d',
      mode = 'markers',
      color=df$fitness)
  return (fig)
}

# Creates a 3D fitness landscape and returns a plot_ly rendering of it.
# Valid colors arguments:
# colors = "PuBuGn"
# colors = colorRampPalette(c("blue", "orange"))(15)
# colors = magma(50, alpha = 1, begin = 0, end = 1, direction = 1) (viridis, plasma, magma, inferno)
plotFitnessLandscape <- function() {
  fitness_x = seq(-1,1, by=0.05)
  fitness_y = seq(-1,1, by=0.05)
  fitness_z = outer(fitness_x,fitness_y,calculateFitnessTwoTrait)
  
  p <- plot_ly(x=fitness_x,
               y=fitness_y,
               z=fitness_z,
               type='surface',
               colors = "PuBuGn",
               opacity=1) %>%
    layout(scene = list(xaxis = list(title = "Trait A"),
                        yaxis = list(title = "Trait B"),
                        zaxis = list(title = "Fitness"),
                        aspectmode='cube'))
  return (p) 
}



# Show a 2D curve of a population's fitness over generations
plotFitness <- function(df) {
  g <- ggplot(data=df, aes(x=gen, y=fitness)) +
    geom_line() +
    geom_point() +
    labs(x = "Generation", y = "Fitness")
  return (g)
}

# Plot Genetic Values for Two Traits
plotHist <- function(pop) {
  gv_a = gv(pop)[,1]
  gv_b = gv(pop)[,2]
  idx = c(1:length(gv_a))
  df <- data.frame(gv_a, gv_b)
  df <- melt(as.data.table(df))
  g <- ggplot(df, aes(x=value, color=variable)) + geom_histogram(binwidth=1, position='identity')
  return (g)
}

