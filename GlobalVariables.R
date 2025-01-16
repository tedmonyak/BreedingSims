n.popSize = 1000
n.subPopSize = 500
n.ne = n.popSize
n.segSites = 1000
n.markers = 1000
n.qtlPerChr = 5
n.margin <- -0.1
n.h2 <- 0.6
n.initTraitVal <- 1
n.var <- 0.1
n.shape <- 1
n.burnInSelProp <- 0.95
n.selProp <- 0.8
n.gens <- 50
n.burnInGens <- 5
n.sims <- 100
n.chr <- 10
n.RILFams <- 200
n.indPerRILFam <- 4
n.nPops <- 2

# QTL Mapping Parameters
n.mappingMethod <- "hk"
n.errorProb <- 0.001
n.step <- 1
n.cores <- 8

# 
addSnpChip <- TRUE
basicPop <- TRUE
saveQtlPlots <- FALSE


getParams <- function() {
  n.df <- data.frame(
    popSize=n.popSize,
    subPopSize=n.subPopSize,
    ne=n.ne,
    segSites=n.segSites,
    markers=n.markers,
    qtlPerChr=n.qtlPerChr,
    margin=n.margin,
    h2=n.h2,
    initTraitVal=n.initTraitVal,
    var=n.var,
    shape=n.shape,
    burnInSelProp=n.burnInSelProp,
    selProp=n.selProp,
    gens=n.gens,
    burnInGens=n.burnInGens,
    sims=n.sims,
    chr=n.chr,
    RILFams=n.RILFams,
    indPerRILFam=n.indPerRILFam,
    nPops=n.nPops,
    mappingMethod=n.mappingMethod,
    errorProb=n.errorProb,
    step=n.step
  )
  (t(n.df))
}
