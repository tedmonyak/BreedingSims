# Fitness Functions
oneTraitFitFunc <- function(x) {
  return (-(x)^2)
}

twoTraitFitFunc <- function(x) {
  res <- -((x[,1])^2) - ((x[,2])^2)
  return (res)
}

calculateFitnessTwoTrait <- function(x,y) {
  res <- -((x)^2) -(y^2)
  return (res)
}

hetLocus <- function(locus) {
  return (length(unique(locus)) > 1)
}

getUniqueQtl <- function(pop) {
  qtlGeno <- pullQtlGeno(pop,1)
  if (pop@nTraits > 1) {
    for (t in 2:pop@nTraits) {
      qtlGeno <- cbind(qtlGeno, pullQtlGeno(pop, trait=t))
    }
  }
  qtlGeno <- qtlGeno[, !duplicated(colnames(qtlGeno))]
}

getEffectSize <- function(locus,
                          id,
                          pop,
                          methodType) {
  pop <- founderPop
  if (!hetLocus(locus)) {
    return (NA)
  }

  if (methodType == "Additive") {
    eff_sizes <- data.frame(rbind(getUniqueQtl(pop), SP$traits[[trait]]@addEff)[pop@nInd+1,])
    colnames(eff_sizes) <- c("eff_size")
    id <- "1_151"
    if (!(id %in% rownames(eff_sizes))) {
      return (NA)
    }
    return(abs(eff_sizes[id,1]))
  }
  else {
    strs = unlist(strsplit(id, "_"))
    chr = strtoi(strs[1])
    site = strtoi(strs[2])
    popSize <- nInd(pop)
    fitPre <- mean(twoTraitFitFunc(gv(pop)))
    pop1 <- editGenome(pop, ind=c(1:popSize), chr=chr, segSites=site, allele=0, simParam=SP)
    pop2 <- editGenome(pop, ind=c(1:popSize), chr=chr, segSites=site, allele=1, simParam=SP)
    return (abs(mean(twoTraitFitFunc(gv(pop2))) - mean(twoTraitFitFunc(gv(pop1)))))
  }
}

# Plots the trait architecture on a genetic basis
# Works for 1 or 2 trait populations
traitArchitecture <- function(pop, methodType="Additive", trait=1) {
  if (methodType == "Additive") {
    # This assumes that qtlGeno and qtlEffList are in the same order
    eff_sizes <- data.frame(rbind(pullQtlGeno(pop, trait), SP$traits[[trait]]@addEff)[pop@nInd+1,])
    colnames(eff_sizes) <- c("eff_size")
    eff_sizes <- rownames_to_column(eff_sizes, "id")
    eff_sizes["eff_size"] <- lapply(eff_sizes["eff_size"], abs)
    return (eff_sizes)
  } else {
    geno <- getUniqueQtl(pop)
    cols <- colnames(geno)
    nLoci <- length(cols)
    popSize <- nrow(geno)
    eff_sizes <- data.frame(id = character(nLoci),
                            eff_size=numeric(nLoci))
    for (l in 1:nLoci) {
      id <- cols[l]
      eff_sizes$id[l] <- id
      locus = geno[,l]
      eff_sizes$eff_size[l] <- getEffectSize(locus,
                                             id,
                                             pop,
                                             methodType)
      
    }
    eff_sizes <- eff_sizes[apply(eff_sizes!=0, 1, all),]
    eff_sizes <- na.omit(eff_sizes)
    eff_sizes <- unique(eff_sizes)
    eff_sizes
    return (eff_sizes)
  }
  
}
