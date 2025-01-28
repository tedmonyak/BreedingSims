# Title: TRAIT ARCHITECTURE
# Author: Ted Monyak
# Description: Contains functions for examining and plotting trait architectures

# Determine whether a population is segregating at a locus
# locus: a list of genotypes from a population
# Returns: true if a locus is heterozygous within a population
hetLocus <- function(locus) {
  return (length(unique(locus)) > 1)
}
# Get the QTL of all traits in a population
# Since QTL can be overlapping between traits, we must filter out the
# duplicated QTL
# Returns: a dataframe of the QTL (where columns are the QTL and rows are the individuals)
getUniqueQtl <- function(pop) {
  # Get the qtl from trait 1
  qtlGeno <- pullQtlGeno(pop,1)
  # Get the qtl from all other traits
  if (pop@nTraits > 1) {
    for (t in 2:pop@nTraits) {
      qtlGeno <- cbind(qtlGeno, pullQtlGeno(pop, trait=t))
    }
  }
  # Remove duplicate QTL
  qtlGeno <- qtlGeno[, !duplicated(colnames(qtlGeno))]
  qtlGeno <- qtlGeno[,order(colnames(qtlGeno))]
  return (qtlGeno)
}

# Gets the absolute value of the effect size for each qtl
# Assumes there are 2 traits with QTL of additive effects
# Calculates effect size as sqrt(e1^2 + e2^2) (where eN is effect size of trait N)
# pop: The population
# Returns a dataframe with the qtl indices as the rows, and one column: 'eff_size'
getQtlEffectSizes <- function(pop) {
  # For each trait:
  # Merge the QTL dataframes with the dataframe of the effect sizes, and filter
  # out rows 1:pop@nInd, which are the genotypes, which we don't care about
  e1 <- data.frame(rbind(pullQtlGeno(pop, 1),
                         SP$traits[[1]]@addEff)[pop@nInd+1,])
  e2 <- data.frame(rbind(pullQtlGeno(pop, 2),
                         SP$traits[[2]]@addEff)[pop@nInd+1,])

  # Outer join the two dataframes
  eff_sizes <- merge(e1, e2, by="row.names", all=TRUE)
  # Change all NA values to zero
  eff_sizes[is.na(eff_sizes)] <- 0
  colnames(eff_sizes) <- c("snp", "eff1", "eff2")
  # Create a 3rd column called eff_size which is sqrt(e1^2 + e2^2) and move 1st
  # column back to the index spot
  eff_sizes <- eff_sizes %>%
    mutate(eff_sizes, eff_size=sqrt(eff1^2 + eff2^2)) %>%
    column_to_rownames(., var="snp")
  # Only store the generated effect sizes from both traits
  eff_sizes <- eff_sizes[3]
  return(eff_sizes)
}

# Calculates the effect size of the QTL
# locus: a list of genotypes from a population at a particular locus
# id: the name of the QTL in the form chr_loc
# pop: the population in question
# methodType: one of 'Additive' (for determining the additive effect) or
# 'Fitness' (for determining the effect of the allele on fitness)
# trait: the index of the trait to get an effect size for (only matters if method
# type is 'Additive')
# Returns: an effect size > 0, or NA if the locus is not segregating in the
# population, or it isn't a QTL
# Before using this function, you should make sure to pass a locus that is from
# a pullQtlGeno() or getUniqueQtl() call, not a pullSegSiteGeno() call,
# since only QTL have effect sizes of > 0
getEffectSize <- function(locus,
                          id,
                          pop,
                          methodType="Additive") {
  # Only segregating alleles have an effect size that is worth measuring
  if (!hetLocus(locus)) {
    return (NA)
  }
  if (methodType == "Additive") {
    eff_sizes <- getQtlEffectSizes(pop)
    # If the QTL we are looking for is not one of the QTL for this trait, return NA
    if (!(id %in% rownames(eff_sizes))) {
      return (NA)
    }
    # Return the absolute value of the effect size for the allele with name 'id'
    return(eff_sizes[id,1])
  }
  else {
    # Determine the chromosome and site id, to use for editing the genome
    strs = unlist(strsplit(id, "_"))
    chr = strtoi(strs[1])
    site = strtoi(strs[2])
    popSize <- nInd(pop)
    # Modify the population so the allele is set to 0 at all loci
    pop0 <- editGenome(pop, ind=c(1:popSize), chr=chr, segSites=site, allele=0, simParam=SP)
    # Modify the population so the allele is set to 1 at all loci
    pop1 <- editGenome(pop, ind=c(1:popSize), chr=chr, segSites=site, allele=1, simParam=SP)
    # The effect size is the absolute value of the difference in mean fitness values
    # between the populations that differ only at the allele in question
    return (abs(mean(twoTraitFitFunc(gv(pop1))) - mean(twoTraitFitFunc(gv(pop0)))))
  }
}

# Determines the genetic architecture of a trait in a population
# pop: the population of interest
# methodType: one of 'Additive' (for determining the additive effect) or
# 'Fitness' (for determining the effect of the allele on fitness)
# trait: the index of the trait to get an effect size for (only matters if method
# type is 'Additive')
# Returns a dataframe with the following columns: 'id', and 'eff_size'
traitArchitecture <- function(pop, methodType="Additive") {
  if (methodType == "Additive") {
    # Determine the heterozygous / segregating alleles
    hetLoci <- as.data.frame(apply(getUniqueQtl(pop), MARGIN=2, FUN=hetLocus))
    colnames(hetLoci) <- "het"
    # Join the het table with the effect size table
    eff_sizes <- merge(hetLoci, getQtlEffectSizes(pop),by="row.names")
    colnames(eff_sizes)[1] <- "id"
    # Filter out any non-segregating alleles
    eff_sizes <- eff_sizes %>% filter(het)
    eff_sizes <- eff_sizes[c(1,3)]
    return (eff_sizes)
  } else {
    # Get all QTL for all traits in a population
    geno <- getUniqueQtl(pop)
    # Get names of qtl
    cols <- colnames(geno)
    nLoci <- length(cols)
    popSize <- nrow(geno)
    eff_sizes <- data.frame(id = character(nLoci),
                            eff_size=numeric(nLoci))
    # Calculate the effect size for each allele
    for (l in 1:nLoci) {
      id <- cols[l]
      eff_sizes$id[l] <- id
      # locus is all of the genotypes in the population at a given locus
      locus = geno[,l]
      eff_sizes$eff_size[l] <- getEffectSize(locus, id, pop, methodType)
    }
    # Filter out all effect sizes of zero
    eff_sizes <- eff_sizes[apply(eff_sizes!=0, 1, all),]
    # Filter out all NA effect sizes
    eff_sizes <- na.omit(eff_sizes)
    return (eff_sizes)
  }
}

# This function will return the additive effect sizes for the QTL for a given trait,
# and the genetic location in morgans
# pop: The population in question
# trait: the index of the trait to get effect sizes for
# Returns: a dataframe with columns "snp" (the qtl id), "eff_size" (the effect size),
# "pos" (the genetic location in morgans), and "chr", the chromosome
getPerTraitQtlEffectSizesAndLocations <- function(pop, trait) {

  # Merge the QTL dataframes with the dataframe of the effect sizes, and filter
  # out rows 1:pop@nInd, which are the genotypes, which we don't care about
  eff_sizes <- data.frame(rbind(pullQtlGeno(pop, trait),
                                SP$traits[[trait]]@addEff)[pop@nInd+1,])
  colnames(eff_sizes) <- c("eff_size")
  # Make effect sizes positive
  eff_sizes["eff_size"] <- lapply(eff_sizes["eff_size"], abs)

  # Get the genetic map, and flatten it so it isn't grouped by chromosome
  genMap <- as.data.frame(unlist(SP$genMap))
  # Remove the chromosome prefix from each id, to get ids in the form chr_loc
  rownames(genMap) <- sub(".*\\.", "", rownames(genMap))
  # Merge effect sizes and locations on the index column (snp id)
  eff_sizes <- merge(eff_sizes, genMap, by="row.names")
  colnames(eff_sizes) <- c("snp", "eff_size", "pos")
  # Determine the chromosome of each qtl, which is stored in the snp id
  eff_sizes <- eff_sizes %>%
    mutate(chr = sub("_.*", "", snp))
  return(eff_sizes)
}

# This function will create a ggplot of the trait architecture
# pop: The population to calculate trait architecture for
# methodType: to use for the calculation of effect size. Either 'Additive', for
# the underlying additive effect size, or 'Fitness', to determine the effect of
# each allele on fitness
# trait: the index of the trait under question
# Returns: a ggplot
plotTraitArchitecture <- function(pop, methodType="Additive", popName="") {
  eff_sizes <- traitArchitecture(pop, methodType)
  # Create a plot with the effect sizes ranked in descending order
  g <- ggplot(data=eff_sizes, aes(x=reorder(id, -eff_size), y=eff_size)) +
    geom_bar(stat="identity") +
    labs(x = "Variant Id", y = "Effect Size", title=paste0(methodType, "_", popName)) +
    theme +
    theme(axis.text.x = element_text(angle = 45,
                                     vjust = 0.5,
                                     hjust=1,
                                     size = 6,
                                     margin = margin(b = 10)),
          axis.text.y = element_text(margin = margin(l=10, r=10)))
  return (g)
}
