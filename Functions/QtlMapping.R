# Title: QTL Mapping
# Author: Ted Monyak
# Description: Contains functions for doing linkage mapping

# Determines the number of significant LOD peaks in a linkage map
# Also plots the linkage map and a PCA plot of the RIL
# RIL: a recombinant inbred line family
# parentA: one of the parents crossed to create the RIL
# parentB: the other parent crossed to create the RIL
# save_dir: a directory in which to save the plots
# Returns: the number of significant LOD peaks
getSigQtl <- function(RIL, parentA, parentB, save_dir) {
  # Create a cross object used in rqtl
  cross <- getCross(RIL, parentA, parentB, "riself")
  # Remove any null markers
  cross <- drop.nullmarkers(cross)

  # Check to see if there are no markers once null markers have been dropped
  if (length(cross$geno) == 0) {
    return (0)
  }
  # TODO: Remove this block?
  # Create a vector of size n.chr to check whether each chromosome has the minimum
  # number of markers required to do linkage mapping. If not, return '0'
  #markersPerChr <- rep(n.minMarkers, times=n.chr)
  #if (any(nmar(cross) < n.minMarkers)) {
  #  return (0)
  #}
  # Get the phenotype names
  phes <- phenames(cross)[1:2]
  # Run QTL mapping
  cross <- calc.genoprob(cross, step=n.step, error.prob=n.errorProb)
  out.hk <- scanone(cross, pheno.col=phes, method=n.mappingMethod,
                    n.cluster=n.cores)
  operm.hk <- scanone(cross, pheno.col=phes, method=n.mappingMethod,
                      n.perm=1000,n.cluster=n.cores)
  # Extract the signficant QTL based on LOD peaks being above the significance threshold
  sigQtl <- pullSigQTL(cross,
                       pheno.col=phes,
                       s1.output=out.hk,
                       perm.output=operm.hk,
                       returnQTLModel=FALSE,
                       alpha=0.05,
                       controlAcrossCol=TRUE)
  if (saveQtlPlots) {
    fname <- file.path(save_dir, "linkagemap.pdf")
    pdf(fname, width=11, height=4)
    par(mar=c(5,5,5,1))
    cols <- c("forestgreen", "gold2")
    plot(out.hk, type="n", ylim=c(0,max(as.matrix(out.hk[,-c(1:2)]))),
         main="", xlab = "Position", ylab= "LOD Score",
         cex.main=2, cex.lab=1.5, cex.axis=1.5)
    for (i in 1:length(phes)) plot(out.hk, add=T, lodcolumn = i, col = cols[i], lwd=3)
    legend(900, 25, legend=c("Trait A", "Trait B"), fill=cols)
    abline(h=summary(operm.hk[1,1]), col="forestgreen", lty = "dashed", lwd=1)
    abline(h=summary(operm.hk[1,2]), col="gold2", lty = "dashed", lwd=1)
    dev.off()
    # PCA adapted from https://github.com/HighlanderLab/jbancic_alphasimr_plants/blob/main/04_Features/simulateGWAS.R
    # Merge the genotypes of the QTL
    geno = rbind(pullSegSiteGeno(RIL), pullSegSiteGeno(parentA), pullSegSiteGeno(parentB))
    PCA  = dudi.pca(df = geno, center = T, scale = F, scannf = F, nf = 5)
    (VAF = 100 * PCA$eig[1:5] / sum(PCA$eig)) # variance explained
    df.PCA = data.frame(
      "Pop" = c(rep("RIL", RIL@nInd), "Parent A", "Parent B"),
      "PC1" = PCA$l1$RS1,
      "PC2" = PCA$l1$RS2)
    
    ggplot(df.PCA, aes(x = PC1, y = PC2)) +
      geom_point(aes(colour = factor(Pop))) +
      ggtitle("Population structure") +
      xlab(paste("Pcomp1: ", round(VAF[1], 2), "%", sep = "")) +
      ylab(paste("Pcomp2: ", round(VAF[2], 2), "%", sep = ""))
    
    fname <- file.path(save_dir, "PCA.pdf")
    ggplot2::ggsave(filename = fname,
                    device = "pdf",
                    width=8,
                    height=8)
  }

  write.table(getParams(), file.path(save_dir, "params.txt"), col.names=FALSE, quote=FALSE, sep=":\t")
  return (nrow(sigQtl))
}
