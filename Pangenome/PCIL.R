library(AlphaSimR)

set.seed(42)

# 4 recurrent parents
n.recurrentParents <- 4
recurrentParentNames <- c("MOTOMARADI", "MACIA", "IRAT204", "CSM-63E")

# 127 donor parents
n.donorParents <- 127

# Sorghum has 10 chromosomes
n.chr <- 10
# Seg Sites per Chromosome
n.segSites <- 1000

n.progeny <- 1

# number of clans
n.clans <- n.recurrentParents

founders = runMacs(
  nInd=n.recurrentParents + n.donorParents,
  nChr=n.chr,
  segSites=n.segSites,
)

SP = SimParam$new(founders)
SP$setTrackPed(TRUE)

# Modify this to be newMapPop() to supply a genetic map and haplotype map
founderPop <- newPop(founders, simParam = SP)

recurrentParents <- newPop(founderPop[1:n.recurrentParents], id=recurrentParentNames)
donorParents <- newPop(founderPop[(n.recurrentParents+1):(n.recurrentParents+n.donorParents)], id=paste0("D", 1:n.donorParents))

recurrentGeno = pullSegSiteGeno(recurrentParents)
donorGeno = pullSegSiteGeno(donorParents)

# Randomize the assignment of donor lines to recurrent parents
# Modify this to assign particular donors to clans
clans <- sample(rep(c(1:n.clans), times=n.donorParents/n.clans))
# Ensure that clan assignment vector is the same size as the number of donor lines
gap <- n.donorParents - length(clans)
if (gap > 0) {
  clans <- c(clans, sample(1:gap))
}

# Create a crossing scheme where each donor is crossed with the recurrent parent to which it is assigned
# Female parents are the recurrent lines
crossPlan <- matrix(c(clans, rep(1:n.donorParents)), nrow=n.donorParents, ncol=2)


# Adjust nProgeny to get more progeny per cross
F1 <- makeCross2(recurrentParents, donorParents, crossPlan, nProgeny=n.progeny)
BC1F1 <- makeCross2(recurrentParents, F1, crossPlan)
BC1F2 <- self(BC1)
BC1F3 <- self(BC1F2)

PCIL_Geno <- pullSegSiteGeno(BC1F3)

# consider making explicit families of each pedigree
pca_recurrent <- prcomp(recurrentGeno)
eigen_recurrent <- pca_recurrent$sdev^2
genetic_distances_recurrent <- dist(recurrentGeno, method = "euclidean")
hclust_recurrent <- hclust(genetic_distances_recurrent, method = "average")

pca_PCIL <- prcomp(PCIL_Geno)
eigen_PCIL <- pca_PCIL$sdev^2
genetic_distances_PCIL <- dist(PCIL_Geno, method = "euclidean")
hclust_PCIL <- hclust(genetic_distances_PCIL, method = "average")

par(mfrow = c(3,2))
par(mar=c(1,1,1,1))
color_palette_pcil <- rainbow(nInd(BC1F3))
color_palette_recurrent <- c("blue", "red", "yellow", "green")
#color_palette <- (c("blue", "red", "orange", "gold", "yellow", "darkorange3", "deeppink3", "khaki3", "lightsalmon", "tan3", "wheat2"))
#colors <- rep(color_palette[=:11], each = n.progeny) # creating vector of colors

plot(eigen_recurrent, type = "b", col = "black",
     xlab = "Principal Component", ylab = "Eigenvalue", 
     main="Variance Explained PC in Recurrent Parents")

plot(eigen_PCIL, type = "b", col = "black",
     xlab = "Principal Component", ylab = "Eigenvalue", 
     main="Variance Explained PC in PCIL")

plot(x=pca_recurrent$x[,1],y=pca_recurrent$x[,2], main="PCA of Recurrent Parents", col = color_palette_recurrent)
plot(x=pca_PCIL$x[,1],y=pca_PCIL$x[,2], main="PCA of PCIL Families", col = color_palette_pcil)

plot(hclust_recurrent, main = "Clustering of Recurrent Parents", labels = FALSE)

#selected
plot(hclust_PCIL, main = "Clustering of PCIL Families", labels = FALSE)
