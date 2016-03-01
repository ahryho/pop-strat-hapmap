# --- Clustering --- #

source("packages.R")

# --- Function for computing clustering accuracy --- #
accuracyCalc <- function(confTbl, startCol = 1)
{
  corr <- 0;
  for(i in startCol:ncol(confTbl))
  {
    corr <- corr + max(confTbl[, i])  
  }
  accuracy <- corr / sum(confTbl)
  return(accuracy)
}

# --- Function for computing plotting the clusters --- #
plot.clusters <- function(data.clustering, hc.groups, nclusters)
{
  par(mfrow = c(1,1))
  plotcluster(data.clustering, hc.groups, main = paste0("Clustering in ", nclusters, " classes"))
  # grid(20, 20, lty = 1, lwd = 1)
  # text(data.clustering[, 1], data.clustering[, 2], col = hc.groups)
  abline(h = 0, v = 0, col = "gray")
  # legend("bottom", paste("Cl", 1:nclusters, sep =" "), pch = 2, col = 1:nc, bty = "n", horiz = T)
  
  return(recordPlot())
}

genofile <- snpgdsOpen(snpgdsExampleFileName()) # Open the GDS file
pop.group <- read.gdsn(index.gdsn(genofile, path = "sample.annot/pop.group"))

set.seed(1000)
snpset <- snpgdsLDpruning(genofile, ld.threshold = 1);
snpset.id <- unlist(snpset)

pca <- snpgdsPCA(genofile, snp.id = snpset.id, num.thread = 3)

numberEigValues <- length(na.omit(pca$varprop))
(pca.var <- pca$varprop[1:numberEigValues] * 100)

par(mfrow = c(1,1))
plot(pca.var, type = "n",  main = "Scree plot",  xlab = "Principal Components", ylab = "Eigenvalues")
lines(x = 1:numberEigValues, pca.var, type="b", pch = 19, col = "red")
axis(1, seq(0,numberEigValues, 1))

# To determine how much PCs we need to take , the inertia will be computed
# SNP correlations between eigenvactors and SNP genotypes
chr <- read.gdsn(index.gdsn(genofile, "snp.chromosome"))
CORR <- snpgdsPCACorr(pca, genofile)

# To compute the inertia of each column and the total inertia the correlation matrix will be transposed
corrMatr <- as.data.frame(t(CORR$snpcorr))
corrMatr <- corrMatr[-which(is.na(corrMatr[, 1])),]

Inertia <- matrix(0, nrow = 1, ncol = ncol(corrMatr))
for (j in 1:ncol(Inertia))
  Inertia[j] <- sum(corrMatr[, j] ^ 2)
totalInertia <- sum(Inertia)
Inertia <- Inertia / totalInertia * 100
cumInertia <- cumsum(Inertia)
pca.eig <- cbind(pca.var, t(Inertia), cumInertia)
colnames(pca.eig) <- c("eigenvalue", "percentage of variance", "cumulative percenatge of variance") 
pca.eig

# According to Kaiser rule we should to choose the first three PCs, the eigenvalues of which are greater that 1,
# But these PCs don't contain even 50% of total inertia, hence for the further analysis there will be taken
# the first 15 PCs that represent 70% of total variance for the first experiment, and three is for the second one

nPCs = 3
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
data.gen <- data.frame(sample.id = pca$sample.id, pop = factor(pop.group)[match(pca$sample.id, sample.id)],
                       PC = scale(pca$eigenvect[, 1:nPCs], center = T, scale = T), stringsAsFactors = FALSE)

number_of_total_observations <- nrow(data.gen)
data.gen$pop <- as.factor(data.gen$pop)
# levels(data.gen$pop)  # 1 - CEU, 2 - HCB, 3 - JRT, 4 - YRI

data.clustering <- data.gen[, -c(1, 2)]
rownames(data.clustering) <-  data.gen[, 1]

# --- 1. Hierarchical clustering --- #

set.seed(678)
distM <- dist(data.clustering, method = "manhattan") # minkowski, euclidean, manhattan

hc <- hclust(distM, method = "ward.D")
par(mfrow = c(1,2))
barplot(hc$height, main = "Height histogram")
plot(hc, hang = -1, labels = data.gen$pop)
nc <- 4 # number of clusters
rect.hclust(hc, k = nc)

hc.groups <- cutree(hc, k = nc)
plot.clusters(data.clustering, hc.groups, nc)

xtable(result <- table(data.gen$pop, hc.groups))
accuracy.hierarchical.3PCs <- accuracyCalc(result) * 100 # 83.15%

cdg <- aggregate(as.data.frame(data.clustering), list(hc.groups), mean)[, -1]

SSB <- sum(rowSums(cdg ^ 2) * as.numeric(table(hc.groups))) # the sum of squares between clusters
SST <- sum(rowSums(data.clustering ^ 2))
(RS.hclust <- 100 * SSB/SST) # 96.33
(CH.hclust <- (SSB / (nc - 1)) / ((SST - SSB)/(nrow(data.clustering) - nc))) # 2403.813

xtable(cbind(accuracy.hierarchical.3PCs, SSB, SST, RS.hclust, CH.hclust))
data.profiling <- cbind(as.factor(hc.groups), data.clustering) 

catdes(data.profiling, 1)

# --- 2. k-means clustering --- #

set.seed(124)
k.means <- kmeans(data.clustering, cdg)

plot.clusters(data.clustering, k.means$cluster, nc)

result <- table(data.gen$pop, k.means$cluster)
accuracy.kmeans.3PCs <- accuracyCalc(result) # 83.15%

SSB <- sum(rowSums(k.means$centers ^ 2) * k.means$size) # the sum of squares between clusters
SSW <- sum(k.means$withinss) # the sum of squares within clusters
(Ib.kmeans <- 100 * SSB / (SSB + SSW)) # 96.32
(CH.kmeans <- (SSB / (nc - 1)) / (SSW/(nrow(data.gen) - nc))) #2403.813

# --- 3. DBScan --- #

# --- Function for estimation the parameter --- #
dbscan.tune <- function(data, target, epsilon, pts)
{
  for (e in epsilon)
    for (p in pts)
    {
      card.dbscan <- dbscan(data, eps = e, MinPts = p)
      res <- table(target, card.dbscan$cluster)
      accuracies[toString(p), toString(e)] <- accuracyCalc(res, 2)
    }
  accuracies
}

set.seed(121)

epsilon <- (1:10) * 10^(-1)
pts <- 1:10

accuracies <- matrix(nrow = length(pts), ncol = length(epsilon))
colnames(accuracies) <- epsilon
rownames(accuracies) <- pts

accuracies <- dbscan.tune(data.clustering, data.gen$pop, epsilon, pts)

xtable(accuracies)

gen.dbscan <- dbscan(data.clustering, eps = 0.7, MinPts = 3)
plot.clusters(data.clustering, gen.dbscan$cluster, 4)

xtable(result <- table(data.gen$pop, gen.dbscan$cluster))
accuracy.dbscan.3PCs <- accuracyCalc(result, 2) * 100# 81.72

cdg <- aggregate(as.data.frame(data.clustering), list(gen.dbscan$cluster), mean)[, -1]

SSB <- sum(rowSums(cdg ^ 2) * as.numeric(table(gen.dbscan$cluster))) # the sum of squares between clusters
SST <- sum(rowSums(data.clustering ^ 2))
(RS.dbscan <- 100 * SSB/SST) # 92.36
CH.dbscan <- (SSB / (nc - 1)) / ((SST - SSB)/(nrow(data.gen) - nc)) # 1108.632

xtable(cbind(accuracy.dbscan.3PCs, SSB, SST, RS.dbscan, CH.dbscan))

# --- 4. Probalistic clustering (EM - clustering) --- #

emc <- Mclust(data.clustering, G = 1:3) 
summary(emc, parameters = T)
cdg <- aggregate(as.data.frame(data.clustering), list(emc$classification), mean)[, -1]

SSB <- sum(rowSums(cdg ^ 2) * as.numeric(table(emc$classification)))
SST <- sum(rowSums(data.clustering ^ 2))
(RS.emc <- 100 * SSB / SST) # 92.11
(CH.emc <- (SSB / (nc - 1)) / ((SST - SSB)/(nrow(data.gen) - nc))) # 1070.72

xtable(result <- table(data.gen$pop, emc$classification))
accuracy.emc.3PCs <- accuracyCalc(result) * 100 # 83.15

xtable(cbind(accuracy.emc.3PCs, SSB, SST, RS.emc, CH.emc))

# --- Consolidation --- #
k.means.emc <- kmeans(data.clustering, centers = cdg)
k.means.emc$size

SSB <- sum(rowSums(k.means.emc$centers ^ 2) * k.means.emc$size)
SSW <- sum(k.means.emc$withinss)
Ib.kmeans.emc <- 100 * SSB / (SSB + SSW) # 96.32
CH.kmeans.emc <- (SSB / (nc - 1)) / (SSW/(nrow(data.gen) - nc)) # 2403.813 # Calinski & Harabasz Index

emc$BIC
plot(emc, cex = 0.4)
