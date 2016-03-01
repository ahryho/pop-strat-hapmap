# Install the package #
source("http://bioconductor.org/biocLite.R")

# Load the R packages: gdsfmt and SNPRelate
if(!require(gdsfmt)) {
  biocLite("gdsfmt"); require(gdsfmt)}

if(!require(SNPRelate)) {
  biocLite("SNPRelate"); require(SNPRelate)}

if(!require(class)) {
  install.packages("class"); require(class)} # knn

genofile <- snpgdsOpen(snpgdsExampleFileName()) # Open the GDS file

# snp.chromosome, an integer or character mapping for each chromosome. 
# Integer: numeric values 1-26, mapped in order from 1-22, 23=X, 24=XY (the pseudoautosomal region), 
# 25=Y, 26=M (the mitochondrial probes), and 0 for probes with unknown positions; it does not allow NA.

get.attr.gdsn(index.gdsn(genofile, "snp.chromosome")) 

#read.gdsn(index.gdsn(genofile, "snp.chromosome"))

pop.group <- read.gdsn(index.gdsn(genofile, path = "sample.annot/pop.group"))
table(pop.group)
# Yoruba in Ibadan, Nigeria (abbreviation: YRI)
# Japanese in Tokyo, Japan (abbreviation: JPT)
# Han Chinese in Beijing, China (abbreviation: CHB)
# CEPH (Utah residents with ancestry from northern and western Europe) (abbreviation: CEU)

# --- Data Analysis --- #

set.seed(1000)
snpset <- snpgdsLDpruning(genofile, ld.threshold = 1)
head(snpset)
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
# snpgdsPCA(genofile, eigen.cnt = 10)
CORR$snpcorr[, 1:8] # 32 x 9088

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
# But these PCs don't contain even 50% of total inertia, hence for further analysis there will be taken
# the first 15 PCs and 20 PCs that represent 70% and 80% of total variance respectively (two experiments).

nPCs = 20
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
data.gen <- data.frame(sample.id = pca$sample.id, pop = factor(pop.group)[match(pca$sample.id, sample.id)],
                          PC = pca$eigenvect[, 1:nPCs], stringsAsFactors = FALSE)

head(data.gen)
summary(data.gen)

# Plot the principal component pairs for the first 3 PCs
lbls <- paste("PC", 1:3, "\n", format(pca.var[1:3], digits = 2), "%", sep="")
pairs(pca$eigenvect[, 1:3], col = data.gen$pop, labels = lbls)

savepar <- par(mfrow = c(4,1), mai=c(0.3, 0.55, 0.1, 0.25))
for (i in 1:4)
  plot(abs(CORR$snpcorr[i,]), ylim = c(0,1), xlab = "", ylab = paste("PC", i), col = chr, pch = "+")

par(savepar)

# g <- read.gdsn(index.gdsn(genofile, "genotype"), start=c(1,1), count = c(-1, -1))
# genotype <- as.data.frame(g)

# --- Classification --- #
# data.gen is a dataset used for analysis. Each observation represents a person, who took apart in experiment
# each column is a PC, which was obtained after PCA and represents a linear combination of
# all SNPs which were determined.

number_of_total_observations <- nrow(data.gen)
data.gen$pop <- as.factor(data.gen$pop)
levels(data.gen$pop) <- 1:4 # 1 - CEU, 2 - HCB, 3 - JRT, 4 - YRI
train.data.gen <- data.gen[1:round(0.67 * number_of_total_observations), -1]
test.data.gen <- data.gen[-(1:round(0.67 * number_of_total_observations)), -1]
  

