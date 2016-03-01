# Install the package #
source("http://bioconductor.org/biocLite.R")

# Load the R packages: gdsfmt and SNPRelate
if(!require(gdsfmt)) {
  biocLite("gdsfmt"); require(gdsfmt)}

if(!require(SNPRelate)) {
  biocLite("SNPRelate"); require(SNPRelate)}

if(!require(class)) {
  install.packages("class"); require(class)} # knn

if(!require(caret)) {
  install.packages("caret"); require(caret)} # knn (k-fold CV), nb

if(!require(e1071)) {
  install.packages("e1071"); require(e1071)} # nb, svm

if(!require("kernlab")) {
  install.packages("kernlab"); require(kernlab)} # ksvm

if(!require("rpart")) {
  install.packages("rpart"); require(rpart)} # decision tree

if(!require(rattle)) {
  install.packages("rattle"); require(rattle)} # for the decision tree plot (fancyPlot)

if(!require(randomForest)) {
  install.packages("randomForest"); require(randomForest)}

if(!require(neuralnet)) {
  install.packages("neuralnet"); require(neuralnet)} # NN

if(!require(ROCR)) {
  install.packages("ROCR"); require(ROCR)} # ROC curve

if(!require(pROC)) {
  install.packages("pROC"); require(pROC)} # ROC curve

if(!require(fpc)) {
  install.packages("fpc"); require(fpc)} # pamk, dbscan

if(!require(mclust)) {
  install.packages("mclust"); require(mclust)} # E-M clustering

if(!require(FactoMineR)) {
  install.packages("FactoMineR"); require(FactoMineR)} 