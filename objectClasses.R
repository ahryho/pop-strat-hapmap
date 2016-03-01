Prediction <- setClass(Class = "Prediction", 
                        slots = c(confusionMatrix.train = "table",
                                  confusionMatrix.test = "table",
                                  train.error = "numeric", 
                                  test.error = "numeric",
                                  Precision.train = "numeric",
                                  Recall.train = "numeric",
                                  Precision.test = "numeric",
                                  Recall.test = "numeric",
                                  auc = "numeric",
                                  predicted = "numeric"
                                  ),
                        prototype = list(train.error = 0,
                                         test.error = 0,
                                         predicted = 0)
                        )

#setMethod(f = "plot.prediction", "Prediction", plot.prediction(.knn))

myknn <- setClass(Class = "myknn", 
                 slots = c(k.estimation.table = "matrix",
                           k = "numeric"),
                 prototype = list(k = 1),
                 contains = "Prediction")

mysvm <- setClass(Class = "mysvm",
                  slots = c(
                    cross.errors = "matrix",
                    nSV = "matrix",
                    cost = "numeric",
                    sigma = "numeric"
                    # model = "svm"
                    ),
                  contains = "Prediction")

nb <- setClass(Class = "nb", # Naive Bayes
         contains = "Prediction"
)

dt <- setClass(Class = "dt", # Decision trees
               slots = c(
                 # model.tree = "rpart",
                 model.cptable = "matrix",
                 plot.cptable = "recordedplot",
                 plot.tree = "recordedplot"
               ),
         contains = "Prediction"
)

nn <- setClass(Class = "nn", # Neural Network
         contains = "Prediction"
)

lda.qda <- setClass(Class = "lda.qda",
                    slot = c( x = "matrix"),
               contains = "Prediction"
)

rf <- setClass(Class = "rf", # Random Forest
               slots = c(
                 plot.model = "recordedplot",
                 plot.importance = "recordedplot",
                 plot.roc.curve = "recordedplot"
               ),
         contains = "Prediction"
)