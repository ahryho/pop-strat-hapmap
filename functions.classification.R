setGeneric(name = "setTrainError", 
           def = function(theObject, error)
           {
             standardGeneric("setTrainError")
           })

setMethod(f = "setTrainError", 
          signature = "Prediction",
          definition = function(theObject, error){
            theObject@train.error <- as.numeric(error)
            return(theObject)
          })

setGeneric(name = "setTestError", 
           def = function(theObject, error)
           {
             standardGeneric("setTestError")
           })

setMethod(f = "setTestError", 
          signature = "Prediction",
          definition = function(theObject, error){
            theObject@test.error <- as.numeric(error)
            return(theObject)
          })

setGeneric(name = "setPrediction", 
           def = function(theObject, prediction)
           {
             standardGeneric("setPrediction")
           })

setMethod(f = "setPrediction", 
          signature = "Prediction",
          definition = function(theObject, prediction){
            theObject@predicted <- as.numeric(prediction)
            return(theObject)
          })

# --- knn methods --- #
setGeneric(name = "knn.k.estimation", 
                      def = function(theObject, train.data.set, cv = "loo")
                        {
                          standardGeneric("knn.k.estimation")
                      })

# --- Function for estimation the "k" parameter --- #

setMethod(f = "knn.k.estimation",
                signature = "myknn",
                definition = function(theObject, train.data.set, cv = "loo")
                  {
                    neighbours <- 1:10
                    errors <- matrix(nrow = length(neighbours), ncol = 2)
                    colnames(errors) <- c("Nr of neighbours", "Error")

                    train.data <- subset(train.data.set, select = - pop)
                    train.class <- train.data.set$pop # "pop" is a target class
  
                    for (K in neighbours)
                    {
                      error <- 0
                      switch(cv,
                             
                        "loo" = { # knn with LOOCV
                          prediction <- knn.cv(train.data, train.class, k = K, prob = F)
                          confusion.matrix <- table(prediction, train.class)
                          error <- 100 - sum(diag(confusion.matrix)) / sum(confusion.matrix) * 100
                          }, 
                        
                        "k-fold" = {
                          n = 10
                          for (population in levels(train.data$pop))
                          {
                            train.data.temp <- train.data[train.data$pop != population,]
                            validation.data <- train.data[train.data$pop == population,]
                            validation.class <- validation.data$pop
                            prediction <- knn(train.data.temp[, -1], validation.data[, -1], train.data.temp$pop, k = K, prob = F) 
                            confusion.matrix <- table(prediction, validation.class)
                            error <- error + 100 - sum(diag(confusion.matrix)) / sum(confusion.matrix) * 100 / length(levels(train.data$pop))
                          }
                        }
                      )
                      errors[K, "Nr of neighbours"] <- neighbours[K]
                      errors[K, "Error"] <- error
                    }
                    theObject@confusionMatrix.train <- confusion.matrix
                    theObject@k.estimation.table <- errors
                    theObject@k <- as.numeric(errors[errors[, "Error"] == (min(errors[, "Error"])), "Nr of neighbours"][1])
                    theObject <- setTrainError(theObject, errors[errors[,"Nr of neighbours"] == theObject@k, "Error"])
                    theObject@Precision.train <- diag(theObject@confusionMatrix.train) / rowSums(theObject@confusionMatrix.train)
                    theObject@Recall.train <- diag(theObject@confusionMatrix.train) / colSums(theObject@confusionMatrix.train)
                return(theObject)
                }
)

# split.data <- function(data.set){
#   data <- subset(data.set, select = - pop)
#   class <- data.set$pop # "pop" is a target class
#   return(list(data, class))
# }

setGeneric(name = "knn.prediction", 
           def = function(theObject, train.data.set, test.data.set)
           {
             standardGeneric("knn.prediction")
           })

setMethod(f = "knn.prediction",
          signature = "myknn",
          definition = function(theObject, train.data.set, test.data.set)
          {
            train.data <- subset(train.data.set, select = - pop)
            train.class <- train.data.set$pop # "pop" is a target class
            
            test.data <- subset(test.data.set, select = - pop)
            test.class <- test.data.set$pop # "pop" is a target class
            
            prediction <- knn(train.data, test.data, train.class, k = theObject@k, prob = F)
            
            confusion.matrix <- table(prediction, test.class)
            error <- 100 - sum(diag(confusion.matrix)) / sum(confusion.matrix) * 100
            
            theObject@confusionMatrix.test <- confusion.matrix
            theObject <- setTestError(theObject, error)
            theObject <- setPrediction(theObject, prediction)
            theObject@Precision.test <- diag(theObject@confusionMatrix.test) / rowSums(theObject@confusionMatrix.test)
            theObject@Recall.test <- diag(theObject@confusionMatrix.test) / colSums(theObject@confusionMatrix.test)
            
            .auc <- multiclass.roc(test.class, as.numeric(prediction))
            theObject@auc <- .auc$auc[1]
            
            return(theObject)
          }
)

# --- Naive Bayes --- #

setGeneric(name = "nb.prediction", 
           def = function(theObject, train.data.set, test.data.set)
           {
             standardGeneric("nb.prediction")
           })

setMethod(f = "nb.prediction",
          signature = "nb",
          definition = function(theObject, train.data.set, test.data.set)
          {
            train.data <- subset(train.data.set, select = - pop)
            train.class <- train.data.set$pop # "pop" is a target class
            
            test.data <- subset(test.data.set, select = - pop)
            test.class <- test.data.set$pop # "pop" is a target class
            
            model <- naiveBayes(pop ~ ., data = train.data.set)
            # --- Training set --- #
            prediction <- predict(model, train.data)
            confusion.matrix <- table(train.class, prediction)
            error <- 100 - sum(diag(confusion.matrix)) / sum(confusion.matrix) * 100
            
            theObject@confusionMatrix.train <- confusion.matrix
            theObject <- setTrainError(theObject, error)
            
            # --- Test set --- #
            prediction <- predict(model, newdata = test.data)
            confusion.matrix <- table(test.class, prediction)
            error <- 100 - sum(diag(confusion.matrix)) / sum(confusion.matrix) * 100
            
            theObject@confusionMatrix.test <- confusion.matrix
            theObject <- setTestError(theObject, error)
            theObject <- setPrediction(theObject, prediction)
            
            .auc <- multiclass.roc(test.class, as.numeric(prediction))
            theObject@auc <- .auc$auc[1]
            
            theObject@Precision.train <- diag(theObject@confusionMatrix.train) / rowSums(theObject@confusionMatrix.train)
            theObject@Recall.train <- diag(theObject@confusionMatrix.train) / colSums(theObject@confusionMatrix.train)
            theObject@Precision.test <- diag(theObject@confusionMatrix.test) / rowSums(theObject@confusionMatrix.test)
            theObject@Recall.test <- diag(theObject@confusionMatrix.test) / colSums(theObject@confusionMatrix.test)
            
            return(theObject)
          }
)

# --- SVM --- #
setGeneric(name = "svm.estimation", 
           def = function(theObject, train.data.set, which.kernel = "RBF")
           {
             standardGeneric("svm.estimation")
           })

setMethod(f = "svm.estimation",
          signature = "mysvm",
          definition = function(theObject, train.data.set, which.kernel = "RBF")
          {
            costs <- 10^(-2:1)
            sigmas <- 1 /(2 * c(0.005, 0.125, 0.5, 1, 4))
            
            cross.errors <- matrix(nrow = length(costs), ncol = 2)
            colnames(cross.errors) <- c("Error", "nSV")
            rownames(cross.errors) <- costs
            
            cross.errors.radial <- matrix(nrow = length(sigmas), ncol = length(costs))
            colnames(cross.errors.radial) <- costs
            rownames(cross.errors.radial) <- sigmas
            
            nSV.radial <- matrix(nrow = length(sigmas), ncol = length(costs))
            colnames(nSV.radial) <- costs
            rownames(nSV.radial) <- sigmas
            
            kCV <- nrow(train.data.set) # the number of samples for cross-validation, nrow(train.data.set) means LOOCV
            for (myC in costs)
              {
                error <- 0
                switch(which.kernel,
                       
                       linear = {
                         model <-  ksvm(pop ~ .,data = train.data.set, kernel = 'polydot', C = myC, cross = kCV)
                         cross.errors[toString(myC), "Error"] <- cross(model) * 100
                         cross.errors[toString(myC), "nSV"] <- nSV(model) * 100
                         },
                       quadratic = {
                         for (g in sigmas)
                         {
                            quad <- polydot(degree = 2, scale = 1, offset = 1)
                            model <-  ksvm(pop ~ .,data = train.data.set, kernel = quad, sigma = g, C = myC, cross = kCV)
                            cross.errors.radial[toString(g), toString(myC)] <- cross(model) * 100
                            nSV.radial[toString(g), toString(myC)] <- nSV(model)
                         }
                        },
                       RBF = {
                         for (g in sigmas)
                         {
                           model <- ksvm(pop ~ .,data = train.data.set, C = myC, sigma = g, cross = kCV) 
                           cross.errors.radial[toString(g), toString(myC)] <- cross(model) * 100
                           nSV.radial[toString(g), toString(myC)] <- nSV(model)
                         }
                        },
                       stop("Enter one of 'linear', 'quad', 'RBF'"))
              }
            if(which.kernel == "linear")
              theObject@cross.errors <- cross.errors
            else 
              {
                theObject@cross.errors <- cross.errors.radial
                theObject@nSV <- nSV.radial
              }
            
            return(theObject)
          }
)

setGeneric(name = "svm.prediction", 
           def = function(theObject, train.data.set, test.data.set, which.kernel)
           {
             standardGeneric("svm.prediction")
           })

setMethod(f = "svm.prediction",
          signature = "mysvm",
          definition = function(theObject, train.data.set, test.data.set, which.kernel)
          {
            train.data <- subset(train.data.set, select = - pop)
            train.class <- train.data.set$pop # "pop" is a target class
            
            test.data <- subset(test.data.set, select = - pop)
            test.class <- test.data.set$pop # "pop" is a target class
            
            switch(which.kernel,
                  linear = {model <- svm(train.data, train.class, type = "C-classification", cost = theObject@cost, kernel = "linear", scale = FALSE)},
                  quadratic = {model <- svm(train.data, train.class, type="C-classification", cost = theObject@cost, kernel="polynomial", degree = 2, coef0 = 1, scale = FALSE)},
                  RBF = {model <- svm(train.data, train.class, type = "C-classification", gamma = theObject@sigma, cost = theObject@cost, kernel = "radial", scale = FALSE)},
                  stop("Enter one of 'linear', 'quadratic', 'RBF'")
            )
            
            # --- Training set --- #
            prediction <- predict(model, train.data)
            confusion.matrix <- table(train.class, prediction)
            error <- 100 - sum(diag(confusion.matrix)) / sum(confusion.matrix) * 100
            
            theObject@confusionMatrix.train <- confusion.matrix
            theObject <- setTrainError(theObject, error)
            
            # --- Test set --- #
            prediction <- predict(model, newdata = test.data)
            confusion.matrix <- table(test.class, prediction)
            error <- 100 - sum(diag(confusion.matrix)) / sum(confusion.matrix) * 100
            
            # theObject@model <- model
            theObject@confusionMatrix.test <- confusion.matrix
            theObject <- setTestError(theObject, error)
            theObject <- setPrediction(theObject, prediction)
            
            theObject@Precision.train <- diag(theObject@confusionMatrix.train) / rowSums(theObject@confusionMatrix.train)
            theObject@Recall.train <- diag(theObject@confusionMatrix.train) / colSums(theObject@confusionMatrix.train)
            theObject@Precision.test <- diag(theObject@confusionMatrix.test) / rowSums(theObject@confusionMatrix.test)
            theObject@Recall.test <- diag(theObject@confusionMatrix.test) / colSums(theObject@confusionMatrix.test)
            
            .auc <- multiclass.roc(test.class, as.numeric(prediction))
            theObject@auc <- .auc$auc[1]
            
            return(theObject)
          }
)

setGeneric(name = "dt.prediction", 
           def = function(theObject, train.data.set, test.data.set)
           {
             standardGeneric("dt.prediction")
           })

setMethod(f = "dt.prediction",
          signature = "dt",
          definition = function(theObject, train.data.set, test.data.set)
          {
            train.data <- subset(train.data.set, select = - pop)
            train.class <- train.data.set$pop # "pop" is a target class
            
            test.data <- subset(test.data.set, select = - pop)
            test.class <- test.data.set$pop # "pop" is a target class
            
            model = rpart(pop ~ ., data = train.data.set, parms = list(split='gini'), 
                               control = rpart.control(cp = 0.001, xval = 10, maxdepth = 15))
            
            # Prunining the tree
            par(mfrow = c(1, 2))
            plot(model$cptable[, 2], model$cptable[, 3], col = "purple", type = "o", xlab = "nsplits", ylab = "error")
            lines(model$cptable[, 2], model$cptable[, 4], col = "blue", type = "o")
            grid(20, 20, lwd = 1)
            legend("topleft", c("R(T) training", "R(T) cv"), col = c("purple", "blue"), lty = 1, cex = 0.8, bty = "n", text.width = 6, seg.len = 0.5)
            
            plotcp(model)
            
            theObject@plot.cptable <- recordPlot()       
            theObject@model.cptable <- model$cptable
            
            alfa <- model$cptable[which.min(model$cptable[, 4]), 1]
            model.tree.prune <- prune(model, cp = alfa)
            
            # theObject@model.tree <- model.tree.prune
            
            par(mfrow = c(1,1))
            fancyRpartPlot(model.tree.prune, cex = 0.9)
            theObject@plot.tree <- recordPlot() 
            
            # --- Training data --- #
            set.seed(8934)
            
            pred <- predict(model.tree.prune, newdata = train.data)
            prediction <- NULL
            prediction[pred[, 1] >= 0.5] = "CEU"
            prediction[pred[, 2] >= 0.5] = "HCB"
            prediction[pred[, 3] >= 0.5] = "JPT"
            prediction[pred[, 4] >= 0.5] = "YRI"
            
            theObject@confusionMatrix.train <- table(train.class, prediction)
            theObject@train.error <- 100 - sum(diag(theObject@confusionMatrix.train)) / sum(theObject@confusionMatrix.train) * 100
            
            # --- Test data --- #
            pred <- predict(model.tree.prune, newdata = test.data)
            prediction <- NULL
            prediction[pred[, 1] >= 0.5] = "CEU"
            prediction[pred[, 2] >= 0.5] = "HCB"
            prediction[pred[, 3] >= 0.5] = "JPT"
            prediction[pred[, 4] >= 0.5] = "YRI"
            
            theObject@confusionMatrix.test <- table(test.class, prediction)
            theObject@test.error <- 100 - sum(diag(theObject@confusionMatrix.test)) / sum(theObject@confusionMatrix.test) * 100
            
            theObject@Precision.train <- diag(theObject@confusionMatrix.train) / rowSums(theObject@confusionMatrix.train)
            theObject@Recall.train <- diag(theObject@confusionMatrix.train) / colSums(theObject@confusionMatrix.train)
            theObject@Precision.test <- diag(theObject@confusionMatrix.test) / rowSums(theObject@confusionMatrix.test)
            theObject@Recall.test <- diag(theObject@confusionMatrix.test) / colSums(theObject@confusionMatrix.test)
            
            return(theObject)
          }
)

setGeneric(name = "rf.prediction", 
           def = function(theObject, train.data.set, test.data.set)
           {
             standardGeneric("rf.prediction")
           })

setMethod(f = "rf.prediction",
          signature = "rf",
          definition = function(theObject, train.data.set, test.data.set)
          {
            train.data <- subset(train.data.set, select = - pop)
            train.class <- train.data.set$pop # "pop" is a target class
            
            test.data <- subset(test.data.set, select = - pop)
            test.class <- test.data.set$pop # "pop" is a target class
            
            model <- randomForest(pop ~ ., train.data.set, xtest = test.data, ytest = test.class)
            print(model)
            
            # --- Training data --- #
            theObject@confusionMatrix.train <- as.table(model$confusion)
            theObject@train.error <- 100 - sum(diag(theObject@confusionMatrix.train)) / sum(theObject@confusionMatrix.train) * 100
            
            # --- Test data --- # 
            theObject@confusionMatrix.test <- as.table(model$test$confusion)
            theObject@predicted <- as.numeric(model$test$predicted)
            theObject@test.error <- 100 - sum(diag(theObject@confusionMatrix.test)) / sum(theObject@confusionMatrix.test) * 100
            
            .auc <- multiclass.roc(test.class, as.numeric(model$test$predicted))
            theObject@auc <- .auc$auc[1]
            
            # --- Plot --- #
            plot(model)
            theObject@plot.model <- recordPlot() 
            
            plot(importance(model), lty = 2, pch = 16)
            grid(20, 20, lty = 1, lwd = 1)
            lines(importance(model))
            text(importance(model), rownames(importance(model)), pos = 4, cex = 0.8)
            theObject@plot.importance <- recordPlot()
            
            # ROC curve
            # the positive prediction on the test samples
#             par(mfrow= c(2, 2))
#             for (i in 1:4)
#             {
#               pred <- prediction(model$test$votes[, i], test.class)
#               plot(performance(pred, 'tpr', 'fpr'), main = "ROC curve")
#               abline(0, 1, col="red")
#             }
#             theObject@plot.roc.curve <- recordPlot()
            
            theObject@Precision.train <- diag(theObject@confusionMatrix.train) / rowSums(theObject@confusionMatrix.train)
            theObject@Recall.train <- diag(theObject@confusionMatrix.train) / colSums(theObject@confusionMatrix.train)
            theObject@Precision.test <- diag(theObject@confusionMatrix.test) / rowSums(theObject@confusionMatrix.test)
            theObject@Recall.test <- diag(theObject@confusionMatrix.test) / colSums(theObject@confusionMatrix.test)
            
            return(theObject)
          }
)

setGeneric(name = "lda.qda.prediction", 
           def = function(theObject, train.data.set, test.data.set, discriminant.analysis)
           {
             standardGeneric("lda.qda.prediction")
           })

setMethod(f = "lda.qda.prediction",
          signature = "lda.qda",
          definition = function(theObject, train.data.set, test.data.set, discriminant.analysis)
          {
            train.data <- subset(train.data.set, select = - pop)
            train.class <- train.data.set$pop # "pop" is a target class
            
            test.data <- subset(test.data.set, select = - pop)
            test.class <- test.data.set$pop # "pop" is a target class
            
            switch(discriminant.analysis,
                   LDA = {model <- lda(train.data, train.class)},
                   QDA = {model <- qda(train.data, train.class)},
                   stop("Enter one of 'LDA', 'QDA'")
            )
            
            model.cv <- update(model, CV = T)
            
            # --- Training data --- #
            theObject@confusionMatrix.train <- table(train.class, model.cv$class)
            theObject@train.error <- 100 - sum(diag(theObject@confusionMatrix.train)) / sum(theObject@confusionMatrix.train) * 100
            
            # --- Test data --- #
            prediction <- predict(model, test.data)
            theObject@confusionMatrix.test <- table(test.class, prediction$class)
            theObject@test.error <- 100 - sum(diag(theObject@confusionMatrix.test)) / sum(theObject@confusionMatrix.test) * 100
            
            theObject <- setPrediction(theObject, prediction$class)
            
            if (discriminant.analysis == "LDA")
              theObject@x <- prediction$x
            
            theObject@Precision.train <- diag(theObject@confusionMatrix.train) / rowSums(theObject@confusionMatrix.train)
            theObject@Recall.train <- diag(theObject@confusionMatrix.train) / colSums(theObject@confusionMatrix.train)
            theObject@Precision.test <- diag(theObject@confusionMatrix.test) / rowSums(theObject@confusionMatrix.test)
            theObject@Recall.test <- diag(theObject@confusionMatrix.test) / colSums(theObject@confusionMatrix.test)
            
            .auc <- multiclass.roc(test.class, as.numeric(prediction$class))
            theObject@auc <- .auc$auc[1]
            
            return(theObject)
          }
)

# --- Plot the prediction --- #

plot.prediction <- function(predicted, data.set, var.x = "PC.1", var.y = "PC.2", legend = "", title.name = "")
{
  colors <- c("red", "green", "blue", "yellow")
  mycolors <- colors[as.numeric(predicted)]
  
  save.par <- par(mfrow= c(1, 1))
  pred.plot <- plot(data.set[, c(var.x, var.y)] * 2, xlab = var.x, ylab = var.y, type = "n")
  grid(20, 20, lty = 1, lwd = 1)
  points(data.set[, c(var.x, var.y)], col = mycolors, pch = "*", cex = 2)
  
  for (i in 1:4)
  {
    m1 <- mean(subset(data.set[, var.x], predicted == i))
    m2 <- mean(subset(data.set[, var.y], predicted == i))
    points(m1, m2, pch = 16, cex = 1,col = i)
  }
    
  legend("bottom", legend,  col = colors, pch = "*", lty = 2, cex = 0.75, horiz = T, bty = "n")
  title(title.name)
  
  return(recordPlot())
}

