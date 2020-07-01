#Code for the publication **Maternal blood lipidomics analyse link critical metabolic pathways associated with 
#severe preeclampsia** authored by Yu Liu, Bing He and Mano Maurya *et al*.

#The code is needed by the main program "lipidomic_main.R" to do machine learning analysis on lipidomic data.

#Created by Yu Liu on Jun 29, 2020 (uiluy@umich.edu, Department of Computational Medicine and Bioinformatics, 
#University of Michigan, MI)

lilikoi.machine_learning11 <- function (PDSmatrix = lilikoimat, 
                                        measurementLabels = lilikoilabels, 
                                        significantPathways = significantPathways, 
                                        selectmod = 'LDA', skippam = FALSE, cvnum = 10, 
                                        randomseed = 7, dividep = 0.8, dividseed = NULL){
  
  library(caret)
  library(pROC)
  library(ggplot2)
  library(gbm)
  library(PRROC)
  
  cancer_df <- data.frame(t((PDSmatrix[significantPathways, 
                                       ])), Label = measurementLabels, check.names = T)
  colnames(cancer_df)[which(names(cancer_df) == "Label")] <- "subtype"
  performance_training <- matrix(rep(0, len = 56), nrow = 8)
  performance_testing <- matrix(rep(0, len = 56), nrow = 8)
  performance <- matrix(rep(0, len = 56), nrow = 8)
  models <- list()
  
  if(is.null(dividseed)){
    set.seed(randomseed + 1)
  }else{
    set.seed(dividseed)
  }
  
  trainIndex <- createDataPartition(cancer_df$subtype, p = dividep, 
                                    list = FALSE, times = 1)
  trainDf <- cancer_df[trainIndex, ]
  testDf <- cancer_df[-trainIndex, ]
  
  library(MLmetrics)
  
  MySummary  <- function(data, lev = NULL, model = NULL){
    a1 <- defaultSummary(data, lev, model)
    b1 <- twoClassSummary(data, lev, model)
    c1 <- prSummary(data, lev, model)
    d1 <- multiClassSummary(data, lev, model)
    out <- c(a1, b1, c1, d1)
    
    return(out)
  }
  
  control <- trainControl(method = "cv", number = cvnum, classProbs = TRUE, 
                          summaryFunction = MySummary)
  
  
  set.seed(randomseed)
  garbage <- utils::capture.output(fit.cart <- train(subtype ~ 
                                                       ., data = trainDf, method = "rpart", trControl = control, 
                                                     metric = "ROC"))
  models[[1]] <- fit.cart
  
  writemetrics <- function(fitres = fit.cart, blank = performance, colidx = 1){
    
    blank[1, colidx] <- max(fitres$results$ROC)
    blank[2, colidx] <- fitres$results$Sens[which.max(fitres$results$ROC)]
    blank[3, colidx] <- fitres$results$Spec[which.max(fitres$results$ROC)]
    blank[4, colidx] <- fitres$results$Accuracy[which.max(fitres$results$ROC)]
    blank[5, colidx] <- fitres$results$Precision[which.max(fitres$results$ROC)]
    blank[6, colidx] <- fitres$results$Recall[which.max(fitres$results$ROC)]
    blank[7, colidx] <- fitres$results$F1[which.max(fitres$results$ROC)]
    blank[8, colidx] <- fitres$results$Balanced_Accuracy[which.max(fitres$results$ROC)]
    
    return(blank)
    
  }
  
  performance <- writemetrics(fitres = fit.cart, blank = performance, colidx = 1)
  
  predictdata <- function(fitres = fit.cart, blank = performance_training, colidx = 1, 
                          dat = trainDf){
    
    cartClasses <- stats::predict(fitres, newdata = dat, type = "prob")
    cartClasses1 <- stats::predict(fitres, newdata = dat)
    cartConfusion <- confusionMatrix(data = cartClasses1, dat$subtype)
    cart.ROC <- roc(predictor = cartClasses$Normal, response = dat$subtype, 
                    levels = rev(levels(dat$subtype)))
    blank[1, colidx] <- as.numeric(cart.ROC$auc)
    blank[2, colidx] <- cartConfusion$byClass[1]
    blank[3, colidx] <- cartConfusion$byClass[2]
    blank[4, colidx] <- cartConfusion$overall[1]
    blank[5, colidx] <- cartConfusion$byClass[5]
    blank[6, colidx] <- cartConfusion$byClass[6]
    blank[7, colidx] <- cartConfusion$byClass[7]
    blank[8, colidx] <- cartConfusion$byClass[11]
    
    return(blank)
    
    
  }
  
  performance_training <- predictdata(fitres = fit.cart, blank = performance_training, colidx = 1, 
                                      dat = trainDf)
  
  performance_testing <- predictdata(fitres = fit.cart, blank = performance_testing, colidx = 1, 
                                     dat = testDf)
  
  
  
  set.seed(randomseed)
  garbage <- suppressWarnings(utils::capture.output(fit.lda <- train(subtype ~ 
                                                                       ., data = trainDf, method = "lda", trControl = control, 
                                                                     metric = "ROC", trace = F)))
  models[[2]] <- fit.lda
  
  performance <- writemetrics(fitres = fit.lda, blank = performance, colidx = 2)
  
  performance_training <- predictdata(fitres = fit.lda, blank = performance_training, colidx = 2, 
                                      dat = trainDf)
  
  performance_testing <- predictdata(fitres = fit.lda, blank = performance_testing, colidx = 2, 
                                     dat = testDf)
  
  
  set.seed(randomseed)
  garbage <- utils::capture.output(fit.svm <- train(subtype ~ 
                                                      ., data = trainDf, method = "svmRadial", trControl = control, 
                                                    metric = "ROC"))
  models[[3]] <- fit.svm
  performance <- writemetrics(fitres = fit.svm, blank = performance, colidx = 3)
  
  performance_training <- predictdata(fitres = fit.svm, blank = performance_training, colidx = 3, 
                                      dat = trainDf)
  
  performance_testing <- predictdata(fitres = fit.svm, blank = performance_testing, colidx = 3, 
                                     dat = testDf)
  
  
  
  set.seed(randomseed)
  garbage <- utils::capture.output(fit.rf <- train(subtype ~ 
                                                     ., data = trainDf, method = "rf", trControl = control, 
                                                   metric = "ROC"))
  models[[4]] <- fit.rf
  
  performance <- writemetrics(fitres = fit.rf, blank = performance, colidx = 4)
  
  performance_training <- predictdata(fitres = fit.rf, blank = performance_training, colidx = 4, 
                                      dat = trainDf)
  
  performance_testing <- predictdata(fitres = fit.rf, blank = performance_testing, colidx = 4, 
                                     dat = testDf)
  
  
  set.seed(randomseed)
  garbage <- suppressWarnings(utils::capture.output(fit.gbm <- train(subtype ~ 
                                                                       ., data = trainDf, method = "gbm", trControl = control, 
                                                                     metric = "ROC")))
  models[[5]] <- fit.gbm
  
  performance <- writemetrics(fitres = fit.gbm, blank = performance, colidx = 5)
  
  performance_training <- predictdata(fitres = fit.gbm, blank = performance_training, colidx = 5, 
                                      dat = trainDf)
  
  performance_testing <- predictdata(fitres = fit.gbm, blank = performance_testing, colidx = 5, 
                                     dat = testDf)
  
  
  if(skippam == FALSE){
    
    set.seed(randomseed)
    garbage <- utils::capture.output(fit.pam <- train(subtype ~ 
                                                        ., data = trainDf, method = "pam", trControl = control, 
                                                      metric = "ROC"))
    models[[6]] <- fit.pam
    
    performance <- writemetrics(fitres = fit.pam, blank = performance, colidx = 6)
    
    performance_training <- predictdata(fitres = fit.pam, blank = performance_training, colidx = 6, 
                                        dat = trainDf)
    
    performance_testing <- predictdata(fitres = fit.pam, blank = performance_testing, colidx = 6, 
                                       dat = testDf)
    
  }
  
  set.seed(randomseed)
  garbage <- suppressWarnings(utils::capture.output(fit.log <- train(subtype ~ 
                                                                       ., data = trainDf, method = "glmnet", trControl = control, 
                                                                     metric = "ROC")))
  models[[7]] <- fit.log
  
  performance <- writemetrics(fitres = fit.log, blank = performance, colidx = 7)
  
  performance_training <- predictdata(fitres = fit.log, blank = performance_training, colidx = 7, 
                                      dat = trainDf)
  
  performance_testing <- predictdata(fitres = fit.log, blank = performance_testing, colidx = 7, 
                                     dat = testDf)
  
  
  p5 <- tryCatch({
    graphics::plot(graphics::plot(varImp(fit.cart, scale = FALSE, 
                                         top = 20), main = "RPART"))
  }, error = function(err){
    NULL
  })
  
  p6 <- graphics::plot(graphics::plot(varImp(fit.lda, scale = FALSE, 
                                             top = 20), main = "LDA"))
  p7 <- graphics::plot(graphics::plot(varImp(fit.svm, scale = FALSE, 
                                             top = 20), main = "SVM"))
  p8 <- graphics::plot(graphics::plot(varImp(fit.rf, scale = FALSE, 
                                             top = 20), main = "RF"))
  p9 <- graphics::plot(graphics::plot(varImp(fit.gbm, scale = FALSE, 
                                             top = 20), main = "GBM"))
  p10 <- graphics::plot(graphics::plot(varImp(fit.pam, scale = FALSE, 
                                              top = 20), main = "PAM"))
  p11 <- graphics::plot(graphics::plot(varImp(fit.log, scale = FALSE, 
                                              top = 20), main = "LOG"))
  
  
  
  getROC <- function(fitres = fit.cart, dat = testDf){
    
    resClasses <- stats::predict(fitres, newdata = dat, type = "prob")
    res.ROC <- roc(predictor = resClasses$Normal, response = dat$subtype, 
                    levels = rev(levels(dat$subtype)))
    
    return(res.ROC)
    
  }
  
  cart.ROC <- getROC(fitres = fit.cart, dat = testDf)
  lda.ROC <- getROC(fitres = fit.lda, dat = testDf)
  svm.ROC <- getROC(fitres = fit.svm, dat = testDf)
  rf.ROC <- getROC(fitres = fit.rf, dat = testDf)
  gbm.ROC <- getROC(fitres = fit.gbm, dat = testDf)
  
  if(skippam == FALSE){
    pam.ROC <- getROC(fitres = fit.pam, dat = testDf)
  }
  
  log.ROC <- getROC(fitres = fit.log, dat = testDf)
  
  
  
  smooth_method <- "binormal"
  graphics::plot(cart.ROC, col = "red", cex.lab = 1.5)
  graphics::par(new = TRUE)
  graphics::plot(lda.ROC, col = "green", cex.lab = 1.5)
  graphics::par(new = TRUE)
  graphics::plot(svm.ROC, col = "black", cex.lab = 1.5)
  graphics::par(new = TRUE)
  graphics::plot(rf.ROC, col = "orange", cex.lab = 1.5)
  graphics::par(new = TRUE)
  graphics::plot(gbm.ROC, col = "blue", cex.lab = 1.5)
  graphics::par(new = TRUE)
  
  if(skippam == FALSE){
    graphics::plot(pam.ROC, col = "hotpink", cex.lab = 1.5)
    graphics::par(new = TRUE)
  }
  
  graphics::plot(log.ROC, col = "lightgoldenrod2", main = "Testing ROC", 
                 cex.lab = 1.5)
  graphics::legend(0.2, 0.4, legend = c("RPART", "LDA", "SVM", 
                                        "RF", "GBM", "PAM", "LOG"), col = c("red", "green", "black", 
                                                                            "orange", "blue", "hotpink", "lightgoldenrod2"), lty = 1:2, 
                   cex = 1)
  
  testingroclist <- list()
  testingroclist$RPART <- cart.ROC
  testingroclist$LDA <- lda.ROC
  testingroclist$SVM <- svm.ROC
  testingroclist$RF <- rf.ROC
  testingroclist$GBM <- gbm.ROC
  testingroclist$PAM <- pam.ROC
  testingroclist$LLOG <- log.ROC
  
  selectroc <- testingroclist[[selectmod]]
  
  
  
  getPRROC <- function(fitres = fit.cart, dat = testDf){
    
    resClasses <- stats::predict(fitres, newdata = dat, type = 'prob')
    press <- resClasses$Cancer
    
    trues <- as.character(dat$subtype)
    truess <- rep(0, length(trues))
    truess[trues == 'Cancer'] <- 1
    
    probj <- pr.curve(scores.class0 = press, weights.class0 = truess, curve = TRUE)
    prauc <- probj$auc.integral
    
    
    plot.scores.AUC <- function(y = truess, y.hat = press, measure = "ppv", x.measure = "tpr", 
                                auc = prauc) {
      # plot ROC curve
      library(ROCR)
      pr <- prediction(y.hat, y)
      prf <- performance(pr, measure = measure, x.measure = x.measure)
      
      plot(prf, main = paste0("Curve (AUC-PR: ", round(auc, 2), ")"))
      
      return(prf)
      
    }
    
    prf <- plot.scores.AUC()
    
    return(list(probj = probj, prplotobj = prf))

    
  }
  
  
  cart.PRROC <- tryCatch({
    getPRROC(fitres = fit.cart, dat = testDf)
  }, error = function(err){
    NULL
  })
  
  lda.PRROC <- getPRROC(fitres = fit.lda, dat = testDf)
  
  svm.PRROC <- tryCatch({
    getPRROC(fitres = fit.svm, dat = testDf)
  }, error = function(err){
    NULL
  })
  
  rf.PRROC <- getPRROC(fitres = fit.rf, dat = testDf)
  
  gbm.PRROC <- tryCatch({
    getPRROC(fitres = fit.gbm, dat = testDf)
  }, error = function(err){
    NULL
  })
   
  if(skippam == FALSE){
    pam.PRROC <- getPRROC(fitres = fit.pam, dat = testDf)
  }
  
  log.PRROC <- getPRROC(fitres = fit.log, dat = testDf)
  
  
  
  testingprroclist <- list()
  testingprroclist$RPART <- cart.PRROC
  testingprroclist$LDA <- lda.PRROC
  testingprroclist$SVM <- svm.PRROC
  testingprroclist$RF <- rf.PRROC
  testingprroclist$GBM <- gbm.PRROC
  testingprroclist$PAM <- pam.PRROC
  testingprroclist$LLOG <- log.PRROC
  
  selectprroc <- testingprroclist[[selectmod]]
  
  
  
  
  makeperformance <- function(mat = performance, datname = 'performance'){
    
    performance_list <- list()
    performance_list[[1]] <- mat
    list_performance <- performance_list
    
    AUC_performance <- lapply(list_performance, function(x) x[1, ])
    
    SENS_performance <- lapply(list_performance, function(x) x[2, ])
    
    SPEC_performance <- lapply(list_performance, function(x) x[3, ])

    F1_performance <- lapply(list_performance, function(x) x[7, ])

    Balanced_accuracy_performance <- lapply(list_performance, function(x) x[8, ])
    
    outputauc <- do.call(rbind, lapply(AUC_performance, matrix, ncol = 7, 
                                       byrow = TRUE))

    outputsens <- do.call(rbind, lapply(SENS_performance, matrix, ncol = 7, 
                                        byrow = TRUE))
    
    outputspec <- do.call(rbind, lapply(SPEC_performance, matrix, ncol = 7, 
                                        byrow = TRUE))

    outputf1 <- do.call(rbind, lapply(F1_performance, matrix, ncol = 7, 
                                      byrow = TRUE))

    outputacc <- do.call(rbind, lapply(Balanced_accuracy_performance, 
                                       matrix, ncol = 7, byrow = TRUE))

    
    AUC_mean <- apply(outputauc, 2, mean)
    SENS_mean <- apply(outputsens, 2, mean)
    SPEC_mean <- apply(outputspec, 2, mean)
    F1_mean <- apply(outputf1, 2, mean)
    Balanced_accuracy_mean <- apply(outputacc, 2, mean)
    
    performance_only <- t(t(rep(datname, 7)))
    
    performance_data <- data.frame(AUC = data.frame(AUC = t((t(AUC_mean)))), 
                                   SENS = data.frame(SENS = t((t(SENS_mean)))), 
                                   SPEC = data.frame(SPEC = t((t(SPEC_mean)))), 
                                   F1 = data.frame(F1 = t((t(F1_mean)))), 
                                   Balanced_accuracy = 
                                     data.frame(Balanced_accuracy = t((t(Balanced_accuracy_mean)))), 
                                   datname = performance_only, 
                                   Algorithm = (rep(t(c("RPART", "LDA", "SVM", "RF", "GBM", "PAM", "LOG")), 1)))
    
    return(performance_data)
  }
  
  performance_data <- makeperformance(mat = performance, datname = 'performance')
  performance_data_train <- makeperformance(mat = performance_training, datname = 'training')
  performance_data_test <- makeperformance(mat = performance_testing, datname = 'testing')
  
  meltdata <- function(performance_data = performance_data, algname = NULL){
    
    library(reshape)
    
    performance_data <- performance_data[order(-performance_data$F1),]
    modorder <- as.character(performance_data$Algorithm)
    melted_performance_data <- suppressMessages(melt(performance_data))

    Algorithm <- NULL
    value <- NULL
    variable <- NULL
    textLabels <- geom_text(aes(x = Algorithm, label = round(value, 2)), 
                            position = position_dodge(width = 1), 
                            vjust = -0.5, size = 2)
    
    melted_performance_data$Algorithm <- factor(melted_performance_data$Algorithm, levels = modorder, 
                                                ordered = TRUE)
    
    if(is.null(algname)){
      algname <- unique(performance_data$datname)
    }
    
    p <- ggplot(data = melted_performance_data, aes(x = Algorithm, y = value, fill = variable)) + 
      geom_bar(stat = "identity", position = position_dodge()) + 
      xlab("") + ylab("") + ggtitle(algname) + 
      theme(plot.title = element_text(hjust = 0.5), 
            axis.text = element_text(size = 15, face = "bold"), 
            axis.title = element_text(size = 14, face = "bold")) + 
      labs(fill = "") + textLabels
    print(p)
    
    return(p)
    
  }
  
  p_performance <- meltdata(performance_data = performance_data)
  
  p_training <- meltdata(performance_data = performance_data_train)
  
  p_testing <- meltdata(performance_data = performance_data_test)
  
  
  best_model <- models[which.max(performance_testing[1, ])]
  method <- (unlist(best_model)[[1]])
  if (method == "glmnet") {
    method <- "log"
  }
  if (method == "svmRadial") {
    method <- "svm"
  }
  
  dd <- subset(performance_data_test, Algorithm == toupper(method))
  p_best <- meltdata(performance_data = dd)
  
  
  mlResults <- list()
  mlResults$models <- models
  
  mlResults$performance <- performance_data
  mlResults$performance_training <- performance_data_train
  mlResults$performance_test <- performance_data_test
  
  
  mlResults$train_inx <- trainIndex
  plots <- list()
  plots$mlResults <- mlResults
  
  plots[[paste0(selectmod, '_testing_ROC')]] <- selectroc
  plots[[paste0(selectmod, '_tesing_PRROC')]] <- selectprroc
  
  plots$p_performance <- p_performance
  plots$p_training <- p_training
  plots$p_testing <- p_testing
  plots$p_best <- p_best
  
  return(plots)
}
