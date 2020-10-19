#Code for the publication **Maternal blood lipidomics analyse link critical metabolic pathways associated with 
#severe preeclampsia** authored by Yu Liu, Bing He and Mano Maurya *et al*.

#The code is needed by the main program "lipidomic_main.R" to do machine learning analysis on lipidomic data.

#Created by Yu Liu on Jun 29, 2020 (uiluy@umich.edu, Department of Computational Medicine and Bioinformatics, 
#University of Michigan, MI)

lilikoi.machine_learning11 <- function (PDSmatrix = lilikoimat, 
                                        measurementLabels = lilikoilabels, 
                                        significantPathways = significantPathways, 
                                        selectmod = 'LDA', skippam = FALSE, cvnum = 10, times = 10,
                                        randomseed = 2000, dividep = 0.8, dividseed = 1998){
  
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
                                    list = FALSE, times = times)
  performance_data_final<-NULL
  performance_data_train_final<-NULL
  performance_data_test_final<-NULL
  for (i in 1:times){
  trainIndexs<-trainIndex[,i]
  trainDf <- cancer_df[trainIndexs, ]
  testDf <- cancer_df[-trainIndexs, ]

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
    cartConfusion <- confusionMatrix(data = cartClasses1, as.factor(dat$subtype))
    cart.ROC <- roc(predictor = cartClasses$Normal, response =  as.factor(dat$subtype), 
                    levels = rev(levels( as.factor(dat$subtype))))
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
  
  # p6 <- graphics::plot(graphics::plot(varImp(fit.lda, scale = FALSE, 
  #                                            top = 20), main = "LDA"))
  # p7 <- graphics::plot(graphics::plot(varImp(fit.svm, scale = FALSE, 
  #                                            top = 20), main = "SVM"))
  # p8 <- graphics::plot(graphics::plot(varImp(fit.rf, scale = FALSE, 
  #                                            top = 20), main = "RF"))
  # p9 <- graphics::plot(graphics::plot(varImp(fit.gbm, scale = FALSE, 
  #                                            top = 20), main = "GBM"))
  # p10 <- graphics::plot(graphics::plot(varImp(fit.pam, scale = FALSE, 
  #                                             top = 20), main = "PAM"))
  # p11 <- graphics::plot(graphics::plot(varImp(fit.log, scale = FALSE, 
  #                                             top = 20), main = "LOG"))
  
  
  
  getROC <- function(fitres = fit.cart, dat = testDf){
    
    resClasses <- stats::predict(fitres, newdata = dat, type = "prob")
    res.ROC <- roc(predictor = resClasses$Normal, response = as.factor(dat$subtype), 
                   levels = rev(levels(as.factor(dat$subtype))))
    
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
  
  
  
  
  makeperformance <- function(mat = performance, datname = 'Performance'){
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
  
  performance_data <- makeperformance(mat = performance, datname = 'Performance')
  performance_data_train <- makeperformance(mat = performance_training, datname = 'Training')
  performance_data_test <- makeperformance(mat = performance_testing, datname = 'Testing')
  
  performance_data$Repeats<-i
  performance_data_train$Repeats<-i
  performance_data_test$Repeats<-i
  
  performance_data_final<-rbind(performance_data_final,performance_data)
  performance_data_train_final<-rbind(performance_data_train_final,performance_data_train)
  performance_data_test_final<-rbind(performance_data_test_final,performance_data_test)
  }  
  performance_data_final[is.na(performance_data_final)]<-0
  performance_data_train_final[is.na(performance_data_train_final)]<-0
  performance_data_test_final[is.na(performance_data_test_final)]<-0
  performance_mean<-function(performance_data_final=performance_data){
    AUC_mean_final<-tapply(performance_data_final$AUC,performance_data_final$Algorithm,mean)
    SENS_mean_final<-tapply(performance_data_final$SENS,performance_data_final$Algorithm,mean)
    SPEC_mean_final<-tapply(performance_data_final$SPEC,performance_data_final$Algorithm,mean)
    F1_mean_final<-tapply(performance_data_final$F1,performance_data_final$Algorithm,mean)
    Balanced_accuracy_mean_final<-tapply(performance_data_final$Balanced_accuracy,performance_data_final$Algorithm,mean)
    performance_data_final_mean <- data.frame(AUC = AUC_mean_final, 
                                   F1 = F1_mean_final, 
                                   Balanced_accuracy = Balanced_accuracy_mean_final, 
                                   datname = unique(performance_data_final$datname), 
                                   Algorithm = names(AUC_mean_final))
    return(performance_data_final_mean)
  }
  performance_selection<-function(performance_data_final=performance_data){
    performance_data_final<-performance_data_final[which(rowSums(performance_data_final==0)==0),]
    AUC_final<-tapply(performance_data_final$AUC,performance_data_final$Repeats,sum)
    SENS_final<-tapply(performance_data_final$SENS,performance_data_final$Repeats,sum)
    SPEC_final<-tapply(performance_data_final$SPEC,performance_data_final$Repeats,sum)
    F1_mean_final<-tapply(performance_data_final$F1,performance_data_final$Repeats,sum)
    Balanced_accuracy_mean_final<-tapply(performance_data_final$Balanced_accuracy,performance_data_final$Repeats,sum)
    performance_scores <- cbind(AUC_final,F1_mean_final,Balanced_accuracy_mean_final)
    performance_score <-apply(performance_scores, 1, sum)
    performance_data_selected<-unique(subset(performance_data_final, Repeats==which.max(performance_score)))
    new.list<-colnames(performance_data_selected)
    new.list<-setdiff(new.list,c("SENS","SPEC","Repeats"))
    performance_data_selected<-performance_data_selected[,new.list]
    return(performance_data_selected)
  }
  if(mean(performance_data_test_final$F1)<0.5){
  performance_data <- performance_mean(performance_data_final)
  performance_data_train <- performance_mean(performance_data_train_final)
  performance_data_test <- performance_mean(performance_data_test_final)
  }else{
    performance_data <- performance_selection(performance_data_final)
    performance_data_train <- performance_selection(performance_data_train_final)
    performance_data_test <- performance_selection(performance_data_test_final)
  }
  performance_sd<-function(performance_data_final=performance_data){
    library(Rmisc)
    AUC_sd <- summarySE(performance_data_final, measurevar = "AUC",groupvars = "Algorithm")$sd
    SENS_sd <- summarySE(performance_data_final, measurevar = "SENS",groupvars = "Algorithm")$sd
    SPEC_sd <- summarySE(performance_data_final, measurevar = "SPEC",groupvars = "Algorithm")$sd
    F1_sd <- summarySE(performance_data_final, measurevar = "F1",groupvars = "Algorithm")$sd
    Balanced_accuracy_sd <- summarySE(performance_data_final, measurevar = "Balanced_accuracy",groupvars = "Algorithm")$sd
    performance_data_final_sd <- data.frame(AUC = AUC_sd, 
                                              F1 = F1_sd, 
                                              Balanced_accuracy = Balanced_accuracy_sd, 
                                              datname = unique(performance_data_final$datname), 
                                              Algorithm = summarySE(performance_data_final, measurevar = "AUC",groupvars = "Algorithm")$Algorithm)
    performance_data_final_sd[is.na(performance_data_final_sd)]<-0
    return(performance_data_final_sd)
  }
  performance_data_sd <- performance_sd(performance_data_final)
  performance_data_train_sd <- performance_sd(performance_data_train_final)
  performance_data_test_sd <- performance_sd(performance_data_test_final) 
  
  meltdata <- function(performance_data = performance_data, performance_data_sd = performance_data_sd, algname = NULL){
    rownames(performance_data_sd) <-performance_data_sd$Algorithm
    performance_data_sd<-performance_data_sd[performance_data$Algorithm,]
    library(reshape)
    
    performance_data <- performance_data[order(-performance_data$F1),]
    modorder <- as.character(performance_data$Algorithm)
    melted_performance_data <- suppressMessages(melt(performance_data))
    melted_performance_data_sd <- suppressMessages(melt(performance_data_sd))
    temp<-1-melted_performance_data$value
    temp.table<-data.frame(temp1=temp,temp2=melted_performance_data_sd$value)
    final.sd<-apply(temp.table,1,min)
    melted_performance_data$SD <- final.sd
    
    
    Algorithm <- NULL
    value <- NULL
    variable <- NULL
    textLabels <- geom_text(aes(x = Algorithm, label = round(value, 2)), 
                            position = position_dodge(width = 0.9), 
                            vjust = -0.5, size = 3)
    
    melted_performance_data$Algorithm <- factor(melted_performance_data$Algorithm, levels = modorder, 
                                                ordered = TRUE)
    
    if(is.null(algname)){
      algname <- unique(performance_data$datname)
    }
    
    p <- ggplot(data = melted_performance_data, aes(x = Algorithm, y = value, fill = variable)) + 
      geom_bar(stat = "identity", position = position_dodge()) + 
      geom_errorbar(aes(ymin = value, ymax = value + SD), 
                    width = 0.2, position = position_dodge(0.9),alpha = 0.3) +
      xlab("") + ylab("") + ggtitle(algname) + 
      theme(plot.title = element_text(hjust = 0.5), 
            axis.text = element_text(size = 15, face = "bold"), 
            axis.title = element_text(size = 14, face = "bold")) + 
      labs(fill = "") + textLabels
    print(p)
    
    return(p)
    
  }
  
  p_performance <- meltdata(performance_data = performance_data, performance_data_sd = performance_data_sd)
  p_training <- meltdata(performance_data = performance_data_train, performance_data_sd = performance_data_train_sd)
  p_testing <- meltdata(performance_data = performance_data_test, performance_data_sd = performance_data_test_sd)
  
  
  dd <- subset(performance_data_test, Algorithm == selectmod)
  dd_sd <- subset(performance_data_test_sd, Algorithm == selectmod)
  p_best <- meltdata(performance_data = dd, performance_data_sd = dd_sd)
  
  
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
