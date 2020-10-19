  setwd("./")
  
  newdat <- readRDS('newlipid_9_13_20.rds')

  newpd <- readRDS('newpd_pdonly.rds')
    
  library(gbm)
  
  library(caret)
  
  lilikoimat <- newdat[-1]
  
  newpd <- newpd[row.names(newdat),]
  newdatpd <- cbind(newdat, newpd)
  
  lilikoimat <- t(lilikoimat)
  
  lilikoilabels <- newdat$Label
  lilikoilabels[lilikoilabels == 'Preeclampsia'] <- 'Cancer'
  lilikoilabels[lilikoilabels == 'Control'] <- 'Normal'
  
  lilikoicands <- row.names(lilikoimat)
  
  
  source('./machine_learning11.R')
  
  
  lilikoires <- lilikoi.machine_learning11(PDSmatrix = lilikoimat, 
                                          measurementLabels = lilikoilabels, 
                                          significantPathways = lilikoicands, 
                                          selectmod = 'RF',
                                          dividep = 0.8, 
                                          dividseed = 1996,
                                          times = 10)
  pdf(file = "Training.pdf",width = 8, height = 8)
  print(lilikoires$p_training)
  dev.off()
  pdf(file = "Testing.pdf",width = 4, height = 8)
  print(lilikoires$p_best)
  dev.off()
