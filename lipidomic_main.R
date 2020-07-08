#Code for the publication **Maternal blood lipidomics analyse link critical metabolic pathways associated with 
#severe preeclampsia** authored by Yu Liu, Bing He and Mano Maurya *et al*.

#The code includes lipidomic data preprocessing, ANOVA, WGCNA, differential lipid identification, lipid pathway 
#mapping, machine learning model construction, etc.

#Created by Yu Liu on Jun 29, 2020 (uiluy@umich.edu, Department of Computational Medicine and Bioinformatics, 
#University of Michigan, MI)

rm(list = ls())

wkdir <- 
  'C:/Users/uiluy/Desktop/projects/preeclampsia_R01/EPIC_preeclampsiaR01_placenta/pipelinepackage/completedata/'

setwd(wkdir)

pd <- readRDS('pdwithsubtype_feature50000_added.rds')
pd$Parity[is.na(pd$Parity)] <- 1



wkdir <- 'C:/Users/uiluy/Desktop/projects/metabolomics/'

doscreen <- TRUE

setwd(wkdir)

datafile <- 'own.txt'

dataval <- read.table(datafile, sep = '\t', header = TRUE, 
                      stringsAsFactors = FALSE, check.names = FALSE)

groupidx <- grepl(pattern = '20190913_S', x = colnames(dataval))
groupidx[2] <- TRUE
dataval <- dataval[,groupidx]
colnames(dataval) <- gsub(pattern = '20190913_', replacement = '', x = colnames(dataval))
colnames(dataval) <- paste0('Sample', colnames(dataval))
colnames(dataval)[1] <- 'compound'

metafile <- 'ownmeta.txt'
meta <- read.table(metafile, sep = '\t', header = TRUE, 
                   stringsAsFactors = FALSE, check.names = FALSE)
meta$id <- as.character(meta$`Sample ID`)
#meta[1,1] <- '0002'
meta$`Sample ID` <- paste0('Sample', meta$`Sample ID`)
meta <- subset(meta, `Sample ID` != 'Sample')
dataval <- dataval[c('compound', meta$`Sample ID`)]

library(impute)

compounds <- dataval$compound

vals <- dataval[2:ncol(dataval)]
vals <- t(t(vals))
row.names(vals) <- 1:nrow(vals)
vals <- impute.knn(vals)
vals <- vals$data
vals <- as.data.frame(vals)
vals$compound <- compounds

mergeblock <- function(block){
  valpart <- block[-ncol(block)]
  valpart <- colMeans(valpart)
  valpart <- as.data.frame(valpart)
  valpart <- t(valpart)
  valpart <- as.data.frame(valpart)
  valpart$compound <- unique(block$compound)
  return(valpart)
  
}

library(plyr)

merged <- ddply(.data = vals, .variables = c('compound'), .fun = mergeblock)
merged <- merged[c('compound', colnames(merged)[1:(ncol(merged)-1)])]

compounds <- merged$compound
vals <- merged[-1]
vals <- t(vals)
colnames(vals) <- compounds
vals <- as.data.frame(vals)

colVars <- function(colval){
  
  avg <- mean(colval)
  colsd <- sum((colval - avg)^2)/(length(colval)-1)
  
  return(colsd)
  
}

colsds <- apply(X = vals, MARGIN = 2, FUN = colVars)
vals <- vals[,colsds != 0]

#write.table(vals, 'ownori.txt', sep = '\t', row.names = TRUE, quote = FALSE)

vals <- log2(vals + 1)

makeplot <- function(orimat = t(groupdata), logtrans = FALSE){
  
  if(logtrans == TRUE){
    orimat <- log2(orimat + 1)
  }
  
  boxdata <- t(orimat)
  boxdata <- as.data.frame(boxdata, stringsAsFactors = FALSE)
  
  library(reshape)
  md <- melt(boxdata)
  names(md) <- c('sample', 'value')
  
  library(ggplot2)
  
  p <- ggplot(md, aes(x = sample, y = value, fill = `sample`))
  print(
    p + geom_boxplot() + 
      xlab('Sample') + ylab('Value') + 
      scale_fill_discrete(guide = FALSE) + 
      theme_bw() + 
      theme(axis.text.x = element_blank()) + 
      ggtitle('Data distribution')
    
  )
  
}

#Standard scale#########
standardnorm <- function(colval){
  
  avg <- mean(colval)
  colsd <- sum((colval - avg)^2)/(length(colval)-1)
  colstandard <- (colval - avg)/sqrt(colsd)
  
  return(colstandard)
}

normvals <- apply(X = vals, MARGIN = 2, FUN = standardnorm)

orivals <- vals
#vals <- normvals


#Quantile normalization##########
library(preprocessCore)

samplenames <- rownames(vals)
metabolitenames <- colnames(vals)
#vals <- normalize.quantiles((t(vals)))
#vals <- t(vals)
rownames(vals) <- samplenames
colnames(vals) <- metabolitenames



#makeplot(orimat = vals, logtrans = FALSE)

vals <- as.data.frame(vals)
vals$samplename <- row.names(vals)
#meta$PRETERM <- as.character(meta$PRETERM)

#meta$PRETERM[meta$PRETERM == '1'] <- 'Case'
#meta$PRETERM[meta$PRETERM == '0'] <- 'Control'

meta <- meta[c('Sample ID', 'Group')]

names(meta) <- c('samplename', 'Label')
newdat <- merge(meta, vals, by = c('samplename'))
row.names(newdat) <- newdat$samplename
newdat <- newdat[-1]
newdat$Label <- factor(x = newdat$Label, 
                       levels = c('Controls', 'Cases'), ordered = TRUE)

#saveRDS(newdat, 'ownmetabolitemeasurements.rds')


#Median fold normalization###########
library(metabolomics)

metinput <- newdat
names(metinput)[1] <- 'Group'
metinput$Group <- as.character(metinput$Group)

metoutput <- Normalise(inputdata = metinput, method = c('median'))
metoutput <- metoutput$output
names(metoutput) <- c('Label', metabolitenames)

metoutput$Label <- factor(x = metoutput$Label, 
                          levels = c('Controls', 'Cases'), ordered = TRUE)
newdat <- metoutput

#saveRDS(metoutput, 'metabolitemeasurements.rds')

#makeplot(orimat = metoutput[-1], logtrans = FALSE)

makeplot(orimat = newdat[-1], logtrans = FALSE)

#saveRDS(newdat, 'ownmetabolitemeasurements.rds')




#Sample mapping#####
samplemappingfile <- 'samplemapping.txt'
samplemapping <- read.table(samplemappingfile, sep = '\t', header = TRUE, stringsAsFactors = FALSE, 
                            check.names = FALSE)

samplemapping$`Sample Name` <- gsub(pattern = 'Garmire_MB_0', replacement = '', x = samplemapping$`Sample Name`)
samplemapping$`Sample Name` <- gsub(pattern = 'Garmire_MB_', replacement = '', x = samplemapping$`Sample Name`)
samplemapping$`Sample Name` <- paste0(samplemapping$Group, '_', samplemapping$`Sample Name`)
samplemapping$`Sample Name` <- gsub(pattern = 'Controls_', replacement = 'Control_', 
                                    x = samplemapping$`Sample Name`)
samplemapping$`Sample Name` <- gsub(pattern = 'Cases_', replacement = 'Preeclampsia_', 
                                    x = samplemapping$`Sample Name`)
samplemapping$`Sample ID` <- paste0('Sample', samplemapping$`Sample ID`)
samplemapping <- subset(samplemapping, `Sample ID` %in% row.names(newdat))
newdat <- newdat[match(row.names(newdat), samplemapping$`Sample ID`),]
row.names(newdat) <- samplemapping$`Sample Name`
newdat$Label <- as.character(newdat$Label)
newdat$Label[newdat$Label == 'Controls'] <- 'Control'
newdat$Label[newdat$Label == 'Cases'] <- 'Preeclampsia'
row.names(pd) <- pd$Sample_Name
pd <- pd[row.names(newdat),]

#pd <- pd[-ncol(pd)]

#write.table(pd, 'pdfile.txt', sep = '\t', quote = FALSE, 
#            row.names = FALSE)
#SOV#########

basic <- c("Sample_Name", "Sample_Group", "Smoker", "Gender", "Ethgroup", "Age", "Parity", "BMI", 
           "GestationalAgeWeeks")
othervars <- names(pd)[-match(unique(c(basic, c('Sample_Group', 'Cluster', 'SGA'))), names(pd))]
pd <- pd[c(basic, othervars)]

#Remove highly correlated variables and nochange variables######

samples <- pd$Sample_Name

sovdat <- newdat

sovdat <- sovdat[-1]
sovdat <- t(sovdat)

sovdat <- sovdat[,samples]


#Organize beta matrix#######
sovdat <- t(sovdat)
sovdat <- as.data.frame(sovdat, stringsAsFactors = FALSE)

row.names(pd) <- pd$Sample_Name
pd <- pd[-1]
names(pd)[1] <- 'Sample_Group'

contvars <- c('Age', 'Parity', 'BMI', 'GestationalAgeWeeks', 'Weight_Gain', 'Mothers_Height', 'Length_of_Stay', 
              'Estimated_Blood_Loss_mL', 'Weight_Baby_g', 'Length_Baby_cm', 
              'One_Min_Score_Baby', 'Five_Min_Score_Baby', 'Head_Circumference_Baby_cm', 
              'Length_of_ROM_Hrs')
contvars <- intersect(contvars, names(pd))



for(i in 1:ncol(pd)){
  varname <- names(pd)[i]
  if(!(varname %in% contvars)){
    pd[,i] <- factor(pd[,i])
  }
}

names1 <- c('Sample_Group', 'Gender', 'Ethgroup', 'Age', 'Parity', 'GestationalAgeWeeks', 
            'Gestational_Diabetes', 'Chronic_Hypertension', 
            'Smoker', 'BMI')

namess <- c('Membrane_Ruptured', 'Abruption', 'Neonatal_Malformations_Baby')





#Calculate F#####
sovdat <- sovdat[row.names(pd),]

#Calculate F#####
Ftab <- data.frame(Sample_Group = numeric(), stringsAsFactors = FALSE)
for(i in 2:ncol(pd)){
  varname <- names(pd)[i]
  Ftab[varname] <- numeric()
}


calF <- function(probe = probecol){
  library(car)
  
  newdata <- pd
  pdnames <- names(newdata)
  newdata$beta <- probe
  
  formstr <- paste0(pdnames, collapse = ' + ')
  formstr <- paste0('beta ~ ', formstr)
  formstr <- as.formula(formstr)
  
  fit <- lm(formstr, data = newdata)
  
  aovfit <- Anova(fit, type = 3, singular.ok = TRUE)
  
  
  F <- aovfit$`F value`
  
  F <- F[2:(length(F)-1)]
  names(F) <- pdnames
  F <- as.data.frame(F, stringsAsFactors = FALSE)
  F <- as.data.frame(t(F))
  row.names(F) <- 1
  
  
  Ftab <- rbind(Ftab, F)
  
  return(Ftab)
}


library(parallel)

sovdatlist <- list()

for(i in 1:ncol(sovdat)){
  sovdatlist[[i]] <- sovdat[,i]
}


Ftab <- mclapply(X = sovdat, FUN = calF, mc.cores = 1)
#Ftab <- apply(X = sovdat[,1:50], MARGIN = 2, FUN = calF)

Ftab <- do.call(rbind, Ftab)


Fmean <- colMeans(Ftab)

Fmean <- Fmean[order(-Fmean)]

Fmean <- data.frame(Factor = names(Fmean), Fstat = as.vector(Fmean), stringsAsFactors = FALSE)

finalvars <- unique(c('Sample_Group', Fmean$Factor[Fmean$Fstat > 1]))


library(ggplot2)

sovplot <- function(restab = MSSmean, clustername = 'Preeclampsia', plottype = 'MSS', 
                    textsize = 20){
  
  resmean <- restab
  samplegroupidx <- match('Sample_Group', resmean$Factor)
  resmean$Factor[samplegroupidx] <- paste0(clustername, '_Control')
  
  if(plottype == 'MSS'){
    ytitle <- 'Mean Square'
    resmean <- resmean[order(-resmean$MSSstat),]
    resmean$Factor <- factor(resmean$Factor, levels = resmean$Factor, ordered = TRUE)
    
    p <- ggplot(data = resmean, mapping = aes(x = Factor, y = MSSstat, fill = Factor))
    print(
      p + geom_bar(stat = 'identity') + 
        ggtitle('Source of Variance (Type 3 Anova)') + 
        ylab(ytitle) + 
        xlab('') + 
        scale_fill_discrete(guide = FALSE) + 
        theme_bw() + 
        theme(axis.text.x = element_text(angle = 60, hjust = 1, size = textsize)) + 
        theme(axis.text.y = element_text(size = textsize))
    )
    
  }else if(plottype == 'pval'){
    ytitle <- '-log2(p-val)'
    resmean <- resmean[order(-resmean$logpval),]
    resmean$Factor <- factor(resmean$Factor, levels = resmean$Factor, ordered = TRUE)
    
    p <- ggplot(data = resmean, mapping = aes(x = Factor, y = logpval, fill = Factor))
    print(
      p + geom_bar(stat = 'identity') + 
        ggtitle('Source of Variance (Type 3 Anova)') + 
        ylab(ytitle) + 
        xlab('') + 
        scale_fill_discrete(guide = FALSE) + 
        geom_hline(yintercept = -log2(0.05), color = 'red', size = 1) + 
        theme_bw() + 
        theme(axis.text.x = element_text(angle = 60, hjust = 1, size = textsize)) + 
        theme(axis.text.y = element_text(size = textsize))
    )
  }else{
    ytitle <- 'F statistic'
    resmean <- resmean[order(-resmean$Fstat),]
    resmean$Factor <- factor(resmean$Factor, levels = resmean$Factor, ordered = TRUE)
    
    p <- ggplot(data = resmean, mapping = aes(x = Factor, y = Fstat, fill = Factor))
    print(
      p + geom_bar(stat = 'identity') + 
        ggtitle('') + 
        ylab(ytitle) + 
        xlab('') + 
        scale_fill_discrete(guide = FALSE) + 
        theme_bw() + 
        theme(axis.text.x = element_text(angle = 60, hjust = 1, size = textsize)) + 
        theme(axis.text.y = element_text(size = textsize))
    )
  }
  
}

sovplot(restab = Fmean, plottype = 'F', textsize = 15)


Fmeantmp <- Fmean
Fmeantmp$Factor <- c('Age', 'Ethgroup', 'BMI', 'ChronicHypertension', 'NeonatalMalformation', 
                     'MembraneRupture', 'Preeclampsia', 'GestationalDiabetes', 'GestationalAgeWeeks', 
                     'Parity', 'PlacentalAbruption', 'Gender', 'Smoker')

Fmeantmp$Factor[Fmeantmp$Factor == 'Gender'] <- 'BabyGender'

sovplot(restab = Fmeantmp, plottype = 'F', textsize = 15)

#correlation matrices############
library(lilikoi)
library(corrplot)

pdval <- pd

for(i in 1:ncol(pdval)){
  pdval[,i] <- as.numeric(pdval[,i])
  
}

pdcor <- cor(pdval)

corrplot(pdcor, method = 'circle', type = 'upper', order = 'hclust', tl.col = 'black', tl.srt = 45, 
         tl.cex = 1, number.cex = 0.3, addrect = 4)





pdvaltmp <- pdval

names(pdvaltmp) <- c('Preeclampsia', 'Smoker', 'BabyGender', 'Ethgroup', 'Age', 'Parity', 'BMI', 
                     'GestationalAgeWeeks', 'GestationalDiabetes', 'ChronicHypertension', 
                     'MembraneRupture', 'PlacentalAbruption', 'NeonatialMalformation')

pdvaltmp <- pdvaltmp[c('Preeclampsia', 'BabyGender', 'Age', 'NeonatialMalformation', 'GestationalAgeWeeks', 
                       'BMI', 'MembraneRupture', 'PlacentalAbruption', 'Ethgroup', 'Parity', 
                       'ChronicHypertension', 'Smoker', 'GestationalDiabetes')]


pdcor <- cor(pdvaltmp)

corrplot(pdcor, method = 'circle', type = 'upper', order = 'original', tl.col = 'black', tl.srt = 45, 
         tl.cex = 1, number.cex = 0.3, addrect = 4)


#Clinical statistic#####
pdc <- pd[grep(pattern = 'Control', x = row.names(pd)),]
pdp <- pd[grep(pattern = 'Preeclampsia', x = row.names(pd)),]

constats <- function(convar = 'Age', datc = pdc, datp = pdp){
  
  names(datc)[grep(pattern = convar, x = names(datc))] <- 'var'
  names(datp)[grep(pattern = convar, x = names(datp))] <- 'var'
  
  meanc <- mean(datc$var)
  meanp <- mean(datp$var)
  sdc <- sd(datc$var)
  sdp <- sd(datp$var)
  
  tres <- t.test(datc$var, datp$var)
  tp <- tres$p.value
  
  res <- data.frame(varname = convar, ctrlmean = meanc, preemean = meanp, ctrlsd = sdc, preesd = sdp, 
                    pval = tp, stringsAsFactors = FALSE)
  return(res)
}

constats('Age')
constats('BMI')
constats('GestationalAgeWeeks')



disstats <- function(disvar = 'Chronic_Hypertension', dat = pd){
  
  names(dat)[grep(pattern = disvar, x = names(dat))] <- 'var'
  
  freqtab <- table(dat$Sample_Group, dat$var)
  
  fres <- fisher.test(freqtab)
  fp <- fres$p.value
  
  res <- list(varname = disvar, freq = freqtab, pval = fp)
  return(res)
}

disstats('Smoker')
disstats('Gender')
disstats('Ethgroup')
disstats('Parity')
disstats('Gestational_Diabetes')
disstats('Chronic_Hypertension')
disstats('Membrane_Ruptured')
disstats('Abruption')
disstats('Neonatal_Malformations_Baby')




#Screen data########
screenedftab <- Ftab[Ftab$Sample_Group > 1,]

screenedfmean <- colMeans(screenedftab)
screenedfmean <- screenedfmean[order(-screenedfmean)]
screenedfmean <- data.frame(Factor = names(screenedfmean), Fstat = as.vector(screenedfmean), 
                            stringsAsFactors = FALSE)
sovplot(restab = screenedfmean, plottype = 'F', textsize = 15)





screenedfmeantmp <- screenedfmean

screenedfmeantmp$Factor <- c('Preeclampsia', 'Age', 'Ethgroup', 'MembraneRupture', 'ChronicHypertension', 
                             'BMI', 'GestationalAgeWeeks', 'NeonatalMalformations', 'GestationalDiabetes', 
                             'Parity', 'BabyGender', 'Smoker', 'PlacentalAbruption')

sovplot(restab = screenedfmeantmp, plottype = 'F', textsize = 15)






screenedmets <- row.names(screenedftab)
screeneddat <- newdat[c('Label', screenedmets)]

#Correlation matrix######
library(pheatmap)

screenedcor <- cor(screeneddat[-1], pdval)
names(screenedcor)[3] <- 'Baby_Gender'
pheatmap(screenedcor, 
         color = colorRampPalette(colors = c('green', 'black', 'red'))(100), 
         show_rownames = FALSE, 
         scale = 'none', fontsize_row = fontsizerow)


newvarorder <- c(8, 11, 12, 7, 4, 6, 5, 13, 1, 3, 10, 2, 9)

pheatmap(screenedcor[,newvarorder], 
         color = colorRampPalette(colors = c('green', 'black', 'red'))(100), 
         show_rownames = FALSE, 
         scale = 'none', fontsize_row = fontsizerow, 
         cluster_cols = FALSE)




screenedcortmp <- screenedcor
colnames(screenedcortmp) <- c('Preeclampsia', 'Smoker', 'BabyGender', 'Ethgroup', 'Age', 'Parity', 'BMI', 
                              'GestationalAgeWeeks', 'GestationalDiabetes', 'ChronicHypertension', 
                              'MembraneRupture', 'PlacentalAbruption', 'NeonatialMalformation')
pheatmap(screenedcortmp, 
         color = colorRampPalette(colors = c('green', 'black', 'red'))(100), 
         show_rownames = FALSE, 
         scale = 'none', fontsize_row = fontsizerow)

pheatmap(screenedcortmp[,newvarorder], 
         color = colorRampPalette(colors = c('green', 'black', 'red'))(100), 
         show_rownames = FALSE, 
         scale = 'none', fontsize_row = fontsizerow, 
         cluster_cols = FALSE)

#Lipid cluster check (3 clusters)###


out <- pheatmap(screenedcortmp, 
                color = colorRampPalette(colors = c('green', 'black', 'red'))(100), 
                show_colnames = TRUE, 
                scale = 'row', show_rownames = FALSE, 
                cluster_rows = TRUE, cluster_cols = TRUE, 
                legend = FALSE, annotation_legend = FALSE)


lipidcolor <- sort(cutree(out$tree_row, k = 3))

mode(lipidcolor) <- 'character'

lipidgroup <- data.frame(Cluster=factor(lipidcolor, labels = paste0('Cluster', seq(3))))

print(
  
  pheatmap(screenedcortmp, annotation_row = lipidgroup, 
           color = colorRampPalette(colors = c('green', 'black', 'red'))(100), 
           show_colnames = TRUE, 
           scale = 'row', show_rownames = FALSE, 
           cluster_rows = TRUE, cluster_cols = TRUE, 
           legend = FALSE, annotation_legend = TRUE)
  
  
)





tcortmp <- t(screenedcortmp)

tout <- pheatmap(tcortmp, 
                 color = colorRampPalette(colors = c('green', 'black', 'red'))(100), 
                 show_colnames = FALSE, 
                 scale = 'column', show_rownames = TRUE, 
                 cluster_rows = TRUE, cluster_cols = TRUE, 
                 legend = FALSE, annotation_legend = FALSE)


lipidcolor <- sort(cutree(tout$tree_col, k = 3))

mode(lipidcolor) <- 'character'

tlipidgroup <- data.frame(Cluster=factor(lipidcolor, labels = paste0('Cluster', seq(3))))



print(
  
  pheatmap(tcortmp, annotation_col = lipidgroup, 
           color = colorRampPalette(colors = c('green', 'black', 'red'))(100), 
           show_colnames = FALSE, 
           scale = 'column', show_rownames = TRUE, 
           cluster_rows = TRUE, cluster_cols = TRUE, 
           legend = TRUE, annotation_legend = TRUE, 
           fontsize_row = 15)
  
)






alldat <- newdat[-1]
lipidcor <- cor(alldat, pdval)

library(pheatmap)

pheatmap(lipidcor, 
         color = colorRampPalette(colors = c('green', 'black', 'red'))(100), 
         show_rownames = FALSE, 
         scale = 'none', fontsize_row = fontsizerow)


lipidcortmp <- lipidcor
colnames(lipidcortmp) <- c('Preeclampsia', 'Smoker', 'BabyGender', 'Ethgroup', 'Age', 'Parity', 'BMI', 
                           'GestationalAgeWeeks', 'GestationalDiabetes', 'ChronicHypertension', 
                           'MembraneRupture', 'PlacentalAbruption', 'NeonatialMalformation')
pheatmap(lipidcortmp, 
         color = colorRampPalette(colors = c('green', 'black', 'red'))(100), 
         show_rownames = FALSE, 
         scale = 'none', fontsize_row = fontsizerow)






library(scales)

lipidspeciecolor <- data.frame(lipid = c('ACar', 'CE', 'Cer-AS', 'Cer_NDS', 'Cer_NS', 
                                         'DAG', 'HexCer-NS', 'LPC', 'LPE', 'OxPC', 
                                         'OxPE', 'PC', 'PE', 'PG', 'PS', 
                                         'SM', 'TAG'), 
                               lipidcolor = c('#F8766D', '#E7851E', '#FF689E', '#D09400', '#B2A100', 
                                              '#89AC00', '#45B500', '#00BC51', '#00C087', '#00C0B2', 
                                              '#00BCD6', '#00B3F2', '#29A3FF', '#9C8DFF', '#D277FF', 
                                              '#F166E8', '#FF61C7'), 
                               stringsAsFactors = FALSE)


lipidgroupenrich <- function(lipidclusteranno = tlipidgroup, lipidspeciecolor = lipidspeciecolor){
  
  getlipidspecies <- function(lipidnames = row.names(lipidclusteranno)){
    
    processed <- gsub(pattern = 'Unknown ', replacement = '', x = lipidnames)
    processed <- gsub(pattern = ' .*$', replacement = '', x = processed)
    processed <- gsub(pattern = '\\(.*$', replacement = '', x = processed)
    
    return(processed)
    
  }
  
  lipidclusteranno$species <- getlipidspecies()
  
  enrichres <- function(targets = sub$species, background = lipidclusteranno$species){
    
    targetssum <- length(targets)
    backsum <- length(background)
    
    uniquetargets <- unique(targets)
    for(i in 1:length(uniquetargets)){
      uniquetarget <- uniquetargets[i]
      a11 <- sum(targets == uniquetarget)
      a12 <- sum(background == uniquetarget)
      a21 <- targetssum - a11
      a22 <- backsum - a12
      mat <- matrix(c(a11, a12, a21, a22), nrow = 2, byrow = TRUE)
      fisheres <- fisher.test(mat)
      fisherp <- fisheres$p.value
      ratio <- (a11/(a11 + a21))/(a12/(a12 + a22))
      if(i == 1){
        fisherps <- fisherp
        ratios <- ratio
      }else{
        fisherps <- c(fisherps, fisherp)
        ratios <- c(ratios, ratio)
      }
    }
    
    res <- data.frame(lipid = uniquetargets, fisherp = fisherps, ratio = ratios, stringsAsFactors = FALSE)
    
    return(res)
    
  }
  
  groupnames <- unique(lipidclusteranno$Cluster)
  
  i <- 1
  
  for(i in 1:length(groupnames)){
    
    groupname <- groupnames[i]
    sub <- lipidclusteranno[lipidclusteranno$Cluster == groupname,]
    
    subspecies <- unique(sub$species)
    
    enrichresult <- enrichres(targets = sub$species, background = lipidclusteranno$species)
    enrichresult$cluster <- groupname
    
    if(i == 1){
      enrichresults <- enrichresult
    }else{
      enrichresults <- rbind(enrichresults, enrichresult)
    }
    
  }
  
  enrichresults <- unique(enrichresults)
  
  plotdat <- subset(enrichresults, fisherp < 0.05)
  plotdat <- merge(plotdat, lipidspeciecolor, by = c('lipid'))
  plotdat <- plotdat[order(plotdat$cluster, -plotdat$fisherp),]
  plotdat$logp <- -log2(plotdat$fisherp)
  
  i <- 1
  for(i in 1:nrow(plotdat)){
    line <- plotdat[i,]
    sub <- subset(lipidclusteranno, Cluster == line$cluster)
    a11 <- sum(sub$species == line$lipid)
    a12 <- sum(lipidclusteranno$species == line$lipid)
    a21 <- nrow(sub) - a11
    a22 <- nrow(lipidclusteranno) - a21
    mat <- matrix(c(a11, a12, a21, a22), nrow = 2, byrow = TRUE)
    
    greaterp <- fisher.test(mat, alternative = 'greater')$p.value
    lessp <- fisher.test(mat, alternative = 'less')$p.value
    
    if(greaterp < lessp){
      line$dir <- 'UP'
    }else{
      line$dir <- 'DN'
    }
    
    if(i == 1){
      lines <- line
    }else{
      lines <- rbind(lines, line)
    }
    
  }
  
  plotdat$logp[lines$dir == 'DN'] <- -plotdat$logp[lines$dir == 'DN']
  
  library(ggplot2)
  
  clusternames <- unique(plotdat$cluster)
  
  i <- 1
  for(i in 1:length(clusternames)){
    clustername <- clusternames[i]
    plotsub <- subset(plotdat, cluster == clustername)
    plotsub$lipidname <- factor(plotsub$lipid, levels = plotsub$lipid, ordered = TRUE)
    
    p <- ggplot(plotsub, aes(x = lipidname, y = logp))
    print(
      p + geom_bar(stat = 'identity', fill = plotsub$lipidcolor) + 
        xlab('') + ylab('') + 
        ggtitle(paste0('Significantly enriched lipids in ', clustername)) + 
        theme_bw() + 
        theme(panel.grid = element_blank()) + 
        theme(axis.text.x = element_text(angle = 90, size = 30)) + 
        theme(axis.text.y = element_text(size = 30)) + 
        theme(plot.title = element_text(size = 25)) + 
        theme(panel.border = element_blank(), axis.line = element_line()) + 
        scale_y_continuous(position = 'right') + 
        coord_flip()
    )
    
    
  }
  
  return(plotdat)
  
  
}


lipidenrichres <- lipidgroupenrich(lipidclusteranno = tlipidgroup, lipidspeciecolor = lipidspeciecolor)



#Correlation between Preeclampsia and Baby gender#####

fullpd <- readRDS('../preeclampsia_R01/EPIC_preeclampsiaR01_placenta/pipelinepackage/completedata/pdwithsubtype_feature50000.rds')

table(pd[c('Sample_Group', 'Gender')])

fisher.test(table(pd[c('Sample_Group', 'Gender')]))$p.value

table(fullpd[c('Sample_Group', 'Gender')])

fisher.test(table(fullpd[c('Sample_Group', 'Gender')]))$p.value















calpval <- function(dat = screeneddat, pddat = pd){
  
  pddat$met <- dat[,2]
  fit <- lm(met ~ ., data = pddat)
  
  aovfit <- Anova(fit, type = 3, singular.ok = TRUE)
  df1s <- data.frame(Factor = row.names(aovfit), Df1 = aovfit$Df, stringsAsFactors = FALSE)
  samplenum <- nrow(dat)
  df1s$Df2 <- samplenum - 1 - df1s$Df1
  
  return(df1s)
  
}

dfs <- calpval()
screenedfmean <- merge(screenedfmean, dfs, by = c('Factor'))
screenedfmean$pval <- pf(screenedfmean$Fstat, screenedfmean$Df1, screenedfmean$Df2, lower.tail = FALSE)

depositnewdat <- newdat

if(doscreen == TRUE){
  newdat <- screeneddat
  finalvars <- unique(c('Sample_Group', screenedfmean$Factor[screenedfmean$Fstat > 1]))
}


#Regression out######

regdat <- newdat[-1]
regpd <- pd

regout <- function(probe = regdat[,1], pddat = regpd, removegest = TRUE){
  
  library(car)
  
  newdata <- pddat
  pdnames <- names(newdata)
  newdata$beta <- probe
  
  if(removegest == TRUE){
    
    pdnames <- pdnames[pdnames != 'GestationalAgeWeeks']
    
  }
  formstr <- paste0(pdnames, collapse = ' + ')
  formstr <- paste0('beta ~ ', formstr)
  formstr <- as.formula(formstr)
  
  fit <- lm(formstr, data = newdata)
  
  groupvals <- newdata$Sample_Group
  groupvals <- as.numeric(groupvals)
  groupvals <- groupvals - 1
  groupcoeff <- fit$coefficients
  groupcoeff <- as.vector(groupcoeff)
  groupcoeff <- groupcoeff[2]
  groupvals[groupvals == 1] <- groupcoeff
  
  resvals <- fit$residuals
  resvals <- as.numeric(resvals)
  
  adjvals <- groupvals - resvals
  
  names(adjvals) <- row.names(pddat)
  adjvals <- data.frame(lipid = adjvals, stringsAsFactors = FALSE)
  
  return(adjvals)
}

library(parallel)

regdatlist <- list()

for(i in 1:ncol(regdat)){
  regdatlist[[i]] <- regout(regdat[,i])
  names(regdatlist)[i] <- colnames(regdat)[i]
  colnames(regdatlist[[i]]) <- names(regdatlist)[i]
}


regoutdat <- do.call(cbind, regdatlist)
grouplabels <- gsub(pattern = '_.*$', replacement = '', x = row.names(regoutdat))

lipidnames <- colnames(regoutdat)
regoutdat$Label <- grouplabels
regoutdat <- regoutdat[c('Label', lipidnames)]



#Limma#####
wogestwk <- FALSE

finalvarstmp <- finalvars

if(wogestwk == TRUE){
  
  finalvars <- finalvars[finalvars != 'GestationalAgeWeeks']
  
}


clusterstr <- 'Preeclampsia to Control'

library(limma)

vars <- paste0('Nonlinear Scaling Preeclampsia + Confounding')

formulastr <- paste0(finalvars, collapse = ' + ')
formulastr <- paste0('~ ', formulastr)
formulastr <- formula(formulastr)

design <- model.matrix(formulastr, data = pd)

normfirst <- newdat[-1]
normfirst <- normfirst[row.names(pd),]
normfirst <- t(normfirst)


fit1 <- lmFit(normfirst, design)
fit2 <- eBayes(fit1)

head(topTable(fit2, coef = 2, n = nrow(fit2)))

topmet <- topTable(fit2, coef = 2, n = nrow(fit2))
topmet <- row.names(topmet[topmet$P.Value < 0.01,])
topdat <- newdat[,c('Label', topmet)]

limmapvals  <- topTable(fit2, coef = 2, n = nrow(fit2))
limmapvals$compound <- row.names(limmapvals)
row.names(limmapvals) <- 1:nrow(limmapvals)
limmapvals <- limmapvals[c('compound')]

#Heatmap#######

drawheatmap <- function(datamat = newdat, title = 'Metabolites', whetherrowname = FALSE, 
                        clustermethod = 'average', 
                        casename = 'Preeclampsia', ctrlname = 'Control', 
                        fontsizerow = 10, 
                        whetherclustergroup = TRUE){
  
  library(pheatmap)
  
  datamat <- datamat[order(datamat$Label),]
  
  datanno <- data.frame(group = datamat$Label, stringsAsFactors = FALSE)
  row.names(datanno) <- row.names(datamat)
  datanno$group <- factor(datanno$group, levels = c(casename, ctrlname), ordered = TRUE)
  datamat <- t(datamat[-1])
  
  library(matrixStats)
  
  rowvar <- rowVars(datamat)
  removerows <- rowvar < 0.01 
  datamat <- datamat[!removerows,]
  
  if(whetherclustergroup == FALSE){
    print(
      pheatmap(datamat, annotation_col = datanno, 
               color = colorRampPalette(colors = c('green', 'black', 'red'))(100), 
               show_colnames = FALSE, 
               main = title, 
               scale = 'row', show_rownames = whetherrowname, 
               cluster_cols = FALSE, fontsize_row = fontsizerow)
    )
  }else{
    print(
      pheatmap(datamat, annotation_col = datanno, 
               color = colorRampPalette(colors = c('green', 'black', 'red'))(100), 
               show_colnames = FALSE, 
               main = title, 
               scale = 'row', show_rownames = whetherrowname, 
               clustering_method = clustermethod, fontsize_row = fontsizerow)
    )
  }
  
  
  
}


#drawheatmap(datamat = topdat, title = 'Top differential lipids', whetherclustergroup = FALSE, 
#            whetherrowname = TRUE)

#drawheatmap()
#drawheatmap(datamat = topdat, title = 'Top metabolites', whetherrowname = TRUE, fontsizerow = 10)

drawheatmap(clustermethod = 'complete')
drawheatmap(datamat = topdat, title = 'Top metabolites', whetherrowname = TRUE, clustermethod = 'complete', fontsizerow = 10)



drawheatmap2 <- function(datamat = newdat, title = 'Metabolites', whetherrowname = FALSE, 
                         clustermethod = 'complete', 
                         casename = 'Preeclampsia', ctrlname = 'Control', 
                         fontsizerow = 10, clusternum = 5){
  
  library(pheatmap)
  
  datanno <- data.frame(group = datamat$Label, stringsAsFactors = FALSE)
  row.names(datanno) <- row.names(datamat)
  datanno$group <- factor(datanno$group, levels = c(casename, ctrlname), ordered = TRUE)
  datamat <- t(datamat[-1])
  
  library(matrixStats)
  
  rowvar <- rowVars(datamat)
  removerows <- rowvar < 0.01 
  datamat <- datamat[!removerows,]
  
  library(pheatmap)
  hcc <- hclust(dist(datamat), method = clustermethod)
  lipidorder <- hcc$order
  lipidcolor <- cutree(hcc, clusternum)
  mode(lipidcolor) <- 'character'
  lipidseq <- row.names(datamat)[lipidorder]
  
  lipidgroup <- data.frame(Cluster=factor(lipidcolor, labels = paste0('Cluster', seq(clusternum))))
  
  hcc2 <- hclust(dist(t(datamat)), method = clustermethod)
  sampleorder <- hcc2$order
  
  sampleseq <- colnames(datamat)[sampleorder]
  
  datmat <- datamat[lipidseq, sampleseq]
  
  print(
    pheatmap(datmat, annotation_col = datanno, annotation_row = lipidgroup, 
             color = colorRampPalette(colors = c('green', 'black', 'red'))(100), 
             show_colnames = FALSE, 
             main = title, 
             scale = 'row', show_rownames = whetherrowname, 
             cluster_rows = FALSE, cluster_cols = FALSE)
  )
  
}

drawheatmap2(clusternum = 5)

#Wilcox comparision######

wilcoxcomp <- function(dat = newdat, pddat = pd, var = 'Gestational_Diabetes', pcutoff = 0.01){
  
  varlevels <- unique(pddat[,var])
  
  if(length(varlevels) != 2){
    return(NULL)
  }
  
  samplelist <- list()
  i <- 1
  
  for(i in 1:length(varlevels)){
    
    varlevel <- varlevels[i]
    samples <- pddat[pddat[,var] == varlevel,]
    samplelist[[i]] <- samples
  }
  
  names(samplelist) <- c('ctrsamples', 'dissamples')
  ctrsamples <- samplelist$ctrsamples
  dissamples <- samplelist$dissamples
  
  ctrgroup <- dat[row.names(ctrsamples),]
  disgroup <- dat[row.names(dissamples),]
  ctrgroup <- ctrgroup[-1]
  disgroup <- disgroup[-1]
  disgroup <- disgroup[,colnames(ctrgroup)]
  ctrgroup <- t(ctrgroup)
  disgroup <- t(disgroup)
  
  i <- 1
  
  for(i in 1:nrow(ctrgroup)){
    
    lipid <- row.names(ctrgroup)[i]
    ctrvals <- ctrgroup[i,]
    disvals <- disgroup[i,]
    
    wilcoxres <- wilcox.test(disvals, ctrvals)
    wilcoxp <- wilcoxres$p.value
    
    ctrmean <- mean(ctrvals)
    dismean <- mean(disvals)
    
    resline <- data.frame(lipid = lipid, dismean = dismean, ctrmean = ctrmean, wilcoxp = wilcoxp, 
                          stringsAsFactors = FALSE)
    
    if(i == 1){
      reslines <- resline
    }else{
      reslines <- rbind(reslines, resline)
    }
    
  }
  
  reslines$wilcoxp.adj <- p.adjust(reslines$wilcoxp, method = 'BH')
  
  wilsiglipids <- subset(reslines, wilcoxp < pcutoff)$lipid
  
  wilcoxdat <- dat[c('Label', wilsiglipids)]
  
  wilcoxdat$Label[row.names(wilcoxdat) %in% row.names(ctrsamples)] <- 'Control'
  wilcoxdat$Label[row.names(wilcoxdat) %in% row.names(dissamples)] <- 'Case'
  
  reslist <- list()
  
  reslines <- reslines[order(reslines$wilcoxp, reslines$wilcoxp.adj),]
  fullreslines <- reslines
  reslines <- subset(reslines, wilcoxp < pcutoff)
  
  reslist[[1]] <- wilcoxdat
  reslist[[2]] <- reslines
  reslist[[3]] <- fullreslines
  
  drawheatmap(datamat = wilcoxdat, title = 'Wilcox Top metabolites', whetherrowname = TRUE, 
              clustermethod = 'complete', fontsizerow = 10, casename = 'Case')
  
  drawheatmap(datamat = wilcoxdat, title = 'Wilcox Top metabolites', whetherrowname = TRUE, fontsizerow = 10, 
              casename = 'Case')
  
  drawheatmap(datamat = wilcoxdat, title = 'Wilcox top differential lipids', whetherrowname = TRUE, 
              fontsizerow = 10, whetherclustergroup = FALSE, casename = 'Case')
  
  names(reslist) <- c('wilcoxdat', 'reslines', 'fullreslines')
  
  return(reslist)
  
  
}

wilcoxdatres <- wilcoxcomp(var = 'Sample_Group')

diabeteswilcoxdatres <- wilcoxcomp(var = 'Gestational_Diabetes')

intersect(colnames(wilcoxdatres$wilcoxdat), colnames(diabeteswilcoxdatres$wilcoxdat))[-1]

pree_diabetes_lipids <- 
  intersect(colnames(wilcoxdatres$wilcoxdat), colnames(diabeteswilcoxdatres$wilcoxdat))[-1]



#Lilkoi Classifiers 1#####
newdattmp <- newdat

library(lilikoi)
library(gbm)

library(caret)

lilikoimat <- newdat[-1]

newpd <- pdval
newpd <- newpd[-1]
newpd <- newpd[-match('GestationalAgeWeeks', colnames(newpd))]
newpd <- newpd[row.names(newdat),]
newdatpd <- cbind(newdat, newpd)


#lilikoimat <- newdatpd[-1]

#lilikoimat <- newpd


lilikoimat <- t(lilikoimat)

lilikoilabels <- newdat$Label
lilikoilabels[lilikoilabels == 'Preeclampsia'] <- 'Cancer'
lilikoilabels[lilikoilabels == 'Control'] <- 'Normal'

lilikoicands <- row.names(lilikoimat)


#lilikoimat <- readRDS('lilikoimat.rds')



#Feature Selection###
library(RWeka)
library(lilikoi)

newdatt <- newdat

#newdatt <- newdatpd

#newdatt <- cbind(newdat[1], newpd)


newdatt$Label <- as.character(newdatt$Label)

newdatt$Label[newdatt$Label == 'Preeclampsia'] <- 'Cancer'
newdatt$Label[newdatt$Label == 'Control'] <- 'Normal'
newdatt$Label <- factor(newdatt$Label, levels = c('Normal', 'Cancer'), ordered = TRUE)





significantPathways <- lilikoi.select_pathways(PDSmatrix = lilikoimat, 
                                               metaboliteMeasurements = newdatt, 
                                               threshold = 0.1, method = 'gain')

significantPathways1 <- significantPathways
#significantPathways <- row.names(lilikoimat)



source('./machine_learning11.R')

lilikoires1 <- lilikoi.machine_learning11(PDSmatrix = lilikoimat, 
                                          measurementLabels = lilikoilabels, 
                                          significantPathways = significantPathways, 
                                          selectmod = 'LDA', 
                                          cvnum = 10, randomseed = 1, 
                                          dividep = 0.8)

logauc1 <- lilikoires1[[3]]


#Lilkoi Classifiers 2#####
newdat <- newdattmp

library(lilikoi)
library(gbm)

library(caret)

lilikoimat <- newdat[-1]

newpd <- pdval
newpd <- newpd[-1]
newpd <- newpd[-match('GestationalAgeWeeks', colnames(newpd))]
newpd <- newpd[row.names(newdat),]
newdatpd <- cbind(newdat, newpd)


lilikoimat <- newdatpd[-1]

#lilikoimat <- newpd


lilikoimat <- t(lilikoimat)

lilikoilabels <- newdat$Label
lilikoilabels[lilikoilabels == 'Preeclampsia'] <- 'Cancer'
lilikoilabels[lilikoilabels == 'Control'] <- 'Normal'

lilikoicands <- row.names(lilikoimat)


#lilikoimat <- readRDS('lilikoimat.rds')



#Feature Selection###
library(RWeka)
library(lilikoi)

newdatt <- newdat

newdatt <- newdatpd

#newdatt <- cbind(newdat[1], newpd)


newdatt$Label <- as.character(newdatt$Label)

newdatt$Label[newdatt$Label == 'Preeclampsia'] <- 'Cancer'
newdatt$Label[newdatt$Label == 'Control'] <- 'Normal'
newdatt$Label <- factor(newdatt$Label, levels = c('Normal', 'Cancer'), ordered = TRUE)





significantPathways <- lilikoi.select_pathways(PDSmatrix = lilikoimat, 
                                               metaboliteMeasurements = newdatt, 
                                               threshold = 0.1, method = 'gain')


#significantPathways <- row.names(lilikoimat)



#source('./machine_learning11.R')

lilikoires1 <- lilikoi.machine_learning11(PDSmatrix = lilikoimat, 
                                          measurementLabels = lilikoilabels, 
                                          significantPathways = significantPathways, 
                                          selectmod = 'LDA', 
                                          cvnum = 10, randomseed = 1, 
                                          dividep = 0.8)

log2auc1 <- lilikoires1[[3]]


#Lilkoi Classifiers 3#####
newdat <- newdattmp

library(lilikoi)
library(gbm)

library(caret)

lilikoimat <- newdat[-1]

newpd <- pdval
newpd <- newpd[-1]
newpd <- newpd[-match('GestationalAgeWeeks', colnames(newpd))]
newpd <- newpd[row.names(newdat),]
newdatpd <- cbind(newdat, newpd)


#lilikoimat <- newdatpd[-1]

lilikoimat <- newpd


lilikoimat <- t(lilikoimat)

lilikoilabels <- newdat$Label
lilikoilabels[lilikoilabels == 'Preeclampsia'] <- 'Cancer'
lilikoilabels[lilikoilabels == 'Control'] <- 'Normal'

lilikoicands <- row.names(lilikoimat)


#lilikoimat <- readRDS('lilikoimat.rds')



#Feature Selection###
library(RWeka)
library(lilikoi)

newdatt <- newdat

#newdatt <- newdatpd

newdatt <- cbind(newdat[1], newpd)


newdatt$Label <- as.character(newdatt$Label)

newdatt$Label[newdatt$Label == 'Preeclampsia'] <- 'Cancer'
newdatt$Label[newdatt$Label == 'Control'] <- 'Normal'
newdatt$Label <- factor(newdatt$Label, levels = c('Normal', 'Cancer'), ordered = TRUE)





significantPathways <- lilikoi.select_pathways(PDSmatrix = lilikoimat, 
                                               metaboliteMeasurements = newdatt, 
                                               threshold = 0.1, method = 'gain')


#significantPathways <- row.names(lilikoimat)



#source('./machine_learning11.R')

lilikoires1 <- lilikoi.machine_learning11(PDSmatrix = lilikoimat, 
                                          measurementLabels = lilikoilabels, 
                                          significantPathways = significantPathways, 
                                          selectmod = 'LDA', 
                                          cvnum = 10, randomseed = 1, 
                                          dividep = 0.8)

log3auc1 <- lilikoires1[[3]]




library(scales)

plot(logauc1$prplotobj, col = hue_pal()(3)[1], cex.lab = 1.5)
par(new = TRUE)
plot(log2auc1$prplotobj, col = hue_pal()(3)[2], cex.lab = 1.5)
par(new = TRUE)
plot(log3auc1$prplotobj, col = hue_pal()(3)[3], cex.lab = 1.5)
legend(0.2, 0.85, 
       legend = c(paste0('lipid AUC-PR ', signif(as.numeric(logauc1$probj$auc.integral), 3)), 
                  paste0('lipid+conf AUC-PR ', signif(as.numeric(log2auc1$probj$auc.integral), 3)), 
                  paste0('conf AUC-PR ', signif(as.numeric(log3auc1$probj$auc.integral), 3))), 
       col = hue_pal()(3), lty = 1:2, cex = 1)




#Overal WGCNA correlation network#####
tmpnewdat <- newdat

newdat <- regoutdat



ctrldat <- subset(newdat, Label == 'Control')
preedat <- subset(newdat, Label != 'Control')

overaldat <- newdat

#Choose the soft-thresholding power#
softval <- NULL

topquantile <- 0.95
internum <- 1000
largevexsize <- 15
largeedgsize <- 15
absolutecut <- 0.4

tag <- 'Preeclampsia to Control'

library(WGCNA)

choosepower <- function(datres = ctrldat, rsqcutline = 0.7){
  
  library(WGCNA)
  
  datExpr <- datres[-1]
  
  powers <- c(c(1:20))
  
  sft <- pickSoftThreshold(datExpr, powerVector = powers, RsquaredCut = rsqcutline, verbose = 5)
  
  
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], 
       xlab = "Soft Threshold (power)", 
       ylab = "Scale Free Topology Model Fit, signed R^2", type = "n", main = paste("Scale independence"))
  
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*(sft$fitIndices[,2]), labels = powers, col = 'red')
  
  abline(h = rsqcutline, col = 'red')
  
  
  
  plot(sft$fitIndices[,1], sft$fitIndices[,5], 
       xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n", 
       main = "Mean connectivity")
  
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels = powers, col = 'red')
  
  return(sft)
  
}

if(is.null(softval)){
  
  ctrlpower <- choosepower(rsqcutline = 0.8)
  
  casepower <- choosepower(datres = preedat, rsqcutline = 0.8)
  
  finalsft <- max(ctrlpower$powerEstimate, casepower$powerEstimate)
  
}else{
  finalsft <- softval
}


makemain <- function(datres = ctrldat, sftpower = finalsft, mergesimilar = TRUE, mergecut = 0.25){
  
  library(WGCNA)
  
  datExpr <- datres[-1]
  probes <- names(datExpr)
  
  #Calculate TOM
  adjacency <- adjacency(datExpr, power = sftpower)
  TOM <- TOMsimilarity(adjacency)
  dimnames(TOM) <- list(probes, probes)
  
  dissTOM <- 1 - TOM
  
  # Call the hierarchical clustering function
  geneTree <- hclust(as.dist(dissTOM), method = "average")
  
  #Module identification using dynamic tree cut
  dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                               deepSplit = 2, pamRespectsDendro = FALSE,
                               minClusterSize = 30)
  
  dynamicColors <- labels2colors(dynamicMods)
  oricolors <- dynamicColors
  finalColors <- dynamicColors
  
  # Calculate eigengenes
  MEList <- moduleEigengenes(datExpr, colors = dynamicColors)
  MEs <- MEList$eigengenes
  
  if(mergesimilar == TRUE){
    
    #Calculate dissimilarity of module eigengenes
    MEDiss <- 1-cor(MEs)
    
    MEDissThres <- mergecut
    
    # Call an automatic merging function
    merge <- mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
    finalColors <- merge$colors
    
  }
  
  plotDendroAndColors(geneTree, finalColors, 
                      c("Tree cut"), 
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
  
  if(mergesimilar == TRUE){
    res <- list(tom = TOM, melist = MEList, mergedmelist = merge, dynamicColors = oricolors, 
                mergedColors = finalColors)
  }else{
    res <- list(tom = TOM, melist = MEList, dynamicColors = oricolors)
  }
  
  return(res)
}

ctrlmain <- makemain()
casemain <- makemain(datres = preedat)


#Own network organization####
ctrlcolors <- ctrlmain$mergedColors
names(ctrlcolors) <- row.names(ctrlmain$tom)

casecolors <- casemain$mergedColors
names(casecolors) <- row.names(casemain$tom)


orgnetwork <- function(TOM = ctrlmain$tom, netdata = ctrldat, colors = ctrlcolors, fixcolor = NULL, 
                       quantilecut = topquantile, abscut = 0.5, 
                       largenodesize = largevexsize, largeedgesize = largeedgsize, nodesize = 5, 
                       edgesize = 2, 
                       disfun = ctrldisfun, intnum = internum, pdfprefix = 'ctrlquant', plot = TRUE, 
                       plotonly = TRUE, calconnec = NULL){
  
  datExpr <- netdata[-1]
  
  probes <- names(datExpr)
  
  modTOM <- TOM
  modProbes <- probes
  
  dimnames(modTOM) <- list(modProbes, modProbes)
  
  if(!is.null(abscut)){
    weightcutoff <- abscut
  }else{
    library(corpcor)
    tomvalues <- sm2vec(modTOM)
    weightcutoff <- as.numeric(quantile(tomvalues, quantilecut))
  }
  
  unicolors <- unique(colors)
  modulecolors <- unicolors[unicolors != 'grey']
  
  i <- 1
  
  modulenodes <- data.frame(nodeName = character(), 
                            nodeColor = character(), 
                            stringsAsFactors = FALSE)
  for(i in 1:length(modulecolors)){
    
    modulecolor <- modulecolors[i]
    modulegenes <- names(colors)[colors == modulecolor]
    subTOM <- modTOM[modulegenes, modulegenes]
    
    if(nrow(subTOM) <= 0){
      next()
    }
    
    subcyt <- exportNetworkToCytoscape(subTOM, 
                                       edgeFile = NULL, nodeFile = NULL, 
                                       weighted = TRUE, 
                                       threshold = weightcutoff, 
                                       nodeNames = row.names(subTOM))
    
    subnodedata <- subcyt$nodeData
    subnodedata <- subnodedata[c(1)]
    if(nrow(subnodedata) == 0){
      next()
    }
    subnodedata$nodeColor <- modulecolor
    subnodedata$nodeName <- as.character(subnodedata$nodeName)
    
    modulenodes <- rbind(modulenodes, subnodedata)
    
    
  }
  
  if(is.null(edgesize)){
    edgeweight <- TRUE
  }else{
    edgeweight <- FALSE
  }
  
  cyt <- exportNetworkToCytoscape(modTOM, 
                                  edgeFile = NULL, 
                                  nodeFile = NULL, 
                                  weighted = edgeweight, 
                                  threshold = weightcutoff, 
                                  nodeNames = modProbes)
  edgedata <- cyt$edgeData
  nodedata <- cyt$nodeData
  edgedata <- edgedata[c(1, 2, 3)]
  nodedata <- nodedata[c(1)]
  
  edgedata$fromNode <- as.character(edgedata$fromNode)
  edgedata$toNode <- as.character(edgedata$toNode)
  nodedata$nodeName <- as.character(nodedata$nodeName)
  
  greynodes <- names(colors)[colors == 'grey']
  
  if(length(greynodes) > 0){
    candnodes <- unique(c(greynodes, modulenodes$nodeName))
    edgedata <- subset(edgedata, (fromNode %in% candnodes) & (toNode %in% candnodes))
    nodedata <- subset(nodedata, nodeName %in% unique(c(edgedata$fromNode, edgedata$toNode)))
  }
  
  
  
  library(igraph)
  
  inet <- graph_from_data_frame(d = edgedata, directed = FALSE, vertices = nodedata)
  
  vertexattr <- vertex.attributes(inet)
  edgeattr <- edge.attributes(inet)
  
  if(is.null(calconnec)){
    
    calconn <- function(edgedat = edgedata){
      dat1 <- edgedat
      dat2 <- edgedat[c('toNode', 'fromNode', 'weight')]
      names(dat1) <- names(dat2) <- c('node1', 'node2', 'weight')
      dat <- rbind(dat1, dat2)
      dat <- unique(dat)
      row.names(dat) <- 1:nrow(dat)
      dat <- dat[c('node1', 'weight')]
      names(dat)[1] <- 'node'
      
      library(plyr)
      
      calsum <- function(block){
        nodename <- unique(block[,1])
        numblock <- block[c(2)]
        numsum <- colSums(numblock)
        resblock <- data.frame(node = nodename, conn = numsum, stringsAsFactors = FALSE)
        row.names(resblock) <- 1:nrow(resblock)
        return(resblock)
        
      }
      
      dat <- dat[order(dat$node),]
      row.names(dat) <- 1:nrow(dat)
      conres <- ddply(.data = dat, .variables = c('node'), .fun = calsum)
      
      return(conres)
    }
    
    nodeconn <- calconn()
    connmax <- max(nodeconn$conn)
    nodeconn$normconn <- nodeconn$conn/connmax
    
    normconnvals <- nodeconn$normconn
    names(normconnvals) <- nodeconn$node
    
  }else{
    nodeconn <- calconnec
    
    normconnvals <- nodeconn$normconn
    names(normconnvals) <- nodeconn$node
  }
  
  
  
  asscolors <- function(wholelabels = colors, nodedat = nodedata){
    
    library(scales)
    
    netlabels <- wholelabels[nodedat$nodeName]
    sigcolors <- names(table(netlabels))[order(-table(netlabels))]
    sigcolors  <- sigcolors[sigcolors != 'grey']
    newcolors <- hue_pal()(length(sigcolors))
    netlabels[netlabels == 'grey'] <- '#FFFFFF'
    
    for(i in 1:length(sigcolors)){
      sigcolor <- sigcolors[i]
      netlabels[netlabels == sigcolor] <- newcolors[i]
      
    }
    
    return(netlabels)
    
  }
  
  if(!is.null(fixcolor)){
    labels <- fixcolor[nodedata$nodeName]
  }else{
    if('grey' %in% unique(colors) & length(unique(colors)) == 1){
      labels <- rep('#FFFFFF', nrow(nodedata))
    }else{
      labels <- asscolors()
    }
  }
  
  colormapping <- colors[names(labels)]
  colormapping <- data.frame(ori = colormapping, new = labels, stringsAsFactors = FALSE)
  colormapping <- unique(colormapping)
  row.names(colormapping) <- 1:nrow(colormapping)
  
  vertexattr$nodeColors <- as.vector(labels)
  vsize <- log10(normconnvals)
  
  esize <- log10(edgeattr$weight)
  
  library(scales)
  vsize <- rescale(vsize, to = c(1, largenodesize))
  
  V(inet)$size <- vsize
  
  esize <- rescale(esize, to = c(1, largeedgsize))
  
  E(inet)$width <- esize
  
  if(!is.null(nodesize)){
    V(inet)$size <- nodesize
  }
  
  if(edgeweight == FALSE){
    E(inet)$width <- edgesize
  }
  
  
  set.seed(7)
  
  if(is.null(disfun)){
    netstyle <- layout_with_fr(inet, coords = NULL, dim = 2, niter = intnum, grid = c("nogrid"))
    
    saveRDS(netstyle, paste0(pdfprefix, '_interation', intnum, '_disfun.rds'))
  }else{
    netstyle <- disfun
  }
  
  if(plot == TRUE){
    
    pdf(paste0(pdfprefix, "wgcna.pdf"), height=10, width=10, useDingbats=FALSE)
    print(
      
      plot(inet, vertex.label = NA, vertex.color = vertexattr$nodeColors, 
           edge.color = '#000000', 
           layout = netstyle)
      
      
    )
    dev.off()
    
  }
  
  
  if(plotonly == FALSE){
    res <- list(inet = inet, nodenormconn = normconnvals, nodeconn = nodeconn, 
                nodecolor = labels, edgeweight = edgeattr$weight, colormap = colormapping)
  }else{
    res <- NULL
  }
  
  return(res)
  
  
}


ctrlquantnetres <- orgnetwork(disfun = NULL, plotonly = FALSE, pdfprefix = 'ctrlquant', edgesize = NULL, 
                              abscut = absolutecut, 
                              
                              plot = FALSE)

preequantnetres <- orgnetwork(TOM = casemain$tom, netdata = preedat, colors = casecolors, 
                              disfun = NULL, pdfprefix = 'preequant', plotonly = FALSE, edgesize = NULL, 
                              abscut = absolutecut, 
                              
                              plot = FALSE)

#Export to Cytoscape####
getclass <- function(nodename = ctrlcyto$nodes$source){
  
  unknowngroupidx <- grep(pattern = 'Unknown ', x = nodename)
  normgroupidx <- setdiff(1:length(nodename), unknowngroupidx)
  
  normgroup <- nodename[normgroupidx]
  unknowngroup <- nodename[unknowngroupidx]
  
  normclasses <- gsub(pattern = ' .*$', replacement = '', x = normgroup)
  unknownclasses <- gsub(pattern = 'Unknown ', replacement = '', x = unknowngroup)
  unknownclasses <- gsub(pattern = ' .*$', replacement = '', x = unknownclasses)
  
  nodeclasses <- data.frame(source = c(normgroup, unknowngroup), 
                            classes = c(normclasses, unknownclasses), 
                            idx = c(normgroupidx, unknowngroupidx), stringsAsFactors = FALSE)
  nodeclasses <- nodeclasses[order(nodeclasses$idx),]
  nodeclasses$labels <- nodeclasses$source
  nodeclasses$labels[unknowngroupidx] <- gsub(pattern = 'Unknown ', replacement = 'u', 
                                              x = nodeclasses$labels[unknowngroupidx])
  nodeclasses$labels[unknowngroupidx] <- gsub(pattern = ';.*$', replacement = '', 
                                              x = nodeclasses$labels[unknowngroupidx])
  nodeclasses <- nodeclasses[order(nodeclasses$idx),]
  nodeclasses <- nodeclasses[-grep(pattern = 'idx', x = names(nodeclasses))]
  
  return(nodeclasses)
  
}

ctrlclass <- getclass(nodename = names(ctrlquantnetres$nodecolor))
preeclass <- getclass(nodename = names(preequantnetres$nodecolor))



classcolors <- data.frame(classes = unique(c(ctrlclass$classes, preeclass$classes)), 
                          stringsAsFactors = FALSE)
classcolors$classcolors <- hue_pal()(nrow(classcolors))

exportcyto <- function(res = ctrlquantnetres, writefile = TRUE, 
                       
                       prefix = 'Control_adj_screen', 
                       classcolor = classcolors){
  
  nodes <- data.frame(source = names(res$nodenormconn), normconn = as.vector(res$nodenormconn), 
                      color = as.vector(res$nodecolor), stringsAsFactors = FALSE)
  edges <- as_edgelist(res$inet)
  edges <- as.data.frame(edges, stringsAsFactors = FALSE)
  names(edges) <- c('source', 'target')
  edges$weight <- res$edgeweight
  
  sourcecolor <- nodes[c('source', 'color')]
  targetcolor <- sourcecolor
  names(targetcolor) <- c('target', 'targetcolor')
  names(sourcecolor) <- c('source', 'sourcecolor')
  
  edges <- merge(edges, sourcecolor, by = 'source')
  edges <- merge(edges, targetcolor, by = 'target')
  edges$finalcolor <- '#FFFFFF'
  edges$finalcolor[edges$sourcecolor == edges$targetcolor] <- 
    edges$sourcecolor[edges$sourcecolor == edges$targetcolor]
  
  nodeclass <- getclass(nodename = nodes$source)
  nodes <- merge(nodes, nodeclass, by = c('source'))
  nodes <- merge(nodes, classcolors, by = c('classes'))
  nodes <- nodes[c('source', 'normconn', 'color', 'labels', 'classes', 'classcolors')]
  
  sourceclass <- nodeclass[c('source', 'labels')]
  names(sourceclass) <- c('source', 'sourcelabel')
  targetclass <- sourceclass
  names(targetclass) <- c('target', 'targetlabel')
  edges <- merge(targetclass, edges, by = c('target'))
  edges <- merge(sourceclass, edges, by = c('source'))
  edges <- edges[c('sourcelabel', 'targetlabel', 'weight', 'finalcolor')]
  names(edges) <- c('source', 'target', 'weight', 'finalcolor')
  
  nodes <- nodes[c('labels', 'normconn', 'source', 'classes', 'classcolors', 'color')]
  names(nodes) <- c('source', 'normconn', 'oriname', 'classes', 'classcolors', 'modulecolor')
  
  cytores <- list(nodes = nodes, edges = edges)
  
  if(writefile == TRUE){
    
    write.table(nodes, paste0(prefix, '_WGCNA_nodes.txt'), sep = '\t', row.names = FALSE, 
                quote = FALSE)
    write.table(edges, paste0(prefix, '_WGCNA_edges.txt'), sep = '\t', row.names = FALSE, 
                quote = FALSE)
  }
  
  return(cytores)
  
}

ctrlcyto <- exportcyto(prefix = 'Control_adj_screen_0.8', writefile = FALSE)
preecyto <- exportcyto(res = preequantnetres, 
                       
                       prefix = 'Preeclampsia_adj_screen_0.8', writefile = FALSE)



ctrloxlipids <- ctrlcyto$nodes[grep(pattern = 'OxP', x = ctrlcyto$nodes$source),]
preeoxlipids <- preecyto$nodes[grep(pattern = 'OxP', x = preecyto$nodes$source),]

ctrltaglipids <- ctrlcyto$nodes[grep(pattern = 'TAG', x = ctrlcyto$nodes$source),]
preetaglipids <- preecyto$nodes[grep(pattern = 'TAG', x = preecyto$nodes$source),]


alloxlipids <- unique(ctrloxlipids$source, preeoxlipids$source)
alltaglipids <- unique(ctrltaglipids$source, preetaglipids$source)





#metabolite class enrichment#####
allclasses <- getclass(nodename = colnames(newdat)[2:ncol(newdat)])

colordic <- c('red', 'cyan')
names(colordic) <- c('#F8766D', '#00BFC4')

moduleenrich <- function(groupnet = preecyto$nodes, backclass = allclasses, 
                         colordict = colordic){
  
  modulenames <- unique(groupnet$modulecolor)
  
  
  classenrich <- function(module = moduleset, back = backclass){
    modulecounts <- table(module$classes)
    backcounts <- table(back$classes)
    modulesum <- sum(modulecounts)
    backsum <- sum(backcounts)
    for(j in 1:length(modulecounts)){
      classname <- names(modulecounts)[j]
      moduleclasscount <- modulecounts[classname]
      backclasscount <- backcounts[classname]
      
      a11 <- moduleclasscount
      a12 <- backclasscount
      a21 <- modulesum - moduleclasscount
      a22 <- backsum - backclasscount
      fishermat <- matrix(c(a11, a12, a21, a22), nrow = 2, byrow = TRUE)
      fisherres <- fisher.test(fishermat, alternative = 'greater')
      fisherp <- fisherres$p.value
      
      if(j == 1){
        fisherps <- fisherp
      }else{
        fisherps <- c(fisherps, fisherp)
      }
      
      
    }
    
    names(fisherps) <- names(modulecounts)
    
    return(fisherps)
    
  }
  
  mod_phe_link <- function(modcolors = moduleset, modulelogps = moduleps, textsize = 50, 
                           modulecolorname = figuremodname){
    
    modclasscolors <- modcolors[c('classes', 'classcolors')]
    modclasscolors <- unique(modclasscolors)
    names(modclasscolors)[1] <- 'classname'
    modulelogps <- merge(modulelogps, modclasscolors, by = c('classname'))
    modulelogps <- modulelogps[order(modulelogps$logp),]
    modulelogps$classname <- factor(modulelogps$classname, levels = modulelogps$classname, 
                                    ordered = TRUE)
    modulelogps$classcolors <- factor(modulelogps$classcolors, levels = modulelogps$classcolors, 
                                      ordered = TRUE)
    row.names(modulelogps) <- 1:nrow(modulelogps)
    
    library(ggplot2)
    
    p <- ggplot(modulelogps, aes(x = classname, y = logp))
    print(
      p + geom_bar(stat = 'identity', fill = modulelogps$classcolors) + 
        xlab('') + ylab('') + 
        ggtitle(paste0('Significantly enriched lipids in module ', modulecolorname)) + 
        theme_bw() + 
        theme(panel.grid = element_blank()) + 
        theme(axis.text.x = element_text(angle = 90, size = textsize)) + 
        theme(axis.text.y = element_text(size = textsize)) + 
        theme(plot.title = element_text(size = 25)) + 
        theme(panel.border = element_blank(), axis.line = element_line()) + 
        scale_y_continuous(position = 'right') + 
        coord_flip()
    )
    
  }
  
  for(i in 1:length(modulenames)){
    modulename <- modulenames[i]
    figuremodname <- as.vector(colordict[modulename])
    
    moduleset <- subset(groupnet, modulecolor == modulename)
    
    moduleps <- classenrich()
    
    moduleps <- data.frame(classname = names(moduleps), fisherp = as.vector(moduleps), 
                           stringsAsFactors = FALSE)
    moduleps$logp <- -log2(moduleps$fisherp)
    moduleps <- moduleps[order(-moduleps$logp),]
    moduleps <- subset(moduleps, fisherp < 0.05)
    
    if(nrow(moduleps) < 1){
      next()
    }
    
    mod_phe_link()
    
    print(moduleps)
  }
  
}

moduleenrich(groupnet = ctrlcyto$nodes)

colordic <- c('red', 'green', 'blue')
names(colordic) <- c('#F8766D', '#00BA38', '#619CFF')

moduleenrich(colordict = colordic)


#Relate modules of the 2 sets###
relatetopquantilemodules <- function(ctrlquantnet = ctrlquantnetres, casequantnet = preequantnetres, 
                                     titlesufix = tag, 
                                     ctrlpdcolornames = 
                                       c('red'), 
                                     preepdcolornames = 
                                       c('red', 'cyan')){
  
  vctrl <- V(ctrlquantnet$inet)
  vcase <- V(casequantnet$inet)
  allgenes <- unique(c(names(vctrl), names(vcase)))
  
  ctrlnodecolors <- ctrlquantnet$nodecolor
  casenodecolors <- casequantnet$nodecolor
  
  ctrlmodulecount <- table(ctrlnodecolors)
  casemodulecount <- table(casenodecolors)
  
  ctrlmodulecount <- ctrlmodulecount[ctrlquantnet$colormap$new]
  casemodulecount <- casemodulecount[casequantnet$colormap$new]
  
  names(ctrlmodulecount) <- ctrlpdcolornames
  names(casemodulecount) <- preepdcolornames
  
  
  ctrluninodecolors <- unique(ctrlnodecolors)
  caseuninodecolors <- unique(casenodecolors)
  
  pTable <- matrix(0, nrow = length(ctrluninodecolors), ncol = length(caseuninodecolors))
  rownames(pTable) <- ctrluninodecolors
  colnames(pTable) <- caseuninodecolors
  countTable <- pTable
  
  for(i in 1:length(ctrluninodecolors)){
    ctrluninodecolor <- ctrluninodecolors[i]
    for(j in 1:length(caseuninodecolors)){
      caseuninodecolor <- caseuninodecolors[j]
      
      ctrlset <- names(ctrlnodecolors)[ctrlnodecolors == ctrluninodecolor]
      caseset <- names(casenodecolors)[casenodecolors == caseuninodecolor]
      unionset <- union(ctrlset, caseset)
      interset <- intersect(ctrlset, caseset)
      
      mat <- c(length(setdiff(allgenes, unionset)), 
               length(setdiff(ctrlset, caseset)), 
               length(setdiff(caseset, ctrlset)), 
               length(interset))
      mat <- matrix(mat, nrow = 2)
      
      fisherp <- fisher.test(mat, alternative = "greater")$p.value
      overlap <- length(interset)
      
      pTable[i, j] <- fisherp
      countTable[i, j] <- overlap
      
    }
    
  }
  
  logpTable <- -log10(pTable)
  #Truncate p values smaller than 10^{-50} to 10^{-50}
  logpTable[is.infinite(logpTable)] <- 1.3*max(logpTable[is.finite(logpTable)]);
  logpTable[logpTable > 50 ] <- 50
  
  # Marginal counts (really module sizes)
  ctrlTotals <- apply(countTable, 1, sum)
  caseTotals <- apply(countTable, 2, sum)
  
  ctrlname <- sub(pattern = '^.*to ', replacement = '', x = titlesufix)
  casename <- sub(pattern = ' to.*$', replacement = '', x = titlesufix)
  
  row.names(logpTable) <- row.names(pTable) <- row.names(countTable) <- ctrlpdcolornames
  colnames(logpTable) <- colnames(pTable) <- colnames(countTable) <- preepdcolornames
  
  ctrlmodulecount <- ctrlmodulecount[order(-as.vector(ctrlmodulecount))]
  casemodulecount <- casemodulecount[order(-as.vector(casemodulecount))]
  
  logpTable <- logpTable[names(ctrlmodulecount), names(casemodulecount)]
  pTable <- pTable[names(ctrlmodulecount), names(casemodulecount)]
  countTable <- countTable[names(ctrlmodulecount), names(casemodulecount)]
  
  logpTable <- matrix(logpTable, nrow = length(names(ctrlmodulecount)), byrow = FALSE)
  pTable <- matrix(pTable, nrow = length(names(ctrlmodulecount)), byrow = FALSE)
  countTable <- matrix(countTable, nrow = length(names(ctrlmodulecount)), byrow = FALSE)
  
  colnames(logpTable) <- colnames(pTable) <- colnames(countTable) <- names(casemodulecount)
  rownames(logpTable) <- rownames(pTable) <- rownames(countTable) <- names(ctrlmodulecount)
  
  pTable[pTable >= 0.05] <- 1
  
  textMat <- paste(countTable, "\n(", signif(pTable, 3), ")", sep = "")
  textMat <- gsub(pattern = '\n(1)', replacement = '', x = textMat, fixed = TRUE)
  
  par(mfrow = c(1,1))
  par(cex = 1.0)
  par(mar = c(8, 12, 2.7, 1) + 0.3)
  
  labeledHeatmap(Matrix = logpTable, 
                 xLabels = paste(" ", colnames(logpTable)), 
                 yLabels = paste(" ", rownames(logpTable)),
                 colorLabels = TRUE,
                 xSymbols = paste0(casename, colnames(logpTable), ": ", as.vector(casemodulecount)),
                 ySymbols = paste0(ctrlname, rownames(logpTable), ": ", as.vector(ctrlmodulecount)),
                 textMatrix = textMat,
                 colors = greenWhiteRed(100)[50:100],
                 main = paste0('Correspondence of modules (', titlesufix, ')'),
                 cex.text = 1.5, cex.lab = 1.0, setStdMargins = FALSE)
  
  
}


relatetopquantilemodules(ctrlpdcolornames = c('red', 'cyan'), 
                         preepdcolornames = c('red', 'green', 'blue'), 
                         titlesufix = 'Case to Control')


#Calculate network summary###
summarizenet <- function(netres = ctrlquantnetres, prefix = ctrlgroupname){
  
  library(igraph)
  
  inet <- netres$inet
  clusters <- netres$nodecolor
  edgeweights <- netres$edgeweight
  
  nnodes <- length(V(inet))
  nedges <- length(E(inet))
  
  dens <- edge_density(graph = inet, loops = FALSE)
  
  numclusters <- factor(x = clusters, levels = names(table(clusters)[order(-table(clusters))]), ordered = TRUE)
  numclusters <- as.numeric(numclusters)
  modul <- modularity(x = inet, membership = numclusters, weights = edgeweights)
  
  res <- data.frame(Group = prefix, Nodes = nnodes, Edges = nedges, Density = dens, Modularity = modul, 
                    stringsAsFactors = FALSE)
  
  #calculte for each module
  edges <- as_edgelist(inet)
  colnames(edges) <- c('source', 'target')
  row.names(edges) <- 1:nrow(edges)
  edges <- as.data.frame(edges, stringsAsFactors = FALSE)
  
  uniqclusters <- unique(clusters)
  
  i <- 1
  
  for(i in 1:length(uniqclusters)){
    uniqcluster <- uniqclusters[i]
    uniqclusternodes <- names(clusters[clusters == uniqcluster])
    uniqclusteredges <- subset(edges, source %in% uniqclusternodes & target %in% uniqclusternodes)
    
    clusterinet <- graph_from_data_frame(d = uniqclusteredges, directed = FALSE)
    
    clusternnodes <- length(uniqclusternodes)
    clusternedges <- nrow(uniqclusteredges)
    
    clusterdens <- edge_density(clusterinet)
    
    clusterres <- data.frame(Group = paste0(prefix, '_', uniqcluster), 
                             Nodes = clusternnodes, 
                             Edges = clusternedges, 
                             Density = clusterdens, stringsAsFactors = FALSE)
    if(i == 1){
      clusterreses <- clusterres
    }else{
      clusterreses <- rbind(clusterreses, clusterres)
    }
    
  }
  
  reslist <- list(overallres = res, clusterres = clusterreses)
  
  #saveRDS(reslist, paste0(prefix, '_netstructure.rds'))
  
  return(reslist)
  
}

ctrlstructure <- summarizenet(prefix = 'Control')
casestructure <- summarizenet(netres = preequantnetres, prefix = 'Preeclampsia')

netstructure <- do.call(rbind, list(ctrlstructure$overallres, casestructure$overallres))
clusterstructure <- do.call(rbind, list(ctrlstructure$clusterres, casestructure$clusterres))

#write.table(netstructure, 'metabolitestructures.txt', sep = '\t', row.names = FALSE, quote = FALSE)
#write.table(clusterstructure, 'metabolitemodstructures.txt', sep = '\t', row.names = FALSE, quote = FALSE)




#Diff connect nodes#####
conndiffcutoff <- 5



ctrlnodes <- ctrlquantnetres$nodeconn$node
preenodes <- preequantnetres$nodeconn$node

ctrlsup <- setdiff(preenodes, ctrlnodes)
preesup <- setdiff(ctrlnodes, preenodes)

ctrlsup <- data.frame(node = ctrlsup, conn = 0, normconn = 0, stringsAsFactors = FALSE)
preesup <- data.frame(node = preesup, conn = 0, normconn = 0, stringsAsFactors = FALSE)

ctrlnodes <- rbind(ctrlquantnetres$nodeconn, ctrlsup)
preenodes <- rbind(preequantnetres$nodeconn, preesup)

ctrlnodes <- ctrlnodes[order(ctrlnodes$node),]
preenodes <- preenodes[order(preenodes$node),]

diffs <- preenodes$conn - ctrlnodes$conn
updiffs <- ctrlnodes$node[diffs > conndiffcutoff]
dndiffs <- ctrlnodes$node[diffs < -1*conndiffcutoff]

ctrlnodes$dircolor <- 'nc'
ctrlnodes$dircolor[ctrlnodes$node %in% updiffs] <- 'up'
ctrlnodes$dircolor[ctrlnodes$node %in% dndiffs] <- 'dn'

preenodes$dircolor <- 'nc'
preenodes$dircolor[preenodes$node %in% updiffs] <- 'up'
preenodes$dircolor[preenodes$node %in% dndiffs] <- 'dn'

names(ctrlnodes) <- names(preenodes) <- c('oriname', 'conn', 'normconn', 'dircolor')

ctrlnodes <- ctrlnodes[c('oriname', 'conn', 'dircolor')]
preenodes <- preenodes[c('oriname', 'conn', 'dircolor')]

ctrlcyto$nodes <- merge(ctrlcyto$nodes, ctrlnodes, by = c('oriname'))
preecyto$nodes <- merge(preecyto$nodes, preenodes, by = c('oriname'))

#write.table(ctrlcyto$nodes, 'Control_adj_screen_0.8_WGCNA_nodes.txt', sep = '\t', quote = FALSE, row.names = FALSE)
#write.table(preecyto$nodes, 'Preeclampsia_adj_screen_0.8_WGCNA_nodes.txt', sep = '\t', quote = FALSE, row.names = FALSE)

#Enrichment analysis#####
direnrich <- function(nodeinfo = preecyto$nodes, backinfo = preenodes, textsize = 50){
  
  modulenames <- unique(nodeinfo$modulecolor)
  backup <- sum(backinfo$dircolor == 'up')
  backdn <- sum(backinfo$dircolor == 'dn')
  
  fisherenrich <- function(a11 = a11, a12 = a12, a21 = a21, a22 = a22){
    
    mat <- matrix(c(a11, a12, a21, a22), nrow = 2, byrow = TRUE)
    
    res <- fisher.test(mat, alternative = 'greater')
    resp <- res$p.value
    
    return(resp)
  }
  
  plotenrich <- function(uppval = upp, dnpval = dnp, textsize = textsize){
    
    library(ggplot2)
    
    modulelogps <- data.frame(classname = c('DN', 'UP'), logp = -log2(c(dnpval, uppval)), 
                              classcolors = c('#00BFC4', '#F8766D'), 
                              stringsAsFactors = FALSE)
    
    p <- ggplot(modulelogps, aes(x = classname, y = logp))
    print(
      p + geom_bar(stat = 'identity', fill = modulelogps$classcolors) + 
        xlab('') + ylab('') + 
        ggtitle(paste0("Change enrichment -log2(Fisher's p-val)")) + 
        theme_bw() + 
        theme(panel.grid = element_blank()) + 
        theme(axis.text.x = element_text(angle = 90, size = textsize)) + 
        theme(axis.text.y = element_text(size = textsize)) + 
        theme(plot.title = element_text(size = 25)) + 
        theme(panel.border = element_blank(), axis.line = element_line()) + 
        scale_y_continuous(position = 'right') + 
        coord_flip()
    )
    
  }
  
  i <- 1
  for(i in 1:length(modulenames)){
    modulename <- modulenames[i]
    moduleset <- subset(nodeinfo, modulecolor == modulename)
    
    upnum <- sum(moduleset$dircolor == 'up')
    dnnum <- sum(moduleset$dircolor == 'dn')
    
    a11 <- upnum
    a12 <- backup
    a21 <- nrow(moduleset) - upnum
    a22 <- nrow(backinfo) - backup
    
    upp <- fisherenrich(a11 = a11, a12 = a12, a21 = a21, a22 = a22)
    
    a11 <- dnnum
    a12 <- backdn
    a21 <- nrow(moduleset) - dnnum
    a22 <- nrow(backinfo) - backdn
    
    dnp <- fisherenrich(a11 = a11, a12 = a12, a21 = a21, a22 = a22)
    
    plotenrich(uppval = upp, dnpval = dnp, textsize = textsize)
    
  }
  
}

direnrich(nodeinfo = ctrlcyto$nodes, backinfo = ctrlnodes)
direnrich()

classenrich <- function(nodeinfo = preecyto$nodes, backinfo = allclasses, textsize = 50){
  
  dirnames <- c('up', 'dn')
  
  classenrich <- function(module = dirset, back = backinfo){
    modulecounts <- table(module$classes)
    backcounts <- table(back$classes)
    modulesum <- sum(modulecounts)
    backsum <- sum(backcounts)
    for(j in 1:length(modulecounts)){
      classname <- names(modulecounts)[j]
      moduleclasscount <- modulecounts[classname]
      backclasscount <- backcounts[classname]
      
      a11 <- moduleclasscount
      a12 <- backclasscount
      a21 <- modulesum - moduleclasscount
      a22 <- backsum - backclasscount
      fishermat <- matrix(c(a11, a12, a21, a22), nrow = 2, byrow = TRUE)
      fisherres <- fisher.test(fishermat, alternative = 'greater')
      fisherp <- fisherres$p.value
      
      if(j == 1){
        fisherps <- fisherp
      }else{
        fisherps <- c(fisherps, fisherp)
      }
      
      
    }
    
    names(fisherps) <- names(modulecounts)
    
    return(fisherps)
    
  }
  
  mod_phe_link <- function(modcolors = dirset, modulelogps = dirps, textsize = textsize, 
                           dirname = dirname){
    
    modclasscolors <- modcolors[c('classes', 'classcolors')]
    modclasscolors <- unique(modclasscolors)
    names(modclasscolors)[1] <- 'classname'
    modulelogps <- merge(modulelogps, modclasscolors, by = c('classname'))
    modulelogps <- modulelogps[order(modulelogps$logp),]
    modulelogps$classname <- factor(modulelogps$classname, levels = modulelogps$classname, 
                                    ordered = TRUE)
    modulelogps$classcolors <- factor(modulelogps$classcolors, levels = modulelogps$classcolors, 
                                      ordered = TRUE)
    row.names(modulelogps) <- 1:nrow(modulelogps)
    
    library(ggplot2)
    
    p <- ggplot(modulelogps, aes(x = classname, y = logp))
    print(
      p + geom_bar(stat = 'identity', fill = modulelogps$classcolors) + 
        xlab('') + ylab('') + 
        ggtitle(paste0('Significantly enriched classes in ', toupper(dirname))) + 
        theme_bw() + 
        theme(panel.grid = element_blank()) + 
        theme(axis.text.x = element_text(angle = 90, size = textsize)) + 
        theme(axis.text.y = element_text(size = textsize)) + 
        theme(plot.title = element_text(size = 25)) + 
        theme(panel.border = element_blank(), axis.line = element_line()) + 
        scale_y_continuous(position = 'right') + 
        coord_flip()
    )
    
  }
  
  dirpslist <- list()
  i <- 1
  for(i in 1:length(dirnames)){
    dirname <- dirnames[i]
    dirset <- subset(nodeinfo, dircolor == dirname)
    
    dirps <- classenrich(module = dirset, back = backinfo)
    
    dirps <- data.frame(classname = names(dirps), fisherp = as.vector(dirps), 
                        stringsAsFactors = FALSE)
    dirps$logp <- -log2(dirps$fisherp)
    dirps <- dirps[order(-dirps$logp),]
    dirps <- subset(dirps, fisherp < 0.05)
    
    if(nrow(dirps) < 1){
      next()
    }
    
    mod_phe_link(modcolors = dirset, modulelogps = dirps, textsize = textsize, dirname = dirname)
    
    dirpslist[[dirname]] <- dirps
    
  }
  
  return(dirpslist)
}

upnodes <- data.frame(oriname = updiffs, dircolor = 'up', stringsAsFactors = FALSE)
dnnodes <- data.frame(oriname = dndiffs, dircolor = 'dn', stringsAsFactors = FALSE)
diffnodes <- rbind(upnodes, dnnodes)

nodesinfo <- rbind(ctrlcyto$nodes, preecyto$nodes)
nodesinfo <- nodesinfo[c('source', 'oriname', 'classes', 'classcolors', 'dircolor')]
nodesinfo <- unique(nodesinfo)
diffnodes <- subset(nodesinfo, dircolor != 'nc')


diffenrichlist <- classenrich(nodeinfo = diffnodes)




#OxPL heatmap##########

diffdat <- newdat[,c('Label', diffnodes$oriname)]

diffdat <- diffdat[order(diffdat$Label),]

names(diffs) <- ctrlnodes$oriname


drawheatmap <- function(datamat = newdat, title = 'Metabolites', whetherrowname = FALSE, 
                        clustermethod = 'average', clustercol = TRUE, clusterrow = TRUE, 
                        casename = 'Preeclampsia', ctrlname = 'Control', 
                        fontsizerow = 10){
  
  library(pheatmap)
  
  datanno <- data.frame(group = datamat$Label, stringsAsFactors = FALSE)
  row.names(datanno) <- row.names(datamat)
  datanno$group <- factor(datanno$group, levels = c(casename, ctrlname), ordered = TRUE)
  datamat <- t(datamat[-1])
  
  library(matrixStats)
  
  rowvar <- rowVars(datamat)
  removerows <- rowvar < 0.01 
  datamat <- datamat[!removerows,]
  
  print(
    pheatmap(datamat, annotation_col = datanno, 
             color = colorRampPalette(colors = c('green', 'black', 'red'))(100), 
             show_colnames = FALSE, 
             main = title, 
             scale = 'row', show_rownames = whetherrowname, 
             cluster_rows = clusterrow, cluster_cols = clustercol, 
             clustering_method = clustermethod, fontsize_row = fontsizerow)
  )
  
}


drawheatmap(datamat = diffdat, title = 'Diff metabolites', whetherrowname = TRUE, fontsizerow = 5)

drawheatmap(datamat = diffdat, title = 'Diff metabolites', whetherrowname = TRUE, fontsizerow = 5, 
            clustercol = FALSE)


drawheatmapconn <- function(datadat = diffdat, title = 'Connectivity Diff Lipids', whetherrowname = TRUE, 
                            clustermethod = 'complete', clustercol = FALSE, clusterrow = FALSE, 
                            casename = 'Preeclampsia', ctrlname = 'Control', 
                            dirinf = diffs, clusterfeature = TRUE, gaps = FALSE, 
                            fontsizerow = 5, clusternum = 3){
  
  library(pheatmap)
  
  dirinfo <- dirinf[colnames(datadat)[-1]]
  dirinfo <- dirinfo[order(-dirinfo)]
  dirinfo <- data.frame(lipids = names(dirinfo), diff = dirinfo, stringsAsFactors = FALSE)
  dirinfo$dir <- 'NC'
  dirinfo$dir[dirinfo$diff > 0] <- 'UP'
  dirinfo$dir[dirinfo$diff < 0] <- 'DN'
  diranno <- dirinfo[c('dir')]
  
  datanno <- data.frame(group = datadat$Label, stringsAsFactors = FALSE)
  row.names(datanno) <- row.names(datadat)
  datanno$group <- factor(datanno$group, levels = c(casename, ctrlname), ordered = TRUE)
  
  datamat <- t(datadat[-1])
  
  library(matrixStats)
  
  rowvar <- rowVars(datamat)
  removerows <- rowvar < 0.01 
  datamat <- datamat[!removerows,]
  
  ctrlmat <- datamat[,row.names(subset(datanno, group == ctrlname))]
  preemat <- datamat[,row.names(subset(datanno, group == casename))]
  
  sampleorder <- function(sampledat = ctrlmat, clustermethod = clustermethod){
    
    library(pheatmap)
    
    hcc <- hclust(dist(t(sampledat)), method = clustermethod)
    sampleorder <- hcc$order
    
    sampleseq <- colnames(sampledat)[sampleorder]
    return(sampleseq)
    
  }
  
  ctrlseq <- sampleorder(clustermethod = clustermethod)
  preeseq <- sampleorder(sampledat = preemat, clustermethod = clustermethod)
  sampleseq <- c(ctrlseq, preeseq)
  
  featurecluster <- function(featuredat = datamat, clustermethod = clustermethod, 
                             clusternum = clusternum){
    
    library(pheatmap)
    
    hcc <- hclust(dist(featuredat), method = clustermethod)
    featureorder <- hcc$order
    featurecolor <- cutree(hcc, clusternum)
    mode(featurecolor) <- 'character'
    featureseq <- row.names(featuredat)[featureorder]
    
    featuregroup <- data.frame(Cluster = factor(featurecolor, labels = paste0('Cluster', seq(clusternum))))
    
    featureres <- list(featureseq = featureseq, featuregroup = featuregroup)
    
    return(featureres)
    
  }
  
  featureseq <- row.names(diranno)
  featuregroup <- diranno
  names(featuregroup) <- 'Cluster'
  featuregroup$Cluster <- factor(featuregroup$Cluster, levels = c('UP', 'DN'), ordered = TRUE)
  
  if(clusterfeature == TRUE){
    
    featureres <- featurecluster(clustermethod = clustermethod, 
                                 clusternum = clusternum)
    featureseq <- featureres$featureseq
    featuregroup <- featureres$featuregroup
    
  }
  
  
  datmat <- datamat[featureseq, sampleseq]
  
  if(gaps == TRUE){
    
    featuregroups <- featuregroup[row.names(datmat),]
    datannos <- datanno[colnames(datmat),]
    
    rowgaps <- as.numeric(featuregroups)
    colgaps <- as.numeric(datannos)
    
    rowgap <- unique(rowgaps)
    colgap <- unique(colgaps)
    
    rowcuts <- c()
    colcuts <- c()
    
    for(i in 1:(length(rowgap) - 1)){
      rownums <- rowgap[1:i]
      rowcutnum <- sum(rowgaps %in% rownums)
      rowcuts <- c(rowcuts, rowcutnum)
      
    }
    
    for(j in 1:(length(colgap) - 1)){
      colnums <- colgap[1:j]
      colcutnum <- sum(colgaps %in% colnums)
      colcuts <- c(colcuts, colcutnum)
    }
    
  }else{
    rowcuts <- NULL
    colcuts <- NULL
  }
  
  print(
    pheatmap(datmat, annotation_col = datanno, annotation_row = featuregroup, 
             color = colorRampPalette(colors = c('green', 'black', 'red'))(100), 
             show_colnames = FALSE, 
             main = title, 
             scale = 'row', show_rownames = whetherrowname, 
             fontsize = fontsizerow, 
             cluster_rows = clusterrow, cluster_cols = clustercol, 
             gaps_row = rowcuts, gaps_col = colcuts, 
             border_color = NA)
  )
  
}


drawheatmapconn(gaps = TRUE, fontsizerow = 5, clusterfeature = FALSE)

drawheatmapconn(gaps = TRUE, fontsizerow = 5, clusterfeature = TRUE, clusternum = 10, clustermethod = 'average')

drawheatmapconn(gaps = TRUE, fontsizerow = 5, clusterfeature = TRUE, clusternum = 10)


drawheatmapconnclassmerge <- function(datadat = diffdat, title = 'Connectivity Diff Lipids', 
                                      whetherrowname = TRUE, 
                                      clustercol = FALSE, 
                                      casename = 'Case', ctrlname = 'Control', 
                                      dirinf = diffs, clusterfeature = TRUE, gaps = TRUE, 
                                      fontsizerow = 5, 
                                      clustermethod = 'complete', 
                                      pddat = pd, addpdbar = FALSE){
  
  library(pheatmap)
  
  colnames(datadat) <- gsub(pattern = '\\(', replacement = ' ', x = colnames(datadat))
  colnames(datadat) <- gsub(pattern = '\\)', replacement = '', x = colnames(datadat))
  
  names(dirinf) <- gsub(pattern = '\\(', replacement = ' ', x = names(dirinf))
  names(dirinf) <- gsub(pattern = '\\)', replacement = '', x = names(dirinf))
  
  dirinfo <- dirinf[colnames(datadat)[-1]]
  dirinfo <- dirinfo[order(-dirinfo)]
  dirinfo <- data.frame(lipids = names(dirinfo), diff = dirinfo, stringsAsFactors = FALSE)
  dirinfo$dir <- 'NC'
  dirinfo$dir[dirinfo$diff > 0] <- 'UP'
  dirinfo$dir[dirinfo$diff < 0] <- 'DN'
  diranno <- dirinfo[c('dir')]
  
  datanno <- data.frame(group = datadat$Label, stringsAsFactors = FALSE)
  row.names(datanno) <- row.names(datadat)
  datanno$group <- factor(datanno$group, levels = c(casename, ctrlname), ordered = TRUE)
  
  datamat <- t(datadat[-1])
  
  library(matrixStats)
  
  rowvar <- rowVars(datamat)
  removerows <- rowvar < 0.01 
  datamat <- datamat[!removerows,]
  
  ctrlmat <- datamat[,row.names(subset(datanno, group == ctrlname))]
  preemat <- datamat[,row.names(subset(datanno, group == casename))]
  
  sampleorder <- function(sampledat = ctrlmat, clustermethod = clustermethod){
    
    library(pheatmap)
    
    hcc <- hclust(dist(t(sampledat)), method = clustermethod)
    sampleorder <- hcc$order
    
    sampleseq <- colnames(sampledat)[sampleorder]
    return(sampleseq)
    
  }
  
  ctrlseq <- sampleorder(clustermethod = clustermethod)
  preeseq <- sampleorder(sampledat = preemat, clustermethod = clustermethod)
  sampleseq <- c(ctrlseq, preeseq)
  
  featureseq <- row.names(diranno)
  featuregroup <- diranno
  names(featuregroup) <- 'Cluster'
  featuregroup$Cluster <- factor(featuregroup$Cluster, levels = c('UP', 'DN'), ordered = TRUE)
  
  datmat <- datamat[featureseq, sampleseq]
  
  if(gaps == TRUE){
    
    featuregroups <- featuregroup[row.names(datmat),]
    datannos <- datanno[colnames(datmat),]
    
    rowgaps <- as.numeric(featuregroups)
    colgaps <- as.numeric(datannos)
    
    rowgap <- unique(rowgaps)
    colgap <- unique(colgaps)
    
    rowcuts <- c()
    colcuts <- c()
    
    for(i in 1:(length(rowgap) - 1)){
      rownums <- rowgap[1:i]
      rowcutnum <- sum(rowgaps %in% rownums)
      rowcuts <- c(rowcuts, rowcutnum)
      
    }
    
    for(j in 1:(length(colgap) - 1)){
      colnums <- colgap[1:j]
      colcutnum <- sum(colgaps %in% colnums)
      colcuts <- c(colcuts, colcutnum)
    }
    
  }else{
    rowcuts <- NULL
    colcuts <- NULL
  }
  
  
  mergeclass <- function(mat = uppart){
    
    lipidnames <- row.names(mat)
    newlipidnames <- gsub(pattern = '\\(', replacement = ' ', x = lipidnames)
    newlipidnames <- gsub(pattern = '\\)', replacement = '', x = lipidnames)
    row.names(mat) <- lipidnames
    lipidnames <- lipidnames[order(lipidnames)]
    
    
    unknowns <- lipidnames[grepl(pattern = 'Unknown ', x = lipidnames)]
    others <- lipidnames[!grepl(pattern = 'Unknown', x = lipidnames)]
    prefixes <- gsub(pattern = ' .*$', replacement = '', x = others)
    prefixorders <- names(table(prefixes)[order(-table(prefixes))])
    
    for(i in 1:length(prefixorders)){
      prefix <- prefixorders[i]
      prefixpart <- others[grepl(pattern = paste0('^', prefix), x = others)]
      if(i == 1){
        prefixparts <- prefixpart
      }else{
        prefixparts <- c(prefixparts, prefixpart)
      }
    }
    
    newlipidnames <- c(prefixparts, unknowns)
    newmat <- mat[newlipidnames,]
    return(newmat)
  }
  
  if(sum(diranno$dir == 'UP') > 0){
    uppart <- datmat[diranno$dir == 'UP',]
    uppart <- mergeclass(uppart)
  }
  
  if(sum(diranno$dir == 'DN') > 0){
    dnpart <- datmat[diranno$dir == 'DN',]
    dnpart <- mergeclass(dnpart)
  }
  
  if(sum(diranno$dir == 'UP') <= 0 | sum(diranno$dir == 'DN') <= 0){
    newmat <- datmat
    newmat <- mergeclass(newmat)
  }else{
    newmat <- rbind(uppart, dnpart)
  }
  
  print(
    pheatmap(newmat, annotation_col = datanno, annotation_row = featuregroup, 
             color = colorRampPalette(colors = c('green', 'black', 'red'))(100), 
             show_colnames = FALSE, 
             #main = title, 
             scale = 'row', show_rownames = whetherrowname, 
             fontsize = 4, 
             cluster_rows = FALSE, cluster_cols = clustercol, 
             gaps_row = rowcuts, gaps_col = colcuts, 
             border_color = NA, 
             cellheight = 4)
  )
  
  if(addpdbar == TRUE){
    
    sampleanno <- pddat
    names(sampleanno) <- c('Preeclampsia', 'Smoker', 'Gender', 'Ethgroup', 'Age', 
                           'Parity', 'BMI', 'GestationalAgeWeeks', 'GestationalDiabetes', 'ChronicHypertension', 
                           'MembraneRupture', 'PlacentalAbruption', 'NeonatalMalformation')
    sampleanno <- sampleanno[c('Preeclampsia', 'GestationalAgeWeeks', 'GestationalDiabetes', 
                               'ChronicHypertension', 'Smoker', 
                               'MembraneRupture', 'Ethgroup', 'NeonatalMalformation', 
                               'PlacentalAbruption', 'Gender', 
                               'Age', 'BMI', 'Parity')]
    selvars <- c('Preeclampsia', 'GestationalAgeWeeks', 'GestationalDiabetes', 'ChronicHypertension', 'Smoker', 
                 'Gender')
    #sampleanno <- sampleanno[,selvars]
    
    print(
      pheatmap(newmat, annotation_col = sampleanno, annotation_row = featuregroup, 
               color = colorRampPalette(colors = c('green', 'black', 'red'))(100), 
               show_colnames = FALSE, 
               main = title, 
               scale = 'row', show_rownames = whetherrowname, 
               fontsize = fontsizerow, 
               cluster_rows = FALSE, cluster_cols = clustercol, 
               gaps_row = rowcuts, gaps_col = colcuts, 
               border_color = NA, 
               annotation_legend = TRUE)
    )
    
    
  }
  
  
}


drawheatmapconnclassmerge(datadat = diffdat, clusterfeature = FALSE, 
                          casename = 'Preeclampsia')







wilcoxdiffs <- subset(wilcoxdatres$reslines, lipid %in% colnames(wilcoxdatres$wilcoxdat))
wilcoxdifflipids <- wilcoxdiffs$lipid

wilcoxdiffs <- wilcoxdiffs$dismean - wilcoxdiffs$ctrmean
names(wilcoxdiffs) <- wilcoxdifflipids

drawheatmapconn(datadat = wilcoxdatres$wilcoxdat, title = 'Wilcox Diff Lipids', 
                dirinf = wilcoxdiffs, clusterfeature = FALSE, gaps = TRUE, 
                casename = 'Case', ctrlname = 'Control')

drawheatmapconnclassmerge(datadat = wilcoxdatres$wilcoxdat, title = 'Wilcox Diff Lipids', 
                          dirinf = wilcoxdiffs, gaps = TRUE, clusterfeature = FALSE, 
                          casename = 'Case', ctrlname = 'Control', fontsizerow = 10)






diabeteswilcoxdiffs <- subset(diabeteswilcoxdatres$reslines, lipid %in% colnames(diabeteswilcoxdatres$wilcoxdat))
diabeteswilcoxdifflipids <- diabeteswilcoxdiffs$lipid

diabeteswilcoxdiffs <- diabeteswilcoxdiffs$dismean - diabeteswilcoxdiffs$ctrmean
names(diabeteswilcoxdiffs) <- diabeteswilcoxdifflipids

drawheatmapconn(datadat = diabeteswilcoxdatres$wilcoxdat, title = 'Wilcox Diff Lipids', 
                dirinf = diabeteswilcoxdiffs, clusterfeature = FALSE, gaps = TRUE, 
                casename = 'Case', ctrlname = 'Control')

drawheatmapconnclassmerge(datadat = diabeteswilcoxdatres$wilcoxdat, title = 'Wilcox Diff Lipids', 
                          dirinf = diabeteswilcoxdiffs, gaps = TRUE, clusterfeature = FALSE, 
                          casename = 'Case', ctrlname = 'Control', fontsizerow = 10)





#Pathway mapping#####

wilcoxtops <- wilcoxdifflipids

#cyanmodule <- names(preequantnetres$nodecolor)[preequantnetres$nodecolor == '#00BFC4']

#Get mz value#####
newdat <- tmpnewdat



posfile <- 'posraw.txt'
negfile <- 'negraw.txt'


getmzrt <- function(posdatfile = posfile, negdatfile = negfile, fulldatfile = 'own.txt', 
                    normdatval = newdat, writefile = FALSE, filename = 'unknownmz.txt'){
  
  posdat <- read.table(posdatfile, sep = '\t', header = TRUE, stringsAsFactors = FALSE, comment.char = '')
  negdat <- read.table(negdatfile, sep = '\t', header = TRUE, stringsAsFactors = FALSE, comment.char = '')
  
  fullval <- read.table(fulldatfile, sep = '\t', header = TRUE, 
                        stringsAsFactors = FALSE, check.names = FALSE)
  unknownval <- fullval[grepl(pattern = 'Unknown', x = fullval$compound_name),]
  unknownval <- unknownval[,c(2, 3, 4, 5)]
  
  posdat <- posdat[,c(1, 2 ,3 ,4, 5, 6)]
  negdat <- negdat[,c(1, 2, 3, 4, 5, 6)]
  totaldat <- rbind(posdat, negdat)
  totaldat <- unique(totaldat)
  
  unknownnames <- gsub(pattern = 'Unknown ', x = unknownval$compound_name, replacement = '')
  #unknownnames <- gsub(pattern = '@.*$', replacement = '', x = unknownnames)
  mzvals <- subset(totaldat, Sample.Name %in% unknownnames)
  
  unknownnames <- paste0('Unknown ', mzvals$Sample.Name)
  rts <- mzvals$Expected.RT
  mzs <- (mzvals$Start.Mass...1 + mzvals$End.Mass...1)/2
  
  normdatval <- normdatval[-1]
  normdatval <- t(normdatval)
  normdatval <- as.data.frame(normdatval)
  normdatval <- normdatval[unknownnames,]
  unknowns <- data.frame(compound = unknownnames, feature = paste0(mzs, '_', rts), stringsAsFactors = FALSE)
  rownames(unknowns) <- unknowns$compound
  unknowns <- unknowns[rownames(normdatval),]
  unknowns <- as.data.frame(unknowns)
  unknownnorms <- cbind(unknowns, normdatval)
  row.names(unknownnorms) <- NULL
  
  if(writefile == TRUE){
    
    write.table(unknownnorms, filename, sep = '\t', row.names = FALSE, quote = FALSE)
  }
  
  return(unknownnorms)
  
}




generatemzfile <- function(posdatfile = posfile, negdatfile = negfile, querynames = bootclinic, 
                           writefile = FALSE, prefix = 'bootclinic'){
  
  posdat <- read.table(posdatfile, sep = '\t', header = TRUE, stringsAsFactors = FALSE, comment.char = '')
  negdat <- read.table(negdatfile, sep = '\t', header = TRUE, stringsAsFactors = FALSE, comment.char = '')
  
  convertnames <- function(orinames = posdat$Name){
    newnames <- gsub(pattern = ';@.*$', replacement = '', x = orinames)
    newnames <- gsub(pattern = '; \\[.*$', replacement = '', x = newnames)
    newnames <- gsub(pattern = ';\\[.*$', replacement = '', x = newnames)
    return(newnames)
  }
  
  posnewnames <- convertnames()
  negnewnames <- convertnames(orinames = negdat$Name)
  posdatnames <- names(posdat)
  negdatnames <- names(negdat)
  posdat$newnames <- posnewnames
  negdat$newnames <- negnewnames
  posdat <- posdat[c('newnames', posdatnames)]
  negdat <- negdat[c('newnames', negdatnames)]
  
  uquerynamesidx <- grep(pattern = 'Unknown', x = querynames)
  if(length(uquerynamesidx) > 0){
    uquerynames <- querynames[uquerynamesidx]
    querynames <- querynames[-uquerynamesidx]
  }
  
  getqueryinfo <- function(query = querynames, dat = posdat){
    
    querydat <- subset(dat, newnames %in% query)
    querydat <- querydat[c('newnames', 'Name', 'Start.Mass...1', 'End.Mass...1')]
    querydat$MassMean <- (querydat$Start.Mass...1 + querydat$End.Mass...1)/2
    querydat$MassDiff <- abs(querydat$Start.Mass...1 - querydat$End.Mass...1)
    return(querydat)
  }
  
  posres <- getqueryinfo()
  negres <- getqueryinfo(dat = negdat)
  pnres <- rbind(posres, negres)
  pnres <- pnres[order(pnres$newnames, pnres$Name),]
  
  if(length(uquerynamesidx) > 0){
    
    newuquerynames <- convertnames(orinames = uquerynames)
    newuquerynames <- gsub(pattern = 'Unknown ', replacement = '', x = newuquerynames)
    uposres <- getqueryinfo(query = newuquerynames)
    unegres <- getqueryinfo(query = newuquerynames, dat = negdat)
    upnres <- rbind(uposres, unegres)
    upnres <- upnres[order(upnres$newnames, upnres$Name),]
    
    uextracts <- data.frame(newnames = character(), Name = character(), Start.Mass...1 = character(), 
                            End.Mass...1 = character(), MassMean = character(), 
                            MassDiff = character(), stringsAsFactors = FALSE)
    for(i in 1:length(uquerynames)){
      uqueryname <- uquerynames[i]
      upat <- sub(pattern = 'Unknown.*;', replacement = '', x = uqueryname)
      upat <- sub(pattern = ' @', replacement = '', x = upat)
      uidx <- grep(pattern = upat, x = upnres$Name, fixed = TRUE)
      if(length(uidx) > 0){
        uextract <- upnres[uidx,]
        uextract$newnames <- uqueryname
        uextracts <- rbind(uextracts, uextract)
      }else{
        next()
      }
    }
    
    if(nrow(uextracts) > 0){
      pnres <- rbind(pnres, uextracts)
    }
    
  }
  
  row.names(pnres) <- 1:nrow(pnres)
  
  if(writefile == TRUE){
    filename <- paste0(prefix, 'mz.txt')
    write.table(x = pnres, file = filename, quote = FALSE, sep = '\t', row.names = FALSE)
    
  }
  
  return(pnres)
}




wilcoxtopinfo <- generatemzfile(querynames = wilcoxtops, writefile = TRUE, prefix = 'wilcoxtops')

#cyanmoduleinfo <- generatemzfile(querynames = cyanmodule, writefile = TRUE, prefix = 'cyanmodule')


#Classify lipids to LIPID MAPS#####
classifylipids <- function(orilipids = lilikoiqueryinfo, write = TRUE, prefix = 'lilikoiquery'){
  
  mzs <- orilipids$MassMean
  adducts <- orilipids$Name
  
  rts <- gsub(pattern = '^.*@', replacement = '', x = adducts)
  
  adducts <- gsub(pattern = '^.*;', replacement = '', x = adducts)
  adducts <- gsub(pattern = '@.*$', replacement = '', x = adducts)
  adducts <- gsub(pattern = '^ \\[', replacement = '[', x = adducts)
  adducts[adducts == '[M+Hac-H]-'] <- '[M+HCOO]-'
  polarities <- rep('Positive', length(adducts))
  polarities[adducts == ''] <- 'Neutral'
  polarities[grep(pattern = '\\]-$', x = adducts)] <- 'Negative'
  
  categories <- orilipids$newnames
  categories <- gsub(pattern = '^Unknown ', replacement = '', x = categories)
  categories <- gsub(pattern = ' .*$', replacement = '', x = categories)
  
  inputlipids <- data.frame(newnames = orilipids$newnames, 
                            Name = orilipids$Name, 
                            rts = rts, 
                            mzs = mzs, 
                            polarities = polarities, 
                            adducts = adducts, 
                            categories = categories, stringsAsFactors = FALSE)
  
  poslipids <- subset(inputlipids, polarities == 'Positive')
  neglipids <- subset(inputlipids, polarities == 'Negative')
  neulipids <- subset(inputlipids, polarities == 'Neutral')
  
  res <- list()
  res$poslipids <- poslipids
  res$neglipids <- neglipids
  res$neulipids <- neulipids
  
  if(write == TRUE){
    
    write.table(poslipids, paste0(prefix, '_poslipids.txt'), sep = '\t', row.names = FALSE, 
                quote = FALSE)
    write.table(neglipids, paste0(prefix, '_neglipids.txt'), sep = '\t', row.names = FALSE, 
                quote = FALSE)
    write.table(neulipids, paste0(prefix, '_neulipids.txt'), sep = '\t', row.names = FALSE, 
                quote = FALSE)
    
  }
  
  return(res)
  
}

wilcoxtopres <- classifylipids(orilipids = wilcoxtopinfo, write = TRUE, prefix = 'wilcoxtops')

#Get qualified positive/negative/neural lipids#####

filterres <- function(lipidres = lilikoiqueryres$poslipids, match1file = 'lilikoiquery_pos_ms_hits.txt', 
                      write = TRUE){
  
  match1res <- read.table(match1file, sep = '\t', header = TRUE, stringsAsFactors = FALSE, row.names = NULL)
  match1res <- subset(match1res, LMSD.Examples == 'Examples')
  match1res <- match1res[-ncol(match1res)]
  names(match1res) <- c('mzs', 'Matched.Mass', 'Delta', 'newnames', 'Formula', 'adducts')
  match1res <- unique(match1res)
  row.names(match1res) <- 1:nrow(match1res)
  match1res$rownames <- row.names(match1res)
  
  exactmatches <- merge(match1res, lipidres, by = c('newnames', 'adducts', 'mzs'))
  exactmatches <- unique(exactmatches)
  
  match1res2 <- subset(match1res, !(rownames %in% exactmatches$rownames))
  match1res2 <- unique(match1res2)
  lipidres2 <- subset(lipidres, !(Name %in% exactmatches$Name))
  lipidres2 <- unique(lipidres2)
  
  possmatches <- merge(match1res2, lipidres2, by = c('mzs', 'adducts'))
  possmatches <- unique(possmatches)
  
  mapsyms <- possmatches$newnames.x
  mssyms <- possmatches$newnames.y
  
  mapcats <- gsub(pattern = ' .*$', replacement = '', x = mapsyms)
  mscats <- gsub(pattern = '^Unknown ', replacement = '', x = mssyms)
  mscats <- gsub(pattern = ' .*$', replacement = '', x = mscats)
  mscats[mscats == 'TAG'] <- 'TG'
  mscats <- gsub(pattern = '^Ox', replacement = '', x = mscats)
  mscats[grep(pattern = '^Cer-.*', x = mscats)] <- 'Cer'
  
  mapcars <- gsub(pattern = '^.* ', replacement = '', x = mapsyms)
  mapcars <- gsub(pattern = ';.*$', replacement = '', x = mapcars)
  
  mscars <- gsub(pattern = ';.*$', replacement = '', x = mssyms)
  mscars <- gsub(pattern = '.* ', replacement = '', x = mscars)
  mscars <- gsub(pattern = '[a-z]', replacement = '', x = mscars)
  
  exactmatches$newnames.y <- exactmatches$newnames
  
  sels <- (mapcats == mscats) & (mapcars == mscars)
  exactmatches2 <- possmatches[sels,]
  names(exactmatches2)[grep(pattern = 'newnames.x', x = names(exactmatches2))] <- 'newnames'
  exactmatches2 <- exactmatches2[names(exactmatches)]
  
  res <- rbind(exactmatches, exactmatches2)
  res <- unique(res)
  
  if(write == TRUE){
    
    prefix <- gsub(pattern = '\\.txt', replacement = '', x = match1file)
    
    write.table(res, paste0(prefix, '_qualified.txt'), sep = '\t', row.names = FALSE, quote = FALSE)
    
  }
  
  
  return(res)
  
  
}




matchres1 <- filterres(lipidres = wilcoxtopres$poslipids, match1file = 'wilcoxtops_pos_ms_hits.txt', 
                       write = TRUE)

matchres1 <- filterres(lipidres = wilcoxtopres$neglipids, match1file = 'wilcoxtops_neg_ms_hits.txt')



#Organize lipid map results#####
getmappedlipids <- function(targetdir = './pos_lipids/', matchnames = unique(matchres1$newnames)){
  
  curdir <- getwd()
  setwd(targetdir)
  lipidfiles <- dir()
  
  if(length(lipidfiles) != length(matchnames)){
    
    return(NULL)
    
  }
  
  lipidfilenums <- gsub(pattern = 'lmsd_examples', replacement = '', x = lipidfiles)
  lipidfilenums <- gsub(pattern = '\\.txt', replacement = '', x = lipidfilenums)
  lipidfilenums <- as.numeric(lipidfilenums)
  lipidfilenums <- lipidfilenums[order(lipidfilenums)]
  lipidfiles <- paste0('lmsd_examples', lipidfilenums, '.txt')
  
  for(i in 1:length(lipidfiles)){
    
    oriname <- matchnames[i]
    
    lipidfile <- lipidfiles[i]
    lipidcont <- read.table(lipidfile, sep = '\t', header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, 
                            fill = TRUE)
    lipidcont <- lipidcont[-ncol(lipidcont)]
    lipidcommonname <- lipidcont$`Common Name`
    lipidcommonname <- unique(lipidcommonname)
    lipidcommonname <- gsub(pattern = '\\[.*\\]$', replacement = '', x = lipidcommonname)
    
    lipidcommonname <- data.frame(oriname = oriname, commonname = lipidcommonname, 
                                  stringsAsFactors = FALSE)
    
    if(i == 1){
      lipidcommonnames <- lipidcommonname
    }else{
      lipidcommonnames <- rbind(lipidcommonnames, lipidcommonname)
    }
    
  }
  
  lipidcommonnames <- unique(lipidcommonnames)
  
  setwd(curdir)
  
  return(lipidcommonnames)
  
}


mappedlipids <- getmappedlipids(targetdir = './wilcoxpos_lipids/')



namemapping <- matchres1[c('newnames', 'newnames.y')]
names(namemapping) <- c('oriname', 'msname')
mappedlipids <- merge(mappedlipids, namemapping, by = c('oriname'))




mappedlipidspos <- mappedlipids

mappedlipidsneg <- mappedlipids

mappedlipids <- rbind(mappedlipidspos, mappedlipidsneg)

mappedlipids <- unique(mappedlipids)




#Use lilikoi#####
library(lilikoi)

metabolitePathwayTable_pos <- lilikoi.metab_to_pathway(mappedlipidspos$commonname, "name")

names(mappedlipidspos)[2] <- 'Query'
metabolitePathwayTable_pos <- merge(mappedlipidspos, metabolitePathwayTable_pos, by = c('Query'))



metabolitePathwayTable_neg <- lilikoi.metab_to_pathway(mappedlipidsneg$commonname, 'name')

names(mappedlipidsneg)[2] <- 'Query'
metabolitePathwayTable_neg <- merge(mappedlipidsneg, metabolitePathwayTable_neg, by = c('Query'))





metabolitePathwayTable_final <- rbind(metabolitePathwayTable_pos, metabolitePathwayTable_neg)
metabolitePathwayTable_final <- unique(metabolitePathwayTable_final)

simpathtab <- metabolitePathwayTable_final[c('Query', 'oriname', 'msname', 
                                             'HMDB', 'PubChem', 'KEGG', 'pathway')]

remainings <- (!is.na(simpathtab$HMDB)) | (!is.na(simpathtab$PubChem)) | (!is.na(simpathtab$KEGG)) | 
  (!is.na(simpathtab$pathway))

remainingtab <- simpathtab[remainings,]



#Wilcox pathways###

hmdbpathes <- data.frame(hmdbids = c(rep('HMDB0011769', 3), 
                                     rep('HMDB0004947', 3), 
                                     rep('HMDB0004956', 3), 
                                     rep('HMDB0007869', 2), 
                                     rep('HMDB0007871', 2), 
                                     rep('HMDB0007921', 4), 
                                     rep('HMDB0007934', 2), 
                                     rep('HMDB0007965', 4), 
                                     rep('HMDB0000564', 4), 
                                     rep('HMDB0007980', 4), 
                                     rep('HMDB0008081', 4), 
                                     rep('HMDB0008012', 4), 
                                     rep('HMDB0008026', 4), 
                                     rep('HMDB0008031', 4), 
                                     rep('HMDB0008040', 4), 
                                     rep('HMDB0008041', 4), 
                                     rep('HMDB0008053', 4), 
                                     rep('HMDB0008072', 4), 
                                     rep('HMDB0008085', 4), 
                                     rep('HMDB0008105', 4), 
                                     rep('HMDB0008118', 4), 
                                     rep('HMDB0008136', 4), 
                                     rep('HMDB0008137', 4), 
                                     rep('HMDB0008150', 4), 
                                     rep('HMDB0008168', 4), 
                                     rep('HMDB0008201', 4), 
                                     rep('HMDB0008276', 4), 
                                     rep('HMDB0008308', 4), 
                                     rep('HMDB0008332', 4), 
                                     rep('HMDB0008340', 4),
                                     rep('HMDB0008364', 4), 
                                     rep('HMDB0008397', 4), 
                                     rep('HMDB0008533', 4), 
                                     rep('HMDB0011202', 4), 
                                     rep('HMDB0008564', 4), 
                                     rep('HMDB0008591', 4), 
                                     rep('HMDB0008595', 4), 
                                     rep('HMDB0008792', 4), 
                                     rep('HMDB0010574', 3), 
                                     rep('HMDB0010587', 3), 
                                     rep('HMDB0010630', 3), 
                                     rep('HMDB0012089', 2), 
                                     rep('HMDB0001348', 7), 
                                     rep('HMDB0005376', 2), 
                                     rep('HMDB0005368', 2), 
                                     rep('HMDB0005381', 2), 
                                     rep('HMDB0005432', 2), 
                                     rep('HMDB0005422', 2), 
                                     rep('HMDB0005395', 2), 
                                     rep('HMDB0005407', 2), 
                                     rep('HMDB0010466', 2), 
                                     rep('HMDB0005450', 2), 
                                     rep('HMDB0005457', 2), 
                                     rep('HMDB0005467', 2), 
                                     rep('HMDB0005468', 2), 
                                     rep('HMDB0011487', 2), 
                                     rep('HMDB0011488', 2), 
                                     rep('HMDB0008909', 4), 
                                     rep('HMDB0011507', 2), 
                                     rep('HMDB0011517', 2), 
                                     rep('HMDB0011518', 2), 
                                     rep('HMDB0009549', 4)
                                     
                                     
), 

hmdbpathes = c('Insulin signaling pathway', 'Apoptosis', 'Lipid metabolism pathway', 
               
               'Insulin signaling pathway', 'Apoptosis', 'Lipid metabolism pathway', 
               
               'Insulin signaling pathway', 'Apoptosis', 'Lipid metabolism pathway', 
               
               'Glycerophospholipid metabolism', 'Lipid metabolism pathway', 
               
               'Glycerophospholipid metabolism', 'Lipid metabolism pathway', 
               
               'phosphatidylcholine biosynthesis', 
               'phosphatidylethanolamine biosynthesis', 
               'Glycerophospholipid metabolism', 
               'Lipid metabolism pathway', 
               
               'Glycerophospholipid metabolism', 'Lipid metabolism pathway', 
               
               'phosphatidylcholine biosynthesis', 
               'phosphatidylethanolamine biosynthesis', 
               'Glycerophospholipid metabolism', 
               'Lipid metabolism pathway', 
               
               'phosphatidylcholine biosynthesis', 
               'phosphatidylethanolamine biosynthesis', 
               'Glycerophospholipid metabolism', 
               'Lipid metabolism pathway', 
               
               'phosphatidylcholine biosynthesis', 
               'phosphatidylethanolamine biosynthesis', 
               'Glycerophospholipid metabolism', 
               'Lipid metabolism pathway', 
               
               'phosphatidylcholine biosynthesis', 
               'phosphatidylethanolamine biosynthesis', 
               'Glycerophospholipid metabolism', 
               'Lipid metabolism pathway', 
               
               'phosphatidylcholine biosynthesis', 
               'phosphatidylethanolamine biosynthesis', 
               'Glycerophospholipid metabolism', 
               'Lipid metabolism pathway', 
               
               'phosphatidylcholine biosynthesis', 
               'phosphatidylethanolamine biosynthesis', 
               'Glycerophospholipid metabolism', 
               'Lipid metabolism pathway', 
               
               'phosphatidylcholine biosynthesis', 
               'phosphatidylethanolamine biosynthesis', 
               'Glycerophospholipid metabolism', 
               'Lipid metabolism pathway', 
               
               'phosphatidylcholine biosynthesis', 
               'phosphatidylethanolamine biosynthesis', 
               'Glycerophospholipid metabolism', 
               'Lipid metabolism pathway', 
               
               'phosphatidylcholine biosynthesis', 
               'phosphatidylethanolamine biosynthesis', 
               'Glycerophospholipid metabolism', 
               'Lipid metabolism pathway', 
               
               'phosphatidylcholine biosynthesis', 
               'phosphatidylethanolamine biosynthesis', 
               'Glycerophospholipid metabolism', 
               'Lipid metabolism pathway', 
               
               'phosphatidylcholine biosynthesis', 
               'phosphatidylethanolamine biosynthesis', 
               'Glycerophospholipid metabolism', 
               'Lipid metabolism pathway', 
               
               'phosphatidylcholine biosynthesis', 
               'phosphatidylethanolamine biosynthesis', 
               'Glycerophospholipid metabolism', 
               'Lipid metabolism pathway', 
               
               'phosphatidylcholine biosynthesis', 
               'phosphatidylethanolamine biosynthesis', 
               'Glycerophospholipid metabolism', 
               'Lipid metabolism pathway', 
               
               'phosphatidylcholine biosynthesis', 
               'phosphatidylethanolamine biosynthesis', 
               'Glycerophospholipid metabolism', 
               'Lipid metabolism pathway', 
               
               'phosphatidylcholine biosynthesis', 
               'phosphatidylethanolamine biosynthesis', 
               'Glycerophospholipid metabolism', 
               'Lipid metabolism pathway', 
               
               'phosphatidylcholine biosynthesis', 
               'phosphatidylethanolamine biosynthesis', 
               'Glycerophospholipid metabolism', 
               'Lipid metabolism pathway', 
               
               'phosphatidylcholine biosynthesis', 
               'phosphatidylethanolamine biosynthesis', 
               'Glycerophospholipid metabolism', 
               'Lipid metabolism pathway', 
               
               'phosphatidylcholine biosynthesis', 
               'phosphatidylethanolamine biosynthesis', 
               'Glycerophospholipid metabolism', 
               'Lipid metabolism pathway', 
               
               'phosphatidylcholine biosynthesis', 
               'phosphatidylethanolamine biosynthesis', 
               'Glycerophospholipid metabolism', 
               'Lipid metabolism pathway', 
               
               'phosphatidylcholine biosynthesis', 
               'phosphatidylethanolamine biosynthesis', 
               'Glycerophospholipid metabolism', 
               'Lipid metabolism pathway', 
               
               'phosphatidylcholine biosynthesis', 
               'phosphatidylethanolamine biosynthesis', 
               'Glycerophospholipid metabolism', 
               'Lipid metabolism pathway', 
               
               'phosphatidylcholine biosynthesis', 
               'phosphatidylethanolamine biosynthesis', 
               'Glycerophospholipid metabolism', 
               'Lipid metabolism pathway', 
               
               'phosphatidylcholine biosynthesis', 
               'phosphatidylethanolamine biosynthesis', 
               'Glycerophospholipid metabolism', 
               'Lipid metabolism pathway', 
               
               'phosphatidylcholine biosynthesis', 
               'phosphatidylethanolamine biosynthesis', 
               'Glycerophospholipid metabolism', 
               'Lipid metabolism pathway', 
               
               'phosphatidylcholine biosynthesis', 
               'phosphatidylethanolamine biosynthesis', 
               'Glycerophospholipid metabolism', 
               'Lipid metabolism pathway', 
               
               'phosphatidylcholine biosynthesis', 
               'phosphatidylethanolamine biosynthesis', 
               'Glycerophospholipid metabolism', 
               'Lipid metabolism pathway', 
               
               'phosphatidylcholine biosynthesis', 
               'phosphatidylethanolamine biosynthesis', 
               'Glycerophospholipid metabolism', 
               'Lipid metabolism pathway', 
               
               'phosphatidylcholine biosynthesis', 
               'phosphatidylethanolamine biosynthesis', 
               'Glycerophospholipid metabolism', 
               'Lipid metabolism pathway', 
               
               'phosphatidylcholine biosynthesis', 
               'phosphatidylethanolamine biosynthesis', 
               'Glycerophospholipid metabolism', 
               'Lipid metabolism pathway', 
               
               'phosphatidylcholine biosynthesis', 
               'phosphatidylethanolamine biosynthesis', 
               'Glycerophospholipid metabolism', 
               'Lipid metabolism pathway',
               
               'phosphatidylcholine biosynthesis', 
               'phosphatidylethanolamine biosynthesis', 
               'Glycerophospholipid metabolism', 
               'Lipid metabolism pathway', 
               
               'Cardiolipin Biosynthesis', 
               'Glycerophospholipid metabolism', 
               'Lipid metabolism pathway', 
               
               'Cardiolipin Biosynthesis', 
               'Glycerophospholipid metabolism', 
               'Lipid metabolism pathway', 
               
               'Cardiolipin Biosynthesis', 
               'Glycerophospholipid metabolism', 
               'Lipid metabolism pathway', 
               
               'Apoptosis', 'Lipid metabolism pathway', 
               
               'Globoid Cell Leukodystrophy', 
               'Metachromatic Leukodystrophy (MLD)', 
               'Sphingolipid Metabolism', 
               'Apoptosis', 
               'Lipid metabolism pathway', 
               'Gaucher Disease', 
               'Fabry Disease', 
               
               'De Novo Triacylglycerol Biosynthesis', 
               'Lipid metabolism pathway', 
               
               'De Novo Triacylglycerol Biosynthesis', 
               'Lipid metabolism pathway', 
               
               'De Novo Triacylglycerol Biosynthesis', 
               'Lipid metabolism pathway', 
               
               'De Novo Triacylglycerol Biosynthesis', 
               'Lipid metabolism pathway', 
               
               'De Novo Triacylglycerol Biosynthesis', 
               'Lipid metabolism pathway', 
               
               'De Novo Triacylglycerol Biosynthesis', 
               'Lipid metabolism pathway', 
               
               'De Novo Triacylglycerol Biosynthesis', 
               'Lipid metabolism pathway', 
               
               'De Novo Triacylglycerol Biosynthesis', 
               'Lipid metabolism pathway', 
               
               'De Novo Triacylglycerol Biosynthesis', 
               'Lipid metabolism pathway', 
               
               'De Novo Triacylglycerol Biosynthesis', 
               'Lipid metabolism pathway', 
               
               'De Novo Triacylglycerol Biosynthesis', 
               'Lipid metabolism pathway', 
               
               'De Novo Triacylglycerol Biosynthesis', 
               'Lipid metabolism pathway', 
               
               'Glycerophospholipid metabolism', 
               'Lipid metabolism pathway', 
               
               'Glycerophospholipid metabolism', 
               'Lipid metabolism pathway', 
               
               'phosphatidylcholine biosynthesis', 
               'phosphatidylethanolamine biosynthesis', 
               'Glycerophospholipid metabolism', 
               'Lipid metabolism pathway', 
               
               'Glycerophospholipid metabolism', 
               'Lipid metabolism pathway', 
               
               'Glycerophospholipid metabolism', 
               'Lipid metabolism pathway', 
               
               'Glycerophospholipid metabolism', 
               'Lipid metabolism pathway', 
               
               'phosphatidylcholine biosynthesis', 
               'phosphatidylethanolamine biosynthesis', 
               'Glycerophospholipid metabolism', 
               'Lipid metabolism pathway'
               
               
), 

stringsAsFactors = FALSE)

hmdbpathes <- unique(hmdbpathes)








pubchempathes <- data.frame(pubchemids = c(rep('5383576', 1), 
                                           rep('5283571', 1), 
                                           rep('129657', 3), 
                                           rep('131150', 3), 
                                           rep('52922294', 2), 
                                           rep('24778654', 3), 
                                           rep('24778679', 2), 
                                           rep('452110', 9), 
                                           rep('24778715', 3), 
                                           rep('24778719', 3), 
                                           rep('52922461', 2), 
                                           rep('53478691', 2), 
                                           rep('3082163', 2), 
                                           rep('52922655', 3), 
                                           rep('52922657', 3), 
                                           rep('52922667', 2), 
                                           rep('53478723', 3), 
                                           rep('53478745', 2), 
                                           rep('24778939', 4), 
                                           rep('53478769', 2), 
                                           rep('53478785', 2), 
                                           rep('52922727', 3), 
                                           rep('52922751', 3), 
                                           rep('52922781', 2), 
                                           rep('52922843', 2), 
                                           rep('24779046', 2), 
                                           rep('24779063', 3), 
                                           rep('52923163', 2), 
                                           rep('52923187', 2), 
                                           rep('53478949', 2), 
                                           rep('52923221', 2), 
                                           rep('52923453', 2), 
                                           rep('53480677', 2), 
                                           rep('53479177', 2), 
                                           rep('52923553', 2), 
                                           rep('52923569', 2), 
                                           rep('53479495', 2), 
                                           rep('52941750', 2), 
                                           rep('52926476', 2), 
                                           rep('52927225', 2), 
                                           rep('5283588', 6), 
                                           rep('9543987', 1), 
                                           rep('9544170', 2), 
                                           rep('9544166', 2), 
                                           rep('9543989', 1), 
                                           rep('9544167', 1), 
                                           rep('16058371', 2), 
                                           rep('9544680', 2), 
                                           rep('25240385', 2), 
                                           rep('9544681', 2), 
                                           rep('9544748', 2), 
                                           rep('25240382', 2), 
                                           rep('9544749', 2), 
                                           rep('52924174', 2), 
                                           rep('42607465', 2), 
                                           rep('52924776', 2)), 
                            
                            pubchempathes = c('Sphingolipids metabolism pathway', 
                                              
                                              'Sphingolipids metabolism pathway', 
                                              
                                              'Phosphatidylcholine Biosynthesis', 
                                              'Phosphatidylcholine/Phosphatidylethanolamine Biosynthesis', 
                                              'Phosphatidylethanolamine Biosynthesis', 
                                              
                                              'Phosphatidylcholine Biosynthesis', 
                                              'Phosphatidylcholine/Phosphatidylethanolamine Biosynthesis', 
                                              'Phosphatidylethanolamine Biosynthesis', 
                                              
                                              'Phosphatidylcholine Biosynthesis', 
                                              'Phosphatidylethanolamine Biosynthesis', 
                                              
                                              'Phosphatidylcholine Biosynthesis', 
                                              'Phosphatidylcholine/Phosphatidylethanolamine Biosynthesis', 
                                              'Phosphatidylethanolamine Biosynthesis', 
                                              
                                              'Phosphatidylcholine Biosynthesis', 
                                              'Phosphatidylethanolamine Biosynthesis', 
                                              
                                              'Phosphatidylcholine Biosynthesis', 
                                              'Phosphatidylcholine/Phosphatidylethanolamine Biosynthesis', 
                                              'Phosphatidylethanolamine Biosynthesis', 
                                              'Acetylcholine Synthesis', 
                                              'Choline Metabolism', 
                                              'Lysolipid Incorporation into ER', 
                                              'Lysolipid Incorporation into Mitochondira', 
                                              'Phospholipid Biosynthesis', 
                                              'Phospholipid Remodeling (phosphatidylcholine)', 
                                              
                                              'Phosphatidylcholine Biosynthesis', 
                                              'Phosphatidylcholine/Phosphatidylethanolamine Biosynthesis', 
                                              'Phosphatidylethanolamine Biosynthesis', 
                                              
                                              'Phosphatidylcholine Biosynthesis', 
                                              'Phosphatidylcholine/Phosphatidylethanolamine Biosynthesis', 
                                              'Phosphatidylethanolamine Biosynthesis', 
                                              
                                              'Phosphatidylcholine Biosynthesis', 
                                              'Phosphatidylethanolamine Biosynthesis', 
                                              
                                              'Phosphatidylcholine Biosynthesis', 
                                              'Phosphatidylethanolamine Biosynthesis', 
                                              
                                              'Phosphatidylcholine Biosynthesis', 
                                              'Phosphatidylethanolamine Biosynthesis', 
                                              
                                              'Phosphatidylcholine Biosynthesis', 
                                              'Phosphatidylcholine/Phosphatidylethanolamine Biosynthesis', 
                                              'Phosphatidylethanolamine Biosynthesis', 
                                              
                                              'Phosphatidylcholine Biosynthesis', 
                                              'Phosphatidylcholine/Phosphatidylethanolamine Biosynthesis', 
                                              'Phosphatidylethanolamine Biosynthesis', 
                                              
                                              'Phosphatidylcholine Biosynthesis', 
                                              'Phosphatidylethanolamine Biosynthesis', 
                                              
                                              'Phosphatidylcholine Biosynthesis', 
                                              'Phosphatidylcholine/Phosphatidylethanolamine Biosynthesis', 
                                              'Phosphatidylethanolamine Biosynthesis', 
                                              
                                              'Phosphatidylcholine Biosynthesis', 
                                              'Phosphatidylethanolamine Biosynthesis', 
                                              
                                              'Phosphatidylcholine Biosynthesis', 
                                              'Phosphatidylcholine/Phosphatidylethanolamine Biosynthesis', 
                                              'Phosphatidylethanolamine Biosynthesis', 
                                              'Phospholipid desaturation', 
                                              
                                              'Phosphatidylcholine Biosynthesis', 
                                              'Phosphatidylethanolamine Biosynthesis', 
                                              
                                              'Phosphatidylcholine Biosynthesis', 
                                              'Phosphatidylethanolamine Biosynthesis', 
                                              
                                              'Phosphatidylcholine Biosynthesis', 
                                              'Phosphatidylethanolamine Biosynthesis', 
                                              'Phospholipid desaturation', 
                                              
                                              'Phosphatidylcholine Biosynthesis', 
                                              'Phosphatidylcholine/Phosphatidylethanolamine Biosynthesis', 
                                              'Phosphatidylethanolamine Biosynthesis', 
                                              
                                              'Phosphatidylcholine Biosynthesis', 
                                              'Phosphatidylethanolamine Biosynthesis', 
                                              
                                              'Phosphatidylcholine Biosynthesis', 
                                              'Phosphatidylethanolamine Biosynthesis', 
                                              
                                              'Phosphatidylcholine Biosynthesis', 
                                              'Phosphatidylethanolamine Biosynthesis', 
                                              
                                              'Phosphatidylcholine Biosynthesis', 
                                              'Phosphatidylethanolamine Biosynthesis', 
                                              'Phospholipid desaturation', 
                                              
                                              'Phosphatidylcholine Biosynthesis', 
                                              'Phosphatidylethanolamine Biosynthesis', 
                                              
                                              'Phosphatidylcholine Biosynthesis', 
                                              'Phosphatidylethanolamine Biosynthesis', 
                                              
                                              'Phosphatidylcholine Biosynthesis', 
                                              'Phosphatidylethanolamine Biosynthesis', 
                                              
                                              'Phosphatidylcholine Biosynthesis', 
                                              'Phosphatidylethanolamine Biosynthesis', 
                                              
                                              'Phosphatidylcholine Biosynthesis', 
                                              'Phosphatidylethanolamine Biosynthesis', 
                                              
                                              'Phosphatidylcholine Biosynthesis', 
                                              'Phosphatidylethanolamine Biosynthesis', 
                                              
                                              'Phosphatidylcholine Biosynthesis', 
                                              'Phosphatidylethanolamine Biosynthesis', 
                                              
                                              'Phosphatidylcholine Biosynthesis', 
                                              'Phosphatidylethanolamine Biosynthesis', 
                                              
                                              'Phosphatidylcholine Biosynthesis', 
                                              'Phosphatidylethanolamine Biosynthesis', 
                                              
                                              'Phosphatidylcholine Biosynthesis', 
                                              'Phosphatidylethanolamine Biosynthesis', 
                                              
                                              'Cardiolipin Biosynthesis', 
                                              'Phospholipid Biosynthesis', 
                                              
                                              'Cardiolipin Biosynthesis', 
                                              'Phospholipid Biosynthesis', 
                                              
                                              'Cardiolipin Biosynthesis', 
                                              'Phospholipid Biosynthesis', 
                                              
                                              'Fabry Disease', 'Gaucher Disease', 
                                              'Globoid Cell Leukodystrophy', 'Krabbe Disease', 
                                              'Metachromatic Leukodystrophy (MLD)', 
                                              'Sphingolipid Metabolism', 
                                              
                                              'De Novo Triacylglycerol Biosynthesis', 
                                              
                                              'De Novo Triacylglycerol Biosynthesis', 
                                              'Triacylglycerol Degradation', 
                                              
                                              'De Novo Triacylglycerol Biosynthesis', 
                                              'Triacylglycerol Degradation', 
                                              
                                              'De Novo Triacylglycerol Biosynthesis', 
                                              
                                              'De Novo Triacylglycerol Biosynthesis', 
                                              
                                              'De Novo Triacylglycerol Biosynthesis', 
                                              'Triacylglycerol Degradation', 
                                              
                                              'De Novo Triacylglycerol Biosynthesis', 
                                              'Triacylglycerol Degradation', 
                                              
                                              'De Novo Triacylglycerol Biosynthesis', 
                                              'Triacylglycerol Degradation', 
                                              
                                              'De Novo Triacylglycerol Biosynthesis', 
                                              'Triacylglycerol Degradation', 
                                              
                                              'De Novo Triacylglycerol Biosynthesis', 
                                              'Triacylglycerol Degradation', 
                                              
                                              'De Novo Triacylglycerol Biosynthesis', 
                                              'Triacylglycerol Degradation', 
                                              
                                              'De Novo Triacylglycerol Biosynthesis', 
                                              'Triacylglycerol Degradation', 
                                              
                                              'Phosphatidylcholine Biosynthesis', 
                                              'Phosphatidylethanolamine Biosynthesis', 
                                              
                                              'Lysolipid Incorporation into ER', 
                                              'Lysolipid Incorporation into Mitochondira', 
                                              
                                              'Phosphatidylcholine Biosynthesis', 
                                              'Phosphatidylethanolamine Biosynthesis'), 
                            
                            stringsAsFactors = FALSE)


pubchempathes <- unique(pubchempathes)




keggpathes <- data.frame(keggids = c(rep('C00195', 9), 
                                     rep('C00157', 8), 
                                     rep('C00550', 3), 
                                     rep('C00422', 8), 
                                     rep('C00350', 8)), 
                         
                         keggpathes = c('Sphingolipid metabolism', 'Metabolic pathways', 
                                        'Sphingolipid signaling pathway', 'Necroptosis', 
                                        'Neurotrophin signaling pathway', 'Adipocytokine signaling pathway', 
                                        'Insulin resistance', 
                                        'AGE-RAGE signaling pathway in diabetic complications', 
                                        'Leishmaniasis', 
                                        
                                        'Glycerophospholipid metabolism', 'Arachidonic acid metabolism', 
                                        'Linoleic acid metabolism', 'alpha-Linolenic acid metabolism', 
                                        'Metabolic pathways', 'Biosynthesis of secondary metabolites', 
                                        'Retrograde endocannabinoid signaling', 'Choline metabolism in cancer', 
                                        
                                        'Sphingolipid metabolism', 
                                        'Sphingolipid signaling pathway', 
                                        'Necroptosis', 
                                        
                                        'Glycerolipid metabolism', 'Metabolic pathways', 
                                        'Thermogenesis', 'Regulation of lipolysis in adipocytes', 
                                        'Insulin resistance', 'Fat digestion and absorption', 
                                        'Vitamin digestion and absorption', 'Cholesterol metabolism', 
                                        
                                        'Glycosylphosphatidylinositol (GPI)-anchor biosynthesis', 
                                        'Glycerophospholipid metabolism', 'Metabolic pathways', 
                                        'Biosynthesis of secondary metabolites', 'Autophagy', 
                                        'Retrograde endocannabinoid signaling', 
                                        'Pathogenic Escherichia coli infection', 
                                        'Kaposi sarcoma-associated herpesvirus infection'
                         ), 
                         
                         stringsAsFactors = FALSE)

keggpathes <- unique(keggpathes)





names(remainingtab) <- c('Query', 'oriname', 'msname', 'hmdbids', 'pubchemids', 'keggids', 'pathway')

sharedpathes_hmdb_pubchem <- intersect(hmdbpathes$hmdbpathes, pubchempathes$pubchempathes)
sharedpathes_hmdb_kegg <- intersect(hmdbpathes$hmdbpathes, keggpathes$keggpathes)
#sharedpathes_pubchem_kegg <- intersect(pubchempathes$pubchempathes, keggpathes$keggpathes)

basicpart <- remainingtab[c('msname', 'hmdbids', 'pubchemids', 'keggids')]
names(basicpart)[1] <- 'oriname'

hmdbpart <- merge(basicpart, hmdbpathes, by = c('hmdbids'))
pubchempart <- merge(basicpart, pubchempathes, by = c('pubchemids'))
keggpart <- merge(basicpart, keggpathes,  by = c('keggids'))

hmdbpart$source <- 'HMDB'
pubchempart$source <- 'PubChem'
keggpart$source <- 'KEGG'

hmdbpart <- hmdbpart[c('oriname', 'hmdbpathes', 'source')]
pubchempart <- pubchempart[c('oriname', 'pubchempathes', 'source')]
keggpart <- keggpart[c('oriname', 'keggpathes', 'source')]

names(hmdbpart)[2] <- names(pubchempart)[2] <- names(keggpart)[2] <- 'pathes'
pathes <- rbind(hmdbpart, pubchempart, keggpart)
pathes$source[pathes$pathes %in% sharedpathes_hmdb_pubchem] <- 'HMDB & PubChem'
pathes$source[pathes$pathes %in% sharedpathes_hmdb_kegg] <- 'HMDB & KEGG'
pathes <- unique(pathes)


#Wilcox results#####

nodedirs <- subset(wilcoxdatres$reslines, lipid %in% wilcoxtops)
nodedirs$logFC <- nodedirs$dismean - nodedirs$ctrmean
nodedirs <- data.frame(oriname = wilcoxtops, logFC = nodedirs$logFC, stringsAsFactors = FALSE)
nodedirs$dirs <- 'DN'
nodedirs$dirs[nodedirs$logFC > 0] <- 'UP'
pathes <- merge(pathes, nodedirs, by = c('oriname'))

pathdirs <- pathes[c('pathes', 'logFC')]

library(plyr)


removelargepathes <- function(pathres = pathes){
  
  largepathlist <- c('Metabolic pathways', 'Lipid metabolism pathway')
  
  pathres <- subset(pathres, !(pathes %in% largepathlist))
  pathres <- unique(pathres)
  return(pathres)
}

pathes <- removelargepathes(pathres = pathes)

mergesimilarpathes <- function(pathres = pathes){
  
  unipathes <- unique(pathres$pathes)
  unipathes <- unipathes[order(unipathes)]
  
  simpleCap <- function(x) {
    s <- strsplit(x, " ")[[1]]
    paste(toupper(substring(s, 1, 1)), substring(s, 2), sep = "", collapse = " ")
  }
  
  newpathes <- sapply(unipathes, simpleCap)
  newpathes <- as.vector(newpathes)
  comppathes <- data.frame(prepath = unipathes, newpath = newpathes, stringsAsFactors = FALSE)
  changepathes <- subset(comppathes, prepath != newpath)
  
  i <- 1
  for(i in 1:length(changepathes$prepath)){
    prename <- changepathes$prepath[i]
    newname <- changepathes$newpath[i]
    pathres$pathes[pathres$pathes == prename] <- newname
  }
  
  pathres <- unique(pathres)
  
  return(pathres)
  
}

pathes <- mergesimilarpathes(pathres = pathes)
pathes$pathes[pathes$pathes == 'Sphingolipids Metabolism Pathway'] <- 'Sphingolipid Metabolism'



pathdirs <- pathes[c('pathes', 'logFC')]

library(plyr)

mergefun <- function(set){
  
  logfcs <- set$logFC
  logfcs <- sum(logfcs)
  return(logfcs)
  
}


pathdirmerge <- ddply(.data = pathdirs, .variables = c('pathes'), .fun = mergefun)
names(pathdirmerge) <- c('pathes', 'pathlogFC')

nodeinfo <- pathes[c('oriname', 'dirs', 'logFC')]
nodeinfo <- unique(nodeinfo)
pathinfo <- pathdirmerge
pathinfo$dirs <- 'DN'
pathinfo$dirs[pathinfo$pathlogFC > 0] <- 'UP'
pathinfo <- unique(pathinfo)
pathinfo <- pathinfo[c('pathes', 'dirs', 'pathlogFC')]
names(pathinfo) <- names(nodeinfo)
cytoscapenodes <- rbind(nodeinfo, pathinfo)
cytoscapenodes <- unique(cytoscapenodes)


cytoscapeedges <- pathes
cytoscapeedges <- cytoscapeedges[c('oriname', 'pathes')]
names(cytoscapeedges) <- c('lipids', 'pathes')
cytoscapeedges <- unique(cytoscapeedges)

cytoscapenodes$source <- 'path'
cytoscapenodes$source[cytoscapenodes$oriname %in% wilcoxtops] <- 'node'

cytoscapenodes$abslogFC <- abs(cytoscapenodes$logFC)

#write.table(cytoscapenodes, 'wilcoxcytonodes_3pathdatabases.txt', sep = '\t', quote = FALSE, row.names = FALSE)
#write.table(cytoscapeedges, 'wilcoxcytoedges_3pathdatabases.txt', sep = '\t', quote = FALSE, row.names = FALSE)






#Focus on 1 lipid specie#####
focuses <- significantPathways1

vars <- c('Sample_Group', 'Gestational_Diabetes')



newdat <- as.data.frame(newdat)

#plotdat <- cbind(newdat[focus], pd[vars])


diabetesanovaplot <- function(subdat = focusdat, xlabname = var, textsize = 20, dir = 'greater'){
  
  focusname <- names(subdat)[1]
  names(subdat)[1] <- 'PC 32:0'
  
  library(ggplot2)
  
  anovres <- Anova(lm(subdat$`PC 32:0`~subdat$Sample_Group + subdat$var), type = 3)
  pvals <- anovres$`Pr(>F)`
  pvs <- signif(pvals, 3)
  
  samplegroupp <- pvs[2]
  varp <- pvs[3]
  
  subdat$group <- NA
  subdat$group[subdat$var == 'N' & subdat$Sample_Group == 'Control'] <- 'noPE_noGD'
  subdat$group[subdat$var == 'N' & subdat$Sample_Group == 'Preeclampsia'] <- 'PE_noGD'
  subdat$group[subdat$var == 'Y' & subdat$Sample_Group == 'Preeclampsia'] <- 'PE_GD'
  subdat$group <- factor(subdat$group, levels = c('noPE_noGD', 'PE_noGD', 'PE_GD'), ordered = TRUE)
  
  uniqgroups <- levels(subdat$group)
  for(i in 1:(length(uniqgroups) - 1)){
    uniqgroup1 <- uniqgroups[i]
    sub1 <- subset(subdat, group == uniqgroup1)
    for(j in (i+1):length(uniqgroups)){
      uniqgroup2 <- uniqgroups[j]
      sub2 <- subset(subdat, group == uniqgroup2)
      
      tp <- t.test(sub1$`PC 32:0`, sub2$`PC 32:0`)
      tp <- tp$p.value
      if(i == 1 & j == 2){
        tps <- tp
      }else{
        tps <- c(tps, tp)
      }
    }
    
  }
  names(tps) <- paste0('tp', seq(1, length(tps)))
  tps <- signif(tps, 3)
  
  
  library(ggpubr)
  
  mycomparisions <- list(c('noPE_noGD', 'PE_noGD'), c('PE_noGD', 'PE_GD'), c('noPE_noGD', 'PE_GD'))
  symnumargs <- list(cutpoints = c(0, 0.01, 0.05, 1), symbols = c("**", "*", ""))
  methodargs <- list(alternative = dir)
  
  p <- ggplot(subdat, aes(x = group, y = `PC 32:0`, fill = group))
  
  p <- p + geom_boxplot(position = 'dodge') + 
    xlab('Group') + ylab(focusname) + 
    scale_fill_manual(values = hue_pal()(length(unique(subdat$group))), name = xlabname) + 
    
    stat_compare_means(comparisons = mycomparisions, symnum.args = symnumargs, method = 't.test', size = 7, 
                       method.args = methodargs)

  p <- p + 
    theme_bw() + 
    theme(axis.text.x = element_text(size = textsize), 
          axis.title.x = element_text(size = textsize), 
          axis.title.y = element_text(size = textsize)) + 
    theme(legend.position = 'none')
  
  
  print(p)
  
  anovaline <- data.frame(preep = pvs[2], varp = pvs[3], stringsAsFactors = FALSE)
  
  return(anovaline)
  
}


library(scales)



for(n in 1:length(focuses)){
  focus <- focuses[n]
  plotdat <- cbind(newdat[focus], pd[vars])
  
  for(i in 2:length(vars)){
    var <- vars[i]
    focusdat <- plotdat[c(focus, var, 'Sample_Group')]
    if(names(focusdat)[2] == 'Neonatal_Malformations_Baby'){
      focusdat$Neonatal_Malformations_Baby <- as.character(focusdat$Neonatal_Malformations_Baby)
      focusdat$Neonatal_Malformations_Baby[focusdat$Neonatal_Malformations_Baby == 'N'] <- 'NO'
      focusdat$Neonatal_Malformations_Baby[focusdat$Neonatal_Malformations_Baby == 'Y'] <- 'YES'
      focusdat$Neonatal_Malformations_Baby <- factor(focusdat$Neonatal_Malformations_Baby, 
                                                     levels = rev(c('NO', 'YES')), ordered = TRUE)
    }
    
    
    if(var == 'Ethgroup'){
      focusdat$Ethgroup <- factor(as.character(focusdat$Ethgroup), 
                                  levels = c('Caucasian', 'Pacific island', 'Asian', 'Latin'), 
                                  ordered = TRUE)
    }
    names(focusdat)[2] <- 'var'
    
    if(var == 'GestationalAgeWeeks'){
      var <- 'Gestational Age (Weeks)'
    }else if(var == 'Ethgroup'){
      var <- 'Ethnicity'
    }else if(var == 'Neonatal_Malformations_Baby'){
      var <- 'Neonatal_Malformations'
    }
    
    if(var != 'Gestational Age (Weeks)'){
      
      anovaresline <- diabetesanovaplot(textsize = 20)
    }else{
      varmed <- median(focusdat$var)
      focusdat$vardis <- c('High')
      focusdat$vardis[focusdat$var < varmed] <- c('Low')
      focusdat$vardis <- factor(focusdat$vardis, levels = c('Low', 'High'), ordered = TRUE)
      focusdat <- focusdat[c(focus, 'vardis', 'Sample_Group')]
      names(focusdat)[2] <- 'var'
      drawanovadisplot(textsize = 20)
    }
    
  }
  
  if(n == 1){
    anovareslines <- anovaresline
  }else{
    anovareslines <- rbind(anovareslines, anovaresline)
  }
}

anovareslines$lipid <- focuses
