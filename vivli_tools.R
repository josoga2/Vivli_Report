#setwd
setwd("/Users/josoga2/Documents/wale_docs/phd/server/tubin/data/vivli/")
#import main data
#vivli.data <- read.csv('2024_05_28 atlas_antibiotics.csv', header = T)
load(file = './interaction.mat.RData')
load(file = './vivli_data.RData')
drug.Class <- read.csv('Drugs_Vivli - Sheet1.csv', header = T)
#fix drug classes
drug.Class$Name <- paste0(drug.Class$Name, '_I') 

#import libs
library(psych)
library(gplots)
library(vcd)
library(caret)
library(sf)
library(ggplot2)
library("rnaturalearth")
library("rnaturalearthdata")
library(ggplot2)


select_I <- function(selection.Var, dataToUse = vivli.data) {
  
  desc_locs <- colnames(dataToUse)[c(3:5, 7:13)]
  abx_I_locs <- colnames(dataToUse)[seq(15,112,2)]
  abx_N_locs <- colnames(dataToUse)[seq(14,111,2)]
  gen_N_locs <- colnames(dataToUse)[seq(112,135,1)]
  bug_G_locs <- c("Species", "Country", abx_I_locs)
  bug_FD_locs <- c(desc_locs, abx_I_locs)
  
  out.dat <- NULL
  if (selection.Var == 'desc') {
    out.dat <- dataToUse[desc_locs]
  }
  if (selection.Var == 'pheno') {
    out.dat <- dataToUse['Phenotype']
  }
  if (selection.Var == 'abx_I') {
    out.dat <-dataToUse[abx_I_locs]
  }
  if (selection.Var == 'abx_N') {
    out.dat <- dataToUse[abx_N_locs]
  }
  if (selection.Var == 'gen_N') {
    out.dat <- dataToUse[gen_N_locs]
  }
  if (selection.Var == 'bug_G') {
    out.dat <- dataToUse[bug_G_locs]
  }
  if (selection.Var == 'bug_FD') {
    out.dat <- dataToUse[bugs_FD_locs]
  }
  
  return(out.dat)
}

core_resI <- select_I(c('abx_I'))

contigencyMapper <- function(drugA, drugB, data.df) {
  curr.data <- data.frame('A' = data.df[,drugA],
                          'B' = data.df[,drugB])
  #print(head(curr.data))
  cont.curr <- xtabs(~A+B, data = curr.data)[-(1:2),-(1:2)]#[c('0','1'), c('0','1')]
  
  return(cont.curr)
  
}

phi(contigencyMapper('Ciprofloxacin_I', 'Levofloxacin_I', core_resI))


bugGeoMapp <- function(bug = 'All', geo = 'All', genHeat = F, dataToUse = select_I('bug_G')) {
  
  bugGeo.base <- dataToUse
  
  bugGeo.curr <- NULL
  
  if (geo == 'All' & bug =='All') {
    
    bugGeo.curr = bugGeo.base
  }else if (geo == 'All' & bug != 'All') {
    bugGeo.curr = subset(bugGeo.base, Species == bug)
  }else if (geo != 'All' & bug == 'All'){
    bugGeo.curr = subset(bugGeo.base, Country == geo)
  }else if(geo != 'All' & bug != 'All'){
    bugGeo.curr = subset(bugGeo.base, Country == geo & Species == bug)
  }
  
  bugGeo.curr <- subset(bugGeo.curr, select = -c(Species, Country))
  print(dim(bugGeo.curr))
  
  
  interaction.mat <- NULL
  drug.order <- c()
  
  for (i in colnames(bugGeo.curr)) {
    drug.order <- c(drug.order, i)
    score.list <- c()
    for (j in colnames(bugGeo.curr)) {
      #print(paste0(i,'_',j))
      mkl <- contigencyMapper(i,j, bugGeo.curr)
      #print(mkl)
      if (sum(dim(mkl)) == 4) {
        curr_pair <- phi(mkl)
        #assocplot(curr_pair)
        score.list <- c(score.list, curr_pair)
      }else{
        curr_pair <- NA
        score.list <- c(score.list, curr_pair)
      }
      
    }
    interaction.mat <- cbind(interaction.mat, score.list)
  }
  
  colnames(interaction.mat) <- drug.order
  rownames(interaction.mat) <- drug.order
  print(dim(interaction.mat))
  #print(interaction.mat)
  #interaction.mat[1,1] <- -1
  
  if (genHeat) {
    heatmap.2(x = as.matrix(interaction.mat),
              col = rev(hcl.colors(100, palette = 'RdBu')),
              Rowv = F, Colv = F,
              scale = 'none', 
              colsep = c( 0:ncol(interaction.mat)),
              rowsep = c( 0:nrow(interaction.mat)),
              sepcolor = 'black',
              keysize = 1, 
              trace = "none",
              key = T,
              dendrogram = "none",
              cexRow = 0.75, cexCol = 0.75,
              na.color = 'black',
              main = paste0(bug, '_', geo))
    
  }
  
  return(interaction.mat)
  print('remember to fix color code')
}



sensiMatrixer <- function(cond1, cond2, df, bug = NULL, geo = NULL) {
  
  interaction.chi <- NULL
  drug.order <- c()
  
  for (i in colnames(df)) {
    drug.order <- c(drug.order, i)
    score.list <- c()
    score.assoc <- c()
    for (j in colnames(df)) {
      
      print(paste0(i,'_',j))
      mkl <- contigencyMapper(i,j, df)
      
      if (sum(dim(mkl)) >= 3 & sum(mkl) > 1) {
        
        #curr_pair <- assocstats(mkl)$cramer
        try(
          
          {if (length(chisq.test(mkl)$residuals[cond1,cond2]) == 1) {
            curr_pair <- assocstats(mkl)$cramer
            curr_assoc <- chisq.test(mkl)$residuals[cond1,cond2]
            
          }else{
            curr_assoc <- NA
          }}, silent = T
        )
        
        
        score.list <- c(score.list, curr_pair)
        score.assoc <- c(score.assoc, curr_assoc)
      }else{
        curr_pair <- NA
        curr_assoc <- NA
        score.list <- c(score.list, curr_pair)
        score.assoc <- c(score.assoc, curr_assoc)
      }
      
    }
    interaction.chi <- cbind(interaction.chi, score.assoc)
  }
  
  colnames(interaction.chi) <- drug.order
  rownames(interaction.chi) <- drug.order
  
  
  #print(dim(interaction.chi))
  
  heatmap.2(x = as.matrix(interaction.chi),
            col = hcl.colors(100, palette = 'Red-Green'),
            Rowv = F, Colv = F,
            scale = 'none',
            colsep = c( 0:ncol(interaction.mat)),
            rowsep = c( 0:nrow(interaction.mat)),
            sepcolor = 'black',
            keysize = 1, #key.xtickfun = key.xtickfun(),
            trace = "none",
            key = T,
            dendrogram = "none",
            cexRow = 0.75, cexCol = 0.75,
            main = "",
            na.color = 'black')
  
  return(interaction.chi)
  
}

#njt <- sensiMatrixer(cond2 = 'Resistant', cond1 = 'Susceptible', df = core_resI)

# max(njt, na.rm = T)
# 
# 
# ##one hot encoded data
# core_bugs_desc <- select_I('bug_FD')
# dums <- dummyVars("~.", data = core_bugs_desc, sep = '_')
# one_hot_data <- data.frame(predict(dums, newdata = core_bugs_desc))
# dim(one_hot_data)
# one_hot_data <- one_hot_data[, -which(names(one_hot_data) %in% colnames(core_bugs_desc) )]
# dim(one_hot_data)

#1:367, species
#368:383, family
#384: 466, country
#467:468, Gender
#469:475 age group
#476:487 speciality
#488:584 source 
#585:588 In/Out patient
#589:593 pehnotype
#others: susceptibility 

#plotting tools
myDensPl <- function(variables, 
                     shadeCol = "#22B258", bordCol = "black",
                     titleR = "", lineWidth = 3,
                     MM = F, axNot = F, 
                     horizontal = T) {
  
  tempDens <- density(variables, na.rm = T)
  
  if (axNot) {
    plot(tempDens, main = titleR, lwd = lineWidth, col = bordCol, xaxt = "n",
         yaxt = "n", bty = "n", xlab = "", ylab = "",
         frame = F)
  }else{
    plot(tempDens, main = titleR, lwd = lineWidth, col = bordCol)
  }
  
  polygon(tempDens, col = shadeCol, border = bordCol)
  
  #plot mean and median
  if (MM == T) {
    abline(v = mean(variables, na.rm = T), col = 'red')
    abline(v = median(variables, na.rm = T), col = 'black')
    #add legend for the mean and median
    legend("topright", legend = c('median', 'mean'), col = c('black', 'red'), lty = 1, lwd = 2, inset = 0.1)
  }
}

myDensLinesAd <- function(variables, 
                          shadeCol = "#22B258", bordCol = "black",
                          titleR = "", lineWidth = 3,
                          MM = F) {
  
  tempDens <- density(variables, na.rm = T)
  
  lines(tempDens, main = titleR, lwd = lineWidth, col = bordCol)
  polygon(tempDens, col = shadeCol, border = bordCol)
  
  #plot mean and median
  if (MM == T) {
    abline(v = mean(variables, na.rm = T), col = 'red')
    abline(v = median(variables, na.rm = T), col = 'black')
    #add legend for the mean and median
    legend("topright", legend = c('median', 'mean'), col = c('black', 'red'), lty = 1, lwd = 2, inset = 0.1)
  }
}

heatmapSubsetter <- function(InputData) {
  
  heatmap.2(x = as.matrix(InputData),
            col = rev(hcl.colors(100, palette = 'Red-Green')),
            Rowv = F, Colv = F,
            scale = 'none', 
            colsep = c( 0:ncol(InputData)),
            rowsep = c( 0:nrow(InputData)),
            sepcolor = 'black',
            keysize = 1, #key.xtickfun = key.xtickfun(),
            trace = "none",
            key = T,
            dendrogram = "none",
            cexRow = 1, cexCol = 1,
            main = "",
            na.color = 'black',
            key.title = '')
}



ViewTrendperYear <- function(drug1, drug2) {
  
  trending <- c()
  for (i in allYearMap) {
    #print(head(i))
    trending <- c(trending, i[drug1, drug2])
  }
  names(trending) <- yearOrder
  #print(trending)
  
  #plot(trending, 'b')
  plot(x = trending,
       type = 'b', 
       ylim = c(-1,1),
       lwd = 2,
       ylab = 'Phi Score',
       xlab = '',
       main = paste0(drug1, ' ', drug2),
       xaxt = 'n')
  
  axis(side = 1, at = 1:length(yearOrder), labels = yearOrder, las = 2)
  
  abline(h = 0, lty = 2)
  
  print(trending)
  
}

plot_heatmap <- function(matrix_data, axis_text_size = 10, label = NULL) {
  # Check if the input is a matrix
  if (!is.matrix(matrix_data)) {
    stop("Input must be a matrix.")
  }
  
  # Convert the matrix to a data frame in long format
  df <- reshape2::melt(matrix_data)
  
  # Create the heatmap using ggplot
  heatmap_plot <- ggplot(data = df, aes(x = Var1, y = Var2, fill = value)) +
    geom_tile(color = 'black') +
    scale_fill_gradient2(low = "#1a899c", mid = 'white', midpoint = 0, high = "#f69848", na.value = 'black', limits = c(-1,1)) +
    labs(x = "", y = "", fill = "Value") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = axis_text_size),
          axis.text.y = element_text(size = axis_text_size))+
    coord_fixed()+
    geom_text(aes(label = value), color = 'black', size = axis_text_size/3)
  
  # Print the heatmap
  print(heatmap_plot)
}


