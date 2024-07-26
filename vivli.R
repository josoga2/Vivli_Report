#setwd
source("/Users/josoga2/Documents/wale_docs/phd/R/vivli_tools.R")

#downsize to core descriptives
core_desc <- select_I('desc')
plot(table(core_desc$Year)/10000, 'l', las =2)


#downsize to resistance class
core_resI <- select_I(c('abx_I'))


#calculate phi values
assocstats(contigencyMapper('Amikacin_I','Ampicillin_I', core_resI))$cramer

interaction.mat <- NULL
drug.order <- c()
res.Count <- c()
suc.Count <- c()

for (i in colnames(core_resI)) {
  drug.order <- c(drug.order, i)
  score.list <- c()
  for (j in colnames(core_resI)) {
    print(paste0(i,'_',j))
    mkl <- contigencyMapper(i,j, core_resI)
    #print(mkl[,1])
    
    if (sum(dim(mkl)) == 4) {
      curr_pair <- phi(mkl)
      
      res.Count <- c(res.Count, mkl[1,1])
      suc.Count <- c(suc.Count, mkl[2,1])
    #  print(curr_pair)
      #assocplot(curr_pair)
      score.list <- c(score.list, curr_pair)
    }else{
      curr_pair <- NA
      score.list <- c(score.list, curr_pair)
    }
    
  }
  interaction.mat <- cbind(interaction.mat, score.list)
}



dim(interaction.mat)
head(interaction.mat)


colnames(interaction.mat) <- drug.order
rownames(interaction.mat) <- drug.order


par(oma = c(10,0,0,7))

heatmap.2(x = as.matrix(interaction.mat),
          col = rev(hcl.colors(100, palette = 'Red-Green')),
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
          main = "Overview of Global Phi-Score",
          na.color = 'black', 
          key.title = '',
          key.xlab = '',
          key.ylab = '')


dev.off()

#figre 2a and b
resSuc.count <- table(unname(unlist(core_resI)))[3:4]
barplot((resSuc.count)/100000, 
        las = 1, 
        main = 'Antibiogram Counts', 
        ylab = 'x 1 million',
        col = c('seagreen', 'brown1'))

countofPos <- sum(interaction.mat > 0.1, na.rm = T)
countofNeg <- sum(interaction.mat < -0.1, na.rm = T)

barplot(c('Negative' = countofNeg, 'Positive' = countofPos), 
        las = 1, 
        main = 'Phi-Score Count', 
        ylab = '',
        col = c('seagreen', 'brown1'))



#save to data
#save(interaction.mat, file = './interaction.mat.RData')
#save(vivli.data, file = './vivli_data.RData')

for (var in rownames(interaction.mat)) {
  print(var, quote = F)
}

##===
#geospatial analysis
bugGeoMapp(geo = 'Australia', genHeat = T)
unique(vivli.data$Country)
unique(vivli.data$Species)

#median centering of the background data
interaction.mat.norm <- interaction.mat - median(as.vector(interaction.mat), na.rm = T)

#compute country level stats
country.name <- c()
medScore <- c()
curr_p.v <- rep(NA, length(unique(vivli.data$Country)))
medDiff <- c()
no_obs <- c()
bgMed <- median(interaction.mat, na.rm = T)
bug.country.split <- NULL
microbe.name <- c()


#pdf(file = "vivli_Country_ Distro.pdf", paper = 'a4r', width = 120, height = 80)
#par(mfrow = c(8,12), cex=0.25, mar=c(2,2,2,2), oma=c(2,2,2,2), no.readonly = T)

n = 1

for (country in unique(vivli.data$Country)) {
  
  print(country)
  
  country.name <- c(country.name, country)
  curr.country <- bugGeoMapp(geo = country, genHeat = F, bug = 'All')
  
  no_obs <- c(no_obs, nrow(curr.country)) #number of observations
  curr.Med <- median(curr.country, na.rm = T) # median
  medScore <- c(medScore, curr.Med) #add to median collection
  
  
  tryCatch(
    expr = {
      curr_p.v[n] <- t.test(as.vector(curr.country), as.vector(interaction.mat), paired = F)$p.value #compute pvalues
      
      #myDensPl(interaction.mat, bordCol = 'grey', shadeCol = t_col('red', percent = 100), lineWidth = 0.75, titleR = cnty)
      #myDensLinesAd(curr.country, bordCol = 'seagreen', shadeCol = t_col('red', percent = 100), lineWidth = 0.75)
    }, error=function(err) NA
  )
  
  medDiff <- c(medDiff, c(curr.Med - bgMed)) 
  
  n = n+1
  
}

viv.Country.dist <- data.frame('country' = country.name,
                               'med_Score' = medScore,
                               'no_obs' = no_obs,
                               'pv' = curr_p.v,
                               'bh' = p.adjust(curr_p.v, method = 'BH'),
                               'med_diff' = medDiff,
                               'bg_Med' = bgMed)

dev.off()
plot(viv.Country.dist$med_diff, 
     -log10(viv.Country.dist$bh), 
     pch = 21,
     bg = 'grey',
     cex = 1,
     xlim = c(-0.2, 0.3), 
     xlab = 'Median Difference', 
     ylab = '-log10(P Value)',
     main = 'By Country')

with(subset(viv.Country.dist, bh <= 0.05 & med_diff <= -0.1),
     text(x=med_diff, y=-log10(pv), label = country , pch = 21, cex= 1, col = 'red', pos=4))

with(subset(viv.Country.dist, bh <= 0.05 & med_diff >= 0.1),
     points(med_diff, -log10(pv), pch = 21, cex= 2, bg = 'blue'))


with(subset(viv.Country.dist, bh <= 0.05 & med_diff <= -0.1),
     points(med_diff, -log10(pv), pch = 21, cex= 2, bg = 'red'))

library(ggrepel)
ggplot(viv.Country.dist, aes(med_diff, -log10(bh), label = country))+
  geom_point(cex = 5, pch = 21, bg = 'grey')+
  lims(x = c(-0.3, 0.3))+
  geom_text_repel(data = subset(viv.Country.dist, bh <= 0.05 & med_diff <= -0.1), color = 'red')+
  geom_text_repel(data = subset(viv.Country.dist, bh <= 0.05 & med_diff >= 0.1), color = 'blue')+
  geom_vline(xintercept = c(0.1, -0.1), lty = 2, lwd = 0.25)+
  geom_hline(yintercept = c(-log10(0.05)), lty = 2, lwd = 0.25)+
  labs(x = 'Z-score (Median Correction)')+
  theme_bw()

bugGeoMapp(geo = 'India', genHeat = T)

bugGeoMapp(bug = 'Escherichia coli', genHeat = T)

par(oma = c(7.5,0,0,7.5))
PENS <- subset(drug.Class, drug.Class$SubClass == 'Quinolones')$Name
plot_heatmap(interaction.mat[PENS,PENS], axis_text_size = 16)


QD_hit <- subset(drug.Class, drug.Class$SubClass == 'Streptogramin')$Name


par(oma = c(7.5,0,0,0))
QD_scores <- interaction.mat["Minocycline_I",]
QD_scores <- QD_scores[!is.na(QD_scores)]
barplot(QD_scores, las = 2, ylim = c(-1,1), ylab = 'Phi-Score', main = QD_hit)

#species level
medScore.spec <- c()
curr_p.v.spec <- rep(NA, length(unique(vivli.data$Species)))
medDiff.spec <- c()
no_obs.spec <- c()
bgMed.spec <- median(interaction.mat, na.rm = T)
microbe.name <- c()


#pdf(file = "vivli_Country_ Distro.pdf", paper = 'a4r', width = 120, height = 80)
#par(mfrow = c(8,12), cex=0.25, mar=c(2,2,2,2), oma=c(2,2,2,2), no.readonly = T)

n = 1

for (species in unique(vivli.data$Species)) {
  
  print(species)
  
  microbe.name <- c(microbe.name, species)
  curr.spec <- bugGeoMapp(genHeat = F, bug = species, geo = 'All')
  
  no_obs.spec <- c(no_obs.spec, nrow(curr.spec)) #number of observations
  curr.Med <- median(curr.spec, na.rm = T) # median
  medScore.spec <- c(medScore.spec, curr.Med) #add to median collection
  
  
  tryCatch(
    expr = {
      curr_p.v.spec[n] <- t.test(as.vector(curr.spec), as.vector(interaction.mat), paired = F)$p.value #compute pvalues
      
      #myDensPl(interaction.mat, bordCol = 'grey', shadeCol = t_col('red', percent = 100), lineWidth = 0.75, titleR = cnty)
      #myDensLinesAd(curr.country, bordCol = 'seagreen', shadeCol = t_col('red', percent = 100), lineWidth = 0.75)
    }, error=function(err) NA
  )
  
  medDiff.spec <- c(medDiff.spec, c(curr.Med - bgMed.spec)) 
  
  n = n+1
  
}

viv.Species.dist <- data.frame('microbe' = microbe.name,
                               'med_Score' = medScore.spec,
                               'no_obs' = no_obs.spec,
                               'pv' = curr_p.v.spec,
                               'med_diff' = medDiff.spec,
                               'bg_Med' = bgMed.spec)

plot(viv.Species.dist$med_diff, 
     -log10(viv.Species.dist$pv), 
     pch = 21,
     bg = 'grey',
     cex = 1,
     ylim = c(0, 20), 
     xlab = 'Median Difference', 
     ylab = '-log10(P Value)',
     main = 'By Species')

#both

country.name <- c()
species.name <- c()
total_count <- c()
all.median <- c()
bgMedian <- median(as.vector(interaction.mat), na.rm = T) #perform median centering
all.PV <- rep(NA, (length(unique(vivli.data$Country)) * length(unique(vivli.data$Species))) )

n= 1

for (country in unique(vivli.data$Country)) {
  
  country.name <- c(country.name, country)
  
  for (species in unique(vivli.data$Species)) {
    species.name <- c(species.name, species)
    
    #print(species)
    #print(country)
    
    curr_bug_geo <- bugGeoMapp(bug = species, geo = country, genHeat = F)
    total_count <- c(total_count, length(unlist(curr_bug_geo)))
    
    all.median <- c(all.median, median(as.vector(curr_bug_geo), na.rm = T))
    
    tryCatch(
      expr = {
        all.PV[n] <- wilcox.test(x = as.vector(curr_bug_geo), y = as.vector(interaction.mat))$p.value
      }, error=function(err) NA
    )
    
    n = n+1
    
  }
}


speciesPerCountry <- data.frame('country.name' = country.name,
                                'species.name' = species.name,
                                'total_count' = total_count,
                                'all.median' = all.median,
                                'bgMedian' = bgMedian,
                                'all.PV' = all.PV,
                                'median_diff' = all.median - bgMedian)

plot_ly(data = speciesPerCountry, type = 'scatter', x = ~median_diff, y = ~all.PV )

subset(speciesPerCountry, median_diff <0.2 & median_diff >0.1) 
bugGeoMapp(geo ='Germany', bug = 'Escherichia coli', genHeat = T)


#does interaction change with year?

yearOrder <- as.integer(str_sort(unique(vivli.data$Year), numeric = T))

allYearMap <- c()
pdf(file = "yearMap.pdf", width = 12, height = 12)
par(oma = c(6.5,0, 0,6.5))
for (year in yearOrder) {
  print(year)
  yearData <- subset(vivli.data, Year == year)
  curr_Year <- bugGeoMapp(dataToUse = select_I('bug_G', dataToUse = yearData),genHeat = T)
  mtext(year)
  
  allYearMap <- c(allYearMap, list(curr_Year))
}

dev.off()

#searching for trends
names(allYearMap) <- yearOrder

ViewTrendperYear('Amoxycillin.clavulanate_I', 'Ampicillin_I')
plot_heatmap(interaction.mat, axis_text_size = 8)

