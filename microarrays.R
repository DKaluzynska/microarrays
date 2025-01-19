#czyszczenie
rm(list = ls())

#ustawienie ?cie?ki
setwd("C:/Users/daria/Desktop/kody/mikromacierze") 

#wczytanie bibliotek
#BiocManager::install("clusterProfiler")

library("AnnotationDbi")
library("clusterProfiler")
library(biomaRt)
library(plotly)
library(ggplot2)
library(marray)
library(pheatmap)
library(dplyr)
library(tidyverse)
library(genefilter)
library(DT)
library(GEOquery)
library(gplots)
ensembl <- useMart('ensembl', dataset = 'hsapiens_gene_ensembl')

#Wczytanie danych z pliku
arr_1 = read.Agilent('US11193903_252800414566_S01_GE2_1100_Jul11_1_1.txt')
arr_2 = read.Agilent('US11193903_252800414566_S01_GE2_1100_Jul11_1_2.txt')
arr_3 = read.Agilent('US11193903_252800414566_S01_GE2_1100_Jul11_1_3.txt')
arr_4 = read.Agilent('US11193903_252800414566_S01_GE2_1100_Jul11_1_4.txt')
arr_5 = read.Agilent('US11193903_252800414566_S01_GE2_1100_Jul11_2_1.txt')
arr_6 = read.Agilent('US11193903_252800414566_S01_GE2_1100_Jul11_2_2.txt')
arr_7 = read.Agilent('US11193903_252800414566_S01_GE2_1100_Jul11_2_3.txt')
arr_8 = read.Agilent('US11193903_252800414566_S01_GE2_1100_Jul11_2_4.txt')

#przedstawienie obraz?w macierzy
#image(arr_1)
#image(arr_2)
#image(arr_3)
#image(arr_4)
#image(arr_5)
#image(arr_6)
#image(arr_7)
#image(arr_8)

#Podzia? na kana?y
arr1_G = maGf(arr_1)
colnames(arr1_G) <- "expression1"
arr1_R = maRf(arr_1)
colnames(arr1_R) <- "expression1"

arr2_G = maGf(arr_2)
colnames(arr2_G) <- "expression2"
arr2_R = maRf(arr_2)
colnames(arr2_R) <- "expression2"

arr3_G = maGf(arr_3)
colnames(arr3_G) <- "expression3"
arr3_R = maRf(arr_3)
colnames(arr3_R) <- "expression3"

arr4_G = maGf(arr_4)
colnames(arr4_G) <- "expression4"
arr4_R = maRf(arr_4)
colnames(arr4_R) <- "expression4"

arr5_G = maGf(arr_5)
colnames(arr5_G) <- "expression5"

arr5_R = maRf(arr_5)
colnames(arr5_R) <- "expression5"


arr6_G = maGf(arr_6)
colnames(arr6_G) <- "expression6"

arr6_R = maRf(arr_6)
colnames(arr6_R) <- "expression6"

arr7_G = maGf(arr_7)
colnames(arr7_G) <- "expression7"
arr7_R = maRf(arr_7)
colnames(arr7_R) <- "expression7"

arr8_G = maGf(arr_8)
colnames(arr8_G) <- "expression8"
arr8_R = maRf(arr_8)
colnames(arr8_R) <- "expression8"


#ekspresja i logarytmowanie
arr_1_logGf = log2(arr1_G)
arr_1_logGf<- as.data.frame(arr_1_logGf)
colnames(arr_1_logGf) <- "expression"

arr_1_logRf = log2(arr1_R)
arr_1_logRf<- as.data.frame(arr_1_logRf)
colnames(arr_1_logRf) <- "expression"


arr_2_logGf = log2(arr2_G)
arr_2_logGf<- as.data.frame(arr_2_logGf)
colnames(arr_2_logGf) <- "expression"

arr_2_logRf = log2(arr2_R)
arr_2_logRf<- as.data.frame(arr_2_logRf)
colnames(arr_2_logRf) <- "expression"

arr_3_logGf = log2(arr3_G)
arr_3_logGf<- as.data.frame(arr_3_logGf)
colnames(arr_3_logGf) <- "expression"

arr_3_logRf = log2(arr3_R)
arr_3_logRf<- as.data.frame(arr_3_logRf)
colnames(arr_3_logRf) <- "expression"

arr_4_logGf = log2(arr4_G)
arr_4_logGf<- as.data.frame(arr_4_logGf)
colnames(arr_4_logGf) <- "expression"

arr_4_logRf = log2(arr4_R)
arr_4_logRf<- as.data.frame(arr_4_logRf)
colnames(arr_4_logRf) <- "expression"

arr_5_logGf = log2(arr5_G)
arr_5_logGf<- as.data.frame(arr_5_logGf)
colnames(arr_5_logGf) <- "expression"

arr_5_logRf = log2(arr5_R)
arr_5_logRf<- as.data.frame(arr_4_logRf)
colnames(arr_4_logRf) <- "expression"

arr_6_logGf = log2(arr6_G)
arr_6_logGf<- as.data.frame(arr_6_logGf)
colnames(arr_6_logGf) <- "expression"

arr_6_logRf = log2(arr6_R)
arr_6_logRf<- as.data.frame(arr_6_logRf)
colnames(arr_6_logRf) <- "expression"

arr_7_logGf = log2(arr7_G)
arr_7_logGf<- as.data.frame(arr_7_logGf)
colnames(arr_7_logGf) <- "expression"

arr_7_logRf = log2(arr7_R)
arr_7_logRf<- as.data.frame(arr_7_logRf)
colnames(arr_7_logRf) <- "expression"

arr_8_logGf = log2(arr8_G)
arr_8_logGf<- as.data.frame(arr_8_logGf)
colnames(arr_8_logGf) <- "expression"

arr_8_logRf = log2(arr8_R)
arr_8_logRf<- as.data.frame(arr_8_logRf)
colnames(arr_8_logRf) <- "expression"


#wykresy pude?kowe
figGf <- plot_ly(data=arr_1_logGf,
                  y=~expression,
                 name = "Arr1_green",
  type = "box")
figGf <- figGf %>% add_boxplot(data = arr_2_logGf,
                               y=~expression,
                               name = "Arr2_green")
figGf <- figGf %>% add_boxplot(data = arr_3_logGf,
                               y=~expression,
                               name = "Arr3_green")
figGf <- figGf %>% add_boxplot(data = arr_4_logGf,
                               y=~expression,
                               name = "Arr4_green")
figGf <- figGf %>% add_boxplot(data = arr_5_logGf,
                               y=~expression,
                               name = "Arr5_green")
figGf <- figGf %>% add_boxplot(data = arr_6_logGf,
                               y=~expression,
                               name = "Arr6_green")
figGf <- figGf %>% add_boxplot(data = arr_7_logGf,
                               y=~expression,
                               name = "Arr7_green")
figGf <- figGf %>% add_boxplot(data = arr_8_logGf,
                               y=~expression,
                               name = "Arr8_green")


figRf <- plot_ly(data=arr_1_logRf,
                 y=~expression,
                 type = "box",
                 name = "Arr1_red")
figRf <- figRf %>% add_boxplot(data = arr_2_logRf,
                               y=~expression,
                               name = "Arr2_red")
figRf <- figRf %>% add_boxplot(data = arr_3_logRf,
                               y=~expression,
                               name = "Arr3_red")
figRf <- figRf %>% add_boxplot(data = arr_4_logRf,
                               y=~expression,
                               name = "Arr4_red")
figRf <- figRf %>% add_boxplot(data = arr_5_logRf,
                               y=~expression,
                               name = "Arr5_red")
figRf <- figRf %>% add_boxplot(data = arr_6_logRf,
                               y=~expression,
                               name = "Arr6_red")
figRf <- figRf %>% add_boxplot(data = arr_7_logRf,
                               y=~expression,
                               name = "Arr7_red")
figRf <- figRf %>% add_boxplot(data = arr_8_logRf,
                               y=~expression,
                               name = "Arr8_red")

#figGf
#figRf

#####################################################################################
#histogramy

figGf1_hist <- plot_ly(data=arr_1_logGf,
                 x=~expression,
                 type = "histogram",
                 name = "Arr1_green")
figGf2_hist <- plot_ly(data=arr_2_logGf,
                      x=~expression,
                      type = "histogram",
                      name = "Arr2_green")

figGf3_hist <- plot_ly(data=arr_3_logGf,
                      x=~expression,
                      type = "histogram",
                      name = "Arr3_green")

figGf4_hist <- plot_ly(data=arr_4_logGf,
                      x=~expression,
                      type = "histogram",
                      name = "Arr4_green")
figGf5_hist <- plot_ly(data=arr_5_logGf,
                      x=~expression,
                      type = "histogram",
                      name = "Arr5_green")
figGf6_hist <- plot_ly(data=arr_6_logGf,
                      x=~expression,
                      type = "histogram",
                      name = "Arr6_green")
figGf7_hist <- plot_ly(data=arr_7_logGf,
                      x=~expression,
                      type = "histogram",
                      name = "Arr7_green")
figGf8_hist <- plot_ly(data=arr_8_logGf,
                      x=~expression,
                      type = "histogram",
                      name = "Arr8_green")


#figGf1_hist
#figGf2_hist
#figGf3_hist
#figGf4_hist
#figGf5_hist
#figGf6_hist
#figGf7_hist
#figGf8_hist

figRf1_hist <- plot_ly(data=arr_1_logRf,
                       x=~expression,
                       type = "histogram",
                       name = "Arr1_green")
figRf2_hist <- plot_ly(data=arr_2_logRf,
                       x=~expression,
                       type = "histogram",
                       name = "Arr2_green")

figRf3_hist <- plot_ly(data=arr_3_logRf,
                       x=~expression,
                       type = "histogram",
                       name = "Arr3_green")

figRf4_hist <- plot_ly(data=arr_4_logRf,
                       x=~expression,
                       type = "histogram",
                       name = "Arr4_green")
figRf5_hist <- plot_ly(data=arr_5_logRf,
                       x=~expression,
                       type = "histogram",
                       name = "Arr5_green")
figRf6_hist <- plot_ly(data=arr_6_logRf,
                       x=~expression,
                       type = "histogram",
                       name = "Arr6_green")
figRf7_hist <- plot_ly(data=arr_7_logRf,
                       x=~expression,
                       type = "histogram",
                       name = "Arr7_green")
figRf8_hist <- plot_ly(data=arr_8_logRf,
                       x=~expression,
                       type = "histogram",
                       name = "Arr8_green")

#figRf1_hist
#figRf2_hist
#figRf3_hist
#figRf4_hist
#figRf5_hist
#figRf6_hist
#figRf7_hist
#figRf8_hist


#Nazwy i opis
nazwy_arr = maGnames(arr_1)@maInfo

geny1 = nazwy_arr$SystematicName


########################################
###G
# Heatmapa kanal zielony

arr_1_logGf<-sapply(arr_1_logGf, as.numeric)

merged_data_G <-as.data.frame(cbind(arr_1_logGf, arr_2_logGf, arr_3_logGf, arr_4_logGf, arr_5_logGf, arr_6_logGf, arr_7_logGf, arr_8_logGf, geny1))
colnames(merged_data_G)[9] <- "g"

filtered_data_G <- distinct(merged_data_G, g, .keep_all = TRUE)
filtered_data_G <- filtered_data_G %>%  
  filter(grepl("A_",`g`)) 
rownames(filtered_data_G)<- filtered_data_G[,9]
filtered_data_G_1<- (filtered_data_G[,1:8])
filtered_data_G_1<- as.data.frame(sapply(filtered_data_G_1, as.numeric))



probemeans <- apply(filtered_data_G[,1:8], 1, mean)
probesd <- apply(filtered_data_G[,1:8], 1, sd)
plot(probemeans, probesd)
q25 <- quantile(probemeans, 0.25,na.rm = T)
whichtosave <- which(probemeans > q25)
q25logdata <- filtered_data_G[,1:8][whichtosave,]
mydata <- q25logdata[apply(q25logdata, 1, IQR) > 1.5, ]
tdata <- aperm(as.matrix(mydata))
pca <- prcomp(tdata, scale=T)
plot(pca$x, type="n")
text(pca$x, rownames(pca$x), cex=0.5)
pearsonCorr <- as.dist(1 - cor(mydata))
hC <- hclust(pearsonCorr)
plot(hC)
heatmap(as.matrix(mydata), col = greenred(100))

genenames<- rownames(mydata)


annot <- getBM(
  attributes = c('agilent_wholegenome',
                 'wikigene_description',
                 'ensembl_gene_id',
                 'entrezgene_id',
                 'gene_biotype',
                 'external_gene_name'),
  filters = 'agilent_wholegenome',
  values = genenames,
  mart = ensembl)



################################################
##R

# heatmapa kanal czerwony
merged_data_R <-as.data.frame(cbind(arr_1_logRf, arr_2_logRf, arr_3_logRf, arr_4_logRf, arr_5_logRf, arr_6_logRf, arr_7_logRf, arr_8_logRf, nazwy_arr$SystematicName))
colnames(merged_data_R)[9] <- "g"

filtered_data_R <- distinct(merged_data_R, g, .keep_all = TRUE)
filtered_data_R <- filtered_data_R %>%  
  filter(grepl("A_",`g`)) 
rownames(filtered_data_R)<- filtered_data_R[,9]
filtered_data_R_1<- (filtered_data_R[,1:8])
filtered_data_R_1<- as.data.frame(sapply(filtered_data_R_1, as.numeric))



selected_data_R <- mydata_R[only_genes[,1],]
selected_data_R <- cbind(selected_data_R, only_genes[,5])

probemeans_r <- apply(selected_data_R[,1:8], 1, mean)
probesd_r <- apply(selected_data_R[,1:8], 1, sd)
plot(probemeans_r, probesd_r)
q25r <- quantile(probemeans_r, 0.25,na.rm = T)
whichtosaver <- which(probemeans_r > q25r)
q25logdatar <- selected_data_R[,1:8][whichtosaver,]
mydata_R <- q25logdatar[apply(q25logdatar, 1, IQR) > 1.5, ]
tdata_R <- aperm(as.matrix(mydata_R))
pca_R <- prcomp(tdata_R, scale=T)
plot(pca_R$x, type="n")
text(pca_R$x, rownames(pca_R$x), cex=0.5)
pearsonCorr_R <- as.dist(1 - cor(mydata_R[]))
hC_R <- hclust(pearsonCorr_R)
plot(hC_R)
heatmap(as.matrix(mydata_R), col = greenred(100))

genenames_R<- rownames(mydata_R)
annot_R <- getBM(
  attributes = c('agilent_wholegenome',
                 'wikigene_description',
                 'ensembl_gene_id',
                 'entrezgene_id',
                 'gene_biotype',
                 'external_gene_name'),
  filters = 'agilent_wholegenome',
  values = genenames_R,
  mart = ensembl)
only_genes <- annot_R %>% filter(!grepl("processed_pseudogene", gene_biotype))
