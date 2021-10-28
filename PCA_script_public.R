#Execute principal component analysis using flow cytometric data - Grace Edmunds March 2021

#clear workspace
rm(list=ls())

#Load Packages 
library(FactoMineR)
library(factoextra)
library(readr)
library("corrplot")

# setworkingdirectory

# setwd <- <mydirectory>

#Import data (the file is provided in supplementary data, you need to save it to your working directory to use it:

cytometry <- read.table("./171122_CIR_csv_Rdata.csv", 
                 header = TRUE,
                 sep = ",")

# Make copy of raw data called cytometry 2 to keep
cytometry2<-cytometry

#Clean data including col names to numbers etc
colnames(cytometry)<-c("Treatment", paste(1:153, sep = ","))

#Make a table of all of the column names with their numbers
CD_combinations<-colnames(cytometry2)

#Check that we have the right number of individuals in each treatment group
table(cytometry$Treatment)

#Run PCA on whole data frame setting column 1 (treatment) as supplem quali variable
res.pca <- PCA(cytometry[,1:154], graph = FALSE, quali.sup = 1)

#Show what results are in the PCA
res.pca

#Get the variables 
var <- get_pca_var(res.pca)

#Make a description to identify variables most significantly associated with each PC
res.desc <- dimdesc(res.pca, axes = c(1,2), proba = 0.05)

# Description of dimension 1
res.desc$Dim.1

#Get the Cos2 values for the variables in each of 5 PCs and make corrplot with text small enough to read
corrplot(var$cos2, is.corr=FALSE, tl.cex = 0.2)

#Get the contribution of the variables to 5 PCs and make corrplot text small 
corrplot(var$contrib, is.corr=FALSE, tl.cex = 0.2) 

#Get the contribution of variables to PC1
fviz_contrib(res.pca, choice = "var", axes = 1)

#This shows they are important to around the first 50 so show these only 
fviz_contrib(res.pca, choice = "var", axes = 1, top = 50)

#make variables plot with top 20 Cos2 values on it
fviz_pca_var(res.pca, select.var = list(cos2 = 20))

#Make biplot PC1 and PC2 with top 20 Cos2 value variables on, variables as arrows, individuals as dots, colour by treatment (which is qualitative sup variable), 95% confidence elipses 
fviz_pca_biplot(res.pca, select.var = list(cos2 = 20), axes = c(1,2), geom = c("point"), geom.var = c("arrow"),  habillage = "Treatment", 
                addEllipses = TRUE, ellipse.type = "confidence")

#Same biplot but with concentration elipses 
fviz_pca_biplot(res.pca, select.var = list(cos2 = 20), axes = c(1,2), geom = c("point"), geom.var = c("arrow"),  habillage = "Treatment", 
                addEllipses = TRUE)

#Make biplot of PC1 and PC2 with top 20 CONTRIBUTING (vs cos2) variables shown labelled so can ID variables but no eplises
fviz_pca_biplot(res.pca, select.var = list(contrib = 20), axes = c(1,2), geom = c("point"), repel = T,  habillage = "Treatment")

#Same with no variable labels 
fviz_pca_biplot(res.pca, select.var = list(contrib = 20), axes = c(1,2), geom = c("point"), geom.var = c("arrow"),  habillage = "Treatment")

#Same with elipses and no variables 
fviz_pca_biplot(res.pca, select.var = list(contrib = 20), axes = c(1,2), geom = c("point"), geom.var = c("arrow"),  habillage = "Treatment", 
                addEllipses = TRUE)

#Same biplot PC1 and PC3 no elipses
fviz_pca_biplot(res.pca, select.var = list(contrib = 20), axes = c(1,3), geom = c("point"), geom.var = c("arrow"),  habillage = "Treatment")

#Same PC1 and PC4
fviz_pca_biplot(res.pca, select.var = list(contrib = 20), axes = c(1,4), geom = c("point"), geom.var = c("arrow"),  habillage = "Treatment") 
               
#Same PC1 and PC5 
fviz_pca_biplot(res.pca, select.var = list(contrib = 20), axes = c(1,5), geom = c("point"), geom.var = c("arrow"),  habillage = "Treatment") 
                
#Same PC 2 and 3
fviz_pca_biplot(res.pca, select.var = list(contrib = 20), axes = c(2,3), geom = c("point"), geom.var = c("arrow"),  habillage = "Treatment") 
            
#Same PC2 and 4
fviz_pca_biplot(res.pca, select.var = list(contrib = 20), axes = c(2,4), geom = c("point"), geom.var = c("arrow"),  habillage = "Treatment") 
              
#Same PC2 and 5
fviz_pca_biplot(res.pca, select.var = list(contrib = 20), axes = c(2,5), geom = c("point"), geom.var = c("arrow"),  habillage = "Treatment") 

#Same PC3 and PC4
fviz_pca_biplot(res.pca, select.var = list(contrib = 20), axes = c(3,4), geom = c("point"), geom.var = c("arrow"),  habillage = "Treatment")

#Same PC3 and PC5 
fviz_pca_biplot(res.pca, select.var = list(contrib = 20), axes = c(3,5), geom = c("point"), geom.var = c("arrow"),  habillage = "Treatment") 
        
#Same PC4 and PC5
fviz_pca_biplot(res.pca, select.var = list(contrib = 20), axes = c(4,5), geom = c("point"), geom.var = c("arrow"),  habillage = "Treatment")

# To work out what these variables are, label the vector arrows far apart
fviz_pca_biplot(res.pca, select.var = list(cos2 = 20), axes = c(1,2), geom = c("point"), repel = T,  habillage = "Treatment")

#Make new copy of Data 
cytometry3 <- cytometry2

#Clean data including column names to numbers and label treatment and volume etc
colnames(cytometry3)<-c("Treatment", "Volume", paste(2:153, sep = ","))

#set volume as quanti sup variable and run pca
res.pca <- PCA(cytometry[,1:154], graph = FALSE, quali.sup = 1, quanti.sup = 2)

#Visualise with volume as habillage 
fviz_pca_biplot(res.pca2, select.var = list(cos2 = 20), axes = c(1,2), geom.var = c("arrow"), geom.ind = c("point"), habillage = "Volume")

