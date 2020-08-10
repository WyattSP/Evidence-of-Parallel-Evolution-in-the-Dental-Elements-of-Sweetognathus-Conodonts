#title: R Code for Evidence of Parallel Evolution in the Dental Elements of Sweetognathus Conodonts
#author: "Wyatt Petryshen"
#date: '2020-07-31'
#output: html_document

#Required packages for landmark selection and procrustes superimposition, and some data visulization later on. Analysis was run using following versions of each package: Momocs (1.2.9), abind (1.4-5), Morpho (2.7), tidyverse (1.3.0), plotly (4.9.1), car (3.0-7). 
library(Momocs); library(abind); library(Morpho); library(tidyverse); library(plotly); library(car)

#First step is to import data into R markdown document, then gather information of ID and denticle positions prior to procrustes superimposition.
#Import of rds dataset, make sure to set proper working directory; this is of the jpeg outlines converted into xy-coordinates.This will be sampled in Momocs prior to procrustes superimposition.
Sweet.Out.Full <- readRDS(file = "Sweetognathus.Outlines.rds")
#Get specimen ID's with denticle position for each outline.
sample.names.full <- attributes(Sweet.Out.Full)
#Get specimen ID's for each set of denticles.
Specimen.ID <- str_sub(sample.names.full$names, start = 4, end = -1)
#Get denticle position for each outline.
Denticle.ID <- str_sub(sample.names.full$names, start = 1, end =2)
#Assign species identity and locality to data
SID <- c("\\bS1\\b", "\\bS2\\b", "\\bS3\\b", "\\bS4\\b", "\\bS5\\b",
         "\\bS6\\b", "\\bS7\\b", "\\bS8\\b", "\\bS9\\b", "\\bS10\\b",
         "\\bS11\\b", "\\bS12\\b", "\\bS13\\b", "\\bS14\\b", "\\bS15\\b",
         "\\bS16\\b", "\\bS17\\b", "\\bS18\\b", "\\bS19\\b", "\\bS20\\b",
         "\\bS21\\b", "\\bS22\\b", "\\bS23\\b", "\\bS24\\b", "\\bS25\\b",
         "\\bS26\\b", "\\bS27\\b", "\\bS28\\b", "\\bS29\\b", "\\bS30\\b",
         "\\bS31\\b", "\\bS32\\b", "\\bS33\\b", "\\bS34\\b", "\\bS35\\b",
         "\\bS36\\b", "\\bS37\\b", "\\bS38\\b", "\\bS39\\b", "\\bS40\\b",
         "\\bS41\\b", "\\bS42\\b", "\\bS43\\b", "\\bS44\\b", "\\bS45\\b",
         "\\bS46\\b", "\\bS47\\b", "\\bS48\\b", "\\bS49\\b", "\\bS50\\b",
         "\\bS51\\b", "\\bS52\\b", "\\bS53\\b", "\\bS54\\b", "\\bS55\\b",
         "\\bS56\\b", "\\bS57\\b", "\\bS58\\b", "\\bS59\\b")

RID <- c("whitei","whitei","whitei","whitei","whitei",
         "whitei","whitei","whitei", "asymmetrica", "clarki", "clarki",
         "clarki", "clarki", "whitei","whitei","whitei", "whitei","whitei","whitei",
         "obliquidentatus",
         "asymmetrica",
         "obliquidentatus",  "obliquidentatus",  "obliquidentatus",  "obliquidentatus",
         "asymmetrica", "asymmetrica",
         "cf. adenticulatus", "aff. binodosus", "aff. binodosus", "aff. binodosus", "aff. binodosus", "aff. binodosus", "cf. adenticulatus",  "cf. adenticulatus", "aff. binodosus",
         "obliquidentatus", "binodosus",  "obliquidentatus", "binodosus",  "obliquidentatus", "binodosus",
         "binodosus","binodosus","obliquidentatus",
         "binodosus","binodosus","binodosus","binodosus","binodosus","binodosus","binodosus",
         "anceps", "binodosus", "anceps", "binodosus", "anceps", "anceps", "anceps")

Locations <- c("Tensleep","Tensleep","Tensleep","Tensleep","Tensleep","Tensleep","Tensleep","Tensleep",
               "Great Bear (17)", "Great Bear (20)","Great Bear (20)","Great Bear (20)","Great Bear (20)",
               "Florence", "Florence","Florence","Florence","Florence","Florence",
               "Russia (4DT)",  "Russia (4DT)",  "Russia (4DT)",  "Russia (4DT)",  "Russia (4DT)",  "Russia (4DT)",
               "Russia (7DT)", "Russia (7DT)",
               "Bolivia",   "Bolivia",   "Bolivia",   "Bolivia",   "Bolivia",   "Bolivia",   "Bolivia",   "Bolivia",   "Bolivia",
               "Melville",  "Melville",  "Melville",  "Melville",  "Melville",
               "Carlin Canyon (703)", "Carlin Canyon (703)", "Carlin Canyon (703)", "Carlin Canyon (703)",
               "Carlin Canyon (699)", "Carlin Canyon (699)", "Carlin Canyon (699)", "Carlin Canyon (699)", "Carlin Canyon (699)", "Carlin Canyon (699)", "Carlin Canyon (699)",
               "Carlin Canyon (706)",  "Carlin Canyon (706)",  "Carlin Canyon (706)",  "Carlin Canyon (706)",  "Carlin Canyon (706)",  "Carlin Canyon (706)",  "Carlin Canyon (706)",
               "Great Bear (20)","Melville",  "Carlin Canyon (703)", "Russia (4DT)", "Great Bear (20)" )

#Assign species name
Species.ID <- Specimen.ID
for(i in 1:length(SID)){
  Species.ID = gsub(SID[i], RID[i], Species.ID)
}
#Assign locality name
Location <- Specimen.ID
for(i in 1:length(SID)){
  Location = gsub(SID[i], Locations[i], Location)
}

#Outlines need to be aligned before we can subset the data, so data must be converted into an outline object with Momocs
#Turn into Out object
Sweet.Data <- Out(Sweet.Out.Full, fac = Specimen.ID)

#Center data
Outline.Data.Centered <- Sweet.Data %>% coo_center() %>% coo_scale() %>% coo_slidedirection(direction = "up")

#Sample equally spaced points around the outline, argument is set to sample 50 points
Sub.Set.Outline <- coo_sample(Outline.Data.Centered, 50)

#Turn aligned specimens into an array for use in Morpho
Outline.array <- abind(Sub.Set.Outline[[1]], along = 3)

#General procrustes analysis in Morpho using sampled outlines.
Outline.gpa <- procSym(Outline.array, scale = 1)

#Dataset for further analysis
plotly.points <- data.frame(Outline.gpa$PCscores[,1:14], factor(Specimen.ID), Outline.gpa$size, Outline.gpa$rho)
Full <- cbind(plotly.points, Specimen.ID, Denticle.ID, Species.ID, Location)

#Input information on colour and groupings for plotting.
pal = c("red","blue","green","purple", "black", "orange")
Species = c("whitei", "asymmetrica","cf. adenticulatus", "aff. binodosus","anceps", "binodosus", "clarki", "obliquidentatus")
Dent1 = c("D1","D2","D3","D4","D5","D6","D7","D8")
Dent = c("D1","D2","D3","D4","D5","D6","D7","D8","D9") #This is all the denticles, but not used as most species only have 8 denticles.

#Plots of principal component data
library(ggplot2)
#Plots of PCA data grouped by denticle position (Supplementary Information figure S4)
plots = function(x){Full %>% group_by(Species.ID) %>%
    filter(Denticle.ID == x) %>% 
    plot_ly(
      x = ~PC1, 
      y = ~PC2,
      color = ~Species.ID,
      symbol = ~Location,
      colors = pal,
      text = ~factor.Specimen.ID.,
      type = "scatter",
      mode = "markers",
      showlegend = F) %>%
    layout(xaxis = list(title = "PC1", 
                        titlefont = F, 
                        showgrid = T), 
           yaxis = list(title = "PC2", 
                        titlefont = F, 
                        showgrid = T),
           title = "Full Data Set PCA")} #Make sure to highlight and run full block at once, plotly will give an error message otherwise.
subfp = lapply(Dent, plots)
p = subplot(subfp,nrows = 3)
print(p)
#Orca was used to save images to SVG files for import into Adobe to create publication figures.
#orca(p, "Full_pca.svg")

#Testing for group differences between species. Results presented in supplementary information table S2 and table S3.
library(conover.test)
library(dunn.test)
#There are a number of different tests I ran on the PC data 

#Permutation test comparing the distance between group means obtanined by random assigment to groups. Groups are based on specimen ID and run through each denticle (1 to)
permu.Sw <- function(x){
  a = Full %>% filter(Denticle.ID == x)
  u = permudist(data = a$PC1 + a$PC2 , groups = a$Species.ID, rounds = 10000)
  print(u)
}

sapply(Dent1, permu.Sw)

#Anova and TukeyHSD comparing group means based on a classifier

#Grouped based on Location
#To cycle through different PC scores used in ANOVA change res.aov = aov(a$[PCs you want to use] ~ a$Location)
#If you want to view ANOVA results, delete TukeyHSD(res.aov) and replace with summary(res.aov)
Location.aov = function(x){
  a = Full %>% filter(Denticle.ID == x)
  res.aov = aov(a$PC1 ~ a$Location)
  TukeyHSD(res.aov)
}

sapply(Dent1, Location.aov)

#Grouped based on Species 
#To cycle through different PC scores used in ANOVA change res.aov = aov(a$[PCs you want to use] ~ a$Species.ID)
#If you want to view ANOVA results, delete TukeyHSD(res.aov2) and replace with summary(res.aov2)
Species.aov = function(x){
  a = Full %>% filter(Denticle.ID == x)
  res.aov2 = aov(a$PC1 ~ a$Species.ID)
  TukeyHSD(res.aov2)
}

sapply(Dent1, Species.aov)

#Two-way Anova testing Species and Location classifier, then Species in combination with location
#Just change the filter criteria for results for difference denticles
a = Full %>% filter(Denticle.ID == "D1")
model = lm(a$PC1 ~ a$Species.ID + a$Location + a$Species.ID:a$Location)
anova(model)
summary(model)

a = Full
lrfit <- glm(PC1 ~ Species.ID + Location, data = a, family = gaussian)
summary(lrfit)
x = anova(lrfit)
plot(x)

# Using Kruskal-Wallis test and Conover-Iman Test; results for Conover-Iman Test are only valid if results from the Kruskal-Wallis test are signigicant. 
#Kruskal-Wallis serves as a nonparametic alternative to one way ANOVA

#My major concern with reporting any of these results is that I could not decide which p.adjustment method to employ. See p.adjustment.methods for documentation. For all of the tests no p-value adjustments were made. Conover.test does require p-value = alpha/2 for the results to be considered significant.

#Conover-Iman Test with Species as the group classifier
cono_test_species = function(x){
  a = Full %>% filter(Denticle.ID == x)
  conover.test(a$PC1, a$Species.ID, method = "none",T, list = T)
}

sapply(Dent1, cono_test_species)

#Conover-Iman Test with Location as the group classifier
cono_test_location = function(x){
  a = Full %>% filter(Denticle.ID == x)
  conover.test(a$PC1, a$Location, method = "none",T, list = T)
}

sapply(Dent1, cono_test_location)

#Finally I ran Dunn's test as it is also a nonparametric test to compare groups means. 
dunn_test_species = function(x){
  a = Full %>% filter(Denticle.ID == x)
  dunn.test(a$PC1, a$Species.ID,"none",T, list = T)
}

sapply(Dent1, dunn_test_species)

#There is a statistical difference in denticle morphology in Species. The results are all similar between the different parametric and nonparametric tests run. 

#Cluster analysis
#Create dataframe of mean values
#This would technically be employing a mean field approach when studying the between-species trait distibution (doi: 10.1016/j.tree.2011.11.014). We would not be able to conclude on intra-species trait distributions based on the cluster results, we would need to look at the original PCA or RW plots.  
library(pvclust)
#Specimen = sprintf("S%d", seq(1:59))

#To change which PCs you are using in the analysis change Full$[which ever PCs you want to use]
mean.pc1 <- function(x,y){
  mean(Full$PC1[Denticle.ID == x & Species.ID == y])
}

mean.pc2 <- function(x,y){
  mean(Full$PC2[Denticle.ID == x & Species.ID == y])
}

pc1.mean <- function(y){
  sapply(Dent1, mean.pc1, y = y)
}

pc2.mean <- function(y){
  sapply(Dent1, mean.pc2, y = y)
}

#PC1 mean
mSPC1 = sapply(Species, pc1.mean)
#PC2 mean
mSPC2 = sapply(Species, pc2.mean)

#cluster analysis for PC1
clust_PC1 = pvclust(mSPC1, method.hclust= "ward.D2", method.dist = "euclidean", iseed = 1, nboot = 1000)
plot(clust_PC1) #Supplementary Information figure 4j.
pvrect(clust_PC1)


#cluster analysis for PC2
clust_PC2 = pvclust(mSPC2, method.hclust= "ward.D2", method.dist = "euclidean", iseed = 1, nboot = 1000)
plot(clust_PC2) #Supplementary Information figure 4k.
pvrect(clust_PC2)

#Selected plots of clusters were exported into SVG files using orca(clust_PC1, "clust_PC1.svg")
#Combination of PC1 and PC2 weighted on variance.
o <- pvclust((mSPC1*0.74)+(mSPC2*0.12), method.hclust= "ward.D2", method.dist = "euclidean", iseed = 1, nboot = 10000)
plot(o)
pvrect(o)
msplot(o)

#Species Metrics

#Total trait distribution
PC1sd = var(Full$PC1)
PC2sd = var(Full$PC2)

OrderedSpecies = c("anceps","binodosus","aff. binodosus","whitei","cf. adenticulatus","asymmetrica","clarki","obliquidentatus")

#I provide examples only for the first PC. If you wish to calculate PC2, change input to Full$PC2

#range
Range.spec = function(x) {
  max(Full$PC1[Species.ID == x]) - min(Full$PC1[Species.ID == x])
}
#variance
var.cono1 = function(x) {
  var(Full$PC1[Species.ID == x])
}
#standard deviation
sd.cono1 = function(x) {
  sqrt(var(Full$PC1[Species.ID == x]))
}
#mean
mean.cono1 = function(x) {
  mean(Full$PC1[Species.ID == x])
}
#trait means to intraspecific trait widths
wfstat = function(x) {
  mean(Full$PC1[Species.ID == x])/sqrt(var(Full$PC1[Species.ID == x]))
}
#absolute values of trait means to intrapspecific trait widths
absmfsw1 = function(x) {
  abs(mean(Full$PC1[Species.ID == x]))/abs(sqrt(var(Full$PC1[Species.ID == x])))
}
#coefficent of variance
coefofvar1 = function(x) {
  sqrt(var(Full$PC1[Species.ID == x]))/mean(Full$PC1[Species.ID == x])
}

#This is cumbersome but change Species.ID for each speices to calculate limiting similarity. Then run the sapply function below. Do this for each speices.
LS = function(x) {
  (mean(Full$PC1[Species.ID == "anceps"])-mean(Full$PC1[Species.ID == x]))/(sqrt(var(Full$PC1[Species.ID == "anceps"]))-sqrt(var(Full$PC1[Species.ID == x])))
}


limit = sapply(OrderedSpecies, LS) #anceps
limit1 = sapply(OrderedSpecies, LS) #binodosus
limit2 = sapply(OrderedSpecies, LS) #aff. binodosus
limit3 = sapply(OrderedSpecies, LS) #whitei
limit4 = sapply(OrderedSpecies, LS) #cf. adenticulatus
limit5 = sapply(OrderedSpecies, LS) #asymmetrica
limit6 = sapply(OrderedSpecies, LS) #clarki
limit7 = sapply(OrderedSpecies, LS) #obliquidentatus

limitsimilarity = data.frame(limit,limit1,limit2,limit3,limit4,limit5,limit6,limit7)

colnames(limitsimilarity) <- c("Sw. anceps","Sw. binodosus","Sw. aff. binodosus","Sw. whitei","Sw. cf. adenticulatus","Sw. asymmetrica","Sw. clarki","Sw. obliquidentatus" )
rownames(limitsimilarity) <- c("Sw. anceps","Sw. binodosus","Sw. aff. binodosus","Sw. whitei","Sw. cf. adenticulatus","Sw. asymmetrica","Sw. clarki","Sw. obliquidentatus" ) 

View(limitsimilarity) #table 2

rangeSW = sapply(OrderedSpecies, Range.spec)
rangeSW = sapply(OrderedSpecies, Range.spec)
rangpc1 = max(Full$PC1) - min(Full$PC1)
commrange = rangeSW/rangpc1
varSw = sapply(OrderedSpecies, var.cono1)
sdSw = sapply(OrderedSpecies, sd.cono1)
mfSw = sapply(OrderedSpecies, wfstat)
meanSW = sapply(OrderedSpecies, mean.cono1)
absmfSw = sapply(OrderedSpecies, absmfsw1)
coefvar = sapply(OrderedSpecies, coefofvar1)
communityvarPC1 = varSw/PC1sd #community-wide variance
traitspace = c("1","2","2","2","3","3","4","4")
IntraSweet <-data.frame(OrderedSpecies,varSw,sdSw,meanSW,mfSw,absmfSw,coefvar,traitspace,communityvarPC1, commrange)
IntraSweet$OrderedSpecies <- factor(IntraSweet$OrderedSpecies, levels = OrderedSpecies)
View(IntraSweet) #table 1

#This was the glm used in the paper
modeltrait = glm(IntraSweet$mfSw~IntraSweet$traitspace)
plot(modeltrait)
summary(modeltrait)
traitaov = aov(modeltrait)
summary(traitaov)

limitdata <- cbind(limitsimilarity, traitspace)
model <- glm(IntraSweet$meanSW ~ IntraSweet$traitspace)
summary(model)

#Evolutionary Trajectory
#You will need to run the above code from line 260 downwards to calculate metrics for PC2 before plotting trajectories.
fulltraj = plot_ly(
  x = IntraSweet$meanSW,
  y = IntraSweet2$meanSW2,
  color = IntraSweet$traitspace,
  type = "scatter",
  mode = "markers",
  showlegend = F,
  text = IntraSweet$OrderedSpecies) %>%
  add_text(textposition = "top right") %>%
  layout(xaxis = list(title = "PC1", 
                      titlefont = F, 
                      showgrid = T), 
         yaxis = list(title = "PC2", 
                      titlefont = F, 
                      showgrid = T))

###orca(fulltraj, file = "fulltrajectory_sweetognathus.svg")


###

#R code for sensitivity testing presented in the supplementary information 
#Import of rds dataset, make sure to set proper working directory
Outline.Data.re <- readRDS(file = "Full_Data_OutlineSense.rds")

Specimen.ID2 = c("1D_38_2","1D_38_3","1D_38","1D_40_2","1D_40_3","1D_40",
                 "2D_38_2","2D_38_3","2D_38","2D_40_2","2D_40_3","2D_40",
                 "3D_38_2","3D_38_3","3D_38","3D_40_2","3D_40_3","3D_40",
                 "4D_38_2","4D_38_3","4D_38","4D_40_2","4D_40_3","4D_40",
                 "5D_38_2","5D_38_3","5D_38","5D_40_2","5D_40_3","5D_40")

Specimen.ID3 <- attributes(Outline.Data.re)

groupID2 <- c("A","A","A","A","A","A","A","A","A","A","A","A","A","A",
              "B","B","B","B","B","B","B","B","B","B","B","B","B","B",
              "C","C","C","C","C","C","C","C","C","C","C","C","C","C",
              "D","D","D","D","D","D","D","D","D","D","D","D","D",
              "E","E","E","E","E","E","E","E","E","E","E","E","E")

Actual.Spec.ID2 <-  c("15","20","31","33_5","33","38_2","38_3","38","39_2","39_3","39","40_2","40_3","40",
                      "15","20","31","33_5","33","38_2","38_3","38","39_2","39_3","39","40_2","40_3","40",
                      "15","20","31","33_5","33","38_2","38_3","38","39_2","39_3","39","40_2","40_3","40",
                      "15","20","33_5","33","38_2","38_3","38","39_2","39_3","39","40_2","40_3","40",
                      "15","20","33_5","33","38_2","38_3","38","39_2","39_3","39","40_2","40_3","40")

Test_Coo <- Out(Outline.Data.re, fac = groupID2)
Test_Coo2 <- Out(Outline.Data.re, fac = Actual.Spec.ID2)

Test_Coo3 <- Out(Outline.Data.re, fac = Specimen.ID3)

Dent_Test2 <- Test_Coo2 %>% coo_center() %>% coo_scale() %>% coo_slidedirection(direction = "up")

pile(Dent_Test2)
panel(Dent_Test2, names = T)

#Semi-landmark analysis of cross-sections

SD_sample <- coo_sample(Dent_Test2, 50)
#Dent_Test belongs above

SD2 <- array(data = c(SD_sample[[1]]$`1D_15`[,1],SD_sample[[1]]$`1D_15`[,2],
                      SD_sample[[1]]$`1D_20`[,1], SD_sample[[1]]$`1D_20`[,2],
                      SD_sample[[1]]$`1D_31`[,1], SD_sample[[1]]$`1D_31`[,2],
                      SD_sample[[1]]$`1D_33_5`[,1], SD_sample[[1]]$`1D_33_5`[,2],
                      SD_sample[[1]]$`1D_33`[,1], SD_sample[[1]]$`1D_33`[,2],
                      SD_sample[[1]]$`1D_38_2`[,1], SD_sample[[1]]$`1D_38_2`[,2],
                      SD_sample[[1]]$`1D_38_3`[,1], SD_sample[[1]]$`1D_38_3`[,2],
                      SD_sample[[1]]$`1D_38`[,1], SD_sample[[1]]$`1D_38`[,2],
                      SD_sample[[1]]$`1D_39_2`[,1], SD_sample[[1]]$`1D_39_2`[,2],
                      SD_sample[[1]]$`1D_39_3`[,1], SD_sample[[1]]$`1D_40_3`[,2],
                      SD_sample[[1]]$`1D_39`[,1], SD_sample[[1]]$`1D_39`[,2],
                      SD_sample[[1]]$`1D_40_2`[,1], SD_sample[[1]]$`1D_40_2`[,2],
                      SD_sample[[1]]$`1D_40_3`[,1], SD_sample[[1]]$`1D_40_3`[,2],
                      SD_sample[[1]]$`1D_40`[,1], SD_sample[[1]]$`1D_40`[,2],
                      
                      SD_sample[[1]]$`2D_15`[,1],SD_sample[[1]]$`2D_15`[,2],
                      SD_sample[[1]]$`2D_20`[,1], SD_sample[[1]]$`2D_20`[,2],
                      SD_sample[[1]]$`2D_31`[,1], SD_sample[[1]]$`2D_31`[,2],
                      SD_sample[[1]]$`2D_33_5`[,1], SD_sample[[1]]$`2D_33_5`[,2],
                      SD_sample[[1]]$`2D_33`[,1], SD_sample[[1]]$`2D_33`[,2],
                      SD_sample[[1]]$`2D_38_2`[,1], SD_sample[[1]]$`2D_38_2`[,2],
                      SD_sample[[1]]$`2D_38_3`[,1], SD_sample[[1]]$`2D_38_3`[,2],
                      SD_sample[[1]]$`2D_38`[,1], SD_sample[[1]]$`2D_38`[,2],
                      SD_sample[[1]]$`2D_39_2`[,1], SD_sample[[1]]$`2D_39_2`[,2],
                      SD_sample[[1]]$`2D_39_3`[,1], SD_sample[[1]]$`2D_40_3`[,2],
                      SD_sample[[1]]$`2D_39`[,1], SD_sample[[1]]$`2D_39`[,2],
                      SD_sample[[1]]$`2D_40_2`[,1], SD_sample[[1]]$`2D_40_2`[,2],
                      SD_sample[[1]]$`2D_40_3`[,1], SD_sample[[1]]$`3D_40_3`[,2],
                      SD_sample[[1]]$`2D_40`[,1], SD_sample[[1]]$`2D_40`[,2],
                      
                      SD_sample[[1]]$`3D_15`[,1],SD_sample[[1]]$`3D_15`[,2],
                      SD_sample[[1]]$`3D_20`[,1], SD_sample[[1]]$`3D_20`[,2],
                      SD_sample[[1]]$`3D_31`[,1], SD_sample[[1]]$`3D_31`[,2],
                      SD_sample[[1]]$`3D_33_5`[,1], SD_sample[[1]]$`3D_33_5`[,2],
                      SD_sample[[1]]$`3D_33`[,1], SD_sample[[1]]$`3D_33`[,2],
                      SD_sample[[1]]$`3D_38_2`[,1], SD_sample[[1]]$`3D_38_2`[,2],
                      SD_sample[[1]]$`3D_38_3`[,1], SD_sample[[1]]$`3D_38_3`[,2],
                      SD_sample[[1]]$`3D_38`[,1], SD_sample[[1]]$`3D_38`[,2],
                      SD_sample[[1]]$`3D_39_2`[,1], SD_sample[[1]]$`3D_39_2`[,2],
                      SD_sample[[1]]$`3D_39_3`[,1], SD_sample[[1]]$`3D_40_3`[,2],
                      SD_sample[[1]]$`3D_39`[,1], SD_sample[[1]]$`3D_39`[,2],
                      SD_sample[[1]]$`3D_40_2`[,1], SD_sample[[1]]$`3D_40_2`[,2],
                      SD_sample[[1]]$`3D_40_3`[,1], SD_sample[[1]]$`3D_38_3`[,2],
                      SD_sample[[1]]$`3D_40`[,1], SD_sample[[1]]$`3D_40`[,2],
                      
                      SD_sample[[1]]$`4D_15`[,1], SD_sample[[1]]$`4D_15`[,2],
                      SD_sample[[1]]$`4D_20`[,1], SD_sample[[1]]$`4D_20`[,2],
                      SD_sample[[1]]$`4D_33_5`[,1], SD_sample[[1]]$`4D_33_5`[,2],
                      SD_sample[[1]]$`4D_33`[,1], SD_sample[[1]]$`4D_33`[,2],
                      SD_sample[[1]]$`4D_38_2`[,1], SD_sample[[1]]$`4D_38_2`[,2],
                      SD_sample[[1]]$`4D_38_3`[,1], SD_sample[[1]]$`4D_38_3`[,2],
                      SD_sample[[1]]$`4D_38`[,1], SD_sample[[1]]$`4D_38`[,2],
                      SD_sample[[1]]$`4D_39_2`[,1], SD_sample[[1]]$`4D_39_2`[,2],
                      SD_sample[[1]]$`4D_39_3`[,1], SD_sample[[1]]$`4D_39_3`[,2],
                      SD_sample[[1]]$`4D_39`[,1], SD_sample[[1]]$`4D_39`[,2],
                      SD_sample[[1]]$`4D_40_2`[,1], SD_sample[[1]]$`4D_40_2`[,2],
                      SD_sample[[1]]$`4D_40_3`[,1], SD_sample[[1]]$`4D_40_3`[,2],
                      SD_sample[[1]]$`4D_40`[,1], SD_sample[[1]]$`4D_40`[,2],
                      
                      SD_sample[[1]]$`5D_15`[,1], SD_sample[[1]]$`5D_15`[,2],
                      SD_sample[[1]]$`5D_20`[,1], SD_sample[[1]]$`5D_20`[,2],
                      SD_sample[[1]]$`5D_33_5`[,1], SD_sample[[1]]$`5D_33_5`[,2],
                      SD_sample[[1]]$`5D_33`[,1], SD_sample[[1]]$`5D_33`[,2],
                      SD_sample[[1]]$`5D_38_2`[,1], SD_sample[[1]]$`5D_38_2`[,2],
                      SD_sample[[1]]$`5D_38_3`[,1], SD_sample[[1]]$`5D_38_3`[,2],
                      SD_sample[[1]]$`5D_38`[,1], SD_sample[[1]]$`5D_38`[,2],
                      SD_sample[[1]]$`5D_39_2`[,1], SD_sample[[1]]$`5D_39_2`[,2],
                      SD_sample[[1]]$`5D_39_3`[,1], SD_sample[[1]]$`5D_39_3`[,2],
                      SD_sample[[1]]$`5D_39`[,1], SD_sample[[1]]$`5D_39`[,2],
                      SD_sample[[1]]$`5D_40_2`[,1], SD_sample[[1]]$`5D_40_2`[,2],
                      SD_sample[[1]]$`5D_40_3`[,1], SD_sample[[1]]$`5D_40_3`[,2],
                      SD_sample[[1]]$`5D_40`[,1], SD_sample[[1]]$`5D_40`[,2]),
             dim = c(50,2,68))

#Procrustes analysis
SD.gpa = procSym(SD2)

#Principle component plots for first and second principle components and first, second, and third principle components.

SD.pca.df <- data.frame(SD.gpa$PCscores[,1],SD.gpa$PCscores[,2], SD.gpa$PCscores[,3], factor(Actual.Spec.ID2))
pal = c("red","blue","green","purple","black","yellow")

SD.pca.df %>% plot_ly(
  x = ~SD.pca.df$SD.gpa.PCscores...1., 
  y = ~SD.pca.df$SD.gpa.PCscores...2.,
  color =  Actual.Spec.ID2,
  colors = pal,
  type = "scatter",
  mode = "markers") %>%
  layout(xaxis = list(title = "PC1", 
                      titlefont = F, 
                      showgrid = T), 
         yaxis = list(title = "PC2", 
                      titlefont = F, 
                      showgrid = T)) 

SD.pca.df %>% plot_ly(
  x = ~SD.pca.df$SD.gpa.PCscores...1., 
  y = ~SD.pca.df$SD.gpa.PCscores...2.,
  z = ~SD.pca.df$SD.gpa.PCscores...3.,
  color = Actual.Spec.ID2,
  colors = pal,
  text = Actual.Spec.ID2,
  hoverinfo = "text",
  type = "scatter3d",
  mode = "markers") %>%
  layout(title = "Outline Landmark Morphospace",
         scene = list(xaxis = list(title = "PC1",
                                   titlefont = F,
                                   showgrid = T),
                      yaxis = list(title = "PC2", 
                                   titlefont = F,
                                   showgrid = T),
                      zaxis = list(title = "PC3",
                                   titlefont = F,
                                   showgrid = T))) %>%
  add_text(textfont =  list( family = "sans serif",
                             size = 14,
                             color = toRGB("grey50")),
           textposition = "topright",
           showlegend = F)

#Thin-spline deformations in morphospace to check for inconsistences. 
ref4 = SD.gpa$mshape
plotRefToTarget(ref4, SD.gpa$orpdata[,,59])
plotRefToTarget(ref4, SD.gpa$orpdata[,,64])
plotRefToTarget(ref4, SD.gpa$orpdata[,,65])
#33_5
plotRefToTarget(ref4, SD.gpa$orpdata[,,32])

#33
plotRefToTarget(ref4, SD.gpa$orpdata[,,5])
#15
plotRefToTarget(ref4, SD.gpa$orpdata[,,6])
#31
plotRefToTarget(ref4, SD.gpa$orpdata[,,56])
#15
plotRefToTarget(ref4, SD.gpa$orpdata[,,1])
#20
plotRefToTarget(ref4, SD.gpa$orpdata[,,16])
#33_5
plotRefToTarget(ref4, SD.gpa$orpdata[,,45])
#38
plotRefToTarget(ref4, SD.gpa$orpdata[,,62])
#38_2
plotRefToTarget(ref4, SD.gpa$orpdata[,,6])
#39_3
plotRefToTarget(ref4, SD.gpa$orpdata[,,38])

#Anova and Tukey's Post Hoc Test for full data set. Null hypothesis is that there is no difference between specimen groups means.

res.aov <- aov(SD.gpa$PCscores[,1] ~ Actual.Spec.ID2)
res.aov2 <- aov(SD.gpa$PCscores[,2] ~ Actual.Spec.ID2)
res.aov3 <- aov(SD.gpa$PCscores[,2] ~ groupID2)
summary(res.aov)
summary(res.aov2)
summary(res.aov3)
a = TukeyHSD(res.aov)
b = TukeyHSD(res.aov2)
c = TukeyHSD(res.aov3)
print(a)
print(b)
print(c)

#Boxplot of first principle component scores grouped by specimen ID
plot_ly(y = ~SD.gpa.3$PCscores[,1], color = ~Specimen.ID3$names, type = "box")  

#Boxplot of first principle component scores grouped by cross-section angle. A = 0˚, B = 2.5˚, C = -2.5˚, D = 5˚, E = -5˚
plot_ly(y = ~SD.gpa$PCscores[,1], color = ~groupID2, type = "box")
