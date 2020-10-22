#title: R Code for Evidence of Parallel Evolution in the Dental Elements of Sweetognathus Conodonts
#author: "Wyatt Petryshen"
#date: '2020-09-24'

#Required packages for landmark selection and Procrustes superimposition, and some data visulization later on. Analysis was run using following versions of each package: Momocs (1.2.9), abind (1.4-5), Morpho (2.7), tidyverse (1.3.0), plotly (4.9.1), car (3.0-7). 
library(Momocs); library(abind); library(tidyverse); library(plotly); library(car); library(ggplot2); library(pvclust); library(geomorph); library(RANN); library(RRPP); library(pbapply)
setwd()
#Import all jpegs into R studio still using Momocs (the program doesn't fuck this up which is nice, and provides a outline 1 pixel thick... for the most part)
#First step is to import data into R markdown document, then gather information of ID and denticle positions prior to procrustes superimposition.
#Import of rds dataset, make sure to set proper working directory; this is of the jpeg outlines converted into xy-coordinates.This will be sampled in Momocs prior to procrustes superimposition.
Sweet.Out.Full <- readRDS(file = "Sweetognathus.Outlines.rds")
#Finds rows with min value of y
# data = closed outline
# y = y value's in data[x,y]
# x = logical vector of min(data[,y])
###
input = Sweet.Out.Full #Just created an extra data set so I don't mess with the original
#Delete duplicate outlines
appended.data <- input[-c(27,93,206)] #removal of double outlines
### This does what it should expect for when outline is more than two pixals thick.
curve.comp <- function(test.data){
  repeat {
    x = apply(test.data, 1, function(y) any(y %in% min(test.data[,2])))
    test.data = test.data[x == "FALSE",]
    if(table(x)["TRUE"] <= 3) {
      repeat {
        x = apply(test.data, 1, function(y) any(y %in% min(test.data[,2])))
        test.data = test.data[x == "FALSE",]
        if(table(x)["TRUE"] <= 3) {return(test.data)} #Might change this to <= 2
        break
      }
    }
  }
}

OpenData <-  pbapply::pblapply(appended.data, curve.comp) #New Open Curves
###
#Get specimen ID's with denticle position for each outline.
sample.names.full <- attributes(OpenData)
#Get specimen ID's for each set of denticles.
Specimen.ID <- str_sub(sample.names.full$names, start = 4, end = -1)
#Get denticle position for each outline.
Denticle.ID <- str_sub(sample.names.full$names, start = 1, end =2)
#Assign species identity and locality to data
{SID <- c("\\bS1\\b", "\\bS2\\b", "\\bS3\\b", "\\bS4\\b", "\\bS5\\b",
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
  
  Lineage <- c("Late Asselian","Late Asselian","Late Asselian","Late Asselian","Late Asselian",
               "Late Asselian","Late Asselian","Late Asselian","Sakmarian","Sakmarian","Sakmarian",
               "Sakmarian","Sakmarian","Late Asselian","Late Asselian","Late Asselian","Late Asselian","Late Asselian","Late Asselian",
               "Late Asselian","Sakmarian","Late Asselian","Late Asselian","Late Asselian","Late Asselian","Sakmarian","Sakmarian",
               "Early Asselian", "Early Asselian","Early Asselian","Early Asselian","Early Asselian","Early Asselian","Early Asselian","Early Asselian","Early Asselian",
               "Late Asselian","Sakmarian", "Late Asselian","Sakmarian","Late Asselian","Sakmarian",
               "Sakmarian","Sakmarian","Late Asselian",
               "Sakmarian","Sakmarian","Sakmarian","Sakmarian","Sakmarian","Sakmarian","Sakmarian",
               "Sakmarian","Sakmarian","Sakmarian","Sakmarian","Sakmarian","Sakmarian","Sakmarian")
  
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
  
  HCTraitGroups <- c("T2","T2","T2","T2","T2",
                     "T2","T2","T2", "T4", "T5", "T5",
                     "T5", "T5", "T2","T2","T2", "T2","T2","T2",
                     "T5","T4", "T5",  "T5",  "T5",  "T5","T4", "T4",
                     "T3", "T2", "T2", "T2", "T2", "T2", "T3",  "T3", "T2",
                     "T5", "T2",  "T5", "T2",  "T5", "T2",
                     "T2","T2","T5","T2","T2","T2","T2","T2","T2","T2",
                     "T1", "T2", "T1", "T2", "T1", "T1", "T1")
  
  Dent = c("D1","D2","D3","D4","D5","D6","D7","D8","D9") #This is all the denticles, but not used as most species only have 8 denticles.
  Species = c("whitei", "asymmetrica","cf. adenticulatus", "aff. binodosus","anceps", "binodosus", "clarki", "obliquidentatus")
  pal = c("red","blue","green","purple", "black", "orange")
}
#Assign species name, locality, and linage
{
  Species.ID <- Specimen.ID
  for(i in 1:length(SID)){
    Species.ID = gsub(SID[i], RID[i], Species.ID)
  }
  #Assign locality name
  Location <- Specimen.ID
  for(i in 1:length(SID)){
    Location = gsub(SID[i], Locations[i], Location)
  }
  
  #Assign by lineage
  LineageGroups <- Specimen.ID
  for(i in 1:length(SID)){
    LineageGroups = gsub(SID[i], Lineage[i], LineageGroups)
  }
  
  #Trait groupings found from HC
  TraitGroups <- Specimen.ID
  for(i in 1:length(SID)){
    TraitGroups = gsub(SID[i], HCTraitGroups[i], TraitGroups)
  }
}
#GPA in geomorph
library(geomorph)
#Digitize Landmarks

#Lets collect semi-landmarks.
#First arrange curves such that {x,y} points along each curve are in order from point nearest to origin {0,0} to largest x-value
### KD-Tree Knn approach
library(RANN) #sort.curve is dependent on the RNNN::nn2() function
RANNsort.curve <- function(x){
  RANNinputdata = x #function input curve
  RANNx = RANNinputdata[RANNinputdata[,2] == min(RANNinputdata[,2]),] #find min values of y and correspounding x coordinate
  RANNx = list(RANNx) #place min values into a list
  RANNy = RANNx[[1]][RANNx[[1]][,1] == min(RANNx[[1]][,1]),] #find start of curve coordinates close to {0,0}
  #Create a new list with this value
  RANNsorted.curve = matrix(RANNy, nrow = 1, ncol = 2) #sorted curve begins with start position closest to {0,0}
  RANNdata.search = RANNinputdata #unsorted curve to search through
  repeat {
    RANNa = tail(RANNsorted.curve,1)#find starting point for sorting
    RANNvaluex = RANNa[1] #set value x to remove
    RANNvaluey = RANNa[2] #set value y to remove
    RANNvalues = c(RANNvaluex,RANNvaluey) #vector position start
    #removes already sorted value
    RANNrow.to.remove = which((RANNdata.search[,1] == RANNvaluex) & (RANNdata.search[,2] == RANNvaluey)) #find indices in unsorted curve for valuex and valuey
    RANNdata.search = RANNdata.search[-RANNrow.to.remove,] #delete coordinates from unsorted row given indices
    #some knn stuff
    RANNvals = nn2(RANNdata.search,query = RANNa, k = min(5, nrow(RANNdata.search)), treetype = "kd") #KD-Tree KNN Algorithm.
    RANN.1 = cbind(t(RANNvals$nn.dists),t(RANNvals$nn.idx)) #combine KDTree-Knn output into a single table
    RANN2.1= as.numeric(min(RANN.1[,1])) #find out line distance
    #placing new values into sorted curve
    RANNvector.index = RANN.1[which(RANN.1[1,] == min(RANN2.1)),2]  #find index for smallest vector
    RANNvector.index = na.omit(RANNvector.index)
    RANNnewrow <-  matrix(RANNdata.search[RANNvector.index[1],], nrow = 1,ncol = 2) #new sorted row to add
    RANNsorted.curve <- rbind(RANNsorted.curve, RANNnewrow) #puts new row into sorted curve list
    
    if(nrow(RANNdata.search) == 2) {
      RANNsorted.curve = RANNsorted.curve[c(nrow(RANNsorted.curve)-4:nrow(RANNsorted.curve)),]
      return(RANNsorted.curve)
      break
    }
  }
}

###Sort curves
sorted.curve <- pbapply::pblapply(OpenData, RANNsort.curve)
#To verify that the curves are sorted I used the recommendations outlined in the Read.Me geomorph file page 30
#Alternatively a line plot using plotly will yield the same results
plot(sorted.curve[[1]], pch =19, cex = 0.1, col = rainbow(nrow(sorted.curve[[1]])))


#Choose Landmark positions
#Now we select semi-landmarks around the curve using geomorph
#calculate.semilandmarks is dependent on the geomorph::digit.curves() function
calculate.semilandmarks <- function(i){
  output = digit.curves(start = i[1,],
                        curve = i,
                        nPoints = 50,
                        closed = F)
  return(output)
}

Semi.Landmarks <- pbapply::pblapply(sorted.curve, calculate.semilandmarks)
#Format that gpagen can use. 
Semi.Landmarks <- abind::abind(Semi.Landmarks, along = 3)
###Procrustes Superimposition
Semi.Landmarks #Semi-landmark coordinates for each specimen
#Defined points that will slide
#sliders <- define.sliders(Semi.Landmarks[,1:2,1],50)
#saveRDS(sliders, file = "sliders.rds")
sliders <- readRDS(file = "sliders.rds")

#This uses a mix of landmarks (start and end point along the outline) and semi-landmarks (50 equally spaced points between start and end point)
S.gpa <- gpagen(
  A = Semi.Landmarks,
  curves = sliders,
  Proj = T,
  print.progress = T
) 

#PCA of the aligned points
#If the function gm.prcomp cannot be located try downloading an older version of geomorph using the hashed code below
                  #install.packages("devtools")
                  #library(devtools)
                  #install_version("geomorph", version = "3.2.0", repos = "http://cran.us.r-project.org")
#If this problem is presistant, download RDS file of saved results with below hashed code
                  #S.pca <– readRDS("S.pca.rds")
                  
S.pca <- gm.prcomp(S.gpa$coords)
summary(S.pca)
plot(S.pca)

#Seperate PCA plots grouped by denticle location
{library(plotly)
  PCA.values <- data.frame(S.pca$x[,1],S.pca$x[,2],Species.ID,Specimen.ID,Location,Denticle.ID, LineageGroups)
  plots = function(x){PCA.values %>% group_by(Species.ID) %>%
      filter(Denticle.ID == x) %>% 
      plot_ly(
        x = ~S.pca.x...1., 
        y = ~S.pca.x...2.,
        color = ~Species.ID,
        symbol = ~Denticle.ID,
        colors = pal,
        text = ~Species.ID,
        type = "scatter",
        mode = "markers",
        showlegend = T) %>%
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
  
  PCA.values %>% 
    plot_ly(
      x = ~S.pca.x...1., 
      y = ~S.pca.x...2.,
      color = ~Species.ID,
      symbol = ~Denticle.ID,
      colors = pal,
      text = ~Species.ID,
      type = "scatter",
      mode = "markers+line",
      showlegend = T) %>%
    layout(xaxis = list(title = "PC1", 
                        titlefont = F, 
                        showgrid = T), 
           yaxis = list(title = "PC2", 
                        titlefont = F, 
                        showgrid = T),
           title = "Full Data Set PCA")
}
#Find Index for Denticles
{
  Denticle_Index = NULL
  Find.Index <- function(x){
    Index_x <- which(PCA.values$Denticle.ID == x)
    Denticle_Index = rbind(Index_x, Denticle_Index)
  }
  Denticle_Index = lapply(Dent, Find.Index)
  names(Denticle_Index) <- paste("Index", Dent, sep = "_")
}
#Find Index for Species 
{
  Species_Index = NULL
  Find.Index.Species <- function(x){
    Index_x <- which(PCA.values$Species.ID == x)
    Species_Index = rbind(Index_x, Species_Index)
  }
  Species_Index = lapply(Species, Find.Index.Species)
  names(Species_Index) <- paste("Index", Species, sep = "_")
}
#Find Index for Species and Denticles
{
  #whitei
  whitei_Denticle_Index = NULL
  whitei.Index.Species.Denticle <- function(x){
    Index_x <- which(PCA.values$Denticle.ID == x & PCA.values$Species.ID== "whitei")
    whitei_Denticle_Index = rbind(Index_x, whitei_Denticle_Index)
  }
  whitei_Denticle_Index <- lapply(Dent, whitei.Index.Species.Denticle)
  names(whitei_Denticle_Index) <- paste("whitei", Dent, sep = "_")
  
  #asymmetrica
  asymmetrica_Denticle_Index = NULL
  asymmetrica.Index.Species.Denticle <- function(x){
    Index_x <- which(PCA.values$Denticle.ID == x & PCA.values$Species.ID== "asymmetrica")
    asymmetrica_Denticle_Index = rbind(Index_x, asymmetrica_Denticle_Index)
  }
  asymmetrica_Denticle_Index <- lapply(Dent, asymmetrica.Index.Species.Denticle)
  names(asymmetrica_Denticle_Index) <- paste("asymmetrica", Dent, sep = "_")
  
  #cf. adenticulatus
  cf.adenticulatus_Denticle_Index = NULL
  cf.adenticulatus.Index.Species.Denticle <- function(x){
    Index_x <- which(PCA.values$Denticle.ID == x & PCA.values$Species.ID == "cf. adenticulatus")
    cf.adenticulatus_Denticle_Index = rbind(Index_x, cf.adenticulatus_Denticle_Index)
  }
  cf.adenticulatus_Denticle_Index <- lapply(Dent, cf.adenticulatus.Index.Species.Denticle)
  names(cf.adenticulatus_Denticle_Index) <- paste("cf. adenticulatus", Dent, sep = "_")
  
  #aff. binodosus
  aff.binodosus_Denticle_Index = NULL
  aff.binodosus.Index.Species.Denticle <- function(x){
    Index_x <- which(PCA.values$Denticle.ID == x & PCA.values$Species.ID== "aff. binodosus")
    aff.binodosus_Denticle_Index = rbind(Index_x, aff.binodosus_Denticle_Index)
  }
  aff.binodosus_Denticle_Index <- lapply(Dent, aff.binodosus.Index.Species.Denticle)
  names(aff.binodosus_Denticle_Index) <- paste("aff. binodosus", Dent, sep = "_")
  
  #anceps
  anceps_Denticle_Index = NULL
  anceps.Index.Species.Denticle <- function(x){
    Index_x <- which(PCA.values$Denticle.ID == x & PCA.values$Species.ID== "anceps")
    anceps_Denticle_Index = rbind(Index_x, anceps_Denticle_Index)
  }
  anceps_Denticle_Index <- lapply(Dent, anceps.Index.Species.Denticle)
  names(anceps_Denticle_Index) <- paste("anceps", Dent, sep = "_")
  
  #binodosus
  binodosus_Denticle_Index = NULL
  binodosus.Index.Species.Denticle <- function(x){
    Index_x <- which(PCA.values$Denticle.ID == x & PCA.values$Species.ID== "binodosus")
    binodosus_Denticle_Index = rbind(Index_x, binodosus_Denticle_Index)
  }
  binodosus_Denticle_Index <- lapply(Dent, binodosus.Index.Species.Denticle)
  names(binodosus_Denticle_Index) <- paste("binodosus", Dent, sep = "_")
  
  #clarki
  clarki_Denticle_Index = NULL
  clarki.Index.Species.Denticle <- function(x){
    Index_x <- which(PCA.values$Denticle.ID == x & PCA.values$Species.ID== "clarki")
    clarki_Denticle_Index = rbind(Index_x, clarki_Denticle_Index)
  }
  clarki_Denticle_Index <- lapply(Dent, clarki.Index.Species.Denticle)
  names(clarki_Denticle_Index) <- paste("clarki", Dent, sep = "_")
  
  #obliquidentatus
  obliquidentatus_Denticle_Index = NULL
  obliquidentatus.Index.Species.Denticle <- function(x){
    Index_x <- which(PCA.values$Denticle.ID == x & PCA.values$Species.ID== "obliquidentatus")
    obliquidentatus_Denticle_Index = rbind(Index_x, obliquidentatus_Denticle_Index)
  }
  obliquidentatus_Denticle_Index <- lapply(Dent, obliquidentatus.Index.Species.Denticle)
  names(obliquidentatus_Denticle_Index) <- paste("obliquidentatus", Dent, sep = "_")
}
#Calculate Mean Shapes from above indexs
{
  #Denticle Mean Shapes
  mean_denticle_shapes = NULL
  Find.Mean <- function(x){
    mean_x <- mshape(S.gpa$coords[,,Denticle_Index[[x]]])
    mean_denticle_shapes = rbind(mean_x, mean_denticle_shapes)
  }
  length.dent <- c(1:9)
  mean_denticle_shapes = lapply(length.dent, Find.Mean)
  names(mean_denticle_shapes) <- paste("Mean", Dent, sep = "_")
  
  #Species Mean Shapes
  mean_species_shapes = NULL
  Find.mean_species_shapes <- function(x){
    mean_x <- mshape(S.gpa$coords[,,Species_Index[[x]]])
    mean_species_shapes = rbind(mean_x, mean_species_shapes)
  }
  length.species <- c(1:8)
  mean_species_shapes = lapply(length.species, Find.mean_species_shapes)
  names(mean_species_shapes) <- paste("Mean", Species, sep = "_")
  
  #whitei mean shapes
  whitei_mean_shapes = NULL
  Find.whitei_mean_shapes <- function(x){
    mean_x <- mshape(S.gpa$coords[,,whitei_Denticle_Index[[x]]])
    whitei_mean_shapes = rbind(mean_x, whitei_mean_shapes)
  }
  whitei_mean_shapes = lapply(length.dent, Find.whitei_mean_shapes)
  names(whitei_mean_shapes) <- paste("whitei_mean", Dent, sep = "_")
  
  #asymmetrica mean shapes
  asymmetrica_mean_shapes = NULL
  Find.asymmetrica_mean_shapes <- function(x){
    mean_x <- mshape(S.gpa$coords[,,asymmetrica_Denticle_Index[[x]]])
    asymmetrica_mean_shapes = rbind(mean_x, asymmetrica_mean_shapes)
  }
  asymmetrica_mean_shapes = lapply(length.dent, Find.asymmetrica_mean_shapes)
  names(asymmetrica_mean_shapes) <- paste("asymmetrica_mean", Dent, sep = "_")
  
  #cf. adenticulatus mean shapes
  cf.adenticulatus_mean_shapes = NULL
  Find.cf.adenticulatus_mean_shapes <- function(x){
    mean_x <- mshape(S.gpa$coords[,,cf.adenticulatus_Denticle_Index[[x]]])
    cf.adenticulatus_mean_shapes = rbind(mean_x, cf.adenticulatus_mean_shapes)
  }
  cf.adenticulatus_mean_shapes = lapply(length.dent, Find.cf.adenticulatus_mean_shapes)
  names(cf.adenticulatus_mean_shapes) <- paste("cf.adenticulatus_mean", Dent, sep = "_")
  
  #aff. binodosus mean shapes
  aff.binodosus_mean_shapes = NULL
  Find.aff.binodosus_mean_shapes <- function(x){
    mean_x <- mshape(S.gpa$coords[,,aff.binodosus_Denticle_Index[[x]]])
    aff.binodosus_mean_shapes = rbind(mean_x, aff.binodosus_mean_shapes)
  }
  aff.binodosus_mean_shapes = lapply(length.dent, Find.aff.binodosus_mean_shapes)
  names(aff.binodosus_mean_shapes) <- paste("aff.binodosus_mean", Dent, sep = "_")
  
  #anceps mean shapes
  anceps_mean_shapes = NULL
  Find.anceps_mean_shapes <- function(x){
    mean_x <- mshape(S.gpa$coords[,,anceps_Denticle_Index[[x]]])
    anceps_mean_shapes = rbind(mean_x, anceps_mean_shapes)
  }
  anceps_mean_shapes = lapply(length.dent, Find.anceps_mean_shapes)
  names(anceps_mean_shapes) <- paste("anceps_mean", Dent, sep = "_")
  
  #binodosus mean shapes
  binodosus_mean_shapes = NULL
  Find.binodosus_mean_shapes <- function(x){
    mean_x <- mshape(S.gpa$coords[,,binodosus_Denticle_Index[[x]]])
    binodosus_mean_shapes = rbind(mean_x, binodosus_mean_shapes)
  }
  binodosus_mean_shapes = lapply(length.dent, Find.binodosus_mean_shapes)
  names(binodosus_mean_shapes) <- paste("binodosus_mean", Dent, sep = "_")
  
  #clarki mean shapes
  clarki_mean_shapes = NULL
  Find.clarki_mean_shapes <- function(x){
    mean_x <- mshape(S.gpa$coords[,,clarki_Denticle_Index[[x]]])
    clarki_mean_shapes = rbind(mean_x, clarki_mean_shapes)
  }
  clarki_mean_shapes = lapply(length.dent, Find.clarki_mean_shapes)
  names(clarki_mean_shapes) <- paste("clarki_mean", Dent, sep = "_")
  
  #obliquidentatus mean shapes
  obliquidentatus_mean_shapes = NULL
  Find.obliquidentatus_mean_shapes <- function(x){
    mean_x <- mshape(S.gpa$coords[,,obliquidentatus_Denticle_Index[[x]]])
    obliquidentatus_mean_shapes = rbind(mean_x, obliquidentatus_mean_shapes)
  }
  obliquidentatus_mean_shapes = lapply(length.dent, Find.obliquidentatus_mean_shapes)
  names(obliquidentatus_mean_shapes) <- paste("obliquidentatus_mean", Dent, sep = "_")
}
#Thin-plate spline deformations
{
  #Denticle 1
  plotRefToTarget(mean_denticle_shapes[[1]], whitei_mean_shapes[[1]], mag = 1)
  plotRefToTarget(mean_denticle_shapes[[1]], asymmetrica_mean_shapes[[1]], mag = 1)
  plotRefToTarget(mean_denticle_shapes[[1]], cf.adenticulatus_mean_shapes[[1]], mag = 1)
  plotRefToTarget(mean_denticle_shapes[[1]], aff.binodosus_mean_shapes[[1]],  mag = 1)
  plotRefToTarget(mean_denticle_shapes[[1]], anceps_mean_shapes[[1]], mag = 1)
  plotRefToTarget(mean_denticle_shapes[[1]], binodosus_mean_shapes[[1]], mag = 1)
  plotRefToTarget(mean_denticle_shapes[[1]], clarki_mean_shapes[[1]], mag = 1)
  plotRefToTarget(mean_denticle_shapes[[1]], obliquidentatus_mean_shapes[[1]], mag = 1)
  #Denticle 2
  plotRefToTarget(mean_denticle_shapes[[2]], whitei_mean_shapes[[2]], mag = 1)
  #plotRefToTarget(mean_denticle_shapes[[2]], asymmetrica_mean_shapes[[2]], mag = 1)
  plotRefToTarget(mean_denticle_shapes[[2]], cf.adenticulatus_mean_shapes[[2]], mag = 1)
  plotRefToTarget(mean_denticle_shapes[[2]], aff.binodosus_mean_shapes[[2]], mag = 1)
  plotRefToTarget(mean_denticle_shapes[[2]], anceps_mean_shapes[[2]], mag = 1)
  plotRefToTarget(mean_denticle_shapes[[2]], binodosus_mean_shapes[[2]], mag = 1)
  plotRefToTarget(mean_denticle_shapes[[2]], clarki_mean_shapes[[2]], mag = 1)
  plotRefToTarget(mean_denticle_shapes[[2]], obliquidentatus_mean_shapes[[2]], mag = 1)
  #Denticle 3
  plotRefToTarget(mean_denticle_shapes[[3]], whitei_mean_shapes[[3]], mag = 1)
  plotRefToTarget(mean_denticle_shapes[[3]], asymmetrica_mean_shapes[[3]], mag = 1)
  plotRefToTarget(mean_denticle_shapes[[3]], cf.adenticulatus_mean_shapes[[3]], mag = 1)
  plotRefToTarget(mean_denticle_shapes[[3]], aff.binodosus_mean_shapes[[3]], mag = 1)
  plotRefToTarget(mean_denticle_shapes[[3]], anceps_mean_shapes[[3]], mag = 1)
  plotRefToTarget(mean_denticle_shapes[[3]], binodosus_mean_shapes[[3]], mag = 1)
  plotRefToTarget(mean_denticle_shapes[[3]], clarki_mean_shapes[[3]], mag = 1)
  plotRefToTarget(mean_denticle_shapes[[3]], obliquidentatus_mean_shapes[[3]], mag = 1)
  #Denticle 4
  plotRefToTarget(mean_denticle_shapes[[4]], whitei_mean_shapes[[4]], mag = 1)
  plotRefToTarget(mean_denticle_shapes[[4]], asymmetrica_mean_shapes[[4]], mag = 1)
  plotRefToTarget(mean_denticle_shapes[[4]], cf.adenticulatus_mean_shapes[[4]], mag = 1)
  plotRefToTarget(mean_denticle_shapes[[4]], aff.binodosus_mean_shapes[[4]], mag = 1)
  plotRefToTarget(mean_denticle_shapes[[4]], anceps_mean_shapes[[4]], mag = 1)
  plotRefToTarget(mean_denticle_shapes[[4]], binodosus_mean_shapes[[4]], mag = 1)
  plotRefToTarget(mean_denticle_shapes[[4]], clarki_mean_shapes[[4]], mag = 1)
  plotRefToTarget(mean_denticle_shapes[[4]], obliquidentatus_mean_shapes[[4]], mag = 1)
  #Denticle 5
  plotRefToTarget(mean_denticle_shapes[[5]], whitei_mean_shapes[[5]],  mag = 1)
  plotRefToTarget(mean_denticle_shapes[[5]], asymmetrica_mean_shapes[[5]], mag = 1)
  plotRefToTarget(mean_denticle_shapes[[5]], cf.adenticulatus_mean_shapes[[5]], mag = 1)
  plotRefToTarget(mean_denticle_shapes[[5]], aff.binodosus_mean_shapes[[5]], mag = 1)
  plotRefToTarget(mean_denticle_shapes[[5]], anceps_mean_shapes[[5]],  mag = 1)
  plotRefToTarget(mean_denticle_shapes[[5]], binodosus_mean_shapes[[5]], mag = 1)
  plotRefToTarget(mean_denticle_shapes[[5]], clarki_mean_shapes[[5]], mag = 1)
  plotRefToTarget(mean_denticle_shapes[[5]], obliquidentatus_mean_shapes[[5]], mag = 1)
  #Denticle 6
  plotRefToTarget(mean_denticle_shapes[[6]], whitei_mean_shapes[[6]], mag = 1)
  plotRefToTarget(mean_denticle_shapes[[6]], asymmetrica_mean_shapes[[6]], mag = 1)
  plotRefToTarget(mean_denticle_shapes[[6]], cf.adenticulatus_mean_shapes[[6]],mag = 1)
  plotRefToTarget(mean_denticle_shapes[[6]], aff.binodosus_mean_shapes[[6]], mag = 1)
  plotRefToTarget(mean_denticle_shapes[[6]], anceps_mean_shapes[[6]], mag = 1)
  plotRefToTarget(mean_denticle_shapes[[6]], binodosus_mean_shapes[[6]], mag = 1)
  plotRefToTarget(mean_denticle_shapes[[6]], clarki_mean_shapes[[6]], mag = 1)
  plotRefToTarget(mean_denticle_shapes[[6]], obliquidentatus_mean_shapes[[6]], mag = 1)
  #Denticle 7
  plotRefToTarget(mean_denticle_shapes[[7]], whitei_mean_shapes[[7]], mag = 1)
  plotRefToTarget(mean_denticle_shapes[[7]], asymmetrica_mean_shapes[[7]], mag = 1)
  plotRefToTarget(mean_denticle_shapes[[7]], cf.adenticulatus_mean_shapes[[7]], mag = 1)
  plotRefToTarget(mean_denticle_shapes[[7]], aff.binodosus_mean_shapes[[7]], mag = 1)
  plotRefToTarget(mean_denticle_shapes[[7]], anceps_mean_shapes[[7]], mag = 1)
  plotRefToTarget(mean_denticle_shapes[[7]], binodosus_mean_shapes[[7]],  mag = 1)
  plotRefToTarget(mean_denticle_shapes[[7]], clarki_mean_shapes[[7]], mag = 1)
  plotRefToTarget(mean_denticle_shapes[[7]], obliquidentatus_mean_shapes[[7]], mag = 1)
  #Denticle 8
  plotRefToTarget(mean_denticle_shapes[[8]], whitei_mean_shapes[[8]], mag = 1)
  plotRefToTarget(mean_denticle_shapes[[8]], asymmetrica_mean_shapes[[8]], mag = 1)
  plotRefToTarget(mean_denticle_shapes[[8]], cf.adenticulatus_mean_shapes[[8]], mag = 1)
  plotRefToTarget(mean_denticle_shapes[[8]], aff.binodosus_mean_shapes[[8]], mag = 1)
  plotRefToTarget(mean_denticle_shapes[[8]], anceps_mean_shapes[[8]], mag = 1)
  plotRefToTarget(mean_denticle_shapes[[8]], binodosus_mean_shapes[[8]], mag = 1)
  plotRefToTarget(mean_denticle_shapes[[8]], clarki_mean_shapes[[8]], mag = 1)
  plotRefToTarget(mean_denticle_shapes[[8]], obliquidentatus_mean_shapes[[8]], mag = 1)
  #Denticle 9
  plotRefToTarget(mean_denticle_shapes[[9]], cf.adenticulatus_mean_shapes[[9]], mag = 1)
  plotRefToTarget(mean_denticle_shapes[[9]], aff.binodosus_mean_shapes[[9]], mag = 1)
  plotRefToTarget(mean_denticle_shapes[[9]], binodosus_mean_shapes[[9]], mag = 1)
  plotRefToTarget(mean_denticle_shapes[[9]], clarki_mean_shapes[[9]], mag = 1)
  plotRefToTarget(mean_denticle_shapes[[9]], obliquidentatus_mean_shapes[[9]], mmag = 1)
}

#PCA Scatter plot with 95% confidence intervals
{confidenceplot <- ggplot(PCA.values, aes(S.pca.x...1.,S.pca.x...2., color = Species.ID), showlegend = T) +
    geom_point() +
    stat_ellipse(type = "norm", linetype = 1) +
    coord_fixed() +
    labs(color = "Species", x = "PC 1", y = "PC 2") +
    theme_minimal() + 
    labs(x = "PC 1: 79.11%", y = "PC 2: 9.08%")+
    scale_x_continuous(breaks = seq(-0.8,0.8,by = 0.2))
  
  orca(confidenceplot, "confidenceplot.svg")}
#Plotreftoshape for PCA axes
{
ds1 <- plotRefToTarget(mshape(S.pca$A[,,]), S.pca$shapes$shapes.PC1$max)
ds2 <- plotRefToTarget(mshape(S.pca$A[,,]), S.pca$shapes$shapes.PC1$min)
ds3 <- plotRefToTarget(mshape(S.pca$A[,,]), S.pca$shapes$shapes.PC2$max)
ds4 <- plotRefToTarget(mshape(S.pca$A[,,]), S.pca$shapes$shapes.PC2$min)
}

#Clustering of species to establish which are most similar
{
mean.pc1 <- function(x,y){
  mean(S.pca$x[,1][Denticle.ID == x & Species.ID == y])
}
mean.pc2 <- function(x,y){
  mean(S.pca$x[,2][Denticle.ID == x & Species.ID == y])
}
mean.pc3 <- function(x,y){
  mean(S.pca$x[,3][Denticle.ID == x & Species.ID == y])
}
pc1.mean <- function(y){
  sapply(Dent[-9], mean.pc1, y = y)
}
pc2.mean <- function(y){
  sapply(Dent[-9], mean.pc2, y = y)
}
pc3.mean <- function(y){
  sapply(Dent[-9], mean.pc3, y = y)
}
mSPC1 = sapply(Species, pc1.mean)
mSPC2 = sapply(Species, pc2.mean)
mSPC3 = sapply(Species, pc3.mean)
y.1 = pvclust(mSPC1*(79.116) + mSPC2*(9.078) + mSPC3*(6.190),method.hclust= "ward.D2", method.dist = "euclidean", iseed = 1, nboot = 9999)
plot(y.1)
}

###Qunatify phenotypic trajectories and comparison between trait spaces
{
#Convert radians to degrees
rad2degree <- function(rad){
  temp = (rad*180)/pi
  return(temp)
}

TrajLabel <- c("S","S","S","S","S",
               "S","S","S", "S", "E", "E",
               "E", "E", "S","S","S", "S","S","S",
               "E",
               "S",
               "E",  "E",  "E",  "E",
               "S", "S",
               "S", "E", "E", "E", "E", "E", "S",  "S", "E",
               "E", "S",  "E", "S",  "E", "S",
               "S","S","E",
               "S","S","S","S","S","S","S",
               "E", "S", "E", "S", "E", "E", "E")

#Assign by lineage
TrajL <- Specimen.ID
for(i in 1:length(SID)){
  TrajL = gsub(SID[i], TrajLabel[i], TrajL)
}

#"Lineages" for trajectory analysis. Split Sakmarian into Early and Late
TrajLineages <- c("Late Asselian","Late Asselian","Late Asselian","Late Asselian","Late Asselian",
                  "Late Asselian","Late Asselian","Late Asselian", "Late Sakmarian", "Late Sakmarian", "Late Sakmarian",
                  "Late Sakmarian", "Late Sakmarian", "Late Asselian","Late Asselian","Late Asselian", "Late Asselian","Late Asselian","Late Asselian",
                  "Late Asselian",
                  "Late Sakmarian",
                  "Late Asselian",  "Late Asselian",  "Late Asselian",  "Late Asselian",
                  "Late Sakmarian", "Late Sakmarian",
                  "Early Asselian", "Early Asselian", "Early Asselian", "Early Asselian", "Early Asselian", "Early Asselian", "Early Asselian",  "Early Asselian", "Early Asselian",
                  "Late Asselian", "Early Sakmarian",  "Late Asselian", "Early Sakmarian",  "Late Asselian", "Early Sakmarian",
                  "Early Sakmarian","Early Sakmarian","Late Asselian",
                  "Early Sakmarian","Early Sakmarian","Early Sakmarian","Early Sakmarian","Early Sakmarian","Early Sakmarian","Early Sakmarian",
                  "Early Sakmarian", "Early Sakmarian", "Early Sakmarian", "Early Sakmarian", "Early Sakmarian", "Early Sakmarian", "Early Sakmarian")

#Assign by lineage
TALineages <- Specimen.ID
for(i in 1:length(SID)){
  TALineages = gsub(SID[i], TrajLineages[i], TALineages)
}

Trajdata <- geomorph.data.frame(coords = two.d.array(S.gpa$coords), Lin = factor(TALineages), EndMembers = factor(TrajL), Species = factor(Species.ID), TraitGroups = as.factor(TraitGroups))
fit <- lm.rrpp(coords ~ Lin * EndMembers, data = Trajdata ,iter = 9999)
reveal.model.designs(fit)
TA <- trajectory.analysis(fit, groups = Trajdata $Lin, traj.pts = Trajdata $EndMembers)
summary(TA, attribute = "MD")
summary(TA, attribute = "TC", angle.type = "deg")

TP <- plot(TA, pch = as.numeric(Trajdata $Lin) + 20, bg = as.numeric(Trajdata $EndMembers) ,cex = 0.7, col = "grey")
add.trajectories(TP, traj.pch = c(21,22,23,24), start.bg = 1, end.bg = 2)
legend("bottomleft", levels(Trajdata $Lin), pch =  c(21,22,23,24), pt.bg = 0.5)

#Figure 5 in main text
TPfig <- plot_ly() %>% layout(yaxis = list(scaleanchor = "x"))
TPfig <- TPfig %>% add_trace(x = ~PC1, y = ~PC2, mode = "markers", type = "scatter", color = Trajdata$EndMembers, symbol = ~Trajdata$Lin, symbols = c('star-triangle-down','square','circle','diamond'), marker = list(size = 3))
TPfig <- TPfig %>% add_trace(x = ~TP[["trajectories"]][["Early Asselian"]][,1], y = ~TP[["trajectories"]][["Early Asselian"]][,2], type = "scatter", mode = "line", name = "Early Asselian", opacity = 1.0,colors = "teal")
TPfig <- TPfig %>% add_trace(x = ~TP[["trajectories"]][["Early Sakmarian"]][,1], y = ~TP[["trajectories"]][["Early Sakmarian"]][,2], type = "scatter", mode = "line", name = "Early Sakmarian", opacity = 1.0,colors = "purple")
TPfig <- TPfig %>% add_trace(x = ~TP[["trajectories"]][["Late Sakmarian"]][,1], y = ~TP[["trajectories"]][["Late Sakmarian"]][,2], type = "scatter", mode = "line", name = "Late Sakmarian", opacity = 1.0,colors = "blue")
TPfig <- TPfig %>% add_trace(x = ~TP[["trajectories"]][["Late Asselian"]][,1], y = ~TP[["trajectories"]][["Late Asselian"]][,2], type = "scatter", mode = "line", name = "Late Asselian", opacity = 1.0,colors = "orange")
TPfig

#Trait space pairwise comparison 
PairMeanfit <- lm.rrpp(coords ~ TraitGroups, data = Trajdata ,iter = 9999)
Meanpair <- pairwise(PairMeanfit, groups = Trajdata$TraitGroups)
Mano_pair <- manova.update(PairMeanfit)
summary(Mano_pair, test = "Pillai")
}

#Denticles 1,2,8,9 will not get results as replications do not occur at each level (insufficient sample size)
#Trajectory Analysis for denticle 1 and comparison between trait spaces
#no results
{
  filter.list <- function(x){
    for(i in x)
      temp = S.gpa$coords[,,i]
    return(temp)
  }
  
  preFilter <-  data.frame(Lin = as.factor(TALineages), EndMembers = as.factor(TrajL), Species = as.factor(Species.ID), Homology = as.factor(Denticle.ID), TraitGroups = as.factor(TraitGroups))
  D1filt.list <- lapply(Denticle_Index$Index_D1, filter.list)
  x1 <- array(unlist(D1filt.list), dim = c(52,2,length(Denticle_Index$Index_D1)))
  dent1attr <- preFilter %>% filter(Homology == "D1")  
  d1traj <- geomorph.data.frame(coords = two.d.array(x1), dent1attr)
  fitD1 <- lm.rrpp(coords ~ Lin * EndMembers, data = d1traj ,iter = 199)
  reveal.model.designs(fitD1)
  TAD1 <- trajectory.analysis(fitD1, groups = d1traj$Lin, traj.pts = d1traj$EndMembers)
  TPD1 <- plot(TAD1, pch = as.numeric(d1traj$Lin) + 20, bg = as.numeric(d1traj$EndMembers) ,cex = 0.7, col = "grey")
  add.trajectories(TPD1, traj.pch = c(21,22,23,24), start.bg = 1, end.bg = 2)
  summary(TAD1, attribute = "MD")
  summary(TAD1, attribute = "TC")
  
  #Trait space pairwise comparison 
  PairD1fit <- lm.rrpp(coords ~ TraitGroups, data = d1traj ,iter = 9999)
  D1pair <- pairwise(PairD1fit, groups = d1traj$TraitGroups)
  Mano_pair1 <- manova.update(PairD1fit)
  summary(Mano_pair1, test = "Pillai")
}
#Trajectory Analysis for denticle 3
#parallel
{
  D3filt.list <- lapply(Denticle_Index$Index_D3, filter.list)
  x3 <- array(unlist(D3filt.list), dim = c(52,2,length(Denticle_Index$Index_D3)))
  dent3attr <- preFilter %>% filter(Homology == "D3")  
  d3traj <- geomorph.data.frame(coords = two.d.array(x3), dent3attr)
  fitD3 <- lm.rrpp(coords ~ Lin * EndMembers, data = d3traj ,iter = 199)
  reveal.model.designs(fitD3)
  TAD3 <- trajectory.analysis(fitD3, groups = d3traj$Lin, traj.pts = d3traj$EndMembers)
  TPD3 <- plot(TAD3, pch = as.numeric(d3traj$Lin) + 20, bg = as.numeric(d3traj$EndMembers) ,cex = 0.7, col = "grey")
  add.trajectories(TPD3, traj.pch = c(21,22,23,24), start.bg = 1, end.bg = 2)
  summary(TAD3, attribute = "MD")
  summary(TAD3, attribute = "TC")
  
  TPfigD3 <- plot_ly() %>% layout(yaxis = list(scaleanchor = "x", title = "PC 2", range = c(-0.3,0.3)), xaxis = list(title = "PC 1", range = c(-0.6,0.4)), showlegend = F, title = "D3")
  TPfigD3 <- TPfigD3 %>% add_trace(x = ~TPD3$pc.points[,1], y = ~TPD3$pc.points[,2], mode = "markers", type = "scatter", color = ~d3traj$EndMembers, symbol = ~d3traj$Lin, symbols = c('star-triangle-down','square','circle','diamond'), marker = list(size = 5.5))
  TPfigD3 <- TPfigD3 %>% add_trace(x = ~TPD3[["trajectories"]][["Early Asselian"]][,1], y = ~TPD3[["trajectories"]][["Early Asselian"]][,2], type = "scatter", mode = "line", name = "Early Asselian", opacity = 1.0,colors = "teal")
  TPfigD3 <- TPfigD3 %>% add_trace(x = ~TPD3[["trajectories"]][["Early Sakmarian"]][,1], y = ~TPD3[["trajectories"]][["Early Sakmarian"]][,2], type = "scatter", mode = "line", name = "Early Sakmarian", opacity = 1.0,colors = "purple")
  TPfigD3 <- TPfigD3 %>% add_trace(x = ~TPD3[["trajectories"]][["Late Sakmarian"]][,1], y = ~TPD3[["trajectories"]][["Late Sakmarian"]][,2], type = "scatter", mode = "line", name = "Late Sakmarian", opacity = 1.0,colors = "blue")
  TPfigD3 <- TPfigD3 %>% add_trace(x = ~TPD3[["trajectories"]][["Late Asselian"]][,1], y = ~TPD3[["trajectories"]][["Late Asselian"]][,2], type = "scatter", mode = "line", name = "Late Asselian", opacity = 1.0,colors = "orange")
  TPfigD3
  orca(TPfigD3, "D3traj.svg")
  
  #Trait space pairwise comparison 
  PairD3fit <- lm.rrpp(coords ~ TraitGroups, data = d3traj ,iter = 9999)
  D3pair <- pairwise(PairD3fit, groups = d3traj$TraitGroups)
  summary(D3pair)
  Mano_pair3 <- manova.update(PairD3fit)
  summary(Mano_pair3, test = "Pillai")
}
#Trajectory Analysis for denticle 4
#paralel
{
  D4filt.list <- lapply(Denticle_Index$Index_D4, filter.list)
  x4 <- array(unlist(D4filt.list), dim = c(52,2,length(Denticle_Index$Index_D4)))
  dent4attr <- preFilter %>% filter(Homology == "D4")  
  d4traj <- geomorph.data.frame(coords = two.d.array(x4), dent4attr)
  fitD4 <- lm.rrpp(coords ~ Lin * EndMembers, data = d4traj ,iter = 199)
  reveal.model.designs(fitD4)
  TAD4 <- trajectory.analysis(fitD4, groups = d4traj$Lin, traj.pts = d4traj$EndMembers)
  TPD4 <- plot(TAD4, pch = as.numeric(d4traj$Lin) + 20, bg = as.numeric(d4traj$EndMembers) ,cex = 0.7, col = "grey")
  add.trajectories(TPD4, traj.pch = c(21,22,23,24), start.bg = 1, end.bg = 2)
  summary(TAD4, attribute = "MD")
  summary(TAD4, attribute = "TC")
  
  TPfigD4 <- plot_ly() %>% layout(yaxis = list(scaleanchor = "x", title = "PC 2", range = c(-0.3,0.3)), xaxis = list(title = "PC 1", range = c(-0.6,0.4)), showlegend = F, title = "D4")
  TPfigD4 <- TPfigD4 %>% add_trace(x = ~TPD4$pc.points[,1], y = ~TPD4$pc.points[,2], mode = "markers", type = "scatter", color = ~d4traj$EndMembers, symbol = ~d4traj$Lin, symbols = c('star-triangle-down','square','circle','diamond'), marker = list(size = 5.5))
  TPfigD4 <- TPfigD4 %>% add_trace(x = ~TPD4[["trajectories"]][["Early Asselian"]][,1], y = ~TPD4[["trajectories"]][["Early Asselian"]][,2], type = "scatter", mode = "line", name = "Early Asselian", opacity = 1.0,colors = "teal")
  TPfigD4 <- TPfigD4 %>% add_trace(x = ~TPD4[["trajectories"]][["Early Sakmarian"]][,1], y = ~TPD4[["trajectories"]][["Early Sakmarian"]][,2], type = "scatter", mode = "line", name = "Early Sakmarian", opacity = 1.0,colors = "purple")
  TPfigD4 <- TPfigD4 %>% add_trace(x = ~TPD4[["trajectories"]][["Late Sakmarian"]][,1], y = ~TPD4[["trajectories"]][["Late Sakmarian"]][,2], type = "scatter", mode = "line", name = "Late Sakmarian", opacity = 1.0,colors = "blue")
  TPfigD4 <- TPfigD4 %>% add_trace(x = ~TPD4[["trajectories"]][["Late Asselian"]][,1], y = ~TPD4[["trajectories"]][["Late Asselian"]][,2], type = "scatter", mode = "line", name = "Late Asselian", opacity = 1.0,colors = "orange")
  TPfigD4
  orca(TPfigD4, "D4traj.svg")
  
  #Trait space pairwise comparison 
  PairD4fit <- lm.rrpp(coords ~ TraitGroups, data = d4traj ,iter = 9999)
  D4pair <- pairwise(PairD4fit, groups = d4traj$TraitGroups)
  summary(D4pair)
  Mano_pair4 <- manova.update(PairD4fit)
  summary(Mano_pair4, test = "Pillai")
}
#Trajectory Analysis for denticle 5
#parallel
{
  preFilter <-  data.frame(Lin = as.factor(TALineages), EndMembers = as.factor(TrajL), Species = as.factor(Species.ID), Homology = as.factor(Denticle.ID))
  dentfilt.list <- lapply(Denticle_Index$Index_D5, filter.list)
  x <- array(unlist(dentfilt.list), dim = c(52,2,57))
  dent5attr <- preFilter %>% filter(Homology == "D5")  
  d5traj <- geomorph.data.frame(coords = two.d.array(x), dent5attr)
  fitD5 <- lm.rrpp(coords ~ Lin * EndMembers, data = d5traj ,iter = 199)
  reveal.model.designs(fitD5)
  TAD5 <- trajectory.analysis(fitD5, groups = d5traj$Lin, traj.pts = d5traj$EndMembers)
  TPD5 <- plot(TAD5, pch = as.numeric(d5traj$Lin) + 20, bg = as.numeric(d5traj$EndMembers) ,cex = 0.7, col = "grey")
  add.trajectories(TPD5, traj.pch = c(21,22,23,24), start.bg = 1, end.bg = 2)
  summary(TAD5, attribute = "MD")
  summary(TAD5, attribute = "TC")
  
  TPfigD5 <- plot_ly() %>% layout(yaxis = list(scaleanchor = "x", title = "PC 2", range = c(-0.3,0.3)), xaxis = list(title = "PC 1", range = c(-0.6,0.4)), showlegend = F, title = "D5")
  TPfigD5 <- TPfigD5 %>% add_trace(x = ~TPD5$pc.points[,1], y = ~TPD5$pc.points[,2], mode = "markers", type = "scatter", color = ~d5traj$EndMembers, symbol = ~d5traj$Lin, symbols = c('star-triangle-down','square','circle','diamond'), marker = list(size = 5.5))
  TPfigD5 <- TPfigD5 %>% add_trace(x = ~TPD5[["trajectories"]][["Early Asselian"]][,1], y = ~TPD4[["trajectories"]][["Early Asselian"]][,2], type = "scatter", mode = "line", name = "Early Asselian", opacity = 1.0,colors = "teal")
  TPfigD5 <- TPfigD5 %>% add_trace(x = ~TPD5[["trajectories"]][["Early Sakmarian"]][,1], y = ~TPD4[["trajectories"]][["Early Sakmarian"]][,2], type = "scatter", mode = "line", name = "Early Sakmarian", opacity = 1.0,colors = "purple")
  TPfigD5 <- TPfigD5 %>% add_trace(x = ~TPD5[["trajectories"]][["Late Sakmarian"]][,1], y = ~TPD4[["trajectories"]][["Late Sakmarian"]][,2], type = "scatter", mode = "line", name = "Late Sakmarian", opacity = 1.0,colors = "blue")
  TPfigD5 <- TPfigD5 %>% add_trace(x = ~TPD5[["trajectories"]][["Late Asselian"]][,1], y = ~TPD4[["trajectories"]][["Late Asselian"]][,2], type = "scatter", mode = "line", name = "Late Asselian", opacity = 1.0,colors = "orange")
  TPfigD5
  orca(TPfigD5, "D5traj.svg")
  
  #Trait space pairwise comparison 
  PairD5fit <- lm.rrpp(coords ~ TraitGroups, data = d5traj ,iter = 9999)
  D5pair <- pairwise(PairD5fit, groups = d5traj$TraitGroups)
  summary(D5pair)
  Mano_pair5 <- manova.update(PairD5fit)
  summary(Mano_pair5, test = "Pillai")
}
#Trajectory Analysis for denticle 6
#paralel
{
  D6filt.list <- lapply(Denticle_Index$Index_D6, filter.list)
  x6 <- array(unlist(D6filt.list), dim = c(52,2,length(Denticle_Index$Index_D6)))
  dent6attr <- preFilter %>% filter(Homology == "D6")  
  d6traj <- geomorph.data.frame(coords = two.d.array(x6), dent6attr)
  fitD6 <- lm.rrpp(coords ~ Lin * EndMembers, data = d6traj ,iter = 199)
  reveal.model.designs(fitD6)
  TAD6 <- trajectory.analysis(fitD6, groups = d6traj$Lin, traj.pts = d6traj$EndMembers)
  TPD6 <- plot(TAD6, pch = as.numeric(d6traj$Lin) + 20, bg = as.numeric(d6traj$EndMembers) ,cex = 0.7, col = "grey")
  add.trajectories(TPD6, traj.pch = c(21,22,23,24), start.bg = 1, end.bg = 2)
  summary(TAD6, attribute = "MD")
  summary(TAD6, attribute = "TC")
  
  TPfigD6 <- plot_ly() %>% layout(yaxis = list(scaleanchor = "x", title = "PC 2", range = c(-0.3,0.3)), xaxis = list(title = "PC 1", range = c(-0.6,0.4)), showlegend = F, title = "D6")
  TPfigD6 <- TPfigD6 %>% add_trace(x = ~TPD6$pc.points[,1], y = ~TPD6$pc.points[,2], mode = "markers", type = "scatter", color = ~d6traj$EndMembers, symbol = ~d6traj$Lin, symbols = c('star-triangle-down','square','circle','diamond'), marker = list(size = 5.5))
  TPfigD6 <- TPfigD6 %>% add_trace(x = ~TPD6[["trajectories"]][["Early Asselian"]][,1], y = ~TPD6[["trajectories"]][["Early Asselian"]][,2], type = "scatter", mode = "line", name = "Early Asselian", opacity = 1.0,colors = "teal")
  TPfigD6 <- TPfigD6 %>% add_trace(x = ~TPD6[["trajectories"]][["Early Sakmarian"]][,1], y = ~TPD6[["trajectories"]][["Early Sakmarian"]][,2], type = "scatter", mode = "line", name = "Early Sakmarian", opacity = 1.0,colors = "purple")
  TPfigD6 <- TPfigD6 %>% add_trace(x = ~TPD6[["trajectories"]][["Late Sakmarian"]][,1], y = ~TPD6[["trajectories"]][["Late Sakmarian"]][,2], type = "scatter", mode = "line", name = "Late Sakmarian", opacity = 1.0,colors = "blue")
  TPfigD6 <- TPfigD6 %>% add_trace(x = ~TPD6[["trajectories"]][["Late Asselian"]][,1], y = ~TPD6[["trajectories"]][["Late Asselian"]][,2], type = "scatter", mode = "line", name = "Late Asselian", opacity = 1.0,colors = "orange")
  TPfigD6
  orca(TPfigD6, "D6traj.svg")
  
  #Trait space pairwise comparison 
  PairD6fit <- lm.rrpp(coords ~ TraitGroups, data = d6traj ,iter = 9999)
  D6pair <- pairwise(PairD6fit, groups = d6traj$TraitGroups)
  summary(D6pair)
  Mano_pair6 <- manova.update(PairD6fit)
  summary(Mano_pair6, test = "Pillai")
}
#Trajectory Analysis for denticle 7
#paralel
{
  D7filt.list <- lapply(Denticle_Index$Index_D7, filter.list)
  x7 <- array(unlist(D7filt.list), dim = c(52,2,length(Denticle_Index$Index_D7)))
  dent7attr <- preFilter %>% filter(Homology == "D7")  
  d7traj <- geomorph.data.frame(coords = two.d.array(x7), dent7attr)
  fitD7 <- lm.rrpp(coords ~ Lin * EndMembers, data = d7traj ,iter = 199)
  reveal.model.designs(fitD7)
  TAD7 <- trajectory.analysis(fitD7, groups = d7traj$Lin, traj.pts = d7traj$EndMembers)
  TPD7 <- plot(TAD7, pch = as.numeric(d7traj$Lin) + 20, bg = as.numeric(d7traj$EndMembers) ,cex = 0.7, col = "grey")
  add.trajectories(TPD7, traj.pch = c(21,22,23,24), start.bg = 1, end.bg = 2)
  summary(TAD7, attribute = "MD")
  summary(TAD7, attribute = "TC")
  
  TPfigD7 <- plot_ly() %>% layout(yaxis = list(scaleanchor = "x", title = "PC 2", range = c(-0.3,0.3)), xaxis = list(title = "PC 1", range = c(-0.6,0.4)), showlegend = F, title = "D7")
  TPfigD7 <- TPfigD7 %>% add_trace(x = ~TPD7$pc.points[,1], y = ~TPD7$pc.points[,2], mode = "markers", type = "scatter", color = ~d7traj$EndMembers, symbol = ~d7traj$Lin, symbols = c('star-triangle-down','square','circle','diamond'), marker = list(size = 5.5))
  TPfigD7 <- TPfigD7 %>% add_trace(x = ~TPD7[["trajectories"]][["Early Asselian"]][,1], y = ~TPD7[["trajectories"]][["Early Asselian"]][,2], type = "scatter", mode = "line", name = "Early Asselian", opacity = 1.0,colors = "teal")
  TPfigD7 <- TPfigD7 %>% add_trace(x = ~TPD7[["trajectories"]][["Early Sakmarian"]][,1], y = ~TPD7[["trajectories"]][["Early Sakmarian"]][,2], type = "scatter", mode = "line", name = "Early Sakmarian", opacity = 1.0,colors = "purple")
  TPfigD7 <- TPfigD7 %>% add_trace(x = ~TPD7[["trajectories"]][["Late Sakmarian"]][,1], y = ~TPD7[["trajectories"]][["Late Sakmarian"]][,2], type = "scatter", mode = "line", name = "Late Sakmarian", opacity = 1.0,colors = "blue")
  TPfigD7 <- TPfigD7 %>% add_trace(x = ~TPD7[["trajectories"]][["Late Asselian"]][,1], y = ~TPD7[["trajectories"]][["Late Asselian"]][,2], type = "scatter", mode = "line", name = "Late Asselian", opacity = 1.0,colors = "orange")
  TPfigD7
  orca(TPfigD7, "D7traj.svg")
  
  #Trait space pairwise comparison 
  PairD7fit <- lm.rrpp(coords ~ TraitGroups, data = d7traj ,iter = 9999)
  D7pair <- pairwise(PairD7fit, groups = d7traj$TraitGroups)
  summary(D7pair)
  Mano_pair7 <- manova.update(PairD7fit)
  summary(Mano_pair7, test = "Pillai")
}

#Looking at Asselian versus Sakmarian Evolution 
###Need start and end points for trajectories
#Don't worry too much about this analysis... I was just trying it out.
{
  SATrajLabel <- c("3","3","3","3","3",
                   "3","3","3", "3", "4", "4",
                   "4", "4", "3","3","3", "3","3","3",
                   "4",
                   "3",
                   "4",  "4",  "4",  "4",
                   "3", "3",
                   "1", "2", "2", "2", "2", "2", "1",  "1", "2",
                   "4", "1","4", "1",  "4", "1",
                   "1","1","4",
                   "1","1","1","1","1","1","1",
                   "2", "1", "2", "1", "2", "2", "2")
  
  #Assign by lineage
  
  SATrajL <- Specimen.ID
  for(i in 1:length(SID)){
    SATrajL = gsub(SID[i], SATrajLabel[i], SATrajL)
  }
  
  #"Lineages" for trajectory analysis. Split Sakmarian into Early and Late
  SATrajLineages <- c("Asselian","Asselian","Asselian","Asselian","Asselian",
                      "Asselian","Asselian","Asselian", "Sakmarian", "Sakmarian", "Sakmarian",
                      "Sakmarian", "Sakmarian", "Asselian","Asselian","Asselian", "Asselian","Asselian","Asselian",
                      "Asselian",
                      "Sakmarian",
                      "Asselian",  "Asselian",  "Asselian",  "Asselian",
                      "Sakmarian", "Sakmarian",
                      "Asselian", "Asselian", "Asselian", "Asselian", "Asselian", "Asselian", "Asselian",  "Asselian", "Asselian",
                      "Asselian", "Sakmarian",  "Asselian", "Sakmarian",  "Asselian", "Sakmarian",
                      "Sakmarian","Sakmarian","Asselian",
                      "Sakmarian","Sakmarian","Sakmarian","Sakmarian","Sakmarian","Sakmarian","Sakmarian",
                      "Sakmarian", "Sakmarian", "Sakmarian", "Sakmarian", "Sakmarian", "Sakmarian", "Sakmarian")
  
  #Assign by lineage
  SATALineages <- Specimen.ID
  for(i in 1:length(SID)){
    SATALineages = gsub(SID[i], SATrajLineages[i], SATALineages)
  }
  
  ###Gonna try trajectory.analysis in RRPP
  SAdata <- geomorph.data.frame(coords = two.d.array(S.gpa$coords), Lin = as.factor(SATALineages), EndMembers = as.factor(SATrajL))
  fitSA <- lm.rrpp(coords ~ Lin * EndMembers, data = SAdata ,iter = 10000)
  reveal.model.designs(fitSA)
  SATA <- trajectory.analysis(fitSA, groups = SAdata$Lin, traj.pts = SAdata$EndMembers)
  summary(SATA, attribute = "MD")
  summary(SATA, attribute = "TC")
  summary(SATA, attribute = "SD")
  
  SATP <- plot(SATA, pch = as.numeric(SAdata$Lin) + 20, bg = as.numeric(SAdata$EndMembers) ,cex = 0.7, col = "grey")
  add.trajectories(SATP, traj.pch = c(21,22), start.bg = 1, end.bg = 4)
  legend("bottomleft", levels(data$Lin), pch =  c(21,22), pt.bg = 0.5)
}

#Density plots for all species
{
#PC 1
d1.df <- PCA.values %>% filter(Denticle.ID == "D1")
d1.line <- data.frame(Species.ID = levels(d1.df$Species.ID),
                      Means = tapply(d1.df$S.pca.x...1., INDEX = d1.df$Species.ID, FUN = mean))
d2.df  <- PCA.values %>% filter(Denticle.ID == "D2")
d2.line <- data.frame(Species.ID = levels(d2.df$Species.ID),
                      Means = tapply(d2.df$S.pca.x...1., INDEX = d2.df$Species.ID, FUN = mean))
d3.df  <- PCA.values %>% filter(Denticle.ID == "D3")
d3.line <- data.frame(Species.ID = levels(d3.df$Species.ID),
                      Means = tapply(d3.df$S.pca.x...1., INDEX = d3.df$Species.ID, FUN = mean))
d4.df  <- PCA.values %>% filter(Denticle.ID == "D4")
d4.line <- data.frame(Species.ID = levels(d4.df$Species.ID),
                      Means = tapply(d4.df$S.pca.x...1., INDEX = d4.df$Species.ID, FUN = mean))
d5.df  <- PCA.values %>% filter(Denticle.ID == "D5")
d5.line <- data.frame(Species.ID = levels(d5.df$Species.ID),
                      Means = tapply(d5.df$S.pca.x...1., INDEX = d5.df$Species.ID, FUN = mean))
d6.df  <- PCA.values %>% filter(Denticle.ID == "D6")
d6.line <- data.frame(Species.ID = levels(d6.df$Species.ID),
                      Means = tapply(d6.df$S.pca.x...1., INDEX = d6.df$Species.ID, FUN = mean))
d7.df  <- PCA.values %>% filter(Denticle.ID == "D7")
d7.line <- data.frame(Species.ID = levels(d7.df$Species.ID),
                      Means = tapply(d7.df$S.pca.x...1., INDEX = d7.df$Species.ID, FUN = mean))
d8.df  <- PCA.values %>% filter(Denticle.ID == "D8")
d8.line <- data.frame(Species.ID = levels(d8.df$Species.ID),
                      Means = tapply(d8.df$S.pca.x...1., INDEX = d8.df$Species.ID, FUN = mean))
d9.df  <- PCA.values %>% filter(Denticle.ID == "D9")
d9.line <- data.frame(Species.ID = levels(d9.df$Species.ID),
                      Means = tapply(d9.df$S.pca.x...1., INDEX = d9.df$Species.ID, FUN = mean))

p = ggplot(data = PCA.values, aes(x = S.pca.x...1.)) +
  stat_density(position = "identity", geom = "line", size = 0.1) +
  xlim(-.65,.65) +
  facet_wrap(~Species.ID) +
  geom_vline(data = d1.line, aes(xintercept = Means), col = 'red', size = .1) + 
  geom_vline(data = d2.line, aes(xintercept = Means), col = 'purple', size = .1) +
  geom_vline(data = d3.line, aes(xintercept = Means), col = 'steelblue4', size = .1) +
  geom_vline(data = d4.line, aes(xintercept = Means), col = 'yellow', size = .1) +
  geom_vline(data = d5.line, aes(xintercept = Means), col = 'wheat', size = .1) +
  geom_vline(data = d6.line, aes(xintercept = Means), col = 'coral', size = .1) +
  geom_vline(data = d7.line, aes(xintercept = Means), col = 'grey1', size = .1) +
  geom_vline(data = d8.line, aes(xintercept = Means), col = 'cyan', size = .1) +
  geom_vline(data = d9.line, aes(xintercept = Means), col = 'orange', size = .1) +
  labs(x = "PC 1", y = "%", color = "Dentilce ID")

orca(p, "TraitplotPCA.svg")

#PC 2
d1.df2 <- PCA.values %>% filter(Denticle.ID == "D1")
d1.line2 <- data.frame(Species.ID = levels(d1.df2$Species.ID),
                      Means = tapply(d1.df2$S.pca.x...2., INDEX = d1.df2$Species.ID, FUN = mean))
d2.df2  <- PCA.values %>% filter(Denticle.ID == "D2")
d2.line2 <- data.frame(Species.ID = levels(d2.df2$Species.ID),
                      Means = tapply(d2.df2$S.pca.x...2., INDEX = d2.df2$Species.ID, FUN = mean))
d3.df2  <- PCA.values %>% filter(Denticle.ID == "D3")
d3.line2 <- data.frame(Species.ID = levels(d3.df2$Species.ID),
                      Means = tapply(d3.df2$S.pca.x...2., INDEX = d3.df2$Species.ID, FUN = mean))
d4.df2  <- PCA.values %>% filter(Denticle.ID == "D4")
d4.line2 <- data.frame(Species.ID = levels(d4.df2$Species.ID),
                      Means = tapply(d4.df2$S.pca.x...2., INDEX = d4.df2$Species.ID, FUN = mean))
d5.df2  <- PCA.values %>% filter(Denticle.ID == "D5")
d5.line2 <- data.frame(Species.ID = levels(d5.df2$Species.ID),
                      Means = tapply(d5.df2$S.pca.x...2., INDEX = d5.df2$Species.ID, FUN = mean))
d6.df2  <- PCA.values %>% filter(Denticle.ID == "D6")
d6.line2 <- data.frame(Species.ID = levels(d6.df2$Species.ID),
                      Means = tapply(d6.df2$S.pca.x...2., INDEX = d6.df2$Species.ID, FUN = mean))
d7.df2  <- PCA.values %>% filter(Denticle.ID == "D7")
d7.line2 <- data.frame(Species.ID = levels(d7.df2$Species.ID),
                      Means = tapply(d7.df2$S.pca.x...2., INDEX = d7.df2$Species.ID, FUN = mean))
d8.df2  <- PCA.values %>% filter(Denticle.ID == "D8")
d8.line2 <- data.frame(Species.ID = levels(d8.df2$Species.ID),
                      Means = tapply(d8.df2$S.pca.x...2., INDEX = d8.df2$Species.ID, FUN = mean))
d9.df2  <- PCA.values %>% filter(Denticle.ID == "D9")
d9.line2 <- data.frame(Species.ID = levels(d9.df2$Species.ID),
                      Means = tapply(d9.df2$S.pca.x...2., INDEX = d9.df2$Species.ID, FUN = mean))

p2 = ggplot(data = PCA.values, aes(x = S.pca.x...2.)) +
  stat_density(position = "identity", geom = "line", size = 0.1) +
  xlim(-.65,.65) +
  facet_wrap(~Species.ID) +
  geom_vline(data = d1.line2, aes(xintercept = Means), col = 'red', size = .1) + 
  geom_vline(data = d2.line2, aes(xintercept = Means), col = 'purple', size = .1) +
  geom_vline(data = d3.line2, aes(xintercept = Means), col = 'steelblue4', size = .1) +
  geom_vline(data = d4.line2, aes(xintercept = Means), col = 'yellow', size = .1) +
  geom_vline(data = d5.line2, aes(xintercept = Means), col = 'wheat', size = .1) +
  geom_vline(data = d6.line2, aes(xintercept = Means), col = 'coral', size = .1) +
  geom_vline(data = d7.line2, aes(xintercept = Means), col = 'grey1', size = .1) +
  geom_vline(data = d8.line2, aes(xintercept = Means), col = 'cyan', size = .1) +
  geom_vline(data = d9.line2, aes(xintercept = Means), col = 'orange', size = .1) +
  labs(x = "PC 1", y = "%", color = "Dentilce ID")

orca(p2, "TraitplotPCA2.svg")
}

#Ecological Metrics on the first two principal components 
{
OrderedSpecies = c("anceps","binodosus","aff. binodosus","whitei","cf. adenticulatus","asymmetrica","clarki","obliquidentatus")
#range
Range.spec = function(x) {
  max(PCA.values$S.pca.x...1.[Species.ID == x]) - min(PCA.values$S.pca.x...1.[Species.ID == x])
}
#variance
var.cono1 = function(x) {
  var(PCA.values$S.pca.x...1.[Species.ID == x])
}
#standard deviation
sd.cono1 = function(x) {
  sqrt(var(PCA.values$S.pca.x...1.[Species.ID == x]))
}
#mean
mean.cono1 = function(x) {
  mean(PCA.values$S.pca.x...1.[Species.ID == x])
}
#trait means to intraspecific trait widths
wfstat = function(x) {
  mean(PCA.values$S.pca.x...1.[Species.ID == x])/sqrt(var(PCA.values$S.pca.x...1.[Species.ID == x]))
}
#absolute values of trait means to intrapspecific trait widths
absmfsw1 = function(x) {
  abs(mean(PCA.values$S.pca.x...1.[Species.ID == x]))/abs(sqrt(var(PCA.values$S.pca.x...1.[Species.ID == x])))
}
#coefficent of variance
coefofvar1 = function(x) {
  sqrt(var(PCA.values$S.pca.x...1.[Species.ID == x]))/mean(PCA.values$S.pca.x...1.[Species.ID == x])
}

#This is cumbersome but change Species.ID for each speices to calculate limiting similarity. Then run the sapply function below. Do this for each speices.
LS.anceps = function(x) {
  (mean(PCA.values$S.pca.x...1.[Species.ID == "anceps"])-mean(PCA.values$S.pca.x...1.[Species.ID == x]))/(sqrt(var(PCA.values$S.pca.x...1.[Species.ID == "anceps"]))-sqrt(var(PCA.values$S.pca.x...1.[Species.ID == x])))
}
LS.binodosus = function(x) {
  (mean(PCA.values$S.pca.x...1.[Species.ID == "binodosus"])-mean(PCA.values$S.pca.x...1.[Species.ID == x]))/(sqrt(var(PCA.values$S.pca.x...1.[Species.ID == "binodosus"]))-sqrt(var(PCA.values$S.pca.x...1.[Species.ID == x])))
}
LS.aff.binodosus = function(x) {
  (mean(PCA.values$S.pca.x...1.[Species.ID == "aff. binodosus"])-mean(PCA.values$S.pca.x...1.[Species.ID == x]))/(sqrt(var(PCA.values$S.pca.x...1.[Species.ID == "aff. binodosus"]))-sqrt(var(PCA.values$S.pca.x...1.[Species.ID == x])))
}
LS.whitei = function(x) {
  (mean(PCA.values$S.pca.x...1.[Species.ID == "whitei"])-mean(PCA.values$S.pca.x...1.[Species.ID == x]))/(sqrt(var(PCA.values$S.pca.x...1.[Species.ID == "whitei"]))-sqrt(var(PCA.values$S.pca.x...1.[Species.ID == x])))
}
LS.cf.adenticulatus = function(x) {
  (mean(PCA.values$S.pca.x...1.[Species.ID == "cf. adenticulatus"])-mean(PCA.values$S.pca.x...1.[Species.ID == x]))/(sqrt(var(PCA.values$S.pca.x...1.[Species.ID == "cf. adenticulatus"]))-sqrt(var(PCA.values$S.pca.x...1.[Species.ID == x])))
}
LS.asymmetrica = function(x) {
  (mean(PCA.values$S.pca.x...1.[Species.ID == "asymmetrica"])-mean(PCA.values$S.pca.x...1.[Species.ID == x]))/(sqrt(var(PCA.values$S.pca.x...1.[Species.ID == "asymmetrica"]))-sqrt(var(PCA.values$S.pca.x...1.[Species.ID == x])))
}
LS.clarki = function(x) {
  (mean(PCA.values$S.pca.x...1.[Species.ID == "clarki"])-mean(PCA.values$S.pca.x...1.[Species.ID == x]))/(sqrt(var(PCA.values$S.pca.x...1.[Species.ID == "clarki"]))-sqrt(var(PCA.values$S.pca.x...1.[Species.ID == x])))
}
LS.obliquidentatus = function(x) {
  (mean(PCA.values$S.pca.x...1.[Species.ID == "obliquidentatus"])-mean(PCA.values$S.pca.x...1.[Species.ID == x]))/(sqrt(var(PCA.values$S.pca.x...1.[Species.ID == "obliquidentatus"]))-sqrt(var(PCA.values$S.pca.x...1.[Species.ID == x])))
}

limit = sapply(OrderedSpecies, LS.anceps) #anceps
limit1 = sapply(OrderedSpecies, LS.binodosus) #binodosus
limit2 = sapply(OrderedSpecies, LS.aff.binodosus) #aff. binodosus
limit3 = sapply(OrderedSpecies, LS.whitei) #whitei
limit4 = sapply(OrderedSpecies, LS.cf.adenticulatus) #cf. adenticulatus
limit5 = sapply(OrderedSpecies, LS.asymmetrica) #asymmetrica
limit6 = sapply(OrderedSpecies, LS.clarki) #clarki
limit7 = sapply(OrderedSpecies, LS.obliquidentatus) #obliquidentatus

limitsimilarity = data.frame(limit,limit1,limit2,limit3,limit4,limit5,limit6,limit7)

colnames(limitsimilarity) <- c("Sw. anceps","Sw. binodosus","Sw. aff. binodosus","Sw. whitei","Sw. cf. adenticulatus","Sw. asymmetrica","Sw. clarki","Sw. obliquidentatus" )
rownames(limitsimilarity) <- c("Sw. anceps","Sw. binodosus","Sw. aff. binodosus","Sw. whitei","Sw. cf. adenticulatus","Sw. asymmetrica","Sw. clarki","Sw. obliquidentatus" ) 

View(limitsimilarity) #table 2

rangeSW = sapply(OrderedSpecies, Range.spec)
rangeSW = sapply(OrderedSpecies, Range.spec)
rangpc1 = max(PCA.values$S.pca.x...1.) - min(PCA.values$S.pca.x...1.)
commrange = rangeSW/rangpc1
varSw = sapply(OrderedSpecies, var.cono1)
sdSw = sapply(OrderedSpecies, sd.cono1)
mfSw = sapply(OrderedSpecies, wfstat)
meanSW = sapply(OrderedSpecies, mean.cono1)
absmfSw = sapply(OrderedSpecies, absmfsw1)
coefvar = sapply(OrderedSpecies, coefofvar1)
communityvarPC1 = varSw/sdSw #community-wide variance
traitspace = c("1","2","2","2","3","4","5","5")
IntraSweet <-data.frame(OrderedSpecies,varSw,sdSw,meanSW,mfSw,absmfSw,coefvar,traitspace,communityvarPC1, commrange)
IntraSweet$OrderedSpecies <- factor(IntraSweet$OrderedSpecies, levels = OrderedSpecies)
View(IntraSweet) #table 1
}
##########
##########
#Sensitivity Qualitative Analysis
#Qualitative assessment of cross sections ability to differentiate morphologies
{Outline.Data.re <- readRDS(file = "Full_Data_OutlineSense.rds")
Actual.Spec.ID2 <-  c("15","20","31","33_5","33","38_2","38_3","38","39_2","39_3","39","40_2","40_3","40",
                      "15","20","31","33_5","33","38_2","38_3","38","39_2","39_3","39","40_2","40_3","40",
                      "15","20","31","33_5","33","38_2","38_3","38","39_2","39_3","39","40_2","40_3","40",
                      "15","20","33_5","33","38_2","38_3","38","39_2","39_3","39","40_2","40_3","40",
                      "15","20","33_5","33","38_2","38_3","38","39_2","39_3","39","40_2","40_3","40")
inputSens = Outline.Data.re #Just created an extra data set so I don't mess with the original
#Delete duplicate outlines
#removal of double outlines
Sens.appended <- inputSens[-c(9,10,11,23,24,38,51,64)] 
#The momocs package incorrectly imports these two cross sections.
### Baseline removal for each specimen
SensData <-  pbapply::pblapply(Sens.appended, curve.comp) #New Open Curves
#Sort points along the curve
#This has been altered from the original function 
#Inserted RANNx = matrix(RANNx, ncol = 2)
#Stop sort at n = 4, last 4 points may not be sorted, probably truncate last 4
RANNsort.curve2 <- function(x){
  RANNinputdata = x #function input curve
  RANNx = RANNinputdata[RANNinputdata[,2] == min(RANNinputdata[,2]),] #find min values of y and correspounding x coordinate
  RANNx = matrix(RANNx,ncol = 2)
  RANNx = list(RANNx) #place min values into a list
  RANNy = RANNx[[1]][RANNx[[1]][,1] == min(RANNx[[1]][,1]),] #find start of curve coordinates close to {0,0}
  #Create a new list with this value
  RANNsorted.curve = matrix(RANNy, nrow = 1, ncol = 2) #sorted curve begins with start position closest to {0,0}
  RANNdata.search = RANNinputdata #unsorted curve to search through
  repeat {
    RANNa = tail(RANNsorted.curve,1)#find starting point for sorting
    RANNvaluex = RANNa[1] #set value x to remove
    RANNvaluey = RANNa[2] #set value y to remove
    RANNvalues = c(RANNvaluex,RANNvaluey) #vector position start
    #removes already sorted value
    RANNrow.to.remove = which((RANNdata.search[,1] == RANNvaluex) & (RANNdata.search[,2] == RANNvaluey)) #find indices in unsorted curve for valuex and valuey
    RANNdata.search = RANNdata.search[-RANNrow.to.remove,] #delete coordinates from unsorted row given indices
    #some knn stuff
    RANNvals = nn2(RANNdata.search,query = RANNa, k = min(5, nrow(RANNdata.search)), treetype = "kd") #KD-Tree KNN Algorithm.
    RANN.1 = cbind(t(RANNvals$nn.dists),t(RANNvals$nn.idx)) #combine KDTree-Knn output into a single table
    RANN2.1= as.numeric(min(RANN.1[,1])) #find out line distance
    #placing new values into sorted curve
    RANNvector.index = RANN.1[which(RANN.1[1,] == min(RANN2.1)),2]  #find index for smallest vector
    RANNvector.index = na.omit(RANNvector.index)
    RANNnewrow <-  matrix(RANNdata.search[RANNvector.index[1],], nrow = 1,ncol = 2) #new sorted row to add
    RANNsorted.curve <- rbind(RANNsorted.curve, RANNnewrow) #puts new row into sorted curve list
    
    if(nrow(RANNdata.search) == 4) {
      RANNsorted.curve = RANNsorted.curve[c(nrow(RANNsorted.curve)-4:nrow(RANNsorted.curve)),]
      return(RANNsorted.curve)
      break
    }
  }
}
Sens.sorted <- pbapply::pblapply(SensData, RANNsort.curve2)

#Choose Landmark positions
#Now we select semi-landmarks around the curve using geomorph
#calculate.semilandmarks is dependent on the geomorph::digit.curves() function
Sens.semis <- pbapply::pblapply(Sens.sorted, calculate.semilandmarks)
#Format that gpagen can use. 
Sens.semis <- abind::abind(Sens.semis, along = 3)
###Procrustes Superimposition
Sens.semis #Semi-landmark coordinates for each specimen
#Defined points that will slide
#Sens.sliders <- define.sliders(Sens.semis[,1:2,1],50)
#saveRDS(Sens.sliders, file = "Sens.sliders.rds")
Sens.sliders <- readRDS(file = "Sens.sliders.rds")

#This uses a mix of landmarks (start and end point along the outline) and semi-landmarks (50 equally spaced points between start and end point)
Sens.gpa <- gpagen(
  A = Sens.semis,
  curves = Sens.sliders,
  Proj = T,
  print.progress = T
)

#PCA of the aligned points
Sens.pca <- gm.prcomp(Sens.gpa$coords)
summary(Sens.pca)
plot(Sens.pca)
plotTangentSpace(Sens.gpa$coords)
#PCA plot of sensitivity outlines
plot_ly(
  x = ~Sens.pca$x[,1], 
  y = ~Sens.pca$x[,2],
  color = ~Actual.Spec.ID2[-c(9,10,11,23,24,38,51,64)],
  colors = pal,
  type = "scatter",
  mode = "markers+line",
  showlegend = T) %>%
  layout(xaxis = list(title = "PC1", 
                      titlefont = F, 
                      showgrid = T), 
         yaxis = list(title = "PC2", 
                      titlefont = F, 
                      showgrid = T),
         title = "Full Data Set PCA")

SensPCA <- ggplot(mapping = aes(Sens.pca$x[,1],Sens.pca$x[,2], color = Actual.Spec.ID2[-c(9,10,11,23,24,38,51,64)]), showlegend = T) +
  geom_point() +
  coord_fixed() +
  labs(color = "Specimen Identification", x = "PC 1", y = "PC 2") +
  theme_minimal() + 
  labs(x = "PC 1: 79.95%", y = "PC 2: 12.67%")+
  scale_x_continuous(breaks = seq(-0.8,0.8,by = 0.2))

orca(SensPCA, "SensPCA.svg")

sensp1max <- plotRefToTarget(mshape(Sens.pca$A),Sens.pca$shapes$shapes.PC1$max)
sensp1min <- plotRefToTarget(mshape(Sens.pca$A),Sens.pca$shapes$shapes.PC1$min)
sensp2max <- plotRefToTarget(mshape(Sens.pca$A),Sens.pca$shapes$shapes.PC2$max)
sensp2min <- plotRefToTarget(mshape(Sens.pca$A),Sens.pca$shapes$shapes.PC2$min)
}

