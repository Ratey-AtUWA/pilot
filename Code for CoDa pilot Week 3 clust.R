source("welcome.R")

#  load packages etc. ####
library(rgr)
library(car)
library(cluster)
library(factoextra)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(flextable)
library(magrittr)

pal11 <- c("black","blue","green4","red2","purple","darkcyan",
          "firebrick","grey","grey40","white","transparent")
palette(pal11)
# --=----=----=----=----=----=----=----=----=----=----=----=----=----=----=--

#  read file and show data ####
cities <- read.csv("cities_Hu_etal_2021.csv", stringsAsFactors = TRUE)
cities$City <- as.character(cities$City)
# --=----=----=----=----=----=----=----=----=----=----=----=----=----=----=--

#  make new Type factor with abbreviated names ####
row.names(cities) <- as.character(cities$City)
cities$sType <- as.character(cities$Type)
cities$sType <- gsub("Compact-Open","CO",cities$sType)
cities$sType <- gsub("Open-Lightweight","OL",cities$sType)
cities$sType <- gsub("Compact","C",cities$sType)
cities$sType <- gsub("Open","O",cities$sType)
cities$sType <- gsub("Industrial","I",cities$sType)
cities$sType <- as.factor(cities$sType)
# --=----=----=----=----=----=----=----=----=----=----=----=----=----=----=--

#  clr-transform data ####
cities_clr <- cities
cities_clr[,c("Compact","Open","Lightweight","Industry")] <- 
  clr(cities_clr[,c("Compact","Open","Lightweight","Industry")])
names(cities);cat("\n");names(cities_clr)
# --=----=----=----=----=----=----=----=----=----=----=----=----=----=----=--

# We are using the same cities land-use dataset from the previous session, from
# Hu *et al*. (2021).

# K-means clustering ####

# K-means clustering is an unsupervised classification method, which is a type
# of machine learning used when you don't know (or don't want to make
# assumptions about) any categories or groups. We do have groupings in our data,
# which are the factors **Type**, **Global**, and **Region**.

# The goal of the K-means clustering algorithm is to find a specified number
# (*K*) groups based on the data. The algorithm first requires an estimate of
# the number of clusters, and there are several ways to do this. The code below
# from the factoextra R package tests different values of *K* and computes the
# 'total within sum of squares' (WSS) based on the distance of observations from
# the 'centroid' (mean) of each cluster, which itself is found by an iterative
# procedure. When the decrease in WSS from *K* to *K*+1 is minimal, the
# algorithm selects that value of *K*.

#  assess optimum clusters for closed and open data ####

require(factoextra)
data0 <- na.omit(cities[,c("Compact","Open","Lightweight","Industry")])
nclus_clos <- fviz_nbclust(data0, kmeans, method = "wss") +
  geom_vline(xintercept = 3, linetype = 2) +
  labs(title="")
#   geom_vline(xintercept = 4, linetype = 2) +
data0 <- na.omit(cities_clr[,c("Compact","Open","Lightweight","Industry")])
nclus_clr <- fviz_nbclust(data0, kmeans, method = "wss") +
  geom_vline(xintercept = 5, linetype = 2) +
  labs(title="")
ggarrange(nclus_clos,nclus_clr,ncol = 2,nrow = 1, 
          labels = c("(a) closed","(b) open (clr)"))
# --=----=----=----=----=----=----=----=----=----=----=----=----=----=----=--

# We have indicated different numbers of clusters for closed and open data,
# based on visual identification of a break in the slope of the WSS vs.
# 'Number of clusters' curve. In the following analyses we will assume the same
# number of clusters (K = 4) for both open and closed data.

## Compute K-means clustering for closed data ####

#  kmeans closed
data0 <- na.omit(cities[,c("sType","Compact","Open","Lightweight","Industry")])
data0[,c("Compact","Open","Lightweight","Industry")] <- 
  scale(data0[,c("Compact","Open","Lightweight","Industry")])
set.seed(123)
cities_clos_kmeans <- kmeans(data0[,2:NCOL(data0)], 4, nstart = 25)
cat("components of output object are:\n")
ls(cities_clos_kmeans)
cat("\nK-means clustering with",length(cities_clos_kmeans$size),
    "clusters of sizes",cities_clos_kmeans$size,"\n\n")
cat("Cluster centers (scaled to z-scores) in K-dimensional space:\n")
cities_clos_kmeans$centers

# --=----=----=----=----=----=----=----=----=----=----=----=----=----=----=--

#  Cities in each K-means cluster from compositionally closed data ####
outtable <- data.frame(Cluster = seq(1,length(cities_clos_kmeans$size),1),
                       Cities = rep("nil",length(cities_clos_kmeans$size)))
for (i in 1:length(cities_clos_kmeans$size)){
  outtable[i,1] <- paste("Cluster",i)
  outtable[i,2] <- paste(names(which(cities_clos_kmeans$cluster==i)), 
                         collapse = " ")}
print(outtable[,2])

# --=----=----=----=----=----=----=----=----=----=----=----=----=----=----=--

# The output object from the kmeans function is a list which contains the
# information we're interested in: the sum-of-squares between and within
# clusters (betweenss, tot.withinss, totss, withinss), the location in *K*
# dimensions of the centers of the clusters (centers), the assignment of each
# observation to a cluster (cluster), the number of observations in each cluster
# (size), and the number of iterations taken to find the solution (iter).

# Applying K-means clustering with 4 clusters to the **closed** cities land-use
# data results in one larger cluster of 25 cities (3) with three smaller
# clusters containing 4-6 cities (1,2, and 4). From the table of cluster
# centers, we get some idea that:

# - Cluster 1 cities have greater proportions of Lightweight land use
# - Cluster 2 cities have greater proportions of Compact land use
# - Cluster 3 cities have similar proportions of all land uses
# - Cluster 4 cities have greater proportions of Open land use

# Interestingly, Lightweight land use does not discriminate between Clusters
# 2-4, since the values of cluster centers 2-4 are the same in the Lightweight
# dimension.

## Compute K-means clustering for open data ####

#  kmeans open
data0 <- na.omit(cities_clr[,c("sType","Compact","Open","Lightweight","Industry")])
data0[,c("Compact","Open","Lightweight","Industry")] <- 
  scale(data0[,c("Compact","Open","Lightweight","Industry")])
set.seed(123)
cities_open_kmeans <- kmeans(data0[,2:NCOL(data0)], 4, nstart = 25)
cat("components of output object are:\n")
ls(cities_open_kmeans)
cat("\nK-means clustering with",length(cities_open_kmeans$size),
    "clusters of sizes",cities_open_kmeans$size,"\n\n")
cat("Cluster centers (scaled to z-scores) in K-dimensional space:\n")
cities_open_kmeans$centers
outtable <- data.frame(Cluster = seq(1,length(cities_open_kmeans$size),1),
                       Cities = rep("nil",length(cities_open_kmeans$size)))

# --=----=----=----=----=----=----=----=----=----=----=----=----=----=----=--

# Cities in each K-means cluster from analysis of clr-transformed data ####
for (i in 1:length(cities_open_kmeans$size)){
  outtable[i,1] <- paste("Cluster",i)
  outtable[i,2] <- paste(names(which(cities_open_kmeans$cluster==i)), 
                         collapse = " ")}
print(outtable)
# --=----=----=----=----=----=----=----=----=----=----=----=----=----=----=--

# Applying K-means clustering with 4 clusters to the **open** cities land-use
# data results in one larger cluster of 25 cities (Cluster 4 - not the same
# cities as the cluster of 25 from analysis of closed data!). There are three
# smaller clusters containing 1-10 cities (1-3). From the table of cluster
# centers, we get might conclude that:

# - The Cluster 1 city has greater Industry, and lower Compact, land use
# - 10 Cluster 2 cities have greater Compact, and lower Open, land uses
# - 4 Cluster 3 cities have greater Lightweight, and lower Industry land uses
# - 25 Cluster 4 cities have somewhat greater Open land use
# 
# An interesting question to ask is whether the clusters are similar to the
# categories we already have in the dataset (Type, Global, and Region). The
# following plot includes labelling to help us see the relationship of K-means
# clusters to the Type category.

# To represent clustering in 3 or more dimensions in a plot, we can use
# principal components to reduce the number of dimensions, but still retain
# information from all the variables. This is done very nicely by the
# fviz_cluster() function from the R package fctoextra* (Kassambara and Mundt,
# 2020). The output from fviz_cluster() is a ggplot, so for more efficient
# presentation and comparison, we save the ggplot2 output to objects, and plot
# these together using the ggarrange() function from the ggpubr R package.

## Plot kmeans clusters for closed and open data

# visualize K-means clusters compared ####
row.names(data0) <- paste0(data0$sType,seq(1,NROW(data0)))
kmeans_viz_clos<- fviz_cluster(cities_clos_kmeans, data = data0[,2:NCOL(data0)],
             palette = c("#800000", "#E7B800", "#FC4E07","purple2"),
             labelsize=10, main = "",
             ellipse.type = "euclid", # Concentration ellipse
             star.plot = TRUE, # Add segments from centroids to items
             repel = F, # if true avoids label overplotting (slow)
             ggtheme = theme_minimal())
kmeans_viz_open<- fviz_cluster(cities_open_kmeans, data = data0[,2:NCOL(data0)],
             palette = c("#800000", "#E7B800", "#FC4E07","purple2"),
             labelsize=10, main = "",
             ellipse.type = "euclid", # Concentration ellipse
             star.plot = TRUE, # Add segments from centroids to items
             repel = F, # if true avoids label overplotting (slow)
             ggtheme = theme_minimal())
ggarrange(kmeans_viz_clos,kmeans_viz_open,ncol = 2,
          labels = c("(a) closed","(b) open (clr)"))
# --=----=----=----=----=----=----=----=----=----=----=----=----=----=----=--

# From the plots we can see that, for both the closed and open cities land use
# data, K-means clustering did not produce clusters that overlap convincingly
# with the city Type category in our data. There is some differentiation; for
# example in both cases there is a cluster (cluster 1 for closed data, cluster 3
# for open data) composed only of Open-Lightweight (OL) cities. It's also clear
# that removing compositional closure in the data makes a difference to both the
# size of the resulting K-means cluster and their composition (in terms of
# cities belonging to each cluster). We haven't checked the relationship of the
# other categories (Global, Region) to the clusters obtained, but some editing
# of the code would answer this question for us!

# --=----=----=----=----=----=----=----=----=----=----=----=----=----=----=--

# K-means using rock data /~\_-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-_/~\ ####
usgs <- read.csv("https://raw.githubusercontent.com/Ratey-AtUWA/compositional_data/main/usgs_ngdb_trimmed.csv", stringsAsFactors = TRUE)
usgs_ppm <- usgs
# FIRST convert percent columns to ppm!!
names(usgs_ppm[,9:20])
usgs_ppm[,9:20] <- usgs_ppm[,9:20]*10000

# now re-order factor levels to have compositionally similar rocks together...
usgs_ppm$Rock <- 
  factor(usgs_ppm$Rock, levels = c("rhyolite","granite","andesite","dacite",
                                   "basalt","gabbro","pyroxenite","peridotite"))
# ...and let's simplify the groups...
usgs_ppm$Lith <- as.character(usgs_ppm$Rock)
usgs_ppm$Lith <- 
  mgsub::mgsub(usgs_ppm$Lith,c("rhyolite","granite","andesite","dacite",
                                "basalt","gabbro","pyroxenite","peridotite"),
                              c("Felsic","Felsic","Intermediate","Intermediate",
                                "Mafic","Mafic","Ultramafic","Ultramafic"))
usgs_ppm$Lith <- as.factor(usgs_ppm$Lith)
str(tapply(usgs_ppm$Rock,usgs_ppm$Lith,droplevels)) # check it

usgs_clr <- usgs_ppm
# using one-line CLR function from Hugo Salinas Matus
usgs_clr[,9:33] <- t(apply(usgs_ppm[,9:33], MARGIN = 1,
                           FUN = function(x){log(x) - mean(log(x))}))

# plot WSS vs. nClust
rockclus_clos <- fviz_nbclust(scale(log10(usgs_ppm[,9:33])), kmeans, 
                              method = "wss") +
  geom_vline(xintercept = 7, linetype = 2) +
  labs(title="")
rockclus_clr <- fviz_nbclust(usgs_clr[,9:33], kmeans, method = "wss") +
  geom_vline(xintercept = 4, linetype = 2) +
  labs(title="")
ggarrange(rockclus_clos,rockclus_clr,ncol = 2,nrow = 1, 
          labels = c("(a) closed","(b) open (clr)"))

# kmeans closed ####
data0 <- na.omit(usgs_ppm[,c(7,9:34)])
data0[,2:26] <- scale(log10(data0[,2:26]))
row.names(data0) <- 
  paste0(substr(data0$Lith,1,1),seq(1,NROW(data0)))
set.seed(123)
usgs_clos_kmeans <- kmeans(data0[,2:26], 4, nstart = 25)
{cat("components of output object are:\n")
  print(ls(usgs_clos_kmeans))
  cat("\nK-means clustering with",length(usgs_clos_kmeans$size),
      "clusters of sizes",usgs_clos_kmeans$size,"\n\n")
  cat("Cluster centers (scaled to z-scores) in K-dimensional space:\n")
  usgs_clos_kmeans$centers}

# kmeans open CLR ####
data0 <- na.omit(usgs_clr[,c(7,9:34)])
row.names(data0) <- 
  paste0(substr(data0$Lith,1,1),seq(1,length(data0$Lith)))
set.seed(123)
usgs_open_kmeans <- kmeans(scale(data0[,2:26]), 4, nstart = 25)
{cat("components of output object are:\n")
  print(ls(usgs_open_kmeans))
  cat("\nK-means clustering with",length(usgs_open_kmeans$size),
      "clusters of sizes",usgs_open_kmeans$size,"\n\n")
  cat("Cluster centers (scaled to z-scores) in K-dimensional space:\n")
  usgs_open_kmeans$centers}

# visualize K-means clusters compared ####
palette(c("black","#80000080","#E7B80080","#FC4E0780","#912CEE80"))
kmeans_viz_clos<- fviz_cluster(usgs_clos_kmeans, data=data0[,2:26],
                          palette = c(2,3,4,5),
                          labelsize=1, pointsize = 2, main = "",
                          ellipse = F, # Concentration ellipse
                          star.plot = F, # if T Add segments centroids to items
                          repel = F, # if T, no label overplotting (slow)
                          ggtheme = theme_minimal())
kmeans_viz_open<- fviz_cluster(usgs_open_kmeans, data = data0[,2:26],
                          palette = c(2,3,4,5),
                          labelsize=1, pointsize = 2, main = "",
                          ellipse = F, # Concentration ellipse
                          star.plot = F, #if T Add segments  centroids to items
                          repel = F, # if true, no label overplotting (slow)
                          ggtheme = theme_minimal())
ggarrange(kmeans_viz_clos,kmeans_viz_open,ncol = 2,
          labels = c("(a) closed","(b) open (clr)")); palette(pal11)

# Rows in each K-means cluster ####
#  from compositionally closed data
outtable <- data.frame(Cluster = paste("Clust",
           rep(seq(1,length(usgs_clos_kmeans$size),1),each=nlevels(data0$Lith)),
           rep(levels(data0$Lith),length(usgs_clos_kmeans$size))),
           count_clos = rep(0,length(usgs_clos_kmeans$size)),
           pct_clos = rep(0,length(usgs_clos_kmeans$size)))
lithNames <- substr(levels(data0$Lith),1,1)
for(j in 1:length(usgs_clos_kmeans$size)){
  for (i in (4*j-3):(4*j)){
    # outtable[i,1] <- paste("Cluster",i)
    outtable[i,"count_clos"] <- length(grep(lithNames[i-((4*j)-4)],
                                            names(usgs_clos_kmeans$cluster[which(usgs_clos_kmeans$cluster==j)])))
    outtable[i,"pct_clos"] <- 100*round(length(grep(lithNames[i-((4*j)-4)],
                               names(usgs_clos_kmeans$cluster[which(usgs_clos_kmeans$cluster==j)])))/as.numeric(table(data0$Lith)[i-((4*j)-4)]),2)
  }
}
print(outtable)

# from open CLR data
outtable2 <- cbind(outtable, 
                   data.frame(count_open=rep(0,length(usgs_open_kmeans$size)),
                              pct_open=rep(0,length(usgs_open_kmeans$size))))
lithNames <- substr(levels(data0$Lith),1,1)
for(j in 1:length(usgs_open_kmeans$size)) {
  for (i in (4*j-3):(4*j)){
    outtable2[i,"count_open"] <- length(grep(lithNames[i-((4*j)-4)], names(usgs_open_kmeans$cluster[which(usgs_open_kmeans$cluster==j)])))
    outtable2[i,"pct_open"] <- 100*round(length(grep(lithNames[i-((4*j)-4)],names(usgs_open_kmeans$cluster[which(usgs_open_kmeans$cluster==j)])))/as.numeric(table(data0$Lith)[i-((4*j)-4)]),2)
  }
}
print(outtable2)

# Hierarchical clustering USGS data _/~\_/~\_/~\_/~\_/~\_/~\_/~\_/~\_/~\_/~ ####
# make dissimilarity matrices
dataHC <- na.omit(usgs_ppm[,c(7,34,9:33)])
dataHC[,3:27] <- scale(log10(dataHC[,3:27]))
row.names(dataHC) <- 
  paste0(substr(dataHC$Lith,1,1),seq(1,NROW(dataHC)))
usgs_clos_diss <- get_dist(dataHC[,3:27], method = "euclidean")
{cat("First 15 rows and columns of closed (scaled log10) distance matrix:\n")
round(as.matrix(usgs_clos_diss)[1:15, 1:15], 1)}

dataHO <- na.omit(usgs_clr[,c(7,34,9:33)])
dataHO[,3:27] <- scale(dataHO[,3:27])
row.names(dataHO) <- 
  paste0(substr(dataHO$Lith,1,1),seq(1,NROW(dataHO)))
usgs_open_diss <- get_dist(dataHO[,3:27], method = "euclidean")
{cat("First 15 rows and columns of closed (scaled log10) distance matrix:\n")
round(as.matrix(usgs_open_diss)[1:15, 1:15], 1)}

# do the H-clustering
usgs_clos_hc <- hclust(usgs_clos_diss, method = "average") # best method
usgs_open_hc <- hclust(usgs_open_diss, method = "average") # best method

# check them both
usgs_clos_coph <- cophenetic(usgs_clos_hc)
cat("Correlation coefficient r =",cor(usgs_clos_diss,usgs_clos_coph),"\n")
usgs_open_coph <- cophenetic(usgs_open_hc)
cat("Correlation coefficient r =",cor(usgs_open_diss,usgs_open_coph),"\n")
cat("\nRule-of-thumb:\n",
 "Cluster tree represents actual distance matrix accurately enough if r>0.75\n")

# plot them both
par(mar=c(5,4,1,1), mfrow=c(2,1), mgp = c(1.7,0.3,0), tcl = 0.25)
plot(usgs_clos_hc, cex=0.8, main=""); mtext("(a) closed", adj = 0.05,cex=1.5)
abline(h = c(9.95,10.3,12.5), col = c("blue2","purple","red3"), lty = c(2,5,1))
# text(3,870, pos=3, labels="Cut for 5 clusters", col = "blue2", offset = 0.2)
text(c(3600,1500,1200),c(9.95,10.3,12.5), pos=3,
     col = c("blue2","purple","red3"), offset = 0.1, 
     labels=c("Cut for 5 clusters","Cut for 4 clusters","Cut for 3 clusters"))
plot(usgs_open_hc, cex=0.8, main=""); mtext("(b) open CLR", adj = 0.05,cex=1.5)
abline(h = c(9.95,10.3,11.5), col = c("blue2","purple","red3"), lty = c(2,5,1))
# text(3,870, pos=3, labels="Cut for 5 clusters", col = "blue2", offset = 0.2)
text(c(2500,1200,800),c(9.95,10.3,11.5), pos=3, col = c("blue2","purple","red3"), 
     offset = 0.2, 
     labels=c("Cut for 5 clusters","Cut for 4 clusters","Cut for 3 clusters"))

# optimum number of clusters
nclus_clos <- fviz_nbclust(dataHC[,3:27], hcut, 
                           method = "silhouette", verbose = F) +
  labs(title="")
nclus_clr <- fviz_nbclust(scale(dataHO[,3:27]), hcut, 
                          method = "silhouette", verbose = F) +
  labs(title="")
ggarrange(nclus_clos,nclus_clr,ncol = 1,nrow = 2, 
          labels = c("(a) closed","(b) open (clr)"))

# We can also assess optimum cluster numbers using the NbClust() function in
# the 'NbClust' package:
  
  \scriptsize

# NbClust package ####
cat("Closed data:\n")
NbClust::NbClust(data = scale(log10(usgs_ppm[,9:33])), 
                 distance = "euclidean", min.nc = 2, max.nc = 15, 
                 method = "complete", index = "silhouette")
cat("\n-------------------------------------------\nOpen CLR data:\n")
NbClust::NbClust(data = scale(usgs_clr[,9:33]), 
                 distance = "euclidean", min.nc = 2, max.nc = 15, 
                 method = "complete", index = "silhouette")

# cut dendrogram closed
usgs_clos_grp <- cutree(usgs_clos_hc, k = 4)
table(usgs_clos_grp)

# cut dendrogram open
usgs_open_grp <- cutree(usgs_open_hc, k = 4)
table(usgs_open_grp)

## Plot dendrograms comparing closed and open data with cuts
#                                     _           _                _  _   
#                                    ( )_        ( )              ( )( )_ 
#     __           ___    __   _ _   | ,_)      _| |   _     ___  |/ | ,_)
#   /'__`\(`\/') /'___) /'__`\( '_`\ | |      /'_` | /'_`\ /' _ `\   | |  
#  (  ___/ >  < ( (___ (  ___/| (_) )| |_    ( (_| |( (_) )| ( ) |   | |_ 
#  `\____)(_/\_)`\____)`\____)| ,__/'`\__)   `\__,_)`\___/'(_) (_)   `\__)
#                             | |                                         
#                             (_)                                         
#      _    _           _                           _                     
#   _ ( )_ ( )         ( )_                        (_ )                   
#  (_)| ,_)|/   ___    | ,_)   _      _        ___  | |    _    _   _   _ 
#  | || |     /',__)   | |   /'_`\  /'_`\    /',__) | |  /'_`\ ( ) ( ) ( )
#  | || |_    \__, \   | |_ ( (_) )( (_) )   \__, \ | | ( (_) )| \_/ \_/ |
#  (_)`\__)   (____/   `\__)`\___/'`\___/'   (____/(___)`\___/'`\___x___/'
#                                                                         
# plot cut dendrogram open
gg_cutden_clos <- 
  fviz_dend(usgs_clos_hc, k = 4, # Cut in five groups
            main = "", cex = 0.7, # label size
            k_colors = c("gray40","red3", "blue2", "purple"),
            color_labels_by_k = TRUE, # color labels by groups
            rect = TRUE, # Add rectangle around groups
            labels_track_height = 1320 # adjust low margin for long labels
            )
gg_cutden_open <- 
  fviz_dend(usgs_open_hc, k = 4, # Cut in five groups
  main = "", cex = 0.7, # label size
  k_colors = c("gray40","red3", "blue2", "purple"),
  color_labels_by_k = TRUE, # color labels by groups
  rect = TRUE, # Add rectangle around groups
  labels_track_height = 10 # adjust low margin for long labels
)
ggarrange(gg_cutden_clos,gg_cutden_open, nrow = 2, 
          labels = c("(a) closed","(b) open (CLR)"))

# ._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._
#
source("thanks.R")
