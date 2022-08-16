library(factoextra) # plotting
library(ggpubr) # plotting
library(magrittr)# pipes
library(corrplot) # plotting

# read cities data ####
cities <- 
  read.csv(paste0("https://raw.githubusercontent.com/Ratey-AtUWA/",
                  "compositional_data/main/cities_Hu_etal_2021.csv"), 
           stringsAsFactors = TRUE)
# Make new Type factor with abbreviated names
row.names(cities) <- as.character(cities$City)
cities$sType <- as.character(cities$Type)
cities$sType <- gsub("Compact-Open","CO",cities$sType)
cities$sType <- gsub("Open-Lightweight","OL",cities$sType)
cities$sType <- gsub("Compact","C",cities$sType)
cities$sType <- gsub("Open","O",cities$sType)
cities$sType <- gsub("Industrial","I",cities$sType)
cities$sType <- as.factor(cities$sType)
# CLR-transform cities data using Hugo's code
cities_clr <- cities
cities_clr[,2:5] <- t(apply(cities[,2:5], MARGIN = 1,
                           FUN = function(x){log(x) - mean(log(x))}))

# compare prcomp() and princomp() functions on cities data ####
data0 <- cities_clr[,c("Compact","Open","Lightweight","Industry")]
pca_cities_open <- prcomp(data0, scale. = TRUE)
{cat("\n--------------------\nVariable weightings (rotations) - open data"," - using prcomp()\n")
print(pca_cities_open$rotation)
cat("...\n\nComponent Variances - CLR-transformed (open) data - using prcomp()\n")
print(round(pca_cities_open$sdev^2,3))
cat("...\n\nProportions of variance explained by each component - using prcomp()",
    "\nCLR-transformed (open) data\n")
print(round(pca_cities_open$sdev^2/sum(pca_cities_open$sdev^2),3))
cat("\nCumulative proportions of variance explained by each component\n")
print(round(cumsum(pca_cities_open$sdev^2/sum(pca_cities_open$sdev^2)),3))
}

princomp_cities_open <- princomp(data0, cor = TRUE)
{cat("\n--------------------\nVariable weightings (rotations) - open data",
     " - using princomp()\n")
  print(princomp_cities_open$loadings)
  cat("...\n\nComponent Variances - CLR-transformed (open) data",
      " - using princomp()\n")
  print(round(princomp_cities_open$sdev^2,3))
  cat("...\n\nProportions of variance explained by each component",
      " - using princomp()\nCLR-transformed (open) data\n")
  print(round(princomp_cities_open$sdev^2/sum(princomp_cities_open$sdev^2),3))
  cat("\nCumulative proportions of variance explained by each component\n")
  print(round(cumsum(princomp_cities_open$sdev^2/sum(princomp_cities_open$sdev^2)),3))
}


# biplots comparing prcomp() and princomp() ####

biplot_pri <- fviz_pca(princomp_cities_open, title = "",
              col.ind = cities_clr$Type, col.var = "black", 
              pch = c(0,1,2,5,6)[cities_clr$Type],
              palette = c("dodgerblue4","sienna4", "red3", "purple", "gold3")) +
              xlim(-4,5) + ylim(-5,4)
biplot_prc <- fviz_pca(pca_cities_open, title = "", col.ind = cities_clr$Type, 
              col.var = "black",  pch = c(0,1,2,5,6)[cities_clr$Type],
              palette = c("dodgerblue4","sienna4", "red3", "purple", "gold3")) +
              xlim(-5,4) + ylim(-5,4)
ggarrange(biplot_prc,biplot_pri, nrow = 1,
          labels = c("(a) prcomp()","(b) princomp()"),
          common.legend = TRUE, legend = "bottom")

# comparing closed and open for geochemical data ####

# read the usgs data ####
usgs <- read.csv("https://raw.githubusercontent.com/Ratey-AtUWA/compositional_data/main/usgs_ngdb_trimmed.csv", stringsAsFactors = TRUE)
usgs_ppm <- usgs
# FIRST convert percent columns to ppm!!
names(usgs_ppm[,9:20])
usgs_ppm[,9:20] <- usgs_ppm[,9:20]*10000

# now re-order factor levels to have compositionally similar rocks together 
usgs_ppm$Rock <- 
  factor(usgs_ppm$Rock, levels = c("rhyolite","granite","andesite","dacite",
                               "basalt","gabbro","pyroxenite","peridotite"))

usgs_clr <- usgs_ppm
# using one-line function from Hugo Salinas Matus -- THANKS!
  #                                ||
  #  _.      _  _ _  ._ _   _   | _||_
  # (_| \/\/(/__>(_) | | | (/_  o \  / 
  #                                \/
usgs_clr[,9:33] <- t(apply(usgs_ppm[,9:33], MARGIN = 1,
                           FUN = function(x){log(x) - mean(log(x))}))
    # FYI
    library(rgr) # to do a comparison
    usgs_clr_rgr <- usgs_ppm
    usgs_clr_rgr[,9:33] <- clr(usgs_ppm[,9:33])
    usgs_clr[1:80,13] - usgs_clr_rgr[1:80,13]

# PCA USGS closed (all %concs converted to ppm)
pca_usgs_clos <- prcomp(log10(usgs_ppm[,9:33]), scale. = TRUE)
{cat("\n--------------------\n",
     "Variable weightings (rotations 1-5) - closed data\n")
print(pca_usgs_clos$rotation[,1:5], digits=3)
cat("\nComponent Variances 1-5 - closed data\n")
print(pca_usgs_clos$sdev[1:5]^2, digits = 3)
cat("\nProportions of variance explained by components 1-5",
    "\nclosed data\n")
print(round(pca_usgs_clos$sdev[1:5]^2/sum(pca_usgs_clos$sdev[1:5]^2),3))
cat("Cumulative proportions of variance explained by components 1-5\n")
print(round(cumsum(pca_usgs_clos$sdev[1:5]^2/
                     sum(pca_usgs_clos$sdev[1:5]^2)),3))
}

# PCA USGS open (all ppm CLR-transformed)
pca_usgs_open <- prcomp(usgs_clr[,9:33], scale. = TRUE)
{cat("\n--------------------\n",
     "Variable weightings (rotations 1-5) - CLR-transformed (open) data\n")
  print(pca_usgs_open$rotation[,1:5], digits = 3)
  cat("\nComponent Variances 1-5 - CLR-transformed (open) data\n")
  print(pca_usgs_open$sdev[1:5]^2, digits = 3)
  cat("\nProportions of variance explained by components 1-5",
      "\nCLR-transformed (open) data\n")
  print(round(pca_usgs_open$sdev[1:5]^2/sum(pca_usgs_open$sdev[1:5]^2),3))
  cat("Cumulative proportions of variance explained by components 1-5\n")
  print(round(cumsum(pca_usgs_open$sdev[1:5]^2/
                       sum(pca_usgs_open$sdev[1:5]^2)),3))
}

par(mfrow=c(1,2), mar=c(4,4,2,2), mgp = c(2,0.6,0), tcl=-0.2)
plot(pca_usgs_clos, ylim=c(0,7), xlab="Component", col = "thistle")
abline(h=1,lty=2, col="blue")
axis(1, at = seq(.6,10)*1.2, labels = seq(1,10), cex.axis = 0.9)
plot(pca_usgs_open, ylim=c(0,7), xlab="Component", col = "moccasin")
abline(h=1, lty=2, col="blue")
axis(1, at = seq(.6,10)*1.2, labels = seq(1,10), cex.axis = 0.9)

#### most important variables in each component ####
# Sort by PC1 -- closed
round(pca_usgs_clos$rotation[rev(order(abs(pca_usgs_clos$rotation[,1]))),1],3)
# Sort by PC1 -- open
round(pca_usgs_open$rotation[rev(order(abs(pca_usgs_open$rotation[,1]))),1],3)

# Sort by PC2 -- closed
round(pca_usgs_clos$rotation[rev(order(abs(pca_usgs_clos$rotation[,2]))),2],3)
# Sort by PC2 -- open
round(pca_usgs_open$rotation[rev(order(abs(pca_usgs_open$rotation[,2]))),2],3)

# ...or using pipes if you don't like lots of nesting... ####
require(magrittr)
sort0 <- pca_usgs_clos$rotation[,1] %>% abs() %>% order() %>% rev()
cat("-=-=- Closed data PC1 -=-=-\n");round(pca_usgs_clos$rotation[sort0,1],3)

sort0 <- pca_usgs_open$rotation[,1] %>% abs() %>% order() %>% rev()
cat("-=-=- Open data PC1 -=-=-\n");round(pca_usgs_open$rotation[sort0,1],3)

sort0 <- pca_usgs_clos$rotation[,2] %>% abs() %>% order() %>% rev()
cat("-=-=- Closed data PC2 -=-=-\n");round(pca_usgs_clos$rotation[sort0,2],3)

sort0 <- pca_usgs_open$rotation[,2] %>% abs() %>% order() %>% rev()
cat("-=-=- Open data PC2 -=-=-\n");round(pca_usgs_open$rotation[sort0,2],3)
rm(sort0)

# compare correlation matrices ####
cor_usgs_clos <- cor(log10(usgs_ppm[,9:33]))
cor_usgs_open <- cor(usgs_clr[,9:33])
require(corrplot)
par(mar=c(1,5,5,3), mfrow=c(1,2), mgp=c(1.5,0.2,0), xpd=T)
corrplot(cor_usgs_clos, title = "\nClosed log10-transformed", tl.col=1, diag=F)
text(par("usr")[2],-1,pos=2,labels="(Negative = red)",col="#B2182B",font=2,cex=1.2)
corrplot(cor_usgs_open, title = "\nOpen CLR-transformed", tl.col=1, diag=F)
text(par("usr")[1],-1,pos=4,labels="(Positive = blue)",col="#1F63A8",font=2,cex=1.2)

# also covariances
cov_usgs_clos <- cov(log10(usgs_ppm[,9:33]))
cov_usgs_open <- cov(usgs_clr[,9:33])
covPal <- colorRampPalette(c("#551A8B", "#A020F0", "#FFFFFF",
                             "#FFA500","#CD8500"))(200)
par(mar=c(1,5,5,3), mfrow=c(1,2), mgp=c(1.5,0.2,0), xpd=T)
corrplot(cov_usgs_clos, title = "\nClosed: log10-transformed", tl.col=1, diag=F,
         is.corr = F, col = covPal)
text(29.4,-1,pos=2,labels="Least covariance",col=covPal[20],font=2,cex=1.2)
corrplot(cov_usgs_open, title = "\nOpen: CLR-transformed", tl.col=1, diag=F,
         is.corr = F, col = covPal)
text(-5,-1,pos=4,labels="Greatest covariance",col=covPal[180],font=2,cex=1.2)

# Finally the biplots ####
# Hint: don't use the 'habillage' option in fviz_pca(). It's veeeeery sloooooow.

fourpair <- c("#B24040","#FF5B5B","#B28559","#E2AA71",
              "#5994B2","#7FD4FF","#9459B2","#D27FFF")
biplot_cl <- fviz_pca(pca_usgs_clos, title = "", geom="point", col.var = 1,
                      pch = rep(c(15,17),4)[usgs_clr$Rock],
                      col.ind = usgs_ppm$Rock, palette = fourpair,
                      alpha.ind = 0.4) + 
  xlim(-7.5,10) + ylim(-10.5,6.5)
biplot_op <- fviz_pca(pca_usgs_open, title = "", geom="point", col.var = 1,
                      pch = rep(c(15,17),4)[usgs_clr$Rock],
                      col.ind = usgs_clr$Rock, palette = fourpair,
                      alpha.ind = 0.4) + 
  xlim(-7.5,10) + ylim(-10.5,6.5)
ggarrange(biplot_cl,biplot_op, nrow = 1,
          labels = c("Closed: log10-transformed", "Open: CLR-transformed"),
          legend = "bottom", common.legend = TRUE)

#
#
#
#
#
#
#
# rm(list=c("data0","fourpair", "biplot_cl","biplot_op","cor_usgs_clos","cor_usgs_open"))