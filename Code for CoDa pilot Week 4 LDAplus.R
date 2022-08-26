# load packages etc
library(rgr)
library(MASS) # essential for LDA !
library(car)
library(factoextra)
library(ggplot2)
library(ggpubr)
library(RcmdrMisc)
library(klaR)
library(TeachingDemos)
library(corrplot)
library(Hmisc)

UWApal <- c("black", "#003087", "#DAAA00", "#8F92C4", "#E5CF7E", 
            "#001D51", "#B7A99F", "#A51890", "#C5003E", "#FDC596", 
            "#AD5D1E","gray40","gray85","#FFFFFF","transparent")
palette(UWApal)

# We are using a curated version of a whole rock major element dataset from
# Hallberg (<https://catalogue.data.wa.gov.au/dataset/hallberg-geochemistry>)

gitpath <- "https://raw.githubusercontent.com/Ratey-AtUWA/pilot/main/"
Hallberg <- read.csv(paste0(gitpath,"Hallberg.csv"), stringsAsFactors = TRUE)

# make new Rock factor with abbreviated names ####
row.names(Hallberg) <- paste0(as.character(Hallberg$Rock),seq(1:NROW(Hallberg)))
Hallberg$sRock <- as.character(Hallberg$Rock)
Hallberg$sRock <- gsub("Basaltic Trachyandesite","C",Hallberg$sRock)
Hallberg$sRock <- gsub("Basaltic Andesite","L",Hallberg$sRock)
Hallberg$sRock <- gsub("Andesite","A",Hallberg$sRock)
Hallberg$sRock <- gsub("Trachybasalt","T",Hallberg$sRock)
Hallberg$sRock <- gsub("Basalt","B",Hallberg$sRock)
Hallberg$sRock <- gsub("Basanite","N",Hallberg$sRock)
Hallberg$sRock <- gsub("Phonotephrite","P",Hallberg$sRock)
Hallberg$sRock <- as.factor(Hallberg$sRock)
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

# Using the *additive log-ratio* transformation to remove closure ####

# The additive logratio (ALR) transformation calculates the logarithm of the
# ratio of variables to one selected variable (which obviously then can not be
# used for further data analysis. The variable used as the denominator is often
# one that behaves "conservatively", such as Al2O3 or TiO2 in geochemistry; see
# the equation below (modified from Grunsky, 2010).
                                                 
# ALR_{i} = log(C_{i}/C_{ref}) for (i = 1, ..., N-1),  where
#
# - ALR_{i} are the ALR-transformed variables 
# - C_{i} are the concentrations of the i elements 
# - C_{ref} is the concentration of the reference variable (denominator)
# - N is the total number of variables.
#
# For multivariate classification methods, Campbell et al. (2009) argue that the
# ALR (or related isometric log-ratio, ILR) transformation is more appropriate
# than the centered log-ratio transformation (CLR) that we have used so far.
                                                 
# alr-transform data ####
Hallberg_alr <- Hallberg
# don't use rgr::alr()
#     Hallberg_alr[,c(11:12,14:24)] <- alr(Hallberg_alr[,11:24], j = 3)
Hallberg_alr[,11:24] <- t(apply(Hallberg[,c(11:24)], MARGIN = 1, 
                                FUN = function(x){log(x) - log(x[3])}))
# then check...
summary(Hallberg_alr[,12:14])
# ... and remove reference variable
Hallberg_alr$Al2O3 <- NULL
names(Hallberg_alr)
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

# LDA -- Linear Discriminant Analysis ####

# NOTE: To implement LDA in R, we use the lda() function in the MASS package
# (Venables and Ripley, 2002). We specify the variables to be considered using a
# formula similar to that used for multiple regression, and we set the prior
# (initial) probabilities of an observation being in a particular category at
# the actual frequencies at which they occur in the data.
                                                 
# We need scale variables to Z-scores, even if we have applied a log-ratio
# transformation (e.g. ALR) to remove compositional closure.
                                                 
## LDA on closed whole rock major element data
                                                 
# We will use LDA to discriminate the rock type, contained in the column 'Rock'
# in the Hallberg dataset. The variables we will use are the major element oxide
# contents (except Al2O3 which was our ALR denominator), SiO2, TiO2, Fe2O2, FeO,
# MnO, MgO, CaO, Na2O, K2O, & P2O5. We're again excluding LOI, CO2 and H2O.

# LDA closed whole rock data ####
data0 <- na.omit(Hallberg)
data0[,c(11:23)] <- scale(log10(data0[,11:23])) # transf just numeric variables
lda_rock_clos <- lda(formula = Rock ~ SiO2 + TiO2 + Fe2O3 + FeO + MnO + 
                       MgO + CaO + Na2O + K2O + P2O5, 
                    data = data0,
                    prior = as.numeric(summary(Hallberg$Rock))/
                      nrow(Hallberg)) 
print(lda_rock_clos)
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

## LDA on open (ALR) whole rock major element data

# LDA open ALR whole rock data ####
data1 <- Hallberg_alr
data1[,11:23] <- scale(data1[,11:23]) # scale just numeric variables
lda_rock_open <- lda(formula = Rock ~ SiO2 + TiO2 + Fe2O3 + FeO + MnO + 
                       MgO + CaO + Na2O + K2O + P2O5, 
                    data = data1,
                    prior = as.numeric(summary(data0$Rock))/nrow(data0)) 
print(lda_rock_open)
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

# The LDA model specified and the prior probabilities are reported back to us in
# the output. The 'Group means' sub-table in the output contains the centroids
# of the categories (8 classes, in our whole rock data) in the n-dimensional
# space defined by the n variables used for classification. The 'Coefficients of
# linear discriminants' sub-table essentially defines the linear discriminant
# functions that separate our categories, since each function (LD1, LD2, etc.)
# is a linear combination of all the variables. Finally, the 'Proportion of
# trace' sub-table gives the _proportions of between-class variance_ that are
# explained by successive discriminant functions (e.g. for the open Hallberg
# rock data, LD1 explains 0.717 (72%) and LD2 explains 0.171 (17%) of variance
# between Rock type categories).

# [NOTE: sometimes we cannot make an LDA model, because the predictors are
# collinear (highly correlated). We may be able to fix this by inspecting a
# correlation matrix for all the predictor variables, and removing one
# variable at a time from correlated pairs, then re-running the LDA procedure.
# By analogy with multiple linear regression, this would mean a
# Pearson's r >= 0.8.]

## Visualising LDA separation ####

### LDA histograms

# We can plot the separation achieved by each linear discriminant (LD) function
# by predicting the classification using the input data, then using the
# **ldahist** function (Venables and Ripley 2002). To see the separation in
# another LDA dimension, we change the subscript in the **predClos\$x[,1]**
# option. Histograms (actually drawn with custom code to plot side-by side) are
# shown in \autoref{ldahistComp}.

# lda hist compare
predClos <- predict(lda_rock_clos, Hallberg[,11:21]) # not CO2, H2O
predOpen <- predict(lda_rock_open, Hallberg[,11:21])
LD1c <- data.frame(Rock=as.character(Hallberg$Rock),LD1=predClos$x[,1])
LD1c$Rock <- factor(LD1c$Rock, levels=levels(Hallberg$Rock))
par(mfcol = c(nlevels(LD1c$Rock),2), mar = c(1,2,1,1), oma = c(1,0,1,0), 
    mgp = c(0.75,0.2,0), tcl=0.15)
for (i in 1:nlevels(LD1c$Rock)){
  with(subset(LD1c, subset=LD1c$Rock==levels(LD1c$Rock)[i]),
       hist(LD1, main = "", breaks = pretty(LD1c$LD1, n=20), col=5,
            xlim = c(min(LD1c$LD1, na.rm=T),max(LD1c$LD1, na.rm=T))))
  box()
  mtext(levels(LD1c$Rock)[i],3,-1.55,adj=0.505, cex = 0.85, font = 2, col=14)
  mtext(levels(LD1c$Rock)[i],3,-1.5, cex = 0.85, font = 2, col = 11)
  if(i==1) mtext("(a) Closed data", 3, 0.5, font=2, cex = 1.4, col = 11)
}

LD1o <- data.frame(Rock=as.character(Hallberg$Rock),LD1=predOpen$x[,1])
LD1o$Rock <- factor(LD1o$Rock, levels=levels(Hallberg$Rock))
for (i in 1:nlevels(LD1o$Rock)){
  with(subset(LD1o, subset=LD1o$Rock==levels(LD1o$Rock)[i]),
       hist(LD1, main = "", breaks = pretty(LD1o$LD1, n=20), col=4,
            xlim = c(min(LD1o$LD1, na.rm=T),max(LD1o$LD1, na.rm=T))))
  box()
  mtext(levels(LD1o$Rock)[i],3,-1.55, adj=0.505, cex = 0.85, font = 2, col=14)
  mtext(levels(LD1o$Rock)[i],3,-1.5, cex = 0.85, font = 2, col = 2)
  if(i==1) mtext("(b) Open data", 3, 0.5, font=2, cex = 1.4, col = 2)
}
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

# The sets of histograms for closed and open data in \autoref{ldahistComp} both
# show some separation of categories, but with overlap. Of course this only
# shows the separation in one dimension, and two or more dimensions may be
# needed to achieve clear separation. We will make plots showing more than one
# LDA dimension later.

### Partition plots

# Another potentially useful way of showing separation of groups in LDA is to
# use a Partition Plot, accessible using the 'partimat()' function from the
# klaR R package (Weihs et al.v 2005).

\scriptsize

# partition plot comparison
require(klaR)
par(mfrow = c(1,2),mar = c(3,3,1,1), mgp= c(1.3,0.2,0), tcl=0.2, font.lab=2)
with(Hallberg, 
     drawparti(Rock, SiO2, TiO2, method="lda",
               image.colors = rainbow(8,s=0.2,v=0.85,end=0.8),
               xlab = expression(bold(paste(SiO[2]," (%)"))), 
               ylab = expression(bold(paste(TiO[2]," (%)"))))
)
mtext("(a)", 3, -1.5, adj=0.05, font=2, cex = 1.2, col=1)
with(Hallberg_alr, 
     drawparti(Rock, SiO2, TiO2, method="lda",
               image.colors = rainbow(8,s=0.2,v=0.85,end=0.8),
               xlab = expression(bold(paste(SiO[2]," (ALR-transformed)"))), 
               ylab = expression(bold(paste(TiO[2]," (ALR-transformed)"))))
)
mtext("(b)", 3, -1.5, adj=0.05, font=2, cex = 1.2)
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

# Partition plots are presented for single pairwise combinations of the
# variables (in this example SiO2 and TiO2) used to make the LDA model. We can
# make different plots by specifying different variables in the 'drawparti()'
# function. Each such plot can be considered to be a different view of the data
# (which of course has multiple dimensions). Colored regions delineate each
# classification area. Any observation that falls within a region is predicted
# to be from a specific category, with apparent mis-classification in a
# different color (but we usually need more than two dimensions for correct
# classification). Each plot also includes the apparent error rate for that view
# of the data.

## Scatter Plots resembling biplots ####

# Scatter-plots showing each variable and observation in linear discriminant
# dimensions, and grouped by category, are useful for visual assessment of how
# well the LDA model separates the observations.

# plot LDA closed ####
par(mfrow = c(2,2), mar = c(3.5,3.5,1,1), mgp = c(1.5,0.2,0), tcl = 0.25,
    lend = "square", ljoin = "mitre", cex.main = 1.1, cex.lab=1.4, font.lab=2)
plot(lda_rock_clos$scaling[,1], lda_rock_clos$scaling[,2],
     xlim = c(-1.5,2), 
     xlab="Linear Discriminant [1]", ylab="Linear Discriminant [2]", 
     main="(a) Variable Coefficients [LD1, LD2]")
abline(v=0,col="grey",lty=2)
abline(h=0,col="grey",lty=2)
text(lda_rock_clos$scaling[,1],lda_rock_clos$scaling[,2],
     labels=names(Hallberg)[c(11,12,14:21)],
     pos=c(2,2,4),cex=1.2,col="blue2",offset=0.2)
mtext("(a)", 3, -1.5, adj = 0.05, cex = 1.2, font = 2)

ldaPred_rock_clos <- predict(lda_rock_clos)
for(i in 1:NROW(lda_rock_clos$scaling)){
  arrows(0,0,lda_rock_clos$scaling[i,1],lda_rock_clos$scaling[i,2],
         length = 0.1, col = 7)
}
legend("bottomright",legend=levels(Hallberg$Rock),
       ncol = 1, bty="o", inset=0.015, col=c(1:3,6,8,9,11,12), 
       title="Rock Type in (b) - (d)", box.col=15,
       pch=c(0:7), pt.lwd = 2, pt.cex = 1.5, cex = 1.4, y.intersp = 0.9)

plot(ldaPred_rock_clos$x[,1], ldaPred_rock_clos$x[,2],
     col=c(1:3,6,8,9,11,12)[Hallberg$Rock],
     pch=c(0:7)[Hallberg$Rock],
     lwd = 2, cex = 1.5, 
     xlab="Linear Discriminant [1]", ylab="Linear Discriminant [2]", 
     main="Predictions for Observations [LD1, LD2]")
abline(v=0,col="grey",lty=2)
abline(h=0,col="grey",lty=2)
mtext("(b)", 3, -1.5, adj = 0.95, cex = 1.2, font = 2)

plot(ldaPred_rock_clos$x[,1], ldaPred_rock_clos$x[,3], 
     col=c(1:3,6,8,9,11,12)[Hallberg$Rock],
     pch=c(0:7)[Hallberg$Rock], 
     lwd = 2, cex = 1.5, 
     xlab="Linear Discriminant [1]", ylab="Linear Discriminant [3]", 
     main="Predictions for Observations [LD1, LD3]")
abline(v=0,col="grey",lty=2)
abline(h=0,col="grey",lty=2)
# text(ldaPred_rock_clos$x[,1], ldaPred_rock_clos$x[,3], 
#      labels=Hallberg$Rock, col=c(1:13)[Hallberg$Rock],
#      pos=1, offset=0.15, cex=0.65)
mtext("(c)", 3, -1.5, adj = 0.05, cex = 1.2, font = 2)

plot(ldaPred_rock_clos$x[,2], ldaPred_rock_clos$x[,4],
     col=c(1:3,6,8,9,11,12)[Hallberg$Rock],
     pch=c(0:7)[Hallberg$Rock], 
     lwd = 2, cex = 1.5, 
     xlab="Linear Discriminant [2]", ylab="Linear Discriminant [4]", 
     main="Predictions for Observations [LD2, LD4]")
abline(v=0,col="grey",lty=2)
abline(h=0,col="grey",lty=2)
# text(ldaPred_rock_clos$x[,2], ldaPred_rock_clos$x[,4], 
#      labels=Hallberg$Rock, col=c(2,4,6,3,12)[Hallberg$Rock],
#      pos=1, offset=0.15, cex=0.65)
mtext("(d)", 3, -1.5, adj = 0.55, cex = 1.2, font = 2)
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

# plot LDA open ####
par(mfrow = c(2,2), mar = c(3.5,3.5,1,1), mgp = c(1.5,0.2,0), tcl = 0.25,
    lend = "square", ljoin = "mitre", cex.main = 1.1, cex.lab=1.4, font.lab=2)
plot(lda_rock_open$scaling[,1], lda_rock_open$scaling[,2], 
     xlim = c(-2.8,0.6), ylim = c(-1.3,1.3), 
     xlab="Linear Discriminant [1]", ylab="Linear Discriminant [2]", 
     main="Variable Coefficients [LD1, LD2]")
abline(v=0,col="grey",lty=2); abline(h=0,col="grey",lty=2)
text(lda_rock_open$scaling[,1], lda_rock_open$scaling[,2], 
     labels=names(Hallberg_alr)[11:21],
     pos = c(1,3,2,4), cex = 1.2, col = "blue2", offset=0.2)
for(i in 1:NROW(lda_rock_open$scaling)){
  arrows(0,0,lda_rock_open$scaling[i,1],lda_rock_open$scaling[i,2],
         length = 0.1, col = 7) }
mtext("(a)", 3, -1.5, adj = 0.95, cex = 1.2, font = 2)
legend("bottomleft", ncol = 1, legend=levels(Hallberg_alr$Rock), 
       col=c(1:3,6,8,9,11,12), pch=c(0:7), pt.lwd = 2,
       bty="o", box.col=15, inset=0.001, x.intersp = 1, 
       title=expression(bold("Rock Type in (b) - (d)")),
       box.lwd=2, pt.cex=1.5, cex=1.4, title.cex = 1.4)

ldaPred_rock_open <- predict(lda_rock_open)

plot(ldaPred_rock_open$x[,1], ldaPred_rock_open$x[,2], 
     col=c(1:3,6,8,9,11,12)[Hallberg_alr$Rock],
     pch=c(0:7)[Hallberg_alr$Rock], lwd=2, 
     cex = 1.5, 
     main="Predictions for Observations [LD1, LD2]", 
     xlab="Linear Discriminant [1]", ylab="Linear Discriminant [2]")
abline(v=0,col="grey",lty=2); abline(h=0,col="grey",lty=2)
# text(ldaPred_rock_open$x[,1], ldaPred_rock_open$x[,2], labels=Hallberg$Rock, 
#      col=c(2,4,6,3,12)[Hallberg$Rock], pos=1, offset=0.15, cex=0.65)
mtext("(b)", 3, -1.5, adj = 0.95, cex = 1.2, font = 2)

plot(ldaPred_rock_open$x[,1], ldaPred_rock_open$x[,3], 
     col=c(1:3,6,8,9,11,12)[Hallberg_alr$Rock],
     pch=c(0:7)[Hallberg_alr$Rock], lwd=2, 
     cex = 1.5, 
     main="Predictions for Observations [LD1, LD3]", 
     xlab="Linear Discriminant [1]", ylab="Linear Discriminant [3]")
abline(v=0,col="grey",lty=2); abline(h=0,col="grey",lty=2)
# text(ldaPred_rock_open$x[,1], ldaPred_rock_open$x[,3], labels=Hallberg$Rock, 
#      col=c(2,4,6,3,12)[Hallberg$Rock], pos=1, offset=0.15, cex=0.65)
mtext("(c)", 3, -1.5, adj = 0.05, cex = 1.2, font = 2)

plot(ldaPred_rock_open$x[,2], ldaPred_rock_open$x[,4], 
     col=c(1:3,6,8,9,11,12)[Hallberg_alr$Rock],
     pch=c(0:7)[Hallberg_alr$Rock], lwd=2, 
     cex = 1.5, 
     main="Predictions for Observations [LD2, LD4]", 
     xlab="Linear Discriminant [2]", ylab="Linear Discriminant [4]")
abline(v=0,col="grey",lty=2); abline(h=0,col="grey",lty=2)
# text(ldaPred_rock_open$x[,2], ldaPred_rock_open$x[,3], labels=Hallberg$Rock, 
#      col=c(2,4,6,3,12)[Hallberg$Rock], pos=1, offset=0.15, cex=0.65)
mtext("(d)", 3, -1.5, adj = 0.95, cex = 1.2, font = 2)
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

# reset par
par(mfrow=c(1,1))
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

# Sometimes (not so much in this example) we see a lot of clustering of the
# predictor variables in variable coefficient plots (a), which may relate to
# spurious relationships between variables caused by compositional closure. This
# clustering is often removed by ALR-transformation.

# Out of the combinations of LDA dimensions selected, the best separation is
# with LD2 vs. LD1, which is not surprising since these dimensions together
# account for about 95% of the between-groups variance -- for both closed and
# open data. There is a lot of apparent overlap, but we can not see the true
# separation in only 2 dimensions. The LDA performed on open data may result in
# slightly clearer separation of samples by Rock category than LDA using closed
# data, but without a multidimensional view this is also hard to be sure about.

# One way that we can get a better idea about the usefulness of our
# classification models is to perform some validation. This involves 'training'
# the model on a subset of our data, and trying to predict the category of a
# different subset using the training model. 

# ├─┴┬┴┬┴┬┴┬┴┬┴┬┴┬┴┬┴┤ ####
# ─────────────────────────────────────────────────────────────────

# Validation of LDA Models for Supervised Classification ####

# make objects for LDA preds, results='hold'}
ldaPred_rock_clos <- predict(lda_rock_clos)
ldaPred_rock_open <- predict(lda_rock_open)

## Inspecting the agreement between actual and predicted categories in LDA 

# To do this easily we just make an R data frame with columns for the actual
# categories (from the original data frame) and the predicted categories (from
# the prediction objects we just made). We add a column telling us if these two
# columns match in each row (which we can see easily, but we use this column to
# calculate a numerical prediction accuracy).

# The code below uses the *head()* function to inspect the first few rows of
# each comparison, but we could easily look at the whole comparison data frames
# using *print()*.

# compare models with reality
closComp <- data.frame(Actual = as.character(Hallberg$Rock),
                       Predicted = as.character(ldaPred_rock_clos$class))
closComp$test <- as.character(Hallberg_alr$Rock) == 
  as.character(ldaPred_rock_clos$class)
k = length(which(closComp$test == TRUE))
cat("Predictions by LDA using closed data:",k,"out of",NROW(Hallberg_alr),
    "=",paste0(round(100*k/NROW(Hallberg_alr),1),"% correct\n"))
head(closComp, n = 10)

openComp <- data.frame(Actual = as.character(Hallberg_alr$Rock),
                       Predicted = as.character(ldaPred_rock_open$class))
openComp$test <- as.character(Hallberg_alr$Rock) == 
  as.character(ldaPred_rock_open$class)
k = length(which(openComp$test == TRUE))
cat("\nPredictions by LDA using open data:",k,"out of",NROW(Hallberg_alr),
    "=",paste0(round(100*k/NROW(Hallberg_alr),1),"% correct\n"))
head(openComp, n = 10)

# For this dataset, it seems as though LDA using either closed open open data is
# not that good at predicting the Rock category for each observation! In the
# output above, Basalt is mis-identified as Dolerite (which might make sense
# given the expected compositional similarity; simplistically, dolerite is a
# more coarse-grained version of basalt). It seems to be more common for Basalt
# to be mis-identified as Metasediment, which may or not make sense depending on
# the composition of the original sediment!
#
# This kind of comparison is not very rigorous, and nor does it address the
# reason we might perform a supervised classification like LDA -- to use data to
# predict *unknown* categories. The ability of LDA to predict unknown categories
# can be addressed by validation procedures, such as the one we investigate
# below.

## Assessment of LDA prediction using a training-validation method

# This can be done a different way, by including the **CV = TRUE** option in the
# **lda()** function, which implements 'leave one out cross validation'. Simply
# stated, this omits one observation at a time, running the LDA on the remaining
# data each time, and predicting the probability of the missed observation being
# in each category (the *posterior probabilities*). [Each observation is
# assigned to the category with the greatest posterior probability.]

# We will use a related method, using 'training' and 'validation' subsets of our
# data. The idea here is that we divide the dataset into two (the code below
# splits the data in half, but other splits could be used, *e.g*. 0.75:0.25).
# The first subset of the data is used to generate an LDA model (*i.e*. a set of
# linear discriminant functions), which are then used to try and predict the
# categories in the other subset. We choose the observations making up each
# subset randomly. Of course, this could give us unreliable results if the
# random selection happens to be somehow unbalanced, so we repeat the
# training-validation process many times to calculate an average prediction
# rate.

# One glitch that can happen is that in selecting a random subset of the data,
# the subset may not include any samples from one or more categories. This
# problem is more likely if our data have some categories having relatively few
# observations. The code below applies an 'error-catching' condition before
# running the training LDA, which is that all categories in a random subset need
# to be populated.

# built-in Leave-One-Out cross-validation ####
        # -- use option CV = TRUE in lda() function

# closed data
lda_rock_clos <- lda(formula = Rock ~ SiO2 + TiO2 + Fe2O3 + FeO + MnO + 
                       MgO + CaO + Na2O + K2O + P2O5, 
                     data = data0,
                     prior = as.numeric(summary(data0$Rock))/
                       nrow(data0), CV = TRUE) 
ls(lda_rock_clos)
loocv <- round(lda_rock_clos$posterior,3)
colnames(loocv)<- substr(colnames(lda_rock_clos$posterior),1,7)
print(loocv)
clipr::write_clip(loocv) # then paste to Excel
# calc proportion correct
ldaComp <- data.frame(Observed=Hallberg$Rock,Predicted=lda_rock_clos$class)
cat("\nProportion of correct classifications",
    signif(100*length(which((ldaComp[,1]==ldaComp[,2])==TRUE))/
             length(ldaComp[,1]),3), "%\n")
# ...so about two thirds correctly classified...

# open data
lda_rock_open <- lda(formula = Rock ~ SiO2 + TiO2 + Fe2O3 + FeO + MnO + 
                       MgO + CaO + Na2O + K2O + P2O5, 
                     data = data1, CV = TRUE,
                     prior = as.numeric(summary(data1$Rock))/nrow(data1)) 
print(lda_rock_open$posterior, digits = 2)
loocv <- round(lda_rock_open$posterior,3)
colnames(loocv)<- substr(colnames(lda_rock_open$posterior),1,7)
print(loocv)
ldaComp <- data.frame(Observed=Hallberg$Rock,Predicted=lda_rock_open$class)
cat("\nProportion of correct classifications",
  signif(100*length(which((ldaComp[,1]==ldaComp[,2])==TRUE))/
          length(ldaComp[,1]),3), "%\n")
# ...only 62% correctly classified !!

# LDA train and predict closed ####

n0 <- 100 # number of iterations
ftrain <- 0.5 # proportion of observations in training set
results <- data.frame(
  Rep = rep(NA, n0),
  matches = rep(NA, n0),
  non_matches = rep(NA, n0),
  success = rep(NA, n0))
train <- sample(1:NROW(Hallberg), round(NROW(Hallberg) * ftrain,0))
# make vector of individual category non-matches
matchByClass <- 
  data.frame(Match1 = rep(0,nlevels(Hallberg$Rock[train]))) 
rownames(matchByClass) <- levels(ldaPred_rock_clos$class)
fMatchXClass <- 
  data.frame(PctMch1 = rep(0,nlevels(Hallberg$Rock[train]))) 
rownames(fMatchXClass) <- levels(ldaPred_rock_clos$class)
# make vector of cumulative category counts in Hallberg[-train] iterations
cc0 <- rep(0,nlevels(Hallberg$Rock))
isOK <- 0 ; i <- 2

for (i in 1:n0) {
  train <- sample(1:NROW(Hallberg), round(NROW(Hallberg) * ftrain,0))
  # set condition requiring all categories to be populated
  if (is.na(match(NA,tapply(Hallberg[train,]$SiO2, 
                            Hallberg[train,]$Rock, sd, na.rm=T))) == TRUE) {
    lda_Rock_train <- lda(formula = Rock ~ SiO2 + TiO2 + Fe2O3 + FeO + 
                            MnO + MgO + CaO + Na2O + K2O + P2O5, 
                          data = Hallberg[train,],
                          prior=as.numeric(summary(Hallberg$Rock[train]))/
                            nrow(Hallberg[train,]))
    ldaPred_rock_clos <- predict(lda_Rock_train, Hallberg[-train,])
    isOK <- isOK + 1
  }
  
  k=0               # number of matches
  m0 <-             # vector of individual category matches 
    as.matrix(rep(0,nlevels(Hallberg$Rock[train]))) 
  rownames(m0) <- levels(Hallberg$Rock)
  m1 <-             # vector of fractional category matches 
    as.matrix(rep(0,nlevels(Hallberg$Rock[train]))) 
  rownames(m1) <- levels(Hallberg$Rock)
  for (jM in 1:NROW(Hallberg[-train,])) {
    for (jS in 1:nlevels(ldaPred_rock_clos$class)) {
      if((ldaPred_rock_clos$class[jM] == levels(ldaPred_rock_clos$class)[jS]) & 
         (Hallberg$Rock[-train][jM] == levels(ldaPred_rock_clos$class)[jS]) ) 
        m0[jS] = m0[jS] + 1
      else  m0[jS] = m0[jS] 
    }
    k = sum(m0)
  }
  cc0 <- cc0 + as.numeric(summary(Hallberg$Rock[-train]))
  m1 <- round(100*m0/as.numeric(summary(Hallberg$Rock[-train])),1)
  matchByClass[,paste0("Match",i)] <- m0
  fMatchXClass[,paste0("PctMch",i)] <- m1
  # output to results data frame: iteration, matches, non-matches, proportion matched
  results[i,] <- c(i, k, NROW(Hallberg[-train,])-k, 
                   signif(k/NROW(Hallberg[-train,]),3))
}
# Output code block
cat(paste("[Based on", n0, "random subsets of",paste0(100*ftrain,"%"),
          "of the dataset to train LDA model\n",
          "      to predict remaining observations]\n"))
cat("Number of obs. in random subsets =",NROW(train),
    " (predicting",NROW(Hallberg)-NROW(train),"samples)\n")
print(numSummary(results[,2:4], statistics=c("mean","sd"))$table)
ns0 <- numSummary(results$success)
t0 <- t.test(results$success)
cat(rep("-\u2013-",24),
    "\nStat. summary for 'success':\nMean = ",round(ns0$table[1],4),
    ", sd = ",round(ns0$table[2],4),
    ", 95% confidence interval = (",
    signif(t0$conf.int[1],3),", ",signif(t0$conf.int[2],4),
    ") (after ",i," reps)\n", sep="")
cat(n0-isOK,"iterations 'failed' due to randomisation missing a category\n\n")
cat("Fraction of matches by category over ALL iterations:\n")
summCats <- data.frame(
  Rock_Type = row.names(matchByClass),
  Total_Matched = rowSums(matchByClass),
  Actual = cc0,
  Percent_Matched = paste0(round(100*(rowSums(matchByClass)/cc0),1),"%"),
  row.names = NULL)
print(summCats)
# tidy up
rm(list = c("n0","ftrain","i","isOK","jS","jM","k","m0","m1","t0","cc0",
            "matchByClass","fMatchXClass","results","train","summCats"))


# The number of iterations makes some difference!
#
# So, it looks like we should run at least 50-100 iterations to get a reasonable
# idea of how well our LDA model performs. More iterations is better, but it
# depends how long you want the run-time to be!


# accuracy vs iterations

par(mar = c(3,3,1,1), mgp = c(1.5,0.2,0), tcl = 0.25, font.lab = 2)
plot(c(10,20,50,100,200,500,1000),
     c(0.5572,0.5986,0.5848,0.5847,0.5837,0.581,0.582), 
     pch=19, cex = 1.2, log = "x", xlim = c(7,1400), ylim = c(0.45,0.65),
     xlab = "Number of iterations", ylab = "Mean accuracy \u00B1 95% CI")
arrows(c(10,20,50,100,200,500,1000),
       c(0.49,0.58,0.572,0.574,0.577,0.576,0.579),
       c(10,20,50,100,200,500,1000),
       c(0.624,0.617,0.598,0.595,0.591,0.586,0.585),
       angle = 90, length = 0.1, code = 3)


# We can run a similar validation process for the ALR-transformed data (we don't
# show the code, as it's effectively identical to the code for validation of LDA
# for closed data, but with references to 'Hallberg' replaced with
# 'Hallberg_alr').


# LDA train and predict open ####

n0 <- 100 # number of iterations
ftrain <- 0.5 # proportion of observations in training set
results <- data.frame(
  Rep = rep(NA, n0),
  matches = rep(NA, n0),
  non_matches = rep(NA, n0),
  success = rep(NA, n0))
train <- sample(1:NROW(Hallberg_alr), round(NROW(Hallberg_alr) * ftrain,0))
# make vector of individual category non-matches
matchByClass <- 
  data.frame(Match1 = rep(0,nlevels(Hallberg_alr$Rock[train]))) 
rownames(matchByClass) <- levels(ldaPred_rock_open$class)
fMatchXClass <- 
  data.frame(PctMch1 = rep(0,nlevels(Hallberg_alr$Rock[train]))) 
rownames(fMatchXClass) <- levels(ldaPred_rock_open$class)
# make vector of cumulative category counts in Hallberg_alr[-train] iterations
cc0 <- rep(0,nlevels(Hallberg_alr$Rock)) 
isOK <- 0 ; i <- 2

for (i in 1:n0) {
  # train <- sample(1:NROW(Hallberg_alr), round(NROW(Hallberg_alr)-5,0))
  train <- sample(1:NROW(Hallberg_alr), round(NROW(Hallberg_alr) * ftrain,0))
  if (is.na(match(NA,tapply(Hallberg_alr[train,]$SiO2, 
                            Hallberg_alr[train,]$Rock, sd, na.rm=T))) == TRUE) {
    lda_Rock_train <- lda(formula = Rock ~ SiO2 + TiO2 + Fe2O3 + FeO + 
                            MnO + MgO + CaO + Na2O + K2O + P2O5, 
                          data = Hallberg_alr[train,],
                          prior=as.numeric(summary(Hallberg_alr$Rock[train]))/
                            nrow(Hallberg_alr[train,]))
    ldaPred_rock_open <- predict(lda_Rock_train, Hallberg_alr[-train,])
    isOK <- isOK + 1
  }
  
  k=0 # number of matches
  m0 <- # vector of individual category matches 
    as.matrix(rep(0,nlevels(Hallberg_alr$Rock[train]))) 
  rownames(m0) <- levels(Hallberg_alr$Rock)
  m1 <- # vector of fractional category matches 
    as.matrix(rep(0,nlevels(Hallberg_alr$Rock[train]))) 
  rownames(m1) <- levels(Hallberg_alr$Rock)
  for (jM in 1:NROW(Hallberg_alr[-train,])) {
    for (jS in 1:nlevels(ldaPred_rock_open$class)) {
      if((ldaPred_rock_open$class[jM] == levels(ldaPred_rock_open$class)[jS]) & 
         (Hallberg_alr$Rock[-train][jM] == levels(ldaPred_rock_open$class)[jS]) ) 
        m0[jS] = m0[jS] + 1
      else  m0[jS] = m0[jS] 
    }
    k = sum(m0)
  }
  cc0 <- cc0 + as.numeric(summary(Hallberg_alr$Rock[-train]))
  m1 <- round(100*m0/as.numeric(summary(Hallberg_alr$Rock[-train])),1)
  matchByClass[,paste0("Match",i)] <- m0
  fMatchXClass[,paste0("PctMch",i)] <- m1
  # output to results data frame: iteration, matches, non-matches, proportion matched
  results[i,] <- c(i, k, NROW(Hallberg_alr[-train,])-k, 
                   signif(k/NROW(Hallberg_alr[-train,]),3))
}

cat(paste("[Based on", n0, "random subsets of",paste0(100*ftrain,"%"),
          "of the dataset to train LDA model\n",
          "      to predict remaining observations]\n"))
cat("Number of obs. in random subsets =",NROW(train),
    " (predicting",NROW(Hallberg_alr)-NROW(train),"samples)\n")
print(numSummary(results[,2:4], statistics=c("mean","sd"))$table)
ns0 <- numSummary(results$success)
t0 <- t.test(results$success)
cat(rep("-\u2013-",24),
    "\nStat. summary for 'success':\nMean = ",round(ns0$table[1],4),
    ", sd = ",round(ns0$table[2],4),
    ", 95% confidence interval = (",
    signif(t0$conf.int[1],3),", ",signif(t0$conf.int[2],4),
    ") (after ",i," reps)\n", sep="")
cat(n0-isOK,"iterations 'failed' due to randomisation missing a category\n\n")
cat("Fraction of matches by category over ALL iterations:\n")
summCats <- data.frame(
  Rock_Type = row.names(matchByClass),
  Total_Matched = rowSums(matchByClass),
  Actual = cc0,
  Percent_Matched = paste0(round(100*(rowSums(matchByClass)/cc0),1),"%"),
  row.names = NULL)
print(summCats)
rm(list = c("n0","ftrain","i","isOK","jS","jM","k","m0","m1","t0","cc0",
            "matchByClass","fMatchXClass","results","train","summCats"))


# If we compare the validation output for the open data with the results for
# closed data, we see very little difference in overall accuracy. Both closed
# and open data generated successful predictions 62.7% of the time, with a
# slightly lower standard deviation for LDA based on open data.
#
# There were small differences (in this validation exercise), between closed and
# open data, in the ability of LDA to predict specific categories. Closed-data
# LDA seemed better at predicting Dolerite, Gabbro, High Mg Basalt, and
# Metasediment. Conversely, LDA using open data seemed better at predicting
# Basalt, Peridotite, Pyroxenite, and Spinel Peridotite.
#
# We may get better prediction accuracy by including different variables in the
# LDA model. The original Hallberg geochemistry data at
# https://catalogue.data.wa.gov.au/dataset/hallberg-geochemistry also contain
# trace element concentrations, which can sometimes provide better
# discrimination than major element concentrations alone. This would be an
# interesting exercise to extend one's skills in data curation and multivariate
# analysis of compositional data.
#
# An issue that we haven't considered yet is whether all the variables we used
# for prediction are necessary. The R package '**klaR**' (Weihs et al. 2005)
# includes the *stepclass()* function, which enables us to refine an LDA model,
# using a procedure similar to stepwise selection of predictors in multiple
# regression.

# remove temp LDA objects
rm(list = c("n0","ftrain","i","e0","jS","jM","k","m0",
            "matchByClass","results","train","loocv"))

# ├─┴┬┴┬┴┬┴┬┴┬┴┬┴┬┴┬┴┤ ####
# ─────────────────────────────────────────────────────────────────

# LDA-PCA -- Linear Discriminant Analysis using Principal Components ####

# We are using a curated version of a whole rock composition dataset from the
# US Geological Survey (<https://mrdata.usgs.gov/ngdb/rock/>). This mostly
# contains data for rock samples from the USA, with a small proportion of data
# from other geographical regions

gitpath<-"https://raw.githubusercontent.com/Ratey-AtUWA/compositional_data/main/"
usgs <- read.csv(paste0(gitpath,"usgs_ngdb_trimmed.csv"), 
                 stringsAsFactors = TRUE)
usgs_ppm <- usgs
usgs_ppm[,9:20] <- t(apply(usgs_ppm[,c(9:20)], MARGIN = 1, 
                           FUN = function(x){x * 1e4}))
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

# First do PCA -- JUST ON CLR-transformed ####

pca_usgs <- prcomp(usgs_clr[,9:33], scale. = TRUE)

### Visualize PCA scree plot (eigenvalues)

par(mar = c(3.,3.5,1.2,1), mfrow = c(1,1), mgp=c(2,0.8,0), font.lab=2)
plot(pca_usgs, main = "", cex.main = 1, col = 4)
axis(1, at=seq(0.8,11.4,l=10), labels=seq(1,10,1), tcl=0.01, col=15, mgp=c(0.5,0.1,0))
mtext("Component", side = 1, line = 1.15, font = 2)
abline(h = 1, col = 11, lty = 2)

# The scree plot for this PCA shows 7 components with
# variances (eigenvalues) greater than 1.

# numerical PCA summary
{pca_usgs$rot[,1:5]
cat("...\n\nComponent Variances - CLR-transformed (open) data\n")
print(pca_usgs$sdev^2)
cat("___\n\nProportions of variance explained by each component",
    "\nCLR-transformed (open) data\n")
print(round(pca_usgs$sdev^2/sum(pca_usgs$sdev^2),3))
cat("___\n\nCumulative proportions of variance explained by each component",
    "\nCLR-transformed (open) data\n")
print(cumsum(round(pca_usgs$sdev^2/sum(pca_usgs$sdev^2),3)))}

# The cumulative variance explained for this PCA shows 8 components before
# cumulative proportion of variance explained > 0.8.

# visualise PC2 vs. PC1
usgsBP <- fviz_pca_biplot(pca_usgs, title = "", geom="point",
                          col.ind = usgs_clr$Lith, 
                          col.var = "black", repel = TRUE,
                          pch = c(19,17,19,17)[usgs_clr$Lith],
                          palette = c("dodgerblue4","sienna","red3","purple"),
                          pointsize = 2.5, alpha.ind = 0.35)
usgsBP + theme_classic()

### Add principal components to data frame

# We will include all principal component scores for each sample in our
# CLR-transformed data frame (code not shown), even though the Kaiser rules
# suggest that at most 7 components contain useful information. We don't have to
# subsequently use them *all* !. [NOTE that we could also perform LDA using
# the information stored in the PCA output object itself.]

# add PCs 1 to 8 to dataframe
usgs_clr <- cbind(usgs_clr, pca_usgs$x[,1:8])
names(usgs_clr)

## LDA on PCA ordination scores of whole rock major element data

# In this exercise we use the first 8 Principal Components as predictor
# variables for the Linear Discriminant Analysis. This is somewhat arbitrary,
# but it's based on the 7-8 suggested by the Kaiser rules.

# LDA-PCA whole rock data ####

# check collinerality of PCs 1-8
round(cor(data0[,3:10]),2)

data0 <- usgs_clr[,c("Rock","Lith",
                     "PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8")]
data0[,2:10] <- scale(data0[,2:10]) # scale just numeric variables

# NOTE - predicting 'Lith' category in this example (not Rock)

lda_pca_usgs <- lda(formula = Lith ~ PC1 + PC2 + PC3 + PC4 + PC5 + 
                       PC6 + PC7 + PC8, data = data0,
                    prior = as.numeric(summary(data0$Lith))/nrow(data0)) 
print(lda_pca_usgs$call) # custom output is tidier for pdf
cat("\nPrior probablilities:\n");print(lda_pca_usgs$prior, digits=2)
cat("\nGroup means:\n");print(lda_pca_usgs$means, digits=2)
{cat("\nCoefficients of linear discriminants:\n")
print(lda_pca_usgs$scaling, digits=3)
cat("\nProportions of between groups variance:\n")
cat(paste0("   LD",seq(1,3)),"\n")
cat(round(lda_pca_usgs$svd^2/sum(lda_pca_usgs$svd^2),4),"\n")}
# histogram showing separation in LD1
pred_lda_pca <- predict(lda_pca_usgs, usgs_clr[,35:42])
ldahist(pred_lda_pca$x[,1], g = usgs_clr$Lith, type = "both")

# pseudo-biplots for LDA on PCs ####
par(mfrow = c(2,2), mar = c(3.5,3.5,1,1), mgp = c(1.5,0.3,0), tcl = 0.25,
    lend = "square", ljoin = "mitre", cex.main = 0.9, font.lab=2)
palette(c("black","#104E8B80", "#A0522D80", "#CD000080", "#A020F080","white"))
plot(lda_pca_usgs$scaling[,1], lda_pca_usgs$scaling[,2], 
     type="n", xlim = c(-1,1),  ylim = c(-0.5,0.5), 
     xlab="Linear Discriminant [1]", ylab="Linear Discriminant [2]", 
     main="Variable Coefficients [LD1, LD2]")
abline(v=0,col="grey",lty=2); abline(h=0,col="grey",lty=2)
for(i in 1:NROW(lda_pca_usgs$scaling)){
  arrows(0,0,lda_pca_usgs$scaling[i,1],lda_pca_usgs$scaling[i,2],
         length = 0.08, col = "grey60") }
text(lda_pca_usgs$scaling[,1]*1.1, lda_pca_usgs$scaling[,2]*1.1, 
     labels=row.names(lda_pca_usgs$scaling[,1:2]),
     cex = 0.9, col = 1, pos=c(1,2,3,4), offset = 0.01)
mtext("(a)", 3, -1.5, adj = 0.05, cex = 1.2, font = 2)
clrz <- matrix(rainbow(20,v=0.8,end=0.8),ncol=10,byrow = T)
attributes(clrz) <- NULL

ldapcaPred_usgs <- predict(lda_pca_usgs)

plot(ldapcaPred_usgs$x[,1], ldapcaPred_usgs$x[,2], 
     cex = 1., lwd=2, col = c(2:5)[usgs_clr$Lith],
     pch = c(15,16,17,18)[usgs_clr$Lith], 
     main = "Predictions for Observations [LD1, LD2] based on PCA", 
     xlab = "Linear Discriminant [1]", ylab = "Linear Discriminant [2]")
abline(v=0,col="grey",lty=2); abline(h=0,col="grey",lty=2)
mtext("(b)", 3, -1.5, adj = 0.05, cex = 1.2, font = 2)

plot(ldapcaPred_usgs$x[,1], ldapcaPred_usgs$x[,3], 
     cex = 1., lwd=2, col = c(2:5)[usgs_clr$Lith],
     pch = c(15,16,17,18)[usgs_clr$Lith], 
     main="Predictions for Observations [LD1, LD3] based on PCA", 
     xlab="Linear Discriminant [1]", ylab="Linear Discriminant [3]")
sf <- 5; abline(v=0,col="grey",lty=2); abline(h=0,col="grey",lty=2)
for(i in 1:NROW(lda_pca_usgs$scaling)){
  arrows(0,0,lda_pca_usgs$scaling[i,1]*sf,lda_pca_usgs$scaling[i,3]*sf,
         length = 0.1, col = 1, lwd = 2) }
shadowtext(lda_pca_usgs$scaling[,1]*sf*1.1, lda_pca_usgs$scaling[,3]*sf*1.1, 
           labels=row.names(lda_pca_usgs$scaling),
           cex = 1.2, col = 1, bg = 6, pos=c(1,2,3,4), offset = 0.01)
mtext("(c)", 3, -1.5, adj = 0.05, cex = 1.2, font = 2)

plot(c(0,1),c(0,1), type="n", bty="n", ann = F, 
     xaxt="n", yaxt="n", xlab="", ylab="")
legend("left", ncol = 1, legend=levels(usgs_clr$Lith), 
       col=c(2:5), pch=c(15,16,17,18), pt.lwd = 2,
       bty="n", box.col="grey90", y.intersp = 1, 
       title="Rock Type in (b) & (c)",
       box.lwd=2, inset=0.02, pt.cex=2.5, cex=2)

# LDA inbuilt cross validation

lda_pca_usgs_cv <- lda(formula = Rock ~ PC1 + PC2 + PC3 + PC4 + PC5 + 
                         PC6 + PC7 + PC8, data = data0,
                       prior = as.numeric(summary(data0$Rock))/nrow(data0),
                       CV = TRUE) 
cat("\nLDA predictions (first several):\n"); head(lda_pca_usgs_cv$class, n=27)
cat("\nPosterior probabilities (first few rows):\n"); head(lda_pca_usgs_cv$posterior)
ldapcaComp <- data.frame(Observed=usgs_clr$Rock,Predicted=lda_pca_usgs_cv$class)
cat("\nComparison of (first several) LDA predictions with reality:\n")
head(ldapcaComp, n = 15)
cat("\nProportion of correct classifications",
  length(which((ldapcaComp[,1]==ldapcaComp[,2])==TRUE))/length(ldapcaComp[,1]),
  "\n")

# So we can correctly classify 60.6% of broader rock types based on composition.

# Stepwise LDA ####

# add PCs 9 to 16 to dataframe
usgs_clr <- cbind(usgs_clr, pca_usgs$x[,9:16])
names(usgs_clr)

# see help(stepclass) for details

lda_pca_step <- stepclass(formula = Rock ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + 
              PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16,
                          data = usgs_clr, method = "lda", 
                          direction = "both", criterion = "CR",
                          improvement = 0.005, output = FALSE)
cat("Final",lda_pca_step$method,"model:",
    as.character(lda_pca_step$formula)[c(2,1,3)],"\n\n")
cat(lda_pca_step$per, signif(lda_pca_step$result.pm[1], 4), sep=":")


# The output from stepwise selection of LDA predictors shows us that we do
# not need to include all the principal components to achieve a
# correctness rate (62.9%) similar to that achieved above (62.7%).

# As is the case with stepwise selection of predictors for multiple linear
# regression, the choice of method options could affect our result. For example,
# you could try changing the *direction*, *criterion*, or *improvement* options
# (run help(package="klaR", topic="stepclass") for details of the choices).
 
# We would re-run our LDA model following stepwise selection, using the formula
# after 'final model:' in the output above.

# final check
data0 <- usgs_clr[,c("Rock","Lith",
                     "PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8",
                     "PC9","PC10","PC11","PC12")]
data0[,3:14] <- scale(data0[,3:14]) # scale just numeric variables

lda_pca_final <- lda(formula = Rock ~ PC1 + PC2 + PC3 + PC4 + PC5 + 
                       PC6 + PC8 + PC9 + PC12, 
                     data = data0,
                     prior = as.numeric(summary(data0$Rock))/
                       nrow(data0), CV = TRUE) 
ls(lda_rock_clos)
loocv <- round(lda_rock_clos$posterior,3)
colnames(loocv)<- substr(colnames(lda_rock_clos$posterior),1,7)
print(loocv)
clipr::write_clip(loocv) # then paste to Excel
# calc proportion correct
ldaComp <- data.frame(Observed=Hallberg$Rock,Predicted=lda_rock_clos$class)
cat("\nProportion of correct classifications",
    signif(100*length(which((ldaComp[,1]==ldaComp[,2])==TRUE))/
             length(ldaComp[,1]),3), "%\n")

# proportion of correct classifications from LOOCV is similar to 
# LDA on PCs 1-8.

# We could go for more detail by looking at which categories are predicted 
# most accurately...
#     (there was code for this somewhere)

# [end code]
