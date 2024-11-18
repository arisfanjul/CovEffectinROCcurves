##%#########################################################################%##
#                                                                             #
####   Comparing ROC and AROC curves (same marker) - A practical example   ####
#                                                                             #
##%#########################################################################%##

#%########################################%#
## Authors: Arís Fanjul Hevia
## Purpose: Test whether the ROC and the AROC curves of the same diagnostic
##          marker are the same or not.
## Data: 11/2024
#%########################################%#

# Packages ------------------------------------------------------------------

# library(refreg)
library(viridis) # for color selection
library(hdrcde) # for estimating the conditional densities
library(ggplot2) #for the graphical representations
library(ggExtra) #for the graphical representations
library(kableExtra) # for presenting the tables

# Loading functions from the functions.R script
source("functions.R")

# Loading the dataset --------------------------------------------------------
data("aegis", package = "refreg") # instalation of package refreg needed

# We will use the following variables:
#. dm: Diabetes mellitus indicator (no, and yes). The healthy/diseased indicator.
#. fpg: the diagnostic marker (fasting plasma glucose levels in mg/dL).
#. age: the covariate that will be considered.

Y = aegis$fpg
YF = aegis$fpg[aegis$dm == "yes"]
YG = aegis$fpg[aegis$dm == "no"]
sampleY = list(YF,YG)

X = aegis$age
XF = aegis$age[aegis$dm == "yes"]
XG = aegis$age[aegis$dm == "no"]
sampleX = list(XF,XG)

# Brief exploratory analysis of the variables --------------------------------

# Sample size for the healthy and the diseased:
table(aegis$dm)

nameG="No-diabetes"
nameF="Diabetes"

# Diagnostic marker
nameY="fpg"
nameYG=paste(nameY,"(",nameG,")")
nameYF=paste(nameY,"(",nameF,")")

Y.sum=rbind(summary(Y), summary(YF), summary(YG))
rownames(Y.sum)=c(nameY,nameYF,nameYG)
kable(Y.sum)

# Covariate
nameX="age"
nameXG=paste(nameX,"(",nameG,")")
nameXF=paste(nameX,"(",nameF,")")

X.sum=rbind(summary(X), summary(XF), summary(XG))
rownames(X.sum)=c(nameX.1,nameXF,nameXG)
kable(X.sum)


# Scatterplot with marginal densities
q = ggplot(aegis) +
  geom_point(aes(x = age, y = fpg, color = dm), alpha = 0.6, shape = 16) +
  scale_color_manual(labels=c(nameG, nameF),
                     values = c("no" = viridis(5)[2], "yes"= viridis(4)[3]))+
  theme_minimal() +
  scale_fill_discrete(labels=c(nameG, nameF))+
  theme(legend.position = "bottom") + 
  labs(x = "age", y = "fpg") 
ggMarginal(q, type = "densigram", groupColour = TRUE, groupFill = TRUE, alpha = 0.3)


# Conditional densities F/G:
plot.densities.FG(XF,YF, XG,YG, x.name=nameX, y.name = nameY, 
                  col1 = viridis(5,alpha=0.4)[2],col2 = viridis(4,alpha=0.4)[3])



# Testing AROC vs ROC for the same marker ------------------------------------

set.seed(1234)
results.1 = Test.ROCAROC(sampleY, sampleX, nboots=200, nr = 2)
results.1$statistics
# T1          T2       TKS
# 1 0.04764318 0.003256029 0.1849691
results.1$pvalues
# pvalue.T1 pvalue.T2 pvalue.TKS
# 1      0.16     0.155      0.695
plot(results.1$p, results.1$ROC, type = "l", col = "orange")
points(results.1$p, results.1$AROC, type = "l", col = "brown")

set.seed(1234)
results.1.1000 = Test.ROCAROC(sampleY, sampleX, nboots=1000, nr = 2)
results.1.1000$statistics
# T1          T2       TKS
# 1 0.04764318 0.003256029 0.1849691
results.1.1000$pvalues
# pvalue.T1 pvalue.T2 pvalue.TKS
# 1     0.184     0.198      0.711

set.seed(1234)
results.2 = Test.ROCAROC(sampleY, sampleX, nboots=200, nr = 3)
results.2$statistics
# T1          T2       TKS
# 1 0.04009048 0.002517042 0.2833548
results.2$pvalues
# pvalue.T1 pvalue.T2 pvalue.TKS
# 1      0.27      0.29      0.115

set.seed(1234)
results.3 = Test.ROCAROC(sampleY, sampleX, nboots=200, nr = 4)
results.3$statistics
# T1           T2       TKS
# 1 0.01677266 0.0008004622 0.2972644
results.3$pvalues
#   pvalue.T1 pvalue.T2 pvalue.TKS
# 1      0.97     0.925      0.245

plot(results.1$p,results.1$ROC,type="l",col = "orange",lwd=2,lty=1,xlab="p",ylab="ROC(p)",xlim=c(0,1),ylim=c(0,1),axes = FALSE, asp = 1)
axis(1)
axis(2)
points(results.1$p,results.1$AROC,type="l", col="darkred",lwd=2, lty=1)
legend("bottomright",legend=c("ROC (pooled)", "AROC"),lwd=2,lty=1, col=c("orange", "darkred"),bty = "n")


plot(results.1$p,results.1$ROC,type="l",col = "orange",lwd=2,lty=1,xlab="p",ylab="ROC(p)",xlim=c(0,1),ylim=c(0,1),axes = FALSE, asp = 1)
axis(1)
axis(2)
points(results.2$p, results.2$ROC, type = "l", col = "orange", lty = 2, lwd = 2)
points(results.3$p, results.3$ROC, type = "l", col = "orange", lty = 3, lwd = 2)
points(results.1$p, results.1$AROC, type = "l", col = "brown", lwd = 2)
points(results.2$p, results.2$AROC, type = "l", col = "brown", lty = 2, lwd = 2)
points(results.3$p, results.3$AROC, type = "l", col = "brown", lty = 3, lwd = 2)
points(results.1$p,results.1$AROC,type="l", col="darkred",lwd=2, lty=1)

legend("bottomright",legend=c("ROC (pooled)", "AROC"),lwd=2,lty=1, col=c("orange", "darkred"),bty = "n")



