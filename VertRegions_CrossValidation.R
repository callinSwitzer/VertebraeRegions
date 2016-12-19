## Callin Switzer
## 16 Dec 2016
## Vertebrae regions, using linear splines or "knots"

#install packages
ipak <- function(pkg){
     new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
     if(length(new.pkg)) install.packages(new.pkg, dependencies = TRUE)
     sapply(pkg, require, character.only = TRUE)
}

packages <- c( "ggplot2", "reshape2", "MASS", "rms")
ipak(packages)


# load data
load("/Users/callinswitzer/Dropbox/dataAnalysisForOthers/VertebrateMorphology_Katrina/alligator.rda")
load("/Users/callinswitzer/Dropbox/dataAnalysisForOthers/VertebrateMorphology_Katrina/mus.rda")
alli <- as.data.frame(alligator)


pcs <- prcomp(alli[, 2:ncol(alli)], center = TRUE, scale. = TRUE)
summary(pcs)
biplot(pcs)
plot(pcs)
str(pcs)

plot(alli$Vertebra, pcs$x[,1])


hcl.flow = hclust(dist(alli[, 2:ncol(alli)], method = 'euclidian'), method = 'ward.D')
plot(hcl.flow)
groups <- cutree(hcl.flow, k=2) # cut tree into k clusters


plot(pcs$x[,1] ~ alli$Vertebra)
aa <- ols(pcs$x[,1]  ~ lsp(alli$Vertebra, c(9, 14), x=TRUE))
lines(predict(aa), x = alli$Vertebra)


abline(m1)
