## Callin Switzer
## 16 Dec 2016
## Vertebrae regions, using linear splines or "knots"

#install packages
ipak <- function(pkg){
     new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
     if(length(new.pkg)) install.packages(new.pkg, dependencies = TRUE)
     sapply(pkg, require, character.only = TRUE)
}

packages <- c( "ggplot2", "reshape2", "MASS", "rms", "viridis")
ipak(packages)


# load data
load("/Users/callinswitzer/Dropbox/dataAnalysisForOthers/VertebrateMorphology_Katrina/alligator.rda")
load("/Users/callinswitzer/Dropbox/dataAnalysisForOthers/VertebrateMorphology_Katrina/mus.rda")
alli <- as.data.frame(mus)


pcs <- prcomp(alli[, 2:ncol(alli)], center = TRUE, scale. = TRUE)
summary(pcs)
biplot(pcs)
plot(pcs)
str(pcs)

plot(alli[,1], pcs$x[,1])


hcl.flow = hclust(dist(alli[, 2:ncol(alli)], method = 'euclidian'), method = 'ward.D')
plot(hcl.flow)
groups <- cutree(hcl.flow, k=2) # cut tree into k clusters


plot(pcs$x[,1] ~ alli[,1])
aa <- ols(pcs$x[,1]  ~ lsp(alli[,1], c(9, 14), x=TRUE))
lines(predict(aa), x = alli[,1])

threeBreak <- function(breakpoints = c(6, 11)){
     # plot the pc1 vs. vertebrae number
     plot(pcs$x[,1]  ~ alli[,1])
     
     # make the models for each of the three regions
     m1 <- lm(pcs$x[1:breakpoints[1],1]  ~ alli[,1][1:breakpoints[1]])
     lines(predict(m1), x = alli[,1][1:breakpoints[1]], col = 'red')
     
     m2 <- lm(pcs$x[(breakpoints[1] + 1):breakpoints[2],1]  ~ alli[,1][(breakpoints[1] + 1):breakpoints[2]])
     lines(predict(m2), x = alli[,1][(breakpoints[1] + 1):breakpoints[2]], col = 'blue')
     
     m3 <- lm(pcs$x[(breakpoints[2] + 1):length(alli[,1]),1]  ~ alli[,1][(breakpoints[2] + 1):length(alli[,1])])
     lines(predict(m3), x = alli[,1][(breakpoints[2] + 1):length(alli[,1])], col = 'green')
     
     
}


# generate random cross-validation groups for each of the vertebrae
alligator <- cbind(alligator,  rep(sample(1:10), length= 22))
colnames(alligator)[ncol(alligator)] <- "cvGroup"

RegressionWBreaks <- function(breakpoints = c(6, 15), animal = alligator, cvGp = 1){
     # compute pca
     pcs <- prcomp(animal[, 2:ncol(animal)], center = TRUE, scale. = TRUE)
     
     dataset1 <- data.frame(vertNum = animal[,1], pc1 = pcs$x[,1], cvGp = animal[, 'cvGroup'])
     dataset_full <- dataset1
     
     # plot the pc1 vs. vertebrae number
     plot(dataset1[, 1:2])
     
     # remove points to later calculate out-of-sample error
     outOfSamplePoints <- data.frame(dataset1[dataset1$cvGp == cvGp,])
     dataset1[dataset1$cvGp == cvGp,1:2] <- NA
     
     # color points that are being used for calculating lines
     points(dataset1[, 1:2],col = 'grey40', pch = 20)
     points(outOfSamplePoints[,1:2], pch = 5)
     legend("bottomright", legend = c("in sample", "out of sample"), pch = c(20, 5), col = c('grey40', 'black'))
     
     # define color pallette for lines
     cols = viridis::inferno(n  = length(breakpoints) + 1, end = 0.8)
     
     # make an empty list for adding the models
     mod <- list()
     predictedVals <- numeric()
     for(ii in 1:(length(breakpoints) + 1)){
          if(ii == 1){
               d2 <- dataset1[1:breakpoints[1], ]
          }
          else if(ii <= length(breakpoints)){
               d2 <- dataset1[(breakpoints[ii - 1] + 1):breakpoints[ii], ]
          }
          else if(ii > length(breakpoints)){
               d2 <- dataset1[(breakpoints[ii-1] + 1):length(animal[,1]),, ]
          }
          else{
               stop("ERROR")
          }
          mod[[ii]] <- lm(pc1 ~ vertNum, data = d2)
          predVals <- predict(mod[[ii]], newdata = data.frame(vertNum = dataset_full[rownames(d2),1]))
          lines(y = predVals, x = dataset_full[rownames(d2),1], col = cols[ii], lwd = 5)
          predictedVals <- c(predictedVals, predVals) # get predicted values
     }
     # get out of sample errors
     newDF <- data.frame(predv =  predictedVals, dataset_full)
     outOfsamp_preds <- newDF[newDF$cvGp == cvGp, ]
     ers <- with(outOfsamp_preds, pc1 - predv )
     # plot out-of-sample errors
     segments(x0 = outOfsamp_preds$vertNum, y0 = outOfsamp_preds$pc1, x1 = outOfsamp_preds$vertNum,
              y1 = outOfsamp_preds$predv, lwd = 3, col= 'grey')
     abline(v = dataset_full$vertNum[breakpoints] + 0.5)
     
     return(list(outOfsampleVertNums = outOfsamp_preds$vertNum, predictionErrors = ers, mean_abs_err = mean(abs(ers)), numBreaks= length(breakpoints), bkpts = breakpoints, crossValGroup = cvGp))
     
}


RegressionWBreaks(breakpoints = c(6, 12), animal = alligator, cvGp = 1)
nrow(mus)


for(ii in 1:10){
     png(paste("~/Desktop/testStack/", formatC(ii, width = 4, flag = 0), ".png", sep = ""))
     aa <- RegressionWBreaks(breakpoints = c(6, 12), animal = alligator, cvGp = ii)
     dev.off()
     print(aa$mean_abs_err)
}



# generate all possible combinations of one two three breakpoints
animal = alligator
Onebp <- 1:(nrow(animal) - 1)
twobp <- data.frame(bp1 = rep(1:length(animal), each = nrow(animal)), bp2 = 1:nrow(animal))
twobp[twobp[,1] >= twobp[,2], ] <- c(NA, NA)
twobp <- twobp[complete.cases(twobp), ]
twobp <- twobp[twobp[,1] < nrow(animal) & twobp[,2] < nrow(animal), ]

threebp <- data.frame(bp1 = rep(1:nrow(animal), each = nrow(animal)^2), 
                      bp2 = rep(1:nrow(animal), each = nrow(animal)), 
                      bp3 = 1:nrow(animal))

threebp <- threebp[threebp[,1] < threebp[, 2] & threebp[,2] < threebp[,3], ]

rwb_vect <- Vectorize(RegressionWBreaks, vectorize.args = 'breakpoints')

RegressionWBreaks(breakpoints = c(6, 16), animal = alligator, cvGp = 1)
nrow(mus)


apply(data.frame(Onebp),MARGIN = 1,  FUN = function(x) RegressionWBreaks(x, animal = alligator))
apply(twobp,MARGIN = 1,  FUN = function(x) RegressionWBreaks(x, alligator))