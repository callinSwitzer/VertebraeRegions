## Callin Switzer
## 16 Dec 2016
## Update 12 Jan
## Vertebrae regions, using linear splines


#############################################################################
#setup
#############################################################################
#install packages
ipak <- function(pkg){
     new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
     if(length(new.pkg)) install.packages(new.pkg, dependencies = TRUE)
     sapply(pkg, require, character.only = TRUE)
}

packages <- c( "ggplot2", "reshape2", "MASS", "rms", "viridis", "rgl")
ipak(packages)


# load data
load("/Users/callinswitzer/Dropbox/dataAnalysisForOthers/VertebrateMorphology_Katrina/alligator.rda")
load("/Users/callinswitzer/Dropbox/dataAnalysisForOthers/VertebrateMorphology_Katrina/mus.rda")

#############################################################################
# initial data inspection
#############################################################################
alli <- as.data.frame(mus)


pcs <- prcomp(alli[, 2:ncol(alli)], center = TRUE, scale. = TRUE)
summary(pcs)
biplot(pcs)
plot(pcs)
str(pcs)

plot(alli[,1], pcs$x[,1])

# compare hierarchical clustering method
hcl.flow = hclust(dist(alli[, 2:ncol(alli)], method = 'euclidian'), method = 'ward.D')
plot(hcl.flow)
groups <- cutree(hcl.flow, k=2) # cut tree into k clusters


# example of splines with that meet at the breakpoints
plot(pcs$x[,1] ~ alli[,1])
aa <- ols(pcs$x[,1]  ~ lsp(alli[,1], c(9, 14), x=TRUE))
lines(predict(aa), x = alli[,1])


# generate random cross-validation groups for each of the vertebrae
alligator <- cbind(alligator,  rep(sample(1:10), length= 22))
mus <- cbind(mus,  rep(sample(1:10), length= length(mus[,1])))

colnames(alligator)[ncol(alligator)] <- "cvGroup"
colnames(mus)[ncol(mus)] <- "cvGroup"
# alligator[, "cvGroup"] <- sample(alligator[, "cvGroup"]) # reshuffle cvGroup


#############################################################################
# function for conducting regression with different breakpoints
#############################################################################
RegressionWBreaks <- function(breakpoints = c(6, 15), animal = alligator, cvGp = 1, ylim = c(-5, 3)){
     # compute pca
     pcs <- prcomp(animal[, 2:(ncol(animal) -1)], center = TRUE, scale. = TRUE)
     
     dataset1 <- data.frame(vertNum = animal[,1], pc1 = pcs$x[,1], cvGp = animal[, 'cvGroup'])
     dataset_full <- dataset1
     
     # plot the pc1 vs. vertebrae number
     plot(dataset1[, 1:2], ylim = ylim)
     
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
          
          OK <- tryCatch({mod[[ii]] <- lm(pc1 ~ vertNum, data = d2)
          predVals <- predict(mod[[ii]], newdata = data.frame(vertNum = dataset_full[rownames(d2),1]))
          lines(y = predVals, x = dataset_full[rownames(d2),1], col = cols[ii], lwd = 5)
          predictedVals <- c(predictedVals, predVals) # get predicted values, but will give 
          # a warning if there is only one vertebra in a group -- rank deficient.
          }, 
          error = function(e) NA)
          if(is.na(OK[1])) return(NA)
     }
     # get out of sample errors
     if(length(predictedVals) > nrow(dataset_full)){
          print("trying to make a break outside of range of vertebrae")
          return(NA)
     }
     else{
          newDF <- data.frame(predv =  predictedVals, dataset_full)
          outOfsamp_preds <- newDF[newDF$cvGp == cvGp, ]
          ers <- with(outOfsamp_preds, pc1 - predv )
          # plot out-of-sample errors
          segments(x0 = outOfsamp_preds$vertNum, y0 = outOfsamp_preds$pc1, x1 = outOfsamp_preds$vertNum,
                   y1 = outOfsamp_preds$predv, lwd = 3, col= 'grey')
          abline(v = dataset_full$vertNum[breakpoints] + 0.5)
          
          return(list(outOfsampleVertNums = outOfsamp_preds$vertNum, 
                      predictionErrors = ers, mean_abs_err = mean(abs(ers)), 
                      numBreaks= length(breakpoints), bkpts = breakpoints, crossValGroup = cvGp))
          }
}


#############################################################################
# function for conducting regression with different breakpoints, without plotting
#############################################################################
RegressionWBreaks_noPlot <- function(breakpoints = c(6, 15), animal = alligator, cvGp = 1, dataset1 = dataset1){
     dataset_full <- dataset1
     
     # plot the pc1 vs. vertebrae number
     # plot(dataset1[, 1:2])
     # 
     # remove points to later calculate out-of-sample error
     outOfSamplePoints <- data.frame(dataset1[dataset1$cvGp == cvGp,])
     dataset1[dataset1$cvGp == cvGp,1:2] <- NA
     
     # color points that are being used for calculating lines
     # points(dataset1[, 1:2],col = 'grey40', pch = 20)
     # points(outOfSamplePoints[,1:2], pch = 5)
     # legend("bottomright", legend = c("in sample", "out of sample"), pch = c(20, 5), col = c('grey40', 'black'))
     
     # define color pallette for lines
     # cols = viridis::inferno(n  = length(breakpoints) + 1, end = 0.8)
     
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
          
          OK <- tryCatch({mod[[ii]] <- lm(pc1 ~ vertNum, data = d2)
          predVals <- predict(mod[[ii]], newdata = data.frame(vertNum = dataset_full[rownames(d2),1]))
          # lines(y = predVals, x = dataset_full[rownames(d2),1], col = cols[ii], lwd = 5)
          predictedVals <- c(predictedVals, predVals) # get predicted values, but will give 
          # a warning if there is only one vertebra in a group -- rank deficient.
          }, 
          error = function(e) NA)
          if(is.na(OK[1])) return(NA)
     }
     # get out of sample errors
     if(length(predictedVals) > nrow(dataset_full)){
          print("trying to make a break outside of range of vertebrae")
          return(NA)
     }
     else{
          newDF <- data.frame(predv =  predictedVals, dataset_full)
          outOfsamp_preds <- newDF[newDF$cvGp == cvGp, ]
          ers <- with(outOfsamp_preds, pc1 - predv )
          # plot out-of-sample errors
          # segments(x0 = outOfsamp_preds$vertNum, y0 = outOfsamp_preds$pc1, x1 = outOfsamp_preds$vertNum,
          #          y1 = outOfsamp_preds$predv, lwd = 3, col= 'grey')
          # abline(v = dataset_full$vertNum[breakpoints] + 0.5)
          
          return(list(outOfsampleVertNums = outOfsamp_preds$vertNum, 
                      predictionErrors = ers, mean_abs_err = mean(abs(ers)), 
                      numBreaks= length(breakpoints), bkpts = breakpoints, crossValGroup = cvGp))
     }
}



# example of using this function
RegressionWBreaks(breakpoints = c(2,3,4), animal = alligator, cvGp = 4)
RegressionWBreaks_noPlot(breakpoints = c(2,7,14), animal = alligator, cvGp = 4)

##################################################################
# generate all possible combinations of one, two, three, or four breakpoints
# but exclude models that have a single vertebra as a "section"
# I can't think of a good way to use CV with a single vertebra
##################################################################

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
diffs <- t(apply(threebp, MARGIN = 1, function(x) diff(as.numeric(x))))
discardfbp <- apply(diffs, MARGIN = 1, function(x) 1 %in% x)
threebp <- threebp[!discardfbp, ]



fourbp <- data.frame(bp1 = rep(1:nrow(animal), each = nrow(animal)^3), 
                      bp2 = rep(1:nrow(animal), each = nrow(animal)^2), 
                      bp3 = rep(1:nrow(animal), each = nrow(animal)), 
                      bp4 = 1:nrow(animal))
fourbp <- fourbp[fourbp[,1] < fourbp[, 2] & fourbp[,2] < fourbp[,3] &  fourbp[,3] < fourbp[,4], ]
# remove single vertebra groups
diffs <- t(apply(fourbp, MARGIN = 1, function(x) diff(as.numeric(x))))
discardfbp <- apply(diffs, MARGIN = 1, function(x) 1 %in% x)
fourbp <- fourbp[!discardfbp, ]



fivebp <- data.frame(bp1 = rep(1:nrow(animal), each = nrow(animal)^4), 
                     bp2 = rep(1:nrow(animal), each = nrow(animal)^3), 
                     bp3 = rep(1:nrow(animal), each = nrow(animal)^2), 
                     bp4 = rep(1:nrow(animal), each = nrow(animal)),
                     bp5 = 1:nrow(animal))
fivebp <- fivebp[fivebp[,1] < fivebp[, 2] & fivebp[,2] < fivebp[,3] &  
                      fivebp[,3] < fivebp[,4] &  fivebp[,4] < fivebp[,5], ]


# remove single vertebra groups
diffs <- t(apply(fivebp, MARGIN = 1, function(x) diff(as.numeric(x))))
discardfbp <- apply(diffs, MARGIN = 1, function(x) 1 %in% x)
fivebp <- fivebp[!discardfbp, ]


## six bp's
sixbp <- data.frame(bp1 = rep(1:nrow(animal), each = nrow(animal)^5), 
                     bp2 = rep(1:nrow(animal), each = nrow(animal)^4), 
                     bp3 = rep(1:nrow(animal), each = nrow(animal)^3), 
                     bp4 = rep(1:nrow(animal), each = nrow(animal)^2),
                     bp5 = rep(1:nrow(animal), each = nrow(animal)),
                     bp6 = 1:nrow(animal))
sixbp <- sixbp[sixbp[,1] < sixbp[, 2] & sixbp[,2] < sixbp[,3] &  
                    sixbp[,3] < sixbp[,4] &  sixbp[,4] < sixbp[,5] &
                sixbp[,5] < sixbp[,6], ]
# remove single vertebra groups
diffs <- t(apply(sixbp, MARGIN = 1, function(x) diff(as.numeric(x))))
discardfbp <- apply(diffs, MARGIN = 1, function(x) 1 %in% x)
sixbp <- sixbp[!discardfbp, ]

#############################################################################
# LOOK AT ALL INSTANCES OF A SINGLE BREAKPOINT
#############################################################################
erDF <- data.frame()
for(kk in 1:length(Onebp)){
     cvError <- numeric()
     for(ii in 1:10){
          # png(paste("~/Desktop/testStack/", formatC(as.numeric(paste0(kk, ii)), width = 4, flag = 0), ".png", sep = ""))
          aa <- RegressionWBreaks(breakpoints = Onebp[kk], animal = alligator, cvGp = ii)
          # dev.off()
          if(!is.na(aa[1])) cvError[ii] <- aa$mean_abs_err
          else cvError[ii] <- NA
          
     }
     print(kk)
     
    paste('error = ', round(mean(cvError, na.rm = TRUE), 4), "breakpt = ", Onebp[kk])
    erDF <- rbind(erDF, c( Onebp[kk], round(mean(cvError, na.rm = TRUE), 4)))
}
colnames(erDF) <- c("breakpt","meanAbsErr")
erDF 

# visualize mean absolute error
plot(erDF)

#############################################################################
## DO THE SAME FOR ALL INSTANCES OF TWO BREAKPOINTS
#############################################################################
animal = alligator
pcs <- prcomp(animal[, 2:(ncol(animal) -1)], center = TRUE, scale. = TRUE)
dataset1 <- data.frame(vertNum = animal[,1], pc1 = pcs$x[,1], cvGp = animal[, 'cvGroup'])
erDF2 <- matrix(ncol = 3, nrow = 0)
system.time({ # about the same speed as a for loop
     for(kk in 1:nrow(twobp)){
          cvError <- sapply(X =  1:10, FUN = function(x){
               aa <- RegressionWBreaks_noPlot(breakpoints = as.numeric(twobp[kk,]), 
                                              animal = alligator, dataset1 = dataset1,  cvGp = x)
               #dev.off()
               if(!is.na(aa[1])) cvError <- aa$mean_abs_err
               else cvError <- NA 
               return(cvError)
          })
          erDF2 <- rbind(erDF2, c( twobp[kk,], round(mean(cvError, na.rm = TRUE), 4)))
          print(kk)
     }
})

colnames(erDF2) <- c("breakpt1", "breakpt2","meanAbsErr")
erDF2 <- data.frame(erDF2)

# interactive 3D visualization
# plot3d(erDF2, type = 'p')


# make a wireframe for 3D visualization
# s = interp(unlist(erDF2$breakpt1), unlist(erDF2$breakpt2), unlist(erDF2$meanAbsErr))
# with(s,wireframe(z,row.values=x,col.values=y, 
#                  main = "Surface elevation data",
#                  drape = TRUE,
#                  colorkey = TRUE,
#                  screen = list(z = 80, x = -80)))


# find min values
erDF2 <- as.data.frame(apply(erDF2, 2, unlist)) # break at 10, error = 0.5402
erDF2[erDF2$meanAbsErr == min(erDF2$meanAbsErr),] # breaks at 5 and 11, error = 0.2154

#############################################################################
## DO THE SAME FOR ALL INSTANCES OF THREE BREAKPOINTS
#############################################################################
erDF3 <- matrix(ncol = 4, nrow = 0)
system.time({
     for(kk in 1:nrow(threebp)){
          cvError <- numeric()
          alligator[,'cvGroup'] <- sample(alligator[,'cvGroup']) # resample cv group each time
          for(ii in 1:10){
               aa <- RegressionWBreaks_noPlot(breakpoints = as.numeric(threebp[kk,]), dataset1= dataset1,
                                              animal = alligator, cvGp = ii)
               if(!is.na(aa[1])) cvError[ii] <- aa$mean_abs_err
               else cvError[ii] <- NA
          }
          
          erDF3 <- rbind(erDF3, c( threebp[kk,], round(mean(cvError, na.rm = TRUE), 4)))
          print(paste(kk, "of", nrow(threebp)))
     } 
})

erDF3 <- as.data.frame(apply(erDF3, 2, unlist))
colnames(erDF3) <- c("breakpt1", "breakpt2", "breakpt3","meanAbsErr")
erDF3 <- data.frame(erDF3)
na.omit(erDF3[erDF3$meanAbsErr == min(erDF3$meanAbsErr, na.rm = TRUE),]) # 5, 11, and 13, error = 0.1795
RegressionWBreaks(breakpoints = c(5, 11, 13), animal = alligator, cvGp = 99)

#############################################################################
## DO THE SAME FOR ALL INSTANCES OF FOUR BREAKPOINTS
#############################################################################
erDF4 <- matrix(ncol = 5, nrow = 0)
system.time({
     for(kk in 1:nrow(fourbp)){
          cvError <- numeric()
          alligator[,'cvGroup'] <- sample(alligator[,'cvGroup']) # resample cv group each time
          for(ii in 1:10){
               aa <- RegressionWBreaks_noPlot(breakpoints = as.numeric(fourbp[kk,]), dataset1= dataset1,
                                              animal = alligator, cvGp = ii)
               if(!is.na(aa[1])) cvError[ii] <- aa$mean_abs_err
               else cvError[ii] <- NA
          }
          
          erDF4 <- rbind(erDF4, c( fourbp[kk,], round(mean(cvError, na.rm = TRUE), 4)))
          print(paste(kk, "of", nrow(fourbp)))
     } 
})
erDF4 <- as.data.frame(apply(erDF4, 2, unlist))
colnames(erDF4) <- c("breakpt1", "breakpt2", "breakpt3", "breakpt4","meanAbsErr")
erDF4 <- data.frame(erDF4)
na.omit(erDF4[erDF4$meanAbsErr == min(erDF4$meanAbsErr, na.rm = TRUE),]) # 5,10,17,20, err = 0.18
RegressionWBreaks(breakpoints = c(5,10,17,20), animal = alligator, cvGp = 99)




#############################################################################
## DO THE SAME FOR ALL INSTANCES OF FIVE BREAKPOINTS, excluding single vert groups
#############################################################################
erDF5 <- matrix(ncol = 6, nrow = 0)
for(kk in 1:nrow(fivebp)){
     cvError <- numeric()
     for(ii in 1:10){
          aa <- RegressionWBreaks_noPlot(breakpoints = as.numeric(fivebp[kk,]), dataset1= dataset1,
                                         animal = alligator, cvGp = ii)
          if(!is.na(aa[1])) cvError[ii] <- aa$mean_abs_err
          else cvError[ii] <- NA
     }
     
     erDF5 <- rbind(erDF5, c( fivebp[kk,], round(mean(cvError, na.rm = TRUE), 4)))
     print(paste(kk, "of", nrow(fivebp)))
}
erDF5 <- as.data.frame(apply(erDF5, 2, unlist))
colnames(erDF5) <- c("breakpt1", "breakpt2", "breakpt3", "breakpt4" , "bp5","meanAbsErr")
erDF5 <- data.frame(erDF5)
na.omit(erDF5[erDF5$meanAbsErr == min(erDF5$meanAbsErr, na.rm = TRUE),]) # 5  ,11 , 13,17,20

RegressionWBreaks(breakpoints = c(5  ,11 , 13,17,20), animal = alligator, cvGp = 1)

#############################################################################
## DO THE SAME FOR ALL INSTANCES OF SIX BREAKPOINTS, excluding single vert groups
#############################################################################
erDF6 <- matrix(ncol = 7, nrow = 0)
for(kk in 1:nrow(sixbp)){
     cvError <- numeric()
     for(ii in 1:10){
          aa <- RegressionWBreaks_noPlot(breakpoints = as.numeric(sixbp[kk,]), dataset1= dataset1,
                                         animal = alligator, cvGp = ii)
          if(!is.na(aa[1])) cvError[ii] <- aa$mean_abs_err
          else cvError[ii] <- NA
     }
     
     erDF6 <- rbind(erDF6, c( sixbp[kk,], round(mean(cvError, na.rm = TRUE), 4)))
     print(paste(kk, "of", nrow(sixbp)))
}
erDF6 <- as.data.frame(apply(erDF6, 2, unlist))
colnames(erDF6) <- c("breakpt1", "breakpt2", "breakpt3", "breakpt4" , "bp5", "bp6","meanAbsErr")
erDF6 <- data.frame(erDF6)
na.omit(erDF6[erDF6$meanAbsErr == min(erDF6$meanAbsErr, na.rm = TRUE),]) # 3,7,11,13,17,20   error= 0.1899

RegressionWBreaks(breakpoints = c(3,7,11,13,17,20), animal = alligator, cvGp = 1)

















err3 <- numeric()
### Test AREA
alligator[, "cvGroup"] <- sample(alligator[, "cvGroup"]) # reshuffle cvGroup
cvError <- numeric()
for(ii in 1:10){
     # png(paste("~/Desktop/testStack/", formatC(as.numeric(paste0(1, kk, ii)), width = 5, flag = 0), 
     # ".png", sep = ""))
     aa <- RegressionWBreaks(breakpoints = c(5  ,11), animal = alligator, cvGp = ii)
     # dev.off()
     if(!is.na(aa[1])) cvError[ii] <- aa$mean_abs_err
     else cvError[ii] <- NA
}

err3 = c(err3, round(mean(cvError, na.rm = TRUE), 4))
plot(err3)
abline(h = mean(err3))

mean(err3)
mean(err4)
mean(err5)
mean(err6)

alligator[, "cvGroup"] <- sample(alligator[, "cvGroup"]) # reshuffle cvGroup
RegressionWBreaks(breakpoints = c(5 ,10,20), animal = alligator, cvGp = 5)

