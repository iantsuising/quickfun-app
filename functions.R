#'@title Quantile Classifier
#'@description Quantile classifier
#'@param obj An FOBJ for classification
#'@param testset Must be a matrix of testset. Each row is an observation of curve. In terms of the prediction to single curve, view \code{pre.quantile}
#'@param testlabel Testlabel, which length should be equal to the number of rows of testset
#'@param weighted Logical. It is to define whether weight is to be used. Default the same to \code{obj$weighted}.
#'@return A list containing the classification result
#'\item{MCR}{The total misclassification}
#'\item{pre}{The predictive result}
#'\item{FP}{False positive rate if only two groups}
#'\item{FN}{False negative rate if only two groups}
#'@examples
#'rm(list = ls())
#'library(refund)
#'data(DTI)
#'X <- DTI$cca
#'y <- DTI$case
#'t <- seq(0, 1, length.out = ncol(X))
#' #not run
#' #datacheck(X, y, t)
#' allData <- cbind(X, y)
#' allData <- allData[which(DTI$visit == 1), ]
#' allData <- na.omit(allData)
#' y <- allData[, ncol(allData)]
#' X <- allData[, -ncol(allData)]
#' t <- seq(0, 1, length.out = ncol(X))
#' datacheck(X, y, t)
#' Index_0 <- which(y == 0)
#' Index_1 <- which(y == 1)
#'  SplitPara <- 0.8 #Split parameter
#'  trainIndex_0 <- sample(Index_0, SplitPara * length(Index_0))
#'  testIndex_0 <- Index_0[-trainIndex_0]
#'  trainIndex_1 <- sample(Index_1, SplitPara * length(Index_1))
#'  testIndex_1 <- numeric()
#'  for(i in Index_1){
#'      if(i %in% trainIndex_1 ==FALSE){
#'         testIndex_1 <- append(testIndex_1, i)
#'        }
#'     }
#'  trainset <- X[c(trainIndex_0, trainIndex_1), ]
#'  trainlabel <- y[c(trainIndex_0, trainIndex_1)]
#'   testset <- X[c(testIndex_0, testIndex_1), ]
#'  testlabel <- y[c(testIndex_0, testIndex_1)]
#'  w1 <- FOBJ(trainset, trainlabel, t)
#'  c1 <- Classif.quantile(w1, testset, testlabel, weighted = FALSE)
#'@export

Classif.quantile <- function(obj, testset, testlabel, weighted = obj$weighted){
  
  if (obj$optns$smooth == 1){
    
    smoothlist <- apply(testset, 1, function(x){  #use apply() to smooth
      spline(obj$origiTime,x,method = "natural", xout = obj$t)
    })
    # obj$t <- smoothlist[[1]]$x #replace original t with smoothed t
    testset <- matrix(0, nrow = nrow(testset), ncol = length(obj$t))
    for (i in 1:nrow(testset)){
      testset[i, ] <- smoothlist[[i]]$y
    } #repalce original X with smoothed X
  }
    
    if(is.vector(testset)){stop("Not a matrix")}
    
    datacheck(testset, testlabel, obj$t)
    # cat(paste0("The weight of FOBJ is", as.character(obj$weighted)))
    
    if(! weighted){
        
        pre <- apply(testset, 1, FUN = function(x){pre.quantile(x, obj, weighted = FALSE)})
        
    }
    
    else{
        pre <- apply(testset, 1, FUN = function(x){pre.quantile(x, obj, weighted = TRUE)})
        
    }
    mcr <- length(which(pre != testlabel))/length(testlabel)
    
    if(obj$K != 2){
        return(list(MCR = mcr, pre = pre))
    }
    
    else{
        fp <- length(which(testlabel != pre & testlabel == obj$y[1]))/length(which(testlabel == obj$y[1]))
        fn <- length(which(testlabel != pre & testlabel == obj$y[2]))/length(which(testlabel == obj$y[2]))
        return(list(MCR = mcr, FP = fp, FN = fn, pre = pre))
    }
}





#This function can be used for data check
#' @title Check Input Data Tyoe
#' @description \code{datacheck()} is to check whether the data to be input is appropriate
#' @param data The input matrix of observed values. Rows represent repititions or samples, and columns represent observed time points
#' @param label The input vector of labels. ith element of this vector should correspond to the label of ith row of data matrix X.
#' @param time An input time grid vector
#' @return NULL is appropriate, otherwise warning
#' @examples
#'  library(refund)
#'  data(DTI)
#'  X <- DTI$cca
#'  y <- DTI$case
#'  t <- seq(0, 1, length.out = ncol(X))
#'   # not run
#'   # datacheck(X, y, t)
#'  allData <- cbind(X, y)
#'  allData <- allData[which(DTI$visit == 1), ]
#'  allData <- na.omit(allData)
#'  y <- allData[, ncol(allData)]
#'  X <- allData[, -ncol(allData)]
#'  t <- seq(0, 1, length.out = ncol(X))
#'  datacheck(X, y, t)
#'@export
datacheck <-function(data, label, time){
    if (!is.vector(label)){
        stop("Input label should be vector \n")
    }
    if (!is.vector(time)){
        stop("Input time should be vector \n")
    }
    if (!is.numeric(time)){
        stop("Input time should be numbers \n")
    }
    if (!is.matrix(data)){
        stop("Input data should be matrix \n")
    }
    if (nrow(data) != length(label)){
        stop("The amount of rows of input data should be the same as the length of input label \n")
    }
    if (ncol(data) != length(time)){
        stop("The amount of cloumns of input data should be the same as the length of input time \n")
    }
    if (sum(is.na(data)) > 0){
        stop("NA is not permitted")
    }
    if (length(unique(label))>=3){
        warning("The outcome may not be reliable if there are more than 2 classes")
    }
    return(NULL)
}



#'@title Find optimal weights based on the quantileloss
#'@description A weight is multiplied to the quantileloss to select the important projection scores.
#'@param weight The initial weight to be optimized
#'@param scorematrixlist A list of quantileloss matrices
#'@param par.sig The parameter in sigmoid function, default is 1000
#'@return It returns a list, which is the output of \code{optim()} function. Details can be seen in \code{\link{optim}}.
#'@import stats

findweights <- function(weight, scorematrixlist, par.sig = 1000){
    K=length(scorematrixlist)
    
    diffmatrix <- vector("list",K*(K-1))
    
    count <- 1
    for (i in 1:length(scorematrixlist)){  #generate list of diffmatrix
        for (j in 1:length(scorematrixlist)){
            if (i!=j){
                diffmatrix[[count]] <- scorematrixlist[[i]][[j]]-scorematrixlist[[i]][[i]]
                count <- count + 1
            }
        }
    }
    
    fn <- function(x){  #generalized function
        res <- 0
        for (i in 1:length(scorematrixlist)){
            temp <- do.call(rbind, diffmatrix[((i-1)*(K-1)+1):((i-1)*(K-1)+K-1)]) %*% x #will return a long vector
            
            sigtemp <- sigmoid(matrix(temp, nrow = (length(temp)/(K-1))), par.sig = par.sig) #divide into matrix with
            product <- numeric()
            
            if(is.vector(sigtemp)){
                product <- sigtemp
            }
            
            else{
                product <- apply(sigtemp, 1, prod)
            }
            
            res=res+sum(product)
            
        }
        return(-res)
    }
    
    #optimlist <- Rsolnp::solnp(weight,fn,eqfun = sum,eqB = 1,control = list(trace = 0))#,tol = 1e-3))
    
    optimlist <- optim(par=weight, fn=fn,lower = 0, upper = 10,method = "L-BFGS-B")
    # optimlist=optim(par=weight,fn=fn,method = "BFGS")
    
    return(optimlist)
}



#ui=rbind(rep(1,length(weight)),rep(-1,length(weight)))
#ci=c(0.1,-0.1)




# optimlist <- solnp(initialweight,fn)





#'@title Functional Objection
#'@description \code{FOBJ()} is used to compute imformation for within FDA framework.
#'@param X The input matrix of observed values. Rows represent repititions or samples, and columns represent observed time points.
#'@param y The input vector of labels. ith element of this vector should correspond to the label of ith row of data matrix X.
#'@param t An input time grid vector.
#'@param optns A list to control computing such as smoothing, threshold and quantile range.
#'@param weighted Whether weight for quantile loss is to be used. \code{FALSE} is default.
#'@return  It returns a class "FOBJ".
#'\item{DataGroupby}{The list group by the label of each group}
#'\item{y}{Vector. Label of each group}
#'\item{t}{Original time points of observation}
#'\item{K}{Numeber. How many groups}
#'\item{phi}{Matrix. The functional principle componants, listed by colums}
#'\item{scores}{List. The projection scores of each group}
#'\item{cutPoints}{Numebr. How many componants}
#'\item{quantlist}{List. The quantile of projection scores of each group}
#'\item{tau}{Vector. Quantile counter}
#'\item{weight}{Matrix. Weight for classification. Default is matrix 1. The number of colums is equal to \code{tau}}
#'\item{optns}{List. The same to the original optns}
#'\item{weighted}{Logical. The same to the parameter "weighted"}
#'\item{propotion}{Vector. Contains the proportion of variance of each FPC}
#'@import fdapace
#'@import stats
#'@examples
#'rm(list = ls())
#'library(refund)
#'data(DTI)
#'X <- DTI$cca
#'y <- DTI$case
#'t <- seq(0, 1, length.out = ncol(X))
#' #not run
#' #datacheck(X, y, t)
#' allData <- cbind(X, y)
#' allData <- allData[which(DTI$visit == 1), ]
#' allData <- na.omit(allData)
#' y <- allData[, ncol(allData)]
#' X <- allData[, -ncol(allData)]
#' t <- seq(0, 1, length.out = ncol(X))
#' w1 <- FOBJ(X, y, t)
#'@export


FOBJ <- function(X, y, t, optns = list(smooth=FALSE, FVEthreshold=99.99, a0=0.02, M0=10, par.sig=1000), weighted = FALSE){
    
    datacheck(X, y, t) #write a sub function to check data format
    
    origiTime <- t #save t
    #Whether smooth or not
    if (optns$smooth == 1){
        smoothlist <- apply(X, 1, function(x){  #use apply() to smooth
            spline(t,x,method = "natural")
        })
        t <- smoothlist[[1]]$x #replace original t with smoothed t
        X <- matrix(0, nrow = nrow(X), ncol = length(t))
        for (i in 1:nrow(X)){
            X[i, ] <- smoothlist[[i]]$y
        } #repalce original X with smoothed X
    }
    #group X according to label into group X0 and X1 as input for pooledcov(). Here is K class
    
    K <- length(unique(y))
    Xlist <- vector("list", K)
    bindX <- cbind(X, y)
    #construct index list to record original sequence of labels
    ind <- vector("list",K)
    for (j in 1:K){
        
        ind[[j]] <- which(bindX[, ncol(bindX)]==(unique(y))[j]) #record row index for each group i
        
    }
    
    names(ind) <- c(unique(y)[1: K])
    
    for (i in 1:K){
        Xlist[[i]] <- bindX[y==(unique(y))[i],]
        Xlist[[i]] <- matrix(as.numeric(Xlist[[i]][, -ncol(bindX)]),nrow(Xlist[[i]]))
    }
    names(Xlist) <- c(unique(y)[1:K]) #assign names to Xlist
    
    # get covariance function
    Covlist=pooledcov(Xlist,t)
    
    #Compute basis function by eigendecomposition
    temp.phi <- getPhi(Covlist$Cov, optns$FVEthreshold, t)
    
    fitphi <- temp.phi$phi
    
    eigenV <- temp.phi$eigenV
    
    cutPoints <- temp.phi$cutPoints
    
    proportion <- temp.phi$proportion
    
    # compute xi
    xilist=vector("list",length(Xlist))
    for (jj in 1:length(Xlist)){
        xilist[[jj]] <- getXi(Xlist[[jj]], fitphi[, 1:temp.phi$cutPoints], t)
    }
    names(xilist) <- c(unique(y)[1:K])
    J=ncol(xilist[[1]])
    
    #compute quantile for each column(projection scores).
    qlist=vector("list",length(Xlist))
    a0=optns$a0
    M0=optns$M0
    theta=seq(a0, 1-a0, (1-2*a0)/(M0-1))
    Ltheta=(1:length(theta))
    
    for (kk in 1:length(Xlist)){ #quantile matrix for each group
        qlist[[kk]]=apply(xilist[[kk]], 2, function(x){ #replace all loops by apply()
            quantile(x,theta)
        })
    }
    names(qlist) <- c(unique(y)[1:K])
    
    #compute weight for quantile classifier
    if(!weighted){
        weight <- rep(1, J *length(theta))
        dim(weight <- c(J, length(theta)))
    }
    
    else{
        weight <- rep(0, J * length(theta))
        
        scorelist <- vector("list",length(Xlist))
        
        for (K in 1:length(Xlist)){
            for (m in 1:length(Xlist)){
                scorelist[[K]][[m]]=lapply(Ltheta,function(x){
                    quantileloss(xilist[[K]],qlist[[m]][x,  ],theta[x])
                })
            }
            
        }
        names(scorelist) <- c(unique(y)[1:K])
        
        #convert scorelist to array list
        arraylist=vector("list",length(Xlist))
        for (k in 1:length(Xlist)){
            for (m in 1:length(Xlist)){
                arraylist[[k]][[m]]=array(as.numeric(unlist(scorelist[[k]][[m]])), dim=c(nrow(xilist[[k]]), J, length(theta)))
            }
            
        }
        names(arraylist) <- c(unique(y)[1:K])
        
        #convert array to long matrix list
        scorematrixlist=vector("list",length(Xlist))
        for (k in 1:length(Xlist)){
            for (m in 1:length(Xlist)){
                scorematrixlist[[k]][[m]]=matrix(arraylist[[k]][[m]], nrow(xilist[[k]]), (J)*length(theta))
            }
        }
        names(scorematrixlist) <- c(unique(y)[1:K])
        
        optwei <- findweights(scorematrixlist = scorematrixlist,weight = weight,par.sig = optns$par.sig)$par
        weight <- optwei
        dim(weight) <- c(J, length(theta))
    }
    
    out <- list(DataGroupby = Xlist, y = unique(y), origiTime = origiTime, t=t,
    K = K, poolCov = Covlist$Cov, singleCov = Covlist$singlecov,
    phi = fitphi, scores = xilist, cutPoints = cutPoints,
    proportion = proportion,
    quantlist = qlist, tau = theta, weight = weight, optns = optns,
    weighted = weighted)
    
    class(out) = "FOBJ"
    
    return(out)
}


#'@title Get functional priciple componants
#'@description FPCA
#'@param Cov Covariance surface
#'@param FVEthreshold FVE
#'@param t Time points
#'@import fdapace
#'@references Yao F, Müller H G, Wang J L. Functional data analysis for sparse longitudinal data[J]. Journal of the American Statistical Association, 2005, 100(470): 577-590.
#'@export
getPhi <- function(Cov, FVEthreshold=99.99, t){
    eig <- eigen(Cov) #return a list with value and vector
    positiveInd <- eig[['values']] >= 0 #why not use dollar sign?
    if (sum(positiveInd) == 0) {
        stop('All eigenvalues are negative. The covariance estimate is incorrect.')
    }
    d <- eig[['values']][positiveInd] #same?
    
    eigenV <- eig[['vectors']][, positiveInd, drop=FALSE] #extract positive eigenvalue column
    
    
    FVE <- cumsum(d) / sum(d) * 100  # cumulative FVE for all available eigenvalues from fitted cov
    
    cutPoints <- min(which(FVE >= FVEthreshold)) # final number of component chosen based on FVE.Only 29 is choosen in demo.
    
    maxK <- cutPoints #maximum K candidate
    
    d <- d[1:maxK]
    
    proportion <- d/sum(d) #Contribution of each PC
    
    eigenV <- eigenV[, 1:maxK, drop=FALSE]
    
    
    # normalization
    muWork = 1:dim(eigenV)[1] #return the row number of eigenV. Since column need to be reduced
    
    
    phi0 <- apply(eigenV, 2, function(x) { #eigen V is matrix. 2 means apply function to column.
        x <- x / sqrt(fdapace::trapzRcpp(t, x^2))# divide each column by a constant. This integration see paper
        if ( 0 <= sum(x*muWork) )
        return(x)
        else
        return(-x)
    })
    
    return(list(phi = phi0, FVE = FVE, cutPoints = cutPoints, d = d, proportion = proportion))
}


#'@title Get projection scores from the estimated eigenfunctions
#'@description After getting estimated eigenfunctions, this function can be used to compute projection scores of observed functional data.
#'@param X.curve The input matrix of observed values. Rows represent repititions or samples, and columns represent observed time points
#'@param phi A matrix containing normalized eigenfunctions, which is always generated by \code{FOBJ}
#'@param t A input time grid vector
#'@return It returns a matrix:
#'\item{Xi}{A matrix containing the projection scores}
#'@import fdapace
#'@export
getXi <- function(X.curve, phi, t){ #This function get the projection score. See the paper algorithm part.
    
    if(! is.numeric(t)) stop("t should be recorded time points!")
    
    
    xi <- NULL
    
    if(is.vector(X.curve)){
        xi <- apply(apply(phi, 2, function(x) {x * X.curve}), 2, function(x){return(trapzRcpp(t, x))})
    }
    
    else{
        for (i in 1:ncol(phi)){
            xi <- cbind(xi ,apply(t(matrix(rep(phi[, i], nrow(X.curve)), length(t))) * X.curve, 1, function(x){return(trapzRcpp(t, x))}))
        } #The Jth column is the series of Jth scores
    }
    return(xi)
}


#'@title Visulization of FOBJ
#'@description This function is aimed to plot the functional data from different aspects
#'@param x An FOBJ, generated by \code{FOBJ}
#'@param J The number of PCs to be plotted. Only positive number is allowed
#'@method plot FOBJ
#'@import stats
#'@import graphics
#'@import grDevices
#'@export


plot.FOBJ <- function(x, J = 4){
    
    if(! 'FOBJ' %in% class(x)){stop("Invalid input class")}
    
    if(x$K > 4) cat("Only the first 4 groups will be displayed")
    
    if(J <= 0) {stop("Only positive number allowed")}
    
    graColors <- c("red", "blue", "black", "darkgreen")
    PointType <- c(3, 17, 5, 0)
    
    yesorno <- readline(prompt = "To continue hit return/Type n to break")
    if(yesorno == "n"){return("End")}
    else{
        
        plot(x$t, x$DataGroupby[[1]][1, ],
        main = "Functional Data Curve",
        xlab = " ", ylab = " ",
        ylim = c(min(unlist(x$DataGroupby)) - 0.2 * abs(min(unlist(x$DataGroupby))),
        max(unlist(x$DataGroupby)) + 0.2 * abs(max(unlist(x$DataGroupby)))),
        type = "l", lty = 1)
        
        for(i in 1:min(x$K, 4)){
            apply(x$DataGroupby[[i]], 1, FUN = function(y){lines(x$t, y, col = graColors[i], lty = i)})
        }
        legend("topleft", legend = x$y[1:min(x$K, 4)], col = graColors[1:min(x$K, 4)], lty = 1:min(x$K, 4))
    } #Original data
    
    yesorno <- readline(prompt = "To continue hit return/Type n to break")
    if(yesorno == "n"){return("End")}
    else{
        plot(x$t, apply(x$DataGroupby[[1]], 2, mean),
        main = "Mean curve of each group",
        xlab = " ", ylab = " ",
        ylim = c(min(unlist(x$DataGroupby)) - 0.2 * abs(min(unlist(x$DataGroupby))),
        max(unlist(x$DataGroupby)) + 0.2 * abs(max(unlist(x$DataGroupby)))),
        type = "l", lty = 1, col = graColors[1])
        for(i in 2:min(x$K, 4)){
            lines(x$t, apply(x$DataGroupby[[i]], 2, mean), col = graColors[i], lty = i)
        }
        legend("topleft", legend = x$y[1:min(x$K, 4)], col = graColors[1:min(x$K, 4)], lty = 1:min(x$K, 4))
    } #Mean curves
    
    
    
    yesorno <- readline(prompt = "To continue hit return/Type n to break")
    if(yesorno == "n"){return("End")}
    else{
        autofit1 <- abs(mean(x$t) / mean(x$poolCov))
        z1 <- autofit1 * x$poolCov
        lim1 <- range(z1)
        len1 <- lim1[2] - lim1[1] + 1
        colorset1 <- terrain.colors(len1, alpha = 0.5)
        color1 <- colorset1[z1 - lim1[1] + 1]
        persp(5 * x$t, 5 * x$t, z1, xlab = "t", ylab =  "s", zlab = " ", col = color1, main = "Pooled Covariance", theta = 45)
        
    } #Pooled covariance surface
    if(x$K == 2){
        yesorno <- readline(prompt = "To continue hit return/Type n to break")
        if(yesorno == "n"){return("End")}
        else{
            autofit2 <- abs(mean(x$t) / mean(x$singleCov[[1]] - x$singleCov[[2]]))
            z2 <- autofit2 * (x$singleCov[[1]] - x$singleCov[[2]])
            lim2 <- range(z2)
            len2 <- lim2[2] - lim2[1] + 1
            colorset2 <- terrain.colors(len2, alpha = 0.5)
            color2 <- colorset2[z2 - lim2[1] + 1]
            persp(5 * x$t, 5 * x$t, z2, xlab = "t", ylab =  "s", zlab = " ", col = color1, main = "Difference Between Two Covariance Functions", theta = 45)
        }
    }
    
    
    
    
    for(i in 1:min(J,x$cutPoints)){
        yesorno <- readline(prompt = "Hit return to continue/Type n to break")
        if(yesorno == "n"){return("End")}
        else{
            plx <- numeric()
            ply <- numeric()
            for(ii in 1:min(4, x$K)){
                plx <- append(plx, density(x$scores[[ii]][, i])$x)
                ply <- append(ply, density(x$scores[[ii]][, i])$y)
            }
            plot(density(x$scores[[1]][, i])$x, density(x$scores[[1]][, i])$y,
            type = "l", lty = 1, col = graColors[1], xlab = " ", ylab = " ",
            xlim = range(plx), ylim = range(ply))
            for(j in 2:min(4, x$K)){
                lines(density(x$scores[[j]][, i])$x, density(x$scores[[j]][, i])$y,
                col = graColors[j], lty = j)
            }
            legend("topleft", legend = x$y[1:min(4, x$K)],
            col = graColors[1:min(4, x$K)], lty = 1:min(4, x$K))
            title(main = paste0("Density of scores on FPC", as.character(i), " ", "(", as.character(round(x$proportion[i], 3)) ,")"))
        }
    } #Density of FPCs
    
    
    for(i in 1:min(J, x$cutPoints)){
        yesorno <- readline(prompt = "Hit return to continue/Type n to break")
        if(yesorno == "n"){return("End")}
        else{
            theta = seq(0.1, 0.9, 0.1)
            ploty <- numeric()
            for(ii in 1:min(4, x$K)){
                ploty <- append(ploty, x$scores[[ii]][, i])
            }
            plot(theta, quantile(x$scores[[1]][, i], theta), col = graColors[1],
            ylim = range(ploty),
            pch = PointType[1],
            xlab = " ", ylab = " ")
            for(j in 2:min(4, x$K)){
                points(theta, quantile(x$scores[[2]][, i], theta), col = graColors[j], pch = PointType[j])
            }
            
            legend("topleft", legend = x$y[1:min(4, x$K)],
            col = graColors[1:min(4, x$K)],
            pch = PointType[1:min(4, x$K)])
            title(main = paste0("Quantile of scores on FPC", as.character(i), " ",
            "(", as.character(round(x$proportion[i], 3)),")"))
        }
    } #Quantile of FPCs
    
    
    yesorno <- readline(prompt = "To continue hit return/Type n to break")
    if(yesorno == "n"){return("End")}
    else{
        toplotx <- 1:min(J, x$cutPoints)
        plotymean <- numeric()
        for(i in 1:min(4, x$K)){
            plotymean <- append(plotymean, apply(x$scores[[i]][, toplotx], 2, mean))
        }
        
        plot(toplotx, apply(x$scores[[1]][, toplotx], 2, mean), type = "l", lty = 1,
        col = graColors[1], xlab = " ", ylab = " ", ylim = range(plotymean),
        main = "The mean of projection scores from each group",
        xaxt = "n")
        for(j in 2:min(4, x$K)){
            lines(toplotx, apply(x$scores[[j]][, toplotx], 2, mean), lty = j, col = graColors[j])
        }
        legend("topright", legend = x$y[1:min(4, x$K)], col = graColors[1:min(4, x$K)], pch = rep(1, min(4, x$K)))
        axis(1, 1:min(J, x$cutPoints))
    }  #Mean of FPCs
    
    yesorno <- readline(prompt = "To continue hit return/Type n to break")
    if(yesorno == "n"){return("End")}
    else{
        plotyvar <- numeric()
        for(i in 1:min(4, x$K)){
            plotyvar <- append(plotyvar, apply(x$scores[[i]][, toplotx], 2, var))
        }
        
        plot(toplotx, apply(x$scores[[1]][, toplotx], 2, var), type = "l", lty = 1,
        col = graColors[1], xlab = " ", ylab = " ", ylim = range(plotyvar),
        main = "The variance of projection scores from each group",
        xaxt = "n")
        for(j in 2:min(4, x$K)){
            lines(toplotx, apply(x$scores[[j]][, toplotx], 2, var), lty = j, col = graColors[j])
        }
        legend("topright", legend = x$y[1:min(4, x$K)], col = graColors[1:min(4, x$K)], pch = rep(1, min(4, x$K)))
        axis(1, 1:min(J, x$cutPoints))
    } #Var of FPCs
    
    for(i in 1:min(J, x$cutPoints)){
        yesorno <- readline(prompt = "Hit return to continue/Type n to break")
        if(yesorno == "n"){return("End")}
        else{
            plot(x$t, x$phi[, i], type = "l", col = "blue",
            xlab = " ", ylab = " ")
            title(main = paste0("FPC", as.character(i), "(", as.character(round(x$proportion[i], 3)) ,")"))
        }
    } #Plot of FPCs
    
}






#This function can return pooled covariance function given the input data
#X is matrix of curve. Each row means a single curve, each column means a obervation time. It should be
#already groupped when pass into the pooledcov()
#y is a vector of labels.
#'@title Pooled Covariance Function
#'@param Xlist A list of gruoped functional data.
#'@param t Time points
#'@import fdapace

pooledcov <- function(Xlist, t){
    
    nx=lapply(Xlist,function(x){
        nrow(x)
    })
    n=sum(unlist(nx))
    
    outlist=vector("list",length(Xlist))
    for (j in 1:length(Xlist)){
        LX0 <- lapply(seq_len(ncol(t(Xlist[[j]]))), function(i) Xlist[[j]][i, ])
        Lt0 <- rep(list(t), nx[[j]])
        outlist[[j]]=FPCA(LX0,Lt0)$fittedCov
    }
    
    
    Cov=0
    for (ii in 1:length(outlist)){
        Cov <- Cov+ (nx[[ii]] - 1) * outlist[[ii]]
        Cov <- Cov / (n - length(outlist) + 1)
    }
    
    
    return(list(Cov=Cov, singlecov=outlist))
}


#'@title Predict Function Using Quantile
#'@description Assign label to a single curve
#'@param X.curve A vector denotes the discrete points of a single curve.
#'@param obj An FOBJ object containing the information for classification.
#'@param weighted Whether multiple weight by \code{FOBJ} is to be used or not. Default the same to \code{obj$weighted}.
#'@import stats
#'@export
pre.quantile <- function(X.curve, obj, weighted = obj$weighted){
    
    if(! is.vector(X.curve)) stop("X is not a single curve")
    
    Xi0 <- getXi(X.curve, obj$phi, obj$t)
    L <- numeric(obj$K)
    
    if(!weighted){
        
        for(i in 1:obj$K){
            L[i] <- sum(sapply(1:length(obj$tau),
            FUN = function(x){quantileloss(Xi0, obj$quantlist[[i]][x, ], obj$tau[x])}))
        }
        
        return(obj$y[which.min(L)])
    }
    
    else{
        for(i in 1:obj$K){
            L[i] <- sum(sapply(1:length(obj$tau),
            FUN = function(x){
                quantileloss(Xi0, obj$quantlist[[i]][x, ], obj$tau[x])
            }) * obj$weight)
        }
        return(obj$y[which.min(L)])
    }
    
}


#This function will compute quantile loss for a single curve(equation in page 7. L)
#Xi is a vector/matrix of projection score value for ith curve, q is empirical quantile, tau is value of tau
#Will return a matrix, column length equal to the number of projection scores.
#tau here should be a single tau.
#q is vector, should not be matrix
#'@title Computing Quantileloss
#'@description Compute quantile loss for a single curve
#'@param Xi A vector/matrix of projection score value for ith curve
#'@param q Emperical quantile
#'@param tau Single tau
#'@references To be confirmed
#'@export
quantileloss <- function(Xi, q, tau){
    if(is.vector(Xi)){
        res <- abs(Xi - q) * ((tau * (Xi > q)) + (1 - tau)*(Xi <= q))
    }
    
    else{
        s.temp <- Xi - matrix(rep(q, nrow(Xi)), nrow(Xi), length(q), byrow = T)
        res <- tau * s.temp - s.temp * (s.temp < 0)
        # res=t(apply(Xi,1,function(x){
        #   res=tau*(x-q)-(x-q)*(x<q) #use apply to repalce loop
        # }))
    }
    return(res)
}


#'@title Sigmoid Function
#'@description Using Sigmoid function to smoothly assign a label
#'@param x vectors or values to be assigned
#'@param par.sig parameter for sigmoid function
#'@return 0 or 1, based on whether your input number or vector is positive or negative.
#'@examples
#' a <- runif(10, -1, 1)
#' b <- sigmoid(a)
#'@export
sigmoid <- function(x, par.sig=1000){
    res = 1/(1+exp(-par.sig*x));
    return(res)
}
