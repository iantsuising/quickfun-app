#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(fdapace)

source('functions.R', local = TRUE)

# Define server logic required to draw a histogram
shinyServer(function(input, output) {

  vals <- reactiveValues()
  observe({
    vals$path1 <- input$file1
    vals$path2 <- input$file2
    vals$M0 <- input$M0
    # vals$run <- input$run
    vals$pre_smooth <- input$pre_smooth
    vals$FVE <- input$FVE

    if (is.null(vals$path1)) {
      set.seed(12345)
      library(refund)
      data(DTI)
      X <- DTI$cca
      y <- DTI$case
      # y[y==0] <- "Normal"
      # y[y==1] <- "MS"
      t <- seq(0, 1, length.out = ncol(X))
      #not run
      #datacheck(X, y, t)
      allData <- cbind(X, y)
      allData <- allData[which(DTI$visit == 1), ]
      allData <- na.omit(allData)
      y <- allData[, ncol(allData)]
      X <- allData[, -ncol(allData)]
      t <- seq(0, 1, length.out = ncol(X))
      datacheck(X, y, t)
      Index_0 <- which(y == 0)
      Index_1 <- which(y == 1)
      SplitPara <- 0.8 #Split parameter
      trainIndex_0 <- sample(Index_0, SplitPara * length(Index_0))
      testIndex_0 <- Index_0[-trainIndex_0]
      trainIndex_1 <- sample(Index_1, SplitPara * length(Index_1))
      testIndex_1 <- numeric()
      for(i in Index_1){
        if(i %in% trainIndex_1 ==FALSE){
          testIndex_1 <- append(testIndex_1, i)
        }
      }
      trainset <- X[c(trainIndex_0, trainIndex_1), ]
      trainlabel <- y[c(trainIndex_0, trainIndex_1)]
      testset <- X[c(testIndex_0, testIndex_1), ]
      testlabel <- y[c(testIndex_0, testIndex_1)]
      vals$x <- FOBJ(trainset, trainlabel, t,
                     optns = list(smooth=vals$pre_smooth, FVEthreshold=vals$FVE, a0=0.02, M0=vals$M0, par.sig=1000))
      vals$c1 <- Classif.quantile(vals$x, testset, testlabel, weighted = FALSE)

      vals$graColors <- c("red", "blue", "black", "darkgreen")
      # vals$PointType <- c(3, 17, 5, 0)
      J <-  4

      vals$testlabel <- testlabel
    } else {
      req(vals$path1)
      req(vals$path2)
      df.trn <- as.matrix(read.csv(vals$path1$datapath))
      df.tst <- as.matrix(read.csv(vals$path2$datapath))
      pp <- ncol(df.trn)
      trainset <- df.trn[, 1:(pp-1)]
      trainlabel <- df.trn[, pp]
      t <- seq(0,1,length.out = (pp-1))
      testset <- df.tst[, 1:(pp-1)]
      testlabel <- df.tst[, pp]
      vals$x <- FOBJ(trainset, trainlabel, t,
                     optns = list(smooth=vals$pre_smooth, FVEthreshold=vals$FVE, a0=0.02, M0=vals$M0, par.sig=1000))
      vals$c1 <- Classif.quantile(vals$x, testset, testlabel, weighted = FALSE)

      vals$graColors <- c("red", "blue", "black", "darkgreen")
      # vals$PointType <- c(3, 17, 5, 0)
      J <-  4

      vals$testlabel <- testlabel
    }
  })


  # observeEvent(input$run,{


  output$text <- renderText({
    paste0("Misclassification rate on testing set: ", round(vals$c1$MCR,4),".")
  })

  output$pretable <- DT::renderDataTable({
    if(TRUE){
      DT::datatable(data = data.frame("Index" = 1:length(vals$testlabel), "True" = (vals$testlabel), "Predicted" = (vals$c1$pre)),
                    options = list(pageLength = 10),
                    rownames = FALSE)
    }})




  output$distPlot <- renderPlot({
    x <- vals$x
    par(mfrow = c(1,2))
    plot(x$t, x$DataGroupby[[1]][1, ],
         main = "Functional Data Curve",
         xlab = " ", ylab = " ",
         ylim = c(min(unlist(x$DataGroupby)) - 0.2 * abs(min(unlist(x$DataGroupby))),
                  max(unlist(x$DataGroupby)) + 0.2 * abs(max(unlist(x$DataGroupby)))),
         type = "l", lty = 1)

    for(i in 1:min(x$K, 4)){
      apply(x$DataGroupby[[i]], 1, FUN = function(y){lines(x$t, y, col = vals$graColors[i], lty = i)})
    }
    legend("topleft", legend = x$y[1:min(x$K, 4)], col = vals$graColors[1:min(x$K, 4)], lty = 1:min(x$K, 4))


    plot(x$t, apply(x$DataGroupby[[1]], 2, mean),
         main = "Mean curve of each group",
         xlab = " ", ylab = " ",
         ylim = c(min(unlist(x$DataGroupby)) - 0.2 * abs(min(unlist(x$DataGroupby))),
                  max(unlist(x$DataGroupby)) + 0.2 * abs(max(unlist(x$DataGroupby)))),
         type = "l", lty = 1, col = vals$graColors[1])
    for(i in 2:min(x$K, 4)){
      lines(x$t, apply(x$DataGroupby[[i]], 2, mean), col = vals$graColors[i], lty = i)
    }
    legend("topleft", legend = x$y[1:min(x$K, 4)], col = vals$graColors[1:min(x$K, 4)], lty = 1:min(x$K, 4))

  })

  output$distCov <- renderPlot({
    x <- vals$x
    # generate bins based on input$bins from ui.R
    # graColors <- c("red", "blue", "black", "darkgreen")
    # PointType <- c(3, 17, 5, 0)
    J <-  4
    if(x$K == 2) {
      par(mfrow = c(1,2))}
    else {
      par(mfrow = c(1,1))}



    autofit1 <- abs(mean(x$t) / mean(x$poolCov))
    z1 <- autofit1 * x$poolCov
    lim1 <- range(z1)
    len1 <- lim1[2] - lim1[1] + 1
    colorset1 <- terrain.colors(len1, alpha = 0.5)
    color1 <- colorset1[z1 - lim1[1] + 1]
    persp(5 * x$t, 5 * x$t, z1, xlab = "t", ylab =  "s", zlab = " ", col = color1, main = "Pooled Covariance", theta = 45)

    if(x$K == 2){
      autofit2 <- abs(mean(x$t) / mean(x$singleCov[[1]] - x$singleCov[[2]]))
      z2 <- autofit2 * (x$singleCov[[1]] - x$singleCov[[2]])
      lim2 <- range(z2)
      len2 <- lim2[2] - lim2[1] + 1
      colorset2 <- terrain.colors(len2, alpha = 0.5)
      color2 <- colorset2[z2 - lim2[1] + 1]
      persp(5 * x$t, 5 * x$t, z2, xlab = "t", ylab =  "s", zlab = " ", col = color1, main = "Difference Between Two Covariance Functions", theta = 45)
    }

  })


  output$distPdf <- renderPlot({
    x <- vals$x
    par(mfrow = c(2,2))
    for(i in 1:min(4,x$cutPoints)){

      plx <- numeric()
      ply <- numeric()
      for(ii in 1:min(4, x$K)){
        plx <- append(plx, density(x$scores[[ii]][, i])$x)
        ply <- append(ply, density(x$scores[[ii]][, i])$y)
      }
      plot(density(x$scores[[1]][, i])$x, density(x$scores[[1]][, i])$y,
           type = "l", lty = 1, col = vals$graColors[1], xlab = " ", ylab = " ",
           xlim = range(plx), ylim = range(ply))
      for(j in 2:min(4, x$K)){
        lines(density(x$scores[[j]][, i])$x, density(x$scores[[j]][, i])$y,
              col = vals$graColors[j], lty = j)
      }
      legend("topleft", legend = x$y[1:min(4, x$K)],
             col = vals$graColors[1:min(4, x$K)], lty = 1:min(4, x$K))
      title(main = paste0("Density of scores on FPC", as.character(i), " ", "(", as.character(round(x$proportion[i], 3)) ,")"))

    }
  })

  output$distMeanVar <- renderPlot({
    x <- vals$x
    J <- 4
    par(mfrow = c(1,2))
    toplotx <- 1:min(J, x$cutPoints)
    plotymean <- numeric()
    for(i in 1:min(4, x$K)){
      plotymean <- append(plotymean, apply(x$scores[[i]][, toplotx], 2, mean))
    }

    plot(toplotx, apply(x$scores[[1]][, toplotx], 2, mean), type = "l", lty = 1,
         col = vals$graColors[1], xlab = " ", ylab = " ", ylim = range(plotymean),
         main = "The mean of projection scores from each group",
         xaxt = "n")
    for(j in 2:min(4, x$K)){
      lines(toplotx, apply(x$scores[[j]][, toplotx], 2, mean), lty = j, col = vals$graColors[j])
    }
    legend("topright", legend = x$y[1:min(4, x$K)], col = vals$graColors[1:min(4, x$K)], pch = rep(1, min(4, x$K)))
    axis(1, 1:min(J, x$cutPoints))

    plotyvar <- numeric()
    for(i in 1:min(4, x$K)){
      plotyvar <- append(plotyvar, apply(x$scores[[i]][, toplotx], 2, var))
    }

    plot(toplotx, apply(x$scores[[1]][, toplotx], 2, var), type = "l", lty = 1,
         col = vals$graColors[1], xlab = " ", ylab = " ", ylim = range(plotyvar),
         main = "The variance of projection scores from each group",
         xaxt = "n")
    for(j in 2:min(4, x$K)){
      lines(toplotx, apply(x$scores[[j]][, toplotx], 2, var), lty = j, col = vals$graColors[j])
    }
    legend("topright", legend = x$y[1:min(4, x$K)], col = vals$graColors[1:min(4, x$K)], pch = rep(1, min(4, x$K)))
    axis(1, 1:min(J, x$cutPoints))


  })

})
# })

