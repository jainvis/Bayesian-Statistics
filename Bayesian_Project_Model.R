# Script specifying the model
# Vishesh Jain, Abheer Sharma, Deepanshu Jain 

#-------------------------------Source File and Packages--------------------------------

source("~/DBDA2E-utilities.r") # Given by Dr. Haydar Demirhan, RMIT University

#------------------------------Function to Generate Chain-------------------------------

genMCMC = function( data , xName="x" , yName="y" , 
                    numSavedSteps=10000 , thinSteps=1 , saveName=NULL  ,
                    runjagsMethod=runjagsMethodDefault , 
                    nChains=nChainsDefault , xPred = xPred) { 
  require(runjags)
  
  #--------------------------------The Data-----------------------------------------
  
  y = data[,yName]/1000
  x = data.matrix(data[,xName])
  # Do some checking that data make sense:
  if ( any( !is.finite(y) ) ) { stop("All y values must be finite.") }
  if ( any( !is.finite(x) ) ) { stop("All x values must be finite.") }
  # for (i in nrow(x)){
  #   if ( ( !is.finite(x[i]) ) ) { stop("All x values must be finite.") }
  # }
  cat("\nCORRELATION MATRIX OF PREDICTORS:\n ")
  show( round(cor(x),3) )
  cat("\n")
  flush.console()
  # Specify the data in a list, for later shipment to JAGS:
  dataList = list(
    x = x ,
    y = y ,
    Nx = dim(x)[2] ,
    Ntotal = dim(x)[1] ,
    xPred = xPred # Data for predictions
  )
  
  #----------------------------------------------The Model-----------------------------------------
  
  modelString = "
  # Standardize the data:
  data {
    ym <- mean(y)
    ysd <- sd(y)
    for ( i in 1:Ntotal ) {
      zy[i] <- ( y[i] - ym ) / ysd
    }
    for ( j in 1:((Nx-1)) ) {
      xm[j] <- mean(x[,j])
      xsd[j] <- sd(x[,j])
      for ( i in 1:Ntotal ) {
        zx[i,j] <- ( x[i,j] - xm[j] ) / xsd[j]
      }
    }

    # Specify the priors for original beta parameters
    # Prior locations to reflect the expert information
    # pm2.5 in thousands taken for the model to execute
    # Sample mean = 0.09861321
    mu0 <- ym # Set to overall mean a priori based on the interpretation of constant term in regression
    mu[1] <- 0.01 # DEWP - Dew Point 
    mu[2] <- 0.01 # TEMP - Temperature
    mu[3] <- 0.01 # PRES - Pressure
    mu[4] <- 0.01 # lws - Cumulated Wind Speed
    mu[5] <- 0.01 # ls - Cumulated Hours of Snow
    mu[6] <- 0.01 # lr - Cumulated Hours of rain
    mu[7] <- 0.01 # cbwd - Combined wind speed

    # Prior variances to reflect the expert information
    # SampleVariance = 0.008473274
    Var0 <- 0.004 # according to sample variance
    Var[1] <- 0.1 # DEWP: No expert knowledge
    Var[2] <- 0.1 # TEMP: No expert knowledge
    Var[3] <- 0.1 # PRES: no expert knowledge
    Var[4] <- 0.1 # cbwd: No expert knowledge
    Var[5] <- 0.1 # lws: No expert knowledge
    Var[6] <- 0.1 # ls: No expert knowledge
    Var[7] <- 0.1 # lr: No expert knowledge
    
    
    # Compute corresponding prior means and variances for the standardised parameters
    muZ[1:(Nx-1)] <-  mu[1:((Nx-1))] * xsd[1:(Nx-1)] / ysd 
    muZ[7] <- mu[7] # PropertType: Categorical variable not standardised
    muZ0 <- (mu0 + (sum( mu[1:(Nx-1)] * xm[1:(Nx-1)] / xsd[1:(Nx-1)] ) + mu[7])*ysd - ym) / ysd 

    # Compute corresponding prior variances and variances for the standardised parameters
    VarZ[1:(Nx-1)] <- Var[1:(Nx-1)] * ( xsd[1:(Nx-1)]/ ysd )^2
    VarZ[7] <- Var[7] # PropertyType: Categorical variable is not standardised
    VarZ0 <- Var0 / (ysd^2)

  }
  # Specify the model for standardized data:
  model {
    for ( i in 1:Ntotal ) {
      zy[i] ~ dnorm( zbeta0 + sum( zbeta[1:(Nx-1)] * zx[i,1:(Nx-1)] ) + zbeta[7]*x[i,7] , 1/zsigma^2 )
    }

    # Priors vague on standardized scale:
    zbeta0 ~ dnorm( muZ0 , 1/VarZ0 )  
    for ( j in 1:Nx ) {
      zbeta[j] ~ dnorm( muZ[j] , 1/VarZ[j] )
    }

    zsigma ~ dgamma(0.001,0.01) # Standardised sigma distribution

    # Transform to original scale:
    beta[1:(Nx-1)] <- ( zbeta[1:(Nx-1)] / xsd[1:(Nx-1)] )*ysd
    beta[7] <- zbeta[7] # As cbwd was not standardised
    beta0 <- zbeta0*ysd  + ym - (sum( zbeta[1:(Nx-1)] * xm[1:(Nx-1)] / xsd[1:(Nx-1)] )+zbeta[7])*ysd
    sigma <- zsigma*ysd
    
    # Compute predictions at every step of the MCMC
    pred[1:2] <- beta0 + beta[1] * xPred[1:2,1] + beta[2] * xPred[1:2,2] + 
                 beta[4] * xPred[1:2,4] + beta[5] * xPred[1:2,5] + beta[6] * xPred[1:2,6] 

    # + beta[3] * xPred[1:2,3] + beta[7] * xPred[1:2,7]

  }
  " # close quote for modelString
  # Write out modelString to a text file
  writeLines( modelString , con="TEMPmodel.txt" )
  #-----------------------------------------------------------------------------

  #-----------------------------------------------------------------------------
  # RUN THE CHAINS
  parameters = c( "beta0" ,  "beta" ,  "sigma", 
                  "zbeta0" , "zbeta" , "zsigma" , "pred" )
  adaptSteps = 2000  # Number of steps to "tune" the samplers
  burnInSteps = 5000
  runJagsOut <- run.jags( method=runjagsMethod ,
                          model="TEMPmodel.txt" , 
                          monitor=parameters , 
                          data=dataList ,  
                          n.chains=nChains ,
                          adapt=adaptSteps ,
                          burnin=burnInSteps , 
                          sample=ceiling(numSavedSteps/nChains) ,
                          thin=thinSteps ,
                          summarise=FALSE ,
                          plots=FALSE )
  codaSamples = as.mcmc.list( runJagsOut )

  if ( !is.null(saveName) ) {
    save( codaSamples , file=paste(saveName,"Mcmc.Rdata",sep="") )
  }
  return( codaSamples )
} # end function

#===============================================================================

smryMCMC = function(  codaSamples , 
                      saveName=NULL) {
  summaryInfo = NULL
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  paramName = colnames(mcmcMat)
  for ( pName in paramName ) {
    summaryInfo = rbind( summaryInfo , summarizePost( mcmcMat[,pName] ) )
  }
  rownames(summaryInfo) = paramName
  # summaryInfo = rbind( summaryInfo , 
  #                      "log10(nu)" = summarizePost( log10(mcmcMat[,"nu"]) ) )
  if ( !is.null(saveName) ) {
    write.csv( summaryInfo , file=paste(saveName,"SummaryInfo.csv",sep="") )
  }
  
  return( summaryInfo)
}

#===============================================================================

plotMCMC = function( codaSamples , data , xName="x" , yName="y" ,
                     showCurve=FALSE ,  pairsPlot=FALSE ,
                     saveName=NULL , saveType="jpg" ) {
  # showCurve is TRUE or FALSE and indicates whether the posterior should
  #   be displayed as a histogram (by default) or by an approximate curve.
  # pairsPlot is TRUE or FALSE and indicates whether scatterplots of pairs
  #   of parameters should be displayed.
  #-----------------------------------------------------------------------------
  y = data[,yName]/1000
  x = data.matrix(data[,xName])
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  chainLength = NROW( mcmcMat )
  zbeta0 = mcmcMat[,"zbeta0"]
  zbeta  = mcmcMat[,grep("^zbeta$|^zbeta\\[",colnames(mcmcMat))]
  if ( ncol(x)==1 ) { zbeta = matrix( zbeta , ncol=1 ) }
  zsigma = mcmcMat[,"zsigma"]
  beta0 = mcmcMat[,"beta0"]
  beta  = mcmcMat[,grep("^beta$|^beta\\[",colnames(mcmcMat))]
  if ( ncol(x)==1 ) { beta = matrix( beta , ncol=1 ) }
  sigma = mcmcMat[,"sigma"]
  pred  = mcmcMat[,grep("^pred$|^pred\\[",colnames(mcmcMat))] # Added by Vishesh
  pred = pred*1000 # To be checked, converting back to the given pm2.5 values

  #-----------------------------------------------------------------------------
  # Compute R^2 for credible parameters:
  YcorX = cor( y , x ) # correlation of y with each x predictor
  Rsq = zbeta %*% matrix( YcorX , ncol=1 )
  #-----------------------------------------------------------------------------
  if ( pairsPlot ) {
    # Plot the parameters pairwise, to see correlations:
    openGraph()
    nPtToPlot = 1000
    plotIdx = floor(seq(1,chainLength,by=chainLength/nPtToPlot))
    panel.cor = function(x, y, digits=2, prefix="", cex.cor, ...) {
      usr = par("usr"); on.exit(par(usr))
      par(usr = c(0, 1, 0, 1))
      r = (cor(x, y))
      txt = format(c(r, 0.123456789), digits=digits)[1]
      txt = paste(prefix, txt, sep="")
      if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
      text(0.5, 0.5, txt, cex=1.25 ) # was cex=cex.cor*r
    }
    pairs( cbind( beta0 , beta , sigma )[plotIdx,] ,
           labels=c( "beta[0]" , 
                     paste0("beta[",1:ncol(beta),"]\n",xName), 
                     expression(sigma) ) , 
           lower.panel=panel.cor , col="skyblue" )
    if ( !is.null(saveName) ) {
      saveGraph( file=paste(saveName,"PostPairs",sep=""), type=saveType)
    }
  }
  #-----------------------------------------------------------------------------
  # Marginal histograms:
  
  decideOpenGraph = function( panelCount , saveName , finished=FALSE , 
                              nRow=2 , nCol=3 ) {
    # If finishing a set:
    if ( finished==TRUE ) {
      if ( !is.null(saveName) ) {
        saveGraph( file=paste0(saveName,ceiling((panelCount-1)/(nRow*nCol))), 
                   type=saveType)
      }
      panelCount = 1 # re-set panelCount
      return(panelCount)
    } else {
    # If this is first panel of a graph:
    if ( ( panelCount %% (nRow*nCol) ) == 1 ) {
      # If previous graph was open, save previous one:
      if ( panelCount>1 & !is.null(saveName) ) {
        saveGraph( file=paste0(saveName,(panelCount%/%(nRow*nCol))), 
                   type=saveType)
      }
      # Open new graph
      openGraph(width=nCol*7.0/3,height=nRow*2.0)
      layout( matrix( 1:(nRow*nCol) , nrow=nRow, byrow=TRUE ) )
      par( mar=c(4,4,2.5,0.5) , mgp=c(2.5,0.7,0) )
    }
    # Increment and return panel count:
    panelCount = panelCount+1
    return(panelCount)
    }
  }
  
  # Original scale:
  panelCount = 1
  panelCount = decideOpenGraph( panelCount , saveName=paste0(saveName,"PostMarg") )
  histInfo = plotPost( beta0 , cex.lab = 1.75 , showCurve=showCurve ,
                       xlab=bquote(beta[0]) , main="Intercept" )
  for ( bIdx in 1:ncol(beta) ) {
    panelCount = decideOpenGraph( panelCount , saveName=paste0(saveName,"PostMarg") )
    histInfo = plotPost( beta[,bIdx] , cex.lab = 1.75 , showCurve=showCurve ,
                         xlab=bquote(beta[.(bIdx)]) , main=xName[bIdx] )
  }
  panelCount = decideOpenGraph( panelCount , saveName=paste0(saveName,"PostMarg") )
  histInfo = plotPost( sigma , cex.lab = 1.75 , showCurve=showCurve ,
                       xlab=bquote(sigma) , main=paste("Scale") )
  panelCount = decideOpenGraph( panelCount , saveName=paste0(saveName,"PostMarg") )
  histInfo = plotPost( Rsq , cex.lab = 1.75 , showCurve=showCurve ,
                       xlab=bquote(R^2) , main=paste("Prop Var Accntd") )
  for ( pIdx in 1:ncol(pred) ) {
  panelCount = decideOpenGraph( panelCount , saveName=paste0(saveName,"PostMarg") )
  histInfo = plotPost( pred[,pIdx] , cex.lab = 1.75 , showCurve=showCurve ,
                       xlab=bquote(pred[.(pIdx)]) , main="Prediction"[pIdx]) 
  } # Modified by Vishesh #, finished=TRUE : Removed
  # Standardized scale:
  panelCount = 1
  panelCount = decideOpenGraph( panelCount , saveName=paste0(saveName,"PostMargZ") )
  histInfo = plotPost( zbeta0 , cex.lab = 1.75 , showCurve=showCurve ,
                       xlab=bquote(z*beta[0]) , main="Intercept" )
  for ( bIdx in 1:(ncol(beta)-1) ) {
    panelCount = decideOpenGraph( panelCount , saveName=paste0(saveName,"PostMargZ") )
    histInfo = plotPost( zbeta[,bIdx] , cex.lab = 1.75 , showCurve=showCurve ,
                         xlab=bquote(z*beta[.(bIdx)]) , main=xName[bIdx] )
  }
  panelCount = decideOpenGraph( panelCount , saveName=paste0(saveName,"PostMargZ") )
  histInfo = plotPost( zsigma , cex.lab = 1.75 , showCurve=showCurve ,
                       xlab=bquote(z*sigma) , main=paste("Scale") )
  panelCount = decideOpenGraph( panelCount , saveName=paste0(saveName,"PostMargZ") )
  histInfo = plotPost( Rsq , cex.lab = 1.75 , showCurve=showCurve ,
                       xlab=bquote(R^2) , main=paste("Prop Var Accntd") )
  panelCount = decideOpenGraph( panelCount , finished=TRUE , saveName=paste0(saveName,"PostMargZ") )
  
  #-----------------------------------------------------------------------------
}
#===============================================================================
