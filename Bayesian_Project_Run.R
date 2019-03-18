# Execution Script
# Vishesh Jain, Abheer Sharma, Deepanshu Jain 

#--------------------------------Data Read---------------------------------

beijingpm = read.csv("~/PRSA_data_2010.1.1-2014.12.31.csv")

#---------------------------Data Cleaning--------------------------------

# Removing missing values as they are only 4.7%

beijingpm <- beijingpm[!(is.na(beijingpm$pm2.5)),]

# Time variables - Year, Month and Hour are not required
# Also, S.No is removed from the dataset

beijingpm <- beijingpm[,6:13]

## Reorder the columns in the dataset

beijingpm <- beijingpm[,c(2,3,4,6,7,8,5,1)]

#----------------------------------------Data Transformation-----------------------------------------

# We will set the levels of the categorical variable as -2,-1,1,2
# This step is required as the variables of linear regression require numerical values

beijingpm$cbwd <- as.character(beijingpm$cbwd)
beijingpm[beijingpm$cbwd=='NW',]$cbwd <- -2
beijingpm[beijingpm$cbwd=='cv',]$cbwd <- -1
beijingpm[beijingpm$cbwd=='NE',]$cbwd <- 1
beijingpm[beijingpm$cbwd=='SE',]$cbwd <- 2

beijingpm$cbwd <- as.factor(beijingpm$cbwd)

#---------------------------------------Specification for MCMC---------------------------------------

# beijingpm = beijingpm[sample(nrow(beijingpm),500),] # Sample was taken to test the chains

yName = "pm2.5" ; xName = c("DEWP",	"TEMP",	"PRES", "Iws", "Is", "Ir","cbwd")
fileNameRoot = "Project"
numSavedSteps = 20000 ; thinSteps=2

graphFileType = "jpg" 

##---------------------------------Packages & Source Reqruired---------------------------------------- 

# Load the relevant model into R's working memory:

source("~/Bayesian_Project_Model.r")

##---------------------------------------Generating Chains------------------------------------------ 

# Generate the MCMC chain:
startTime = proc.time()
xPred = matrix(c(-11,20,-5,15,1025,1000,1.57,2.5,1,0,0,0,-2,1), nrow = 2, ncol = 7)
colnames(xPred) = c("beta[1]","beta[2]","beta[3]","beta[4]","beta[5]","beta[6]","beta[7]")
mcmcCoda = genMCMC( data=beijingpm , xName=xName , yName=yName , 
                    numSavedSteps=numSavedSteps , thinSteps=thinSteps , 
                    saveName=fileNameRoot , xPred = xPred )

stopTime = proc.time()
duration = stopTime - startTime
show(duration) # To display the duration of the process

##-----------------------------------Display Diagnostics--------------------------------------

# Display diagnostics of chain, for specified parameters:
parameterNames = varnames(mcmcCoda) # get all parameter names
for ( parName in parameterNames ) {
  diagMCMC( codaObject=mcmcCoda , parName=parName , 
            saveName=fileNameRoot , saveType=graphFileType )
}

##-------------------------------------Summary Statistics-------------------------------------- 

# Get summary statistics of chain:

summaryInfo = smryMCMC( mcmcCoda , 
                        saveName=fileNameRoot  )
show(summaryInfo)

##--------------------------------Display posterior information---------------------------------

plotMCMC( mcmcCoda , data=beijingpm , xName=xName , yName=yName , 
          pairsPlot=TRUE , showCurve=FALSE ,
          saveName=fileNameRoot , saveType="jpg" )

##################################--------END---------##########################################

