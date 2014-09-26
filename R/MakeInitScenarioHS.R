spec <- 'HS'
propfemale <- 0.5

                 
# set threshold value of total population size for inclusion of demographic stochasticity
threshold <- 500

# iPCoD PROTOCOL STEP 2

# pmean <- 4568 # population size value from the IAMMAWG MU report should be used for the MU being modelled
pmean <- round(pmean*propfemale)






Surv <- rep(0, 18)

# iPCoD PROTOCOL STEP 2
# INPUT DEMOGRAPHIC RATES FROM HARWOOD & KING (2014) HERE

  Surv[c(1, 7, 13)] <- c(pupSurv, juvSurv, adSurv)

# age1 = age at which a calf becomes independent from its mother, default is 1 
  age1 <- 1

# age2  = age at which a female give birth to her first calf
  age2 <- 4

  Fert <- rep(0, 6); Fert[1] <- Fertility/2.0

# determine initial stable age structure for population from Leslie matrix

  L                  <- array(0, dim = c(10, 10))
  L[1, age2]          <- Fert[1] * Surv[7]
  L[1, (age2 + 1):10] <- Fert[1] * Surv[13]
  index3             <- array(c(2:(age1+1), 1:age1), dim = c( age1 ,2) )
  L[index3]            <- Surv[1]
  mature             <- ifelse(age2 == 9, 9, age2 + 1)
  index1             <- array(c((mature + 1):10, (mature):9), dim = c((10 - mature), 2))
  L[index1]          <- Surv[13]
  L[10, 10]          <- Surv[13]
  index2             <- array(c((age1+2):(age2 + 1), (age1+1):age2), dim = c((age2 - age1), 2))
  L[index2]          <- Surv[7]
  
  ev              <- eigen(L)
  #symmetric = TRUE removed because is causes weird values!
  # ev$val[1] contains growth rate of the undisturbed population
  
  evCmplx         <- ev # in case we need this later on
# age_structure holds the proportion of animals in each age class for a stable age structure  
  age_structure   <- as.numeric(ev$vec[1:10, 1] / (sum(ev$vec[1:10, 1])))


# set % level of environmental variation for calf survival, juvenile survival and fertility from expert elicitation results
EnvStoch <- c(30, 30, 25)
# calculate standard deviations that will result in lower 99% CL matching environmental variation
SDs <- EnvStoch/(3*100)
mu <- Surv[c(1,7)]
mu[3] <- Fertility
# calculate parameters of Beta distribution to match these values for mu and SD
a <- 0
b <- 0
  for (ibeta in 1:3){
    b[ibeta] <- mu[ibeta]*(1-mu[ibeta])^2/SDs[ibeta]^2 + mu[ibeta] - 1
    a[ibeta] <- mu[ibeta]*b[ibeta]/(1-mu[ibeta])
    }


  
# iPCoD PROTOCOL STEP 3

# INPUT PROPORTIONS OF POPULATION IN VULNERABLE SUB-POPULATIONS (DEFAULT: entire population in one sub-population)
 
#   vulnmean <- c(1.0) 

  nvulnmean <- length(vulnmean)
# newvulnmean includes proportion of animals in undisturbed remainder of population, if there is one!
  newvulnmean <- vulnmean
  if(sum(vulnmean)== 1){newvulnmean <- newvulnmean} else {newvulnmean[nvulnmean+1] <- 1 - sum(newvulnmean)}

# SET pile_years TO ZERO IF THERE IS NO PILING
# pile_years <- 1

if (pile_years > 0) {

# iPCoD PROTOCOL STEP 4

# read csv file with schedule of piling activities. YOU MAY NEED TO CHANGE THE FILE NAME HERE

pile <- read.csv(file = '../data/MultPilingOpsMultYears.csv', header = TRUE) ## XXX strip out date column

# removes day labels from first column of csv file
widx <- which(colnames(pile) %in% c('Date', 'DayOfYear'))
if(length(widx) > 0){pile <- pile[,-widx]}
# In case we want to get it programmatically, this works:
library(stringr)
yvec   <- str_match(colnames(pile), 'Operation')
npiles <- length(yvec[!is.na(yvec)])
pilesx <- which(str_match(colnames(pile), 'Operation') %in% 'Operation')
pile  <- pile[, pilesx]

# template for  MultPilingOpsMultYears.csv should ensure that number of rows is an exact multiple of 365
pile_years <- nrow(pile) / 365 
# yearvec indicates year number for each day in pile
yearvec <- rep(1:pile_years, each = 365)


# iPCoD PROTOCOL STEP 5

# input number of piling operations to be modelled (3 in this case)
# pilesx1 <- 3

if(pilesx1 != ncol(pile)){stop('Number of Piling Operations do not match')}
# vulnpile is a matrix indicating which columns of pile are to be combined to provide piling information for each vulnerable sub-population
vulnpile <- matrix(0, nrow = nvulnmean, ncol = pilesx1)

# iPCoD PROTOCOL STEP 6

# indicate which operations will affect each vulnerable sub-population
# in this case there is one vulnerable sub-population that is only affected by operations 1 & 2

vulnpile[1, ] <- c(1, 1, 1)

# repeat this for each sub-population. 

#if(length(vulnpile) != ncol(pile)){stop('Number of Piling Operations do not match')}##### CHANGE THIS FROM length(vulnpile) to ncol(vulnpile)
if(ncol(vulnpile) != ncol(pile)){stop('Number of Piling Operations do not match')}
# indicate which operations will affect vulnerable sub-population 2 etc. - not needed in this case 



# "seasons" determines whether or not the number of animals that are likely to be disturbed each day (NDt) 
# and the number that may experience PTS (NPt) vary by season. seasons = 1, no seasonal variation is the default value; 
# if seasons > 1, remember that the simulated year starts in June, for all species except GS, when it starts in October!
# if seasons = 4, 4 there are seasons in a year. summer = June, July, August; autumn = Sept, Oct, Nov; winter = Dec, Jan, Feb; spring = March, April, May                          
# if seasons = 2 (ie just summer and winter), we will actually require 3 breaks, assuming "summer" = May - October (ie months when water is warmest)

seasons <- 1
inputindex <- seasons
if(seasons == 2) {
  seasons <- seasons + 1
  breaks1 <- c(152,182,31)
} else {
    breaks1 <- c(91,92,91,91)
}

if(seasons == 1){s_index <- c(365)} else {s_index <- breaks1}



# these matrices will hold the number of animals predicted to experience PTS and disturbance during one day for each piling operation
# each row contains the values for a particular season
# each column contains the values for a particular piling operation
daily_NPt <- daily_NDt <- matrix(0, nrow = seasons, ncol = npiles)

# iPCoD PROTOCOL STEP 6
# input number of animals predicted to experience significant disturbance on one day of piling for each operation
# the first 1:inputindex values in the string give the number of animals that will be affected in each season
# for piling operation 1, the next (inputindex+1):(2*inputindex) values give the number of animals
#  that will be affected in each season for piling operation 2, and so on

daily_NDt[1:inputindex,] <- numDT

# now do the same for the number of animals predicted to experience PTS on one day of piling

daily_NPt[1:inputindex,] <- numPT


if (inputindex == 2) {
  daily_NPt[3,] <- daily_NPt[1,]
  daily_NDt[3,] <- daily_NDt[1,]
} 
  
daily_NDt <- daily_NDt*propfemale
daily_NPt <- daily_NPt*propfemale


pile <- data.frame(pile, pbdays = rowSums(pile), pvec = rowSums(pile))
# pile$pbdays indicates number of piling events on a particular day, pile$pvec is a flag indicating whether or not piling has occurred
pile$pvec[pile$pvec > 1] <- 1 

NDt <- matrix(0,nrow(pile), npiles, byrow = TRUE)
NPt <- matrix(0,nrow(pile), npiles, byrow = TRUE)
                   
 for (j in 1:npiles){
    NDt[,j]<- rep(rep(c(daily_NDt[1:seasons,j]),c(s_index[1:seasons])),pile_years)
    NPt[,j]<- rep(rep(c(daily_NPt[1:seasons,j]),c(s_index[1:seasons])),pile_years)    
    }



#creates a column for each vulnerable sub-population indicating the days on which there is piling that will affect it
for(i in 1:nrow(vulnpile)){
  pvec <- as.matrix(pile[, pilesx]) %*% as.matrix(vulnpile[i, ])
  pvec[pvec > 1] <- 1
  colnames(pvec) <- paste('vuln', i, 'pvec', sep = '')
  pile <- cbind(pile, pvec)              
}


# iPCoD PROTOCOL STEP 7

# input number of days of "residual" disturbance. DEFAULT IS 1

# days <- 1

# iPCoD PROTOCOL STEP 8

# decide if PTS can occur on any day (default) or only on the first occasion that an individual is disturbed by 
# change Day1 to TRUE if you want animals to be only vulnerable to PTS on the first day they are disturbed

# Day1 = TRUE # Flag for the different PTS model 

}

# iPCoD PROTOCOL STEP 9

# input number of animals predicted to be killed each year as a result of collisions with tidal energy arrays


# NCollisions <- 0

NCollisions <- NCollisions*propfemale                                                             


#number of years for simulation
# years <- 25


