# config.R - to allow users to set all relevant parameters
# in one location. This will be sourced first during the run
# and parameters will be updated from these inputs

# set species code (BND = bottlenose dolphin, GS = grey seal, HP = harbour porpoise, HS = harbour seal, MW = minke whale)
spec <- 'HS'
# set threshold size for demographic stochasticity
threshold <- 500
# set population size
pmean <- 4658
Surv <- rep(0,18)
# set calf/pup survival
Surv[1] <- 0.6
# set juvenile survival
Surv[7] <- 0.82
# set adult survival
Surv[13] <- 0.85
# set fecundity value
Fertility <- 0.95
# set age at which a calf or pup becomes independent of its mother
age1 <- 1
#set age at which an average female gives birth to her first calf
age2 <- 4
# input proportion of animals in each of the vulnerable sub-population(s), default is that the entire population is vulnerable
vulnmean <- c(1)

nvulnmean <- length(vulnmean)
newvulnmean <- vulnmean
if(sum(vulnmean)== 1){newvulnmean <- newvulnmean} else {newvulnmean[nvulnmean+1] <- 1 - sum(newvulnmean)}

# set number of years on which piling will occur. Set this to zero if there is no piling
pile_years <- 1
# specify name of csv file that contains information on days on which piling will occur between the quotation marks
# read csv file with schedule of piling activities. 
pile <- read.csv(file = '../data/MultPilingOpsMultYears.csv', header = TRUE) ## XXX strip out date column
piling.file <- "../data/MultPilingOpsMultYears.csv"

# input number of piling operations to be modelled
pilesx1 <- 1
vulnpile <- matrix(0, nrow = nvulnmean, ncol = pilesx1)

# indicate which operations will affect each vulnerable sub-population
# in this example, there is one vulnerable sub-population that is affected by operations 1 & 2
vulnpile[1, ] <- c(1, 1)

# "seasons" determines whether or not the number of animals that are likely to be disturbed each day (NDt) 
# and the number that may experience PTS (NPt) vary by season. seasons = 1, no seasonal variation is the default value; 
# if seasons > 1, remember that the simulated year starts in June, for all species except GS, when it starts in October!
# if seasons = 4, 4 there are seasons in a year. summer = June, July, August; autumn = Sept, Oct, Nov; winter = Dec, Jan, Feb; spring = March, April, May                          
# if seasons = 2 (ie just summer and winter), we will actually require 3 breaks, assuming "summer" = May - October (ie months when water is warmest)

seasons <- 1
# I THINK WE MAY HAVE TO GIVE UP ON EXPLICITLY MODELLING THESE VARIATIONS WITHIN A PILING OPERATION. IT’S SIMPLER TO DIVIDE EACH PILING OPERATION INTO SEPARATE SUB-OPERATIONS FOR EACH SEASON.

# input number of animals predicted to experience significant disturbance on one day of piling for each operation
# in this example, there are two operations each of which disturbs 622 animals
daily_NDt <- c(622,622)

# now do the same for the number of animals predicted to experience PTS on one day of piling
daily_NPt <- c(0, 0)

# input number of days of "residual" disturbance. DEFAULT IS 1
days <- 1

# decide if PTS can occur on any day (default) or only on the first occasion that an individual is disturbed by 
# change Day1 to TRUE if you want animals to be only vulnerable to PTS on the first day they are disturbed
Day1 <- FALSE

# load(file = paste(spec, piling.file, sep = ‘’))
