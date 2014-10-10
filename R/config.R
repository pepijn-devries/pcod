# config.R - to allow users to set all relevant parameters
# in one location. This will be sourced first during the run
# and parameters will be updated from these inputs

# Species
spec <- 'HS'
spec <- tolower(spec)
# Management Unit

# Population Size - value from the IAMMAWG MU report should be used for the MU being modelled
pmean  <- 4568

# proportion of population that are females
propfemale <- 0.5

# Pup Survival Rate
pupSurv <- 0.6

# Juvenile Survival Rate
juvSurv <- 0.82

# Adult Survival Rate
adSurv <- 0.85

# Age at first birth

# Fertility Rate
Fertility <- 0.95

# Vulnerable population - input proportions of population in vulnerable sub-populations. Default is
# entire population in one sub-population is vulnerable
vulnmean <- c(1.0)

# Number of piling years - default is 1. Set to 0 if no piling.
pile_years <- 1

# Input the number of piling operations to be modelled
pilesx1 <- 3

# daily_NDt - stores estimates of the number of animals of the species being modelled that may experience disturbance
# that is likley to impair an individual's ability to survive, breed, reproduce, or raise young, or that is likely to 
# result in that individual being displaced from an area for a longer period than normal during one day for each operation
# 
# For example, assume there are three piling operations, and there is a fixed number of animals that will experience PTS
# 30 in the first piling operation, 60 in the second, and 80 in the third
# Then that would look like: numDT <- c(30, 60, 80)
#
# default is 60 animals at each operation. These must be integer values
numDT <- c(60, 60, 60)
if (!isTRUE(all(numDT == floor(numDT)))) stop("'numDT' must only contain integer values")

# daily_NPt - stores estimates of the number of animals of the species being modelled that may experience 
# a permanent shift in the threshold for hearing (PTS) during one day for each operation
# default is 2 animals at each operation
numPT <- c(2, 2, 2)
if (!isTRUE(all(numPT == floor(numPT)))) stop("'numPT' must only contain integer values")

# Residual Days of Disturbance
days <- 1

# Determine if PTS can occur on any day, or only on the first occasion that an individual is disturbed
# Change Day1 to TRUE if you want animals to be vulnerable to PTS only on the first day that they
# are disturbed
Day1 <- FALSE

# Collisions - the number of animals predicted to be killed each year as a result of collisions with tidal
# energy arrays
NCollisions <- 10

# years - the number of years for the simulation; default is 25
years <- 25

# Set the number of times you want to run the simulation. 
# Higher values equals less uncertainty, but longer run times
# default is 500
nboot <- 5#00

# proportion vulnerable - Allowing a subset of disturbed animals to experience residual disturbance