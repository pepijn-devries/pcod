#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ stickEvals ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~ create curves from parameters using a broken stick
#~ curves start at 1, break at some X, decrease to some other greater X, then plateau i.e. a piecewise linear, decreasing "sigmoidal" shape
#~ CRD 31/01/2013
#~ Args:
#~ evals: the points you intend to evaluate. In this application thought to be 0-365 days, however anything non-negative should give a "sensible" result
#~ parameters: numeric, length 3. 
#         First parameter is the X breakpoint where Y becomes <1. Second is the X breakpoint where Y plateaus. 
#         Third is the plateau in 0-1. Note the nature of the data collection means the input here is the reduction from 1, not the value of Y itself

stickEvals<- function(evals, parameters){
  # in case they are coming in as a list or the like
  parameters<- as.numeric(parameters)
  
  # convert from the reduction in Y, to Y (which is 0-1 bounded as a probability)
  parameters[3]<- (1-parameters[3])
  
  # create a vector containing just the plateau/lower value of Y
  outputEvals<- rep(parameters[3], length(evals))
  
  # if to left of first breakpoint, Y=1 i.e. no reduction in Y
  outputEvals<- ifelse(evals<parameters[1], 1, outputEvals)
  
  # if between the break points, then straightline transition from 1 to plateau at the second break point
  outputEvals<- ifelse(evals>=parameters[1]&evals<parameters[2], 1-(evals-parameters[1])*(1-parameters[3])/((parameters[2]-parameters[1])), outputEvals)
  
  # outout the predicted Y values only. To plot etc this corresponds to the input "evals"
  list(x=evals, predictedY=outputEvals)
  
}

      
f.experts.survival_and_fertility <- function(spec){
# dd3 contains resampled values for the number of days required to cause moderate and severe disturbance for adults, 
# dd2 contains the equivalent values for juveniles, 
# dd1 contains the equivalent values for dependents (pups or calves)
 
  dd3 <- matureDays[sample(nrow(matureDays), 1), ] # samples rows at random
  dd3[1]<-1
  dd3[2:3] <- round(dd3[2:3], 0)
  disturb_post3 <- stickEvals(0:365, dd3[c(2, 3, 1)])
  
  if (spec == 'MW') {
    dd2 <- c(1,365,365) } else  {# samples rows at random 
      dd2 <- juveDays[sample(nrow(juveDays), 1), ] 
      }  
      
  dd2[2:3] <- round(dd2[2:3], 0)
  disturb_post2 <- stickEvals(0:365, dd2[c(2, 3, 1)])
    
  if (spec == 'GS'| spec == 'MW')  {
    dd1 <- dd2  } else {
    dd1 <- dependentDays[sample(nrow(dependentDays), 1), ] 
    } # samples rows at random   
  
  dd1[2:3] <- round(dd1[2:3], 0)
  disturb_post1 <- stickEvals(0:365, dd1[c(2, 3, 1)])
  
 
  # choose demographic effects of PTS from posteriors
  # 1. survival as a result of PTS at 2-10 kHz Mature Females
  MFsurv_2_10_probs <- output_mature$Survival_2_10kHz$density
  MFsurv_2_10_X <- output_mature$Survival_2_10kHz$X
  sampleSize <- 1
  MFsurv_2_10_randomDraw <- sample(MFsurv_2_10_X, sampleSize, MFsurv_2_10_probs, replace = T)
  
  # 2. fecundity as a result of PTS at 2-10 kHz
  MFrepro_2_10_probs <- output_mature$Reproduction_2_10kHz$density
  MFrepro_2_10_X <- output_mature$Reproduction_2_10kHz$X
  MFrepro_2_10_randomDraw <- sample(MFrepro_2_10_X, sampleSize, MFrepro_2_10_probs, replace = T)
  
  # 3. Juveniles survival as a result of PTS at 2-10 kHz
  Jsurv_2_10_probs <- output_juvenile$Survival_2_10kHz$density
  Jsurv_2_10_X <- output_juvenile$Survival_2_10kHz$X
  Jsurv_2_10_randomDraw <- sample(Jsurv_2_10_X, sampleSize, Jsurv_2_10_probs, replace = T)
  
  # Dependents survival as a result of PTS at 2-10 kHz
  Dsurv_2_10_probs <- output_dependent$Survival_2_10kHz$density
  Dsurv_2_10_X <- output_dependent$Survival_2_10kHz$X
  Dsurv_2_10_randomDraw <- sample(Dsurv_2_10_X, sampleSize, Dsurv_2_10_probs, replace = T)
  
  
  
  #  Populate Surv[8:12] with the survival multipliers for each juvenile disturbance class
  # we assume juvenile minke whales are not affected by disturbance
  if (spec == 'MW') {
      Surv[8] <- 1
      Surv[9] <- 1  } else {
        Surv[8]  <- ((disturb_post2$predictedY[1] - disturb_post2$predictedY[365])/2) + disturb_post2$predictedY[365]
        Surv[9]  <- disturb_post2$predictedY[365]
        }
  Surv[10] <- 1-Jsurv_2_10_randomDraw
  Surv[11] <- Surv[8] * Surv[10]
  Surv[12] <- Surv[9] * Surv[10]

 
    # Populate Surv[1:6] with the survival multipliers for each dependent (pup/calf) disturbance class
    # we assume disturbance and PTS has the same effect on grey seal pups as it does on juveniles
  if (spec == 'GS') {
    Surv[2:6] <- Surv[8:12]} else { 
    # we assume minke whale calves are not affected by disturbance   
      if (spec == 'MW') {
        Surv[2] <- 1
        Surv[3] <- 1  } else {
          Surv[2] <- ((disturb_post3$predictedY[1] - disturb_post3$predictedY[365])/2) + disturb_post3$predictedY[365]
          Surv[3] <- disturb_post3$predictedY[365] 
          }
    Surv[4] <- 1-Dsurv_2_10_randomDraw
    Surv[5] <- Surv[2] * Surv[4]
    Surv[6] <- Surv[3] * Surv[4]
    }
  
  #  fill Surv[(14:18] with the survival multipliers for each adult disturbance class
  Surv[14] <- 1.0
  Surv[15] <- 1.0
  Surv[16] <- 1-MFsurv_2_10_randomDraw
  Surv[17] <- Surv[14] * Surv[16]
  Surv[18] <- Surv[15] * Surv[16]
  
  #  fill Fert[2:6] with the fertility multipliers values for each adult disturbance class
  Fert[2] <- 0.5
  Fert[3] <- 0
  Fert[4] <- 1-MFrepro_2_10_randomDraw
  Fert[5] <- (1-MFrepro_2_10_randomDraw) * 0.5
  Fert[6] <- 0
  
  # Return them
  list(Surv = Surv, Fert = Fert, dd1 = dd1, dd2 = dd2, dd3 = dd3)
}
      

f.current.survival.values<- function(dist_ages, age1, age2, a, b){
  AnnSurv <- rep(1, 60)
  #determines survival values for current year, accounting for environmental stochasticity
  calfsurv    <- rbeta(1,a[1],b[1])
  for(i1 in 1:age1){
    AnnSurv[dist_ages[i1]] <- calfsurv
    AnnSurv[(dist_ages[i1]+1):(dist_ages[i1]+5)] <- calfsurv*Surv[2:6]
  }
  
  juvsurv    <- rbeta(1,a[2],b[2])
  for(i2 in (age1+1):age2){
    AnnSurv[dist_ages[i2]] <- juvsurv
    AnnSurv[(dist_ages[i2]+1):(dist_ages[i2]+5)] <- juvsurv*Surv[8:12]
  }
  
  adsurv     <- Surv[13]
  for (i3 in (age2+1):10) {
    AnnSurv[dist_ages[i3]] <- adsurv
    AnnSurv[(dist_ages[i3]+1):(dist_ages[i3]+5)] <- adsurv*Surv[14:18]
  }
  
  
  return(AnnSurv)
}



f.disturbance_and_PTS <- function(pile, vuln, days, p_disturb, p_pts){ 
  
#   f.disturbance_and_PTS(pile[yrs, max(pilesx) + 2  + i], vuln[i], days, p_disturb[yrs, i], p_pts[yrs, i])
#   pile <- pile[yrs, max(pilesx) + 2  + i]
#   vuln <- vuln[i]
#   days <- days
#   p_disturb <- p_disturb[yrs, i]
#   p_pts <- p_pts[yrs, i]
  
  
  dist_year      <- matrix(pile, vuln, 365, byrow = TRUE)
  pts_year       <- matrix(0, vuln, 365, byrow = TRUE)
  dist_year1     <- apply(dist_year, 1, function(x){
    
    x       <- as.numeric(x)
    idx     <- which(x == 1)
    x[idx]  <- rbinom(length(idx), 1, p_disturb[idx])
    x[!idx] <- 0
    return(x)
    
  })
  
  dist_year1     <- t(dist_year1)
  
  dist_year2 <- apply(dist_year1, 1, function(x){
    x        <- as.numeric(x)
    patrn    <- rep(1, days + 1)
    occur1   <- function(patrn, x) 
    { 
      m <- length(patrn) 
      n <- length(x) 
      candidate <- seq.int(length=n-m+1) 
      for (i in seq.int(length=m)) { 
        candidate <- candidate[patrn[i] == x[candidate + i - 1]] 
      } 
      candidate 
    } 
    
    idx <- occur1(patrn, x) + 1
    
    for(i in idx){
      if(days == 1){
        if(x[i] == 1 & x[(i - 1)] == 1){x[i] <- 0}
      } else {
        if(x[i] == 1 & x[((i + 1) - days)] == 1){x[i] <- 0; x[i + 1] <- 0}
      }
    }
    return(x)
  })
  dist_year2     <- t(dist_year2)
  pidx           <- which(dist_year2 == 1)
  pts_year[pidx] <- rbinom(length(pidx), 1, p_pts)
  dist_year      <- cbind(dist_year2, rep(0, nrow(dist_year2)), rep(0, nrow(dist_year2)))
  
  list(dist_year = dist_year, pts_year = pts_year)
}



f.disturbance_and_PTS.Model2 <- function(pile, vuln, p_disturb, p_pts){ 
  # Here PTS is only possible on the first day of disturbance.
  #pile[yrs, max(pilesx) + 2  + i], vuln[i], p_disturb[yrs, i], p_pts[yrs, i] # for testing
#   pile <- pile[yrs, max(pilesx) + 2  + i]
#   vuln <- vuln[i] 
#   p_disturb <- p_disturb[yrs, i]
#   p_pts <- p_pts[yrs, i]
  
  dist_year      <- matrix(pile, vuln, 365, byrow = TRUE)
  pts_year       <- matrix(0, vuln, 365, byrow = TRUE)
  dist_year1     <- apply(dist_year, 1, function(x){
    
    x       <- as.numeric(x)                         
    idx     <- which(x == 1)
    x[idx]  <- rbinom(length(idx), 1, p_disturb[idx])      
    x[!idx] <- 0
    return(x)
    
  })
  
  dist_year1     <- t(dist_year1)
  
  
  idx            <- which(dist_year1 == 1, arr.ind = TRUE)
  idxn           <- idx[order(idx[, 'row'], idx[, 'col']),]
  didx           <- duplicated(idxn[,'row'])
  idxn           <- idxn[!didx,]
  pts_year[idxn] <- rbinom(nrow(idxn), 1, p_pts)
  dist_year      <- cbind(dist_year1, rep(0, nrow(dist_year1)), rep(0, nrow(dist_year1)))
  
  list(dist_year = dist_year, pts_year = pts_year)
}



f.disturbance.classes <- function(dist_year, pts_year, N, dist_ages, days, dds){
  
  # the following lines account for residual days of disturbance
  flag <- length(which(dist_year > 0))
  if(flag != 0){
    didx <- which(dist_year == 1, arr.ind = TRUE)
    rdx <- data.frame(row = didx[, 1], fromcol = didx[, 2], tocol = didx[, 2] + days)
    colnames(rdx) <- c('row', 'fromcol', 'tocol')
    for(ii in 1:nrow(rdx)){
      dist_year[rdx[ii, 1],rdx[ii, 'fromcol']:rdx[ii, 'tocol']] <- 1
     }
    }
  dist_year[,366] <-  rowSums(dist_year[, 1:365], na.rm = TRUE)
  dist_year[,367] <-  ifelse(rowSums(pts_year[, 1:365], na.rm = TRUE) > 0, 1, 0)

  individuals <- sum(N[1:60])
  dist_sum <- matrix(0, nrow = individuals, ncol = 2)
    # dist_sum[,1] will contain total number of days of disturbance experienced by each individual in N
  # dist_sum[,2] will flag whether or not an individual experienced PTS

  if (nrow(dist_year) == individuals) {
    dist_sum[, 1] <- dist_year[,366]
    dist_sum[, 2] <- dist_year[,367]} else

 # if "individuals" > vuln, replicate rows of dist_year until dist_sum is full
      {times <- individuals%/%nrow(dist_year)
      size <- times*nrow(dist_year)
      remainder <- individuals - size
      if (times > 0) {
      dist_sum[1:size, 1] <- rep(dist_year[,366],times)
      dist_sum[1:size, 2] <- rep(dist_year[,367],times)}
        if (remainder > 0) {
        dist_sum[(size+1):individuals,1] <- dist_year[1:remainder, 366]
        dist_sum[(size+1):individuals,2] <- dist_year[1:remainder, 367]
        }
      }


# create array that will be used to assign the individuals in dist_sum into the different age classes. Filled with 999s because temp[,?,] can vary in length


  temp <- array (999, dim = c(10,sum(N[1:60]), 2))


  index1 <- 1

  for (l in 1:10){
    total <- N[dist_ages[l]] + N[(dist_ages[l]+3)]
    if(total != 0){
    index2 <- index1 + total - 1
    temp[l, 1:total, 1:2] <- dist_sum[index1:index2, 1:2]
    index1 <- index2 + 1

     begin <- 1
     
# check to deal with situations where no animals have PTS already, ie N[dist_ages[l]+3] = 0
    if(N[dist_ages[l]+3] != 0){
    begin <- N[(dist_ages[l]+3)]+1
    tempPTS <- temp[l,1:N[(dist_ages[l]+3)],1]

   # determine whether animals that already have PTS are disturbed
    N[dist_ages[l]+3] <- length(tempPTS[tempPTS >= 0 & tempPTS < dds[l,2]])
    N[dist_ages[l]+4] <- length(tempPTS[tempPTS >= dds[l,2] & tempPTS < dds[l,3]])
    N[dist_ages[l]+5] <- length(tempPTS[tempPTS >= dds[l,3]& tempPTS < 999])
    }
    if (N[dist_ages[l]] != 0) {
    temp_no_prior_PTS <- matrix(data = temp[l,begin:total,1:2], nrow = N[dist_ages[l]], ncol = 2)
    temp_no_prior_PTS <- matrix(data = temp_no_prior_PTS[order(temp_no_prior_PTS[,2]),], nrow = N[dist_ages[l]], ncol = 2)
    no_newPTS <- length(which(temp_no_prior_PTS[,2] == 0))


 if(no_newPTS == N[dist_ages[l]]) {newPTS <- 0} else {
      newPTS <- temp_no_prior_PTS[(no_newPTS+1):N[dist_ages[l]],1]

# determine whether animals that have PTS are disturbed
      N[dist_ages[l]+3] <- N[dist_ages[l]+3] + length(newPTS[newPTS >= 0 & newPTS < dds[l,2]])
      N[dist_ages[l]+4] <- N[dist_ages[l]+4] + length(newPTS[newPTS >= dds[l,2] & newPTS < dds[l,3]])
      N[dist_ages[l]+5] <- N[dist_ages[l]+5] + length(newPTS[newPTS >= dds[l,3] & newPTS < 999])
      }

# determine whether animals that do not have PTS are disturbed
    tempnoPTS <- temp_no_prior_PTS[1:no_newPTS,1]
    N[dist_ages[l]] <- length(tempnoPTS[tempnoPTS < dds[l,2]])
    N[dist_ages[l]+1] <- length(tempnoPTS[tempnoPTS >= dds[l,2] & tempnoPTS < dds[l,3]])
    N[dist_ages[l]+2] <- length(tempnoPTS[tempnoPTS >= dds[l,3] & tempnoPTS < 999])
    }
    } # end if for cases where total != 0
   
    }

  list(N = N)
  
}


f.next.year <- function(N1, N2, dist_ages, pmean, AnnFert, threshold, age1, age2){
#       tmp                <- next.year(Ndistvuln[,k , t], Ndistvuln[,k,t+1], dist_ages, pmean, AnnFert, threshold)

# apply survival rates to all age and disturbance classes
#  N1 <- round(N1*AnnSurv,0)
  N2temp <- rep(c(0), 60)
  if (pmean > threshold) {N2temp <- round(N1*AnnSurv,0)} else {
    for (irand in 1:60) {
        N2temp[irand] <- rbinom(1, N1[irand], AnnSurv[irand])}
        }

  N2 <- rep(c(0),60)

  
# determine pup production in next year
  tempAd <- 0
  iclass <- 0  
  pups <- 0
  if (pmean > threshold)
  {  for (iFert1 in 0:5){
      iclass[1:(10-age2+1)] <- dist_ages[age2:10] + iFert1
      tempAd <- sum(N2temp[iclass])
      pups <- round(tempAd*AnnFert[iFert1+1], 0)
      N2[1] <- N2[1] + pups
      }
  }    else

  {   for (iFert2 in 0:5){
      iclass[1:(10-age2+1)] <- dist_ages[age2:10] + iFert2
      tempAd <- sum(N2temp[iclass])
      pups <- rbinom(1,tempAd,AnnFert[iFert2+1])
      N2[1] <- N2[1] + pups
      }
  } 

# combine all animals that have experienced PTS into one disturbance class, for each age, and do the same for those 
#   that didn't experience PTS
  for (k1 in dist_ages){
    N2temp[k1] <- sum(N2temp[k1:(k1+2)])
    N2temp[(k1+3)] <- sum(N2temp[(k1+3):(k1+5)])
    }

# increase the ages of animals in age classes 1-9 by one year    
  for (k2 in 1:9){
    N2[dist_ages[k2+1]] <- N2temp[dist_ages[k2]]
    N2[(dist_ages[k2+1]+3)] <- N2temp[(dist_ages[k2]+3)]
    }
# carry over the animals that are in age class 10    
  N2[dist_ages[10]] <- N2[dist_ages[10]] + N2temp[dist_ages[10]]
  N2[(dist_ages[10]+3)] <- N2[(dist_ages[10]+3)] + N2temp[(dist_ages[10]+3)]   

  
  list(N1 = N1, N2 = N2)
}

