rm(list = ls())

# iPCoD PROTOCOL STEP 10
# Set the name of the MakeInit file you want to use

source('PCoDFunctionsFinal.R')
source('config.R')
source('MakeInitScenarioHS.R')


# matrix to store values that will be retained from all simulations, now extended to 5 values from original 4
  dat.out <- array(NA, dim = c(years, 5, nboot)) 

# NdistSum is the total population size for the disturbed population each year, NNotDistSum is the equivalent for the matching undisturbed population
  colnames(dat.out) <- c('NdistSum', 'NnotdistSum','Nnotdist-Ndist', '% decline Ndist', '% decline Nnotdist')
        
# store row/element subscripts for undisturbed age classes 
  dist_ages <- seq(1, 55, 6)

# matrices to store numbers affected in different ways in a more accessible format than Ndist for last iteration of bootstrap.   
  Ndist_out_1 <- Ndist_out_2 <- array(0, dim = c(10, 7, pile_years))
  colnames(Ndist_out_1) <- colnames(Ndist_out_2) <- c('age', 'not disturbed','mod dist', 'severe dist', 'PTS only','PTS+mod dist','PTS+severe dist')

  for(boot_number in 1:nboot){  

  Ndist <- matrix(0, nrow = 60, ncol = years)  
  Nnotdist <- matrix(0, nrow = 10, ncol = years)
  NColl_byage <- rep(0,60)
  sdNDt <- rep(0,nvulnmean)
# create array to hold age structures of each vulnerable sub-population and undisturbed remainder  
  Ndistvuln <- array(0, dim = c(60, length(newvulnmean),(pile_years + 1))) 

# call f.experts.survival&fertility to sample bootstrap values from the expert elicitation
  tmp  <- f.experts.survival_and_fertility(spec) 
  Surv <- tmp$Surv
  Fert <- tmp$Fert

# Number of days of disturbance required for moderate and severe disturbance of dependent (pups/calves) (dd1), juveniles (dd2) and adults (dd3)
  dd1  <- tmp$dd1 
  dd2  <- tmp$dd2
  dd3  <- tmp$dd3

# create matrix with dd's for each age class, to be used in function f.disturbance.classes
  dds <- matrix(0, nrow = 10, ncol = 3)
  for (di in 1:3){
    dds[,di]<- rep(c(dd1[di],dd2[di],dd3[di]),c(age1,(age2-age1),(10-age2)))
  }

# set up identical age structures for disturbed and undisturbed components of population  
  Nnotdist[1:10, 1]   <- round(pmean * age_structure[1:10], 0)
  Ndist[dist_ages, 1] <- Nnotdist[1:10,1]
  
# determine age structures for each vulnerable sub-population, so that their combined age structure is the same as Ndist
  for (ivuln in dist_ages){
    Ndistvuln[ivuln,1:length(newvulnmean),1] <- rmultinom(1,Ndist[ivuln,1],newvulnmean)
    }
 
# vector to hold fertility values for current year
  AnnFert<-rep(c(0,6))
  
 # incorporate uncertainty into NDt, NPt and NCollisions, 95% confidence limits are mean +/- 50%
  sd <- exp(rnorm(1, 0, 0.01))
  sdNDt[1:nvulnmean]  <- sd
  NColl_var <- round(NCollisions*sd)

# skip everything from here to line 186 if there is no piling
  if (pile_years > 0){
  
  NDtvar <- round(NDt * sdNDt)
  NPtvar <- round(NPt*sdNDt)
  
  
  # check for zeros and convert to something small
  NDtvar <- ifelse (NDtvar > 0, NDtvar, 0.000001)
  NPtvar <- ifelse (NPtvar > 0, NPtvar, 0.000001)
  
  
  
  p_disturb  <- matrix(0, nrow = pile_years * 365, ncol = nvulnmean)
  colnames(p_disturb) <- paste("sub-pop", 1:nvulnmean, sep = '')
  p_pts <- p_disturb
  
  
  for (j in 1:(nvulnmean)) { 
    for(k in 1:nrow(p_disturb))  {
      p_disturb[k, j] <- ifelse( sum((pile[k, pilesx] * vulnpile[j, pilesx] * NDtvar[k, j]) / (vulnmean[j] * pmean)) > 1.0, 1.0, 
                                 sum((pile[k, pilesx] * vulnpile[j, pilesx] * NDtvar[k, j]) / (vulnmean[j] * pmean))) 
      
      p_pts[k, j]     <- sum(pile[k, pilesx] * vulnpile[j, pilesx] * NPtvar[k, j]) / sum(pile[k, pilesx] * vulnpile[j, pilesx] * NDtvar[k, j])
    }
  }
  
  p_pts[is.na(p_pts)] <- 0  
  
   for(t in 1:pile_years){ 

#  recalculate vuln so that it matches with rescaled sub-population age structures 
#  Link to PopDyn text A 
    if(dim(Ndistvuln)[2] == 1){vuln <-sum(Ndistvuln[,,t])}else{vuln <-colSums(Ndistvuln[,,t])}
    wvec <- which(vuln[1:nvulnmean] > 1000)
    vuln[wvec] <- 1000   
    
# call f.current.survival.values to determine survival values for each disturbance category in current year
    AnnSurv   <- f.current.survival.values(dist_ages, age1, age2, a, b)
  
 
# determine fertility rates for current year                            
    AnnFert[1] <- (rbeta(1, a[3], b[3]))/2
    AnnFert[2:6] <- Fert[2:6]*AnnFert[1] 



for(i in 1:nvulnmean){ # start the loop over the multiple vulnerable sub-pops; 
      if(Day1){
# call f.disturbance&PTS.Model2 (animals only get PTS the first time they are disturbed) if Day1 is TRUE to determine exposure history for each vulnerable sub-population
        yrs <- which(yearvec == t)
        tmp <- f.disturbance_and_PTS.Model2(pile[yrs, max(pilesx) + 2  + i], vuln[i], p_disturb[yrs, i], p_pts[yrs, i]) 
      } else {
# otherwise call f.disturbance_and_PTS
        yrs <- which(yearvec == t)
        tmp       <- f.disturbance_and_PTS(pile[yrs, max(pilesx) + 2  + i], vuln[i], days, p_disturb[yrs, i], p_pts[yrs, i])
      }
      dist_year <- tmp$dist_year
      pts_year  <- tmp$pts_year  
      vulnpop <- i

# assign individuals in sub-population to different disturbance classes
      tmp       <- f.disturbance.classes(dist_year, pts_year, Ndistvuln[,i,t], dist_ages, days, dds)
      Ndistvuln[,i,t]    <- tmp$N  
      
    } # end i loop over vulnerable sub-populations  

# combine all disturbed sub-populations and store values before any mortality is applied
  if(dim(Ndistvuln)[2]== 1){Ntemp <- Ndistvuln[1:60,,t]} else {Ntemp <- rowSums(Ndistvuln[1:60,,t])}
  
  Ndist_out_1[,2:7,t] <- matrix(Ntemp, nrow = 10, ncol = 6, byrow = TRUE)
    
# apply survival rates to all age-classes in Nnotdist
  Nnotdisttemp <- matrix(0, nrow = 10, ncol = 25)
    if (pmean > threshold) {Nnotdisttemp<- round(Nnotdist[1:10,t]*AnnSurv[dist_ages],0)} else {
          for (irand in 1:10) {
               Nnotdisttemp[irand] <- rbinom(1,Nnotdist[irand,t],AnnSurv[dist_ages[irand]])}
               }
# determine pup/calf production in year t+1 for undisturbed population
    if (pmean > threshold) {Nnotdist[1,t+1] <- round(sum(Nnotdisttemp[age2:10])*AnnFert[1], 0)} else {
      Nnotdist[1,t+1] <- rbinom(1, sum(Nnotdisttemp[age2:10]),AnnFert[1])
               }   

# advance ages of undisturbed animals by one year
   Nnotdist[2:10,t+1] <- Nnotdisttemp[1:9]
   Nnotdist[10,t+1] <- Nnotdisttemp[10] + Nnotdist[10,t+1]
    
# call f.next.year to determine numbers in each age-class of the vulnerable sub-populations in next year
    for (k in 1:length(newvulnmean)){
      tmp                <- f.next.year(Ndistvuln[,k , t], Ndistvuln[,k,t+1], dist_ages, pmean, AnnFert, threshold, age1, age2)
      Ndistvuln[,k,t+1]  <- tmp$N2
    }
    
# combine all vulnerable sub-populations into Ndist    
    if(dim(Ndistvuln)[2]== 1){Ntemp2 <- Ndistvuln[1:60,,t+1]} else {Ntemp2 <- rowSums(Ndistvuln[1:60,,t+1])}
    Ndist[,t+1] <- Ntemp2
    Ndist_out_2[,2:7, t] <- matrix(Ntemp2, nrow = 10, ncol = 6, byrow = TRUE)

# distribute deaths due to NCollisions randomly across age classes
  if (NCollisions != 0) {
  NColl_byage[dist_ages] <-   rmultinom(1, NColl_var , age_structure[1:10])
# account for mortality from other causes
  NColl_byage <- round(NColl_byage*sqrt(AnnSurv),0)
    #remove deaths by collision
  Ndist[,t+1] <- Ndist[,t+1] - NColl_byage
  } 
    
# save output. dat.out[,1,] is the total size of the disturbed population
# dat.out[,2,] is the size of the undisturbed population, dat.out[,3,] is the difference between the two populations
# dat.out[,4,] is the annual change in the size of the disturbed population as a percentage
# dat.out[,5,] is the annual change in the size of the undisturbed population as a percentage
  dat.out[t,1, boot_number] <- round(sum(Ndist[1:60,t])/propfemale, 0)
  dat.out[t,2, boot_number] <- round(sum(Nnotdist[1:10,t])/propfemale, 0)
  dat.out[t,3,boot_number] <- dat.out[t,1,boot_number] - dat.out[t,2,boot_number]
  dat.out[t,4,boot_number] <- 100*((dat.out[t,1,boot_number]/dat.out[1,1,boot_number])^(1/t)-1)
  dat.out[t,5,boot_number] <- 100*((dat.out[t,2,boot_number]/dat.out[1,2,boot_number])^(1/t)-1)
  Ndist_out_1[,1,t] <- Ndist_out_2[,1,t] <- c(0:9)
                        
  }# End t loop over Piling Years
 
}
  for(t in (pile_years + 1):(years-1)){

 # call f.current.survival.values to determine survival values for each disturbance category in current year
    AnnSurv   <- f.current.survival.values(dist_ages, age1, age2, a, b)
  
 # determine fertility rates for current year
    AnnFert[1] <- (rbeta(1, a[3], b[3]))/2
    AnnFert[2:6] <- Fert[2:6]*AnnFert[1] 

# apply survival rates to all age-classes in Nnotdist
     if (pmean > threshold) {Nnotdisttemp<- round(Nnotdist[1:10,t]*AnnSurv[dist_ages],0)} else {
          for (irand in 1:10) {
               Nnotdisttemp[irand] <- rbinom(1,Nnotdist[irand,t],AnnSurv[dist_ages[irand]])}
               }
    
# determine pup production in year t+1
      if (pmean > threshold) {Nnotdist[1,t+1] <- round(sum(Nnotdisttemp[age2:10])*AnnFert[1], 0)} else {
              Nnotdist[1,t+1] <- rbinom(1, sum(Nnotdisttemp[age2:10]),AnnFert[1])
               }
   
# advance ages by one year
    Nnotdist[2:10,t+1] <- Nnotdisttemp[1:9]
    Nnotdist[10,t+1] <- Nnotdisttemp[10] + Nnotdist[10,t+1]
    
# call next.year to determine numbers in each age-class in next year                                                                                                                 
    tmp         <- f.next.year(Ndist[, t], Ndist[, t + 1], dist_ages, pmean, AnnFert, threshold, age1, age2) 

    Ndist[,t+1] <- tmp$N2

# distribute deaths due to collisions across age classes
if (NCollisions != 0) {
  NColl_byage[dist_ages] <-   rmultinom(1, NColl_var , age_structure[1:10])
# account for mortality from other causes
    NColl_byage <- round(NColl_byage*sqrt(AnnSurv))
#remove deaths by collision
    Ndist[,t+1] <- Ndist[,t+1] - NColl_byage
  } 
 
# save output. dat.out[,1,] is the total size of the disturbed population
# dat.out[,2,] is the size of the undisturbed population, dat.out[,3,] is the difference between the two populations
# dat.out[,4,] is the annual change in the size of the disturbed population as a percentage
# dat.out[,5,] is the annual change in the size of the undisturbed population as a percentage
       
   dat.out[t,1, boot_number] <- round(sum(Ndist[1:60,t])/propfemale, 0)
   dat.out[t,2, boot_number] <- round(sum(Nnotdist[1:10,t])/propfemale, 0)
   dat.out[t,3,boot_number] <- dat.out[t,1,boot_number] - dat.out[t,2,boot_number]
  dat.out[t,4,boot_number] <- 100*((dat.out[t,1,boot_number]/dat.out[1,1,boot_number])^(1/t)-1)
  dat.out[t,5,boot_number] <- 100*((dat.out[t,2,boot_number]/dat.out[1,2,boot_number])^(1/t)-1)


  
  }# end loop over post-Piling years
  

boot_number  
} # end the bootstrap

#dat.out
#Ndist_out_1
#Ndist_out_2

# save number of animals in each disturbance class, before and after mortality is applied, for the years piling occurs for the last bootstrap
# save modifications to survival and fertility associated with PTS and disturbance for last bootstrap
# save summary statistics

# iPCoD PROTOCOL STEP 11

# Decide on a name for the file that will hold the output data
# Type this before _OutputFile in the next line

  save(Ndist_out_1, Ndist_out_2, dds, Surv, Fert, dat.out, file = paste(spec, '_OutputFile.rdata', sep = ''))