rm(list = ls())

par(mfrow=c(1,1))
load('HS_OutputFile_Concurrent.rdata')

# check number of bootstraps and number of years in simulation
dimensions <- dim(dat.out)
numsim <- dimensions[3]
number_of_years <- dimensions[1] - 1

# input year in which piling begins, if it's not the first year
start_year <- 1

# determine years for which a plot is required
ref_years <- seq(start_year,number_of_years,6)

# we actually require the first plot to be for the first year after piling starts
ref_years[1] <- ref_years[1] + 1
index <- 1

# create array to store risk values for dist and notdist populations
risk_dist <- risk_notdist <- Add_Risk <- array (NA, dim = c(length(ref_years),5))
colnames(risk_dist) <- colnames(risk_notdist) <- colnames(Add_Risk) <- c('year','   risk 1% decline', '  risk 2% decline', '  risk 5% decline', '  median % change')

# calculate risks and plot risk curves for each reference year

for (i in ref_years){

#  par(mfrow=c(1,1))
# create a cumulative response curve for probability of % change in disturbed population
  reversed <- round((-1*dat.out[i,4,1:numsim]),4)
  diff<-sort(unique(reversed))                                  
  perc<-table(reversed)
  z<-perc/numsim
  pb<-as.vector(z)
  pb2<-cumsum(rev(pb))
  
years_since_start <- i - start_year 

# these lines have been revised with simpler code
ONE= length(which(dat.out[i,4,] <= -1.0))/numsim
TWO= length(which(dat.out[i,4,] <= -2.0))/numsim
FIVE=length(which(dat.out[i,4,] <= -5.0))/numsim
MEDIANs<-round(median(dat.out[i,4,]), digits=4)

risk_dist[index,1] <- i-1
risk_dist[index,2] <- ONE
risk_dist[index,3] <- TWO
risk_dist[index,4] <- FIVE
risk_dist[index,5] <- MEDIANs
MEAN <- -1*mean(dat.out[i,4,])

title<-(paste("Risks",years_since_start,"years after start", sep=" "))
plot(rev(diff),pb2,type="l", ylab="Cumulative Probability",xlab="Percent Decline",yaxt="n",cex.lab=1.5,cex.axis=1.5,lwd=2,xlim=c(0,max(rev(diff))+0.5),main=title)
axis(2, at=c(0,0.25,0.5,0.75,1),cex.axis=1.5)
maxlim<-max(rev(diff))  
segments(x0=-5, y0=pb2[rev(diff)==(rev(diff)[which(abs(rev(diff)-1)== min(abs(rev(diff)-1)))])], x1=1, y1=pb2[rev(diff)==(rev(diff)[which(abs(rev(diff)-1)== min(abs(rev(diff)-1)))])],col="darkblue",lty=2,lwd=2)
segments(x0=1, y0=-5, x1=1, y1=pb2[rev(diff)==(rev(diff)[which(abs(rev(diff)-1)== min(abs(rev(diff)-1)))])],col="darkblue",lty=2,lwd=2)
 
segments(x0=-5, y0=pb2[rev(diff)==(rev(diff)[which(abs(rev(diff)-2)== min(abs(rev(diff)-2)))])], x1=2, y1=pb2[rev(diff)==(rev(diff)[which(abs(rev(diff)-2)== min(abs(rev(diff)-2)))])],col="red",lty=2,lwd=2)
segments(x0=2, y0=-5, x1=2, y1=pb2[rev(diff)==(rev(diff)[which(abs(rev(diff)-2)== min(abs(rev(diff)-2)))])],col="red",lty=2,lwd=2)
   
segments(x0=-5, y0=pb2[rev(diff)==(rev(diff)[which(abs(rev(diff)-5)== min(abs(rev(diff)-5)))])], x1=5, y1=pb2[rev(diff)==(rev(diff)[which(abs(rev(diff)-5)== min(abs(rev(diff)-5)))])],col="darkgreen",lty=2,lwd=2)
segments(x0=5, y0=-5, x1=5, y1=pb2[rev(diff)==(rev(diff)[which(abs(rev(diff)-5)== min(abs(rev(diff)-5)))])],col="darkgreen",lty=2,lwd=2)
   
segments(y0=0.5, x0=-5, y1=0.5, x1= MEAN,col="red",lty=2,lwd=5)
segments(x0=MEAN, y0=-0.5, x1=MEAN, y1=0.5,col="red",lty=2,lwd=5)
   

legend(maxlim-(0.6*maxlim),0.98,legend=c("Prob. 1% decline"),bty="n",cex=1,lty=2,lwd=2,col="darkblue")
legend(maxlim-(0.6*maxlim),0.88,legend=c("Prob. 2% decline"),bty="n",col="red",cex=1,lty=2,lwd=2)
legend(maxlim-(0.6*maxlim),0.78,legend=c("Prob. 5% decline"),bty="n",col="darkgreen",cex=1,lty=2,lwd=2)
legend(maxlim-(0.6*maxlim),0.68,legend=c("Median % decline"), bty="n",col="red",cex=1,lty=2,lwd=5)

#save plot as file
#png(filename=(paste("HP_S1_RiskCurve_", i, "_years.png")), type="cairo", units="in", width=8, height=6.4, pointsize=12, res=96)
#dev.off()


#legend(maxlim-(0.2*maxlim),0.98,legend=ONE,bty="n")
#legend(maxlim-(0.2*maxlim),0.88,legend=TWO,bty="n")
#legend(maxlim-(0.2*maxlim),0.78,legend=FIVE,bty="n")
#legend(maxlim-(0.15*maxlim),0.68,legend=MEDIANs,bty="n")
#legend(maxlim-(0.07*maxlim),0.68,legend="%",bty="n")

#calculate risks for undisturbed populations

ONE_notdist = length(which(dat.out[i,5,] <= -1.0))/numsim
TWO_notdist = length(which(dat.out[i,5,] <= -2.0))/numsim
FIVE_notdist =length(which(dat.out[i,5,] <= -5.0))/numsim
MEDIANs_notdist <- round(median(dat.out[i,5,]), digits=4)

risk_notdist[index,1] <- i-1
risk_notdist[index,2] <- ONE_notdist
risk_notdist[index,3] <- TWO_notdist
risk_notdist[index,4] <- FIVE_notdist
risk_notdist[index,5] <- MEDIANs_notdist

# Calculate additional risks due to disturbance
Add_Risk[index,1] <- i-1

index <- index + 1
}

Add_Risk[,2:4] <- risk_dist[,2:4] - risk_notdist[,2:4]

# Save risk results as .csv file - REMEMBER TO CHANGE FILE NAME IN LINE 119!!!

### make a table of probs ###
MainHeadings<-cbind('Predicted Cumulative Probabilities of Population Decline',' ',' ',' ',' ')
HeadingDist1<-cbind('Disturbed Population:','','','','')
HeadingDist2<-cbind('Year','Prob. 1% decline','Prob. 2% decline','Prob. 5% decline','Median % decline')
Spacer<-cbind(' ',' ',' ',' ',' ')
HeadingNonDist1<-cbind('Undisturbed Population:','','','','')
HeadingNonDist2<-cbind('Year','Prob. 1% decline','Prob. 2% decline','Prob. 5% decline','Median % decline')
HeadingAddRisk1<-cbind('Additional Risk from Operations:','','','','')
HeadingAddRisk2<-cbind('Year','Prob. 1% decline','Prob. 2% decline','Prob. 5% decline','Median % decline')


MainOutput<-rbind(MainHeadings,Spacer,HeadingDist1,HeadingDist2,risk_dist,Spacer,HeadingNonDist1,HeadingNonDist2,risk_notdist,Spacer,HeadingAddRisk1,HeadingAddRisk2,Add_Risk)

write.table(MainOutput,col.names=F,row.names=F, sep=",", file = 'HS_Sequential_Risk.csv')


