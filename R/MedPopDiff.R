rm(list = ls())
library(ggplot2)
library(matrixStats)
library(grid)

  
load('HS_OutputFile_Concurrent.rdata')

# check number of bootstraps and number of years in simulation
dimensions <- dim(dat.out)
nboot <- dimensions[3]
number_of_years <- dimensions[1] - 1

# input year in which piling begins, if it's not the first year
start_year <- 1

# calculate first year and number of years to display

years <- number_of_years - start_year + 1

# calculate values for "years post-piling axis"
time<-c(0:(years-1))

difference<-matrix(0, nrow=years, ncol=nboot)
limits <- matrix(0, nrow=years,ncol=2)

# determine median difference between Ndist and Nnotdist in each year
difference<- dat.out[start_year:number_of_years,3,1:nboot]
diff_med<-rowMedians(difference)
# diff_mean <- rowMeans(difference)
# df=data.frame(x=time,y=diff_med)

# determine 80% "confidence" region
for (t in 1:years){limits[t,] <- quantile(difference[t,], c(0.10,0.90))} 
df = data.frame(x = time, y = diff_med, lcl = limits[,1], ucl = limits[,2])
# calculate range for y-axis values
lower <- min(limits[,1]) - 100
upper <- max(limits[,2]) + 100

# plot annual median difference between Nnotdist and Ndist
z <- ggplot(aes(x, y), data=df) + 
  geom_line(method="loess") + 
  geom_line(aes(x, lcl), lty = 'dotted') +
  geom_line(aes(x, ucl), lty = 'dotted') +
  scale_y_continuous(limits = c(lower,upper)) +
  labs(x = "Year", y = "Median Population Difference") + 
  theme(axis.text=element_text(size=12),axis.title=element_text(size=18,face="bold")) + 
  theme(panel.margin = unit(1, "lines")) # just to space 'em out and avoid the 25 obscuring the 0
# png(file = 'MedPopDiffPlot_HS_Concurrent.png', width = 1200, height = 800)

print(z)

 
# calculate median annual % change in Ndist and 80% confidence region
Ndist_dec <- matrix(0, nrow=years, ncol=nboot)
Ndist_dec <-dat.out[start_year:number_of_years,4,1:nboot]
Ndist_dec_med<-rowMedians(Ndist_dec)
# Ndist_dec_mean <- rowMeans(Ndist_dec)
# df=data.frame(x=time,y=Ndist_dec_med)

for (t in 1:years){limits[t,] <- quantile(Ndist_dec[t,], c(0.10,0.90))}
lcl <- limits[,1]
ucl <- limits[,2]

df = data.frame(x = time, y = Ndist_dec_med, lcl = limits[,1], ucl = limits[,2])
# calculate range for y-axis values
lower <- min(limits[,1]) - 1
upper <- max(limits[,2]) + 1

# plot % annual change in Ndist
z <- ggplot(aes(x, y), data=df) +
  geom_line(method="loess") + 
  geom_line(aes(x, lcl), lty = 'dotted') +
  geom_line(aes(x, ucl), lty = 'dotted') +
  scale_y_continuous(limits = c(lower,upper)) +
  labs(x = "Year", y = "Annual % Change in Population Size") + 
  theme(axis.text=element_text(size=12),axis.title=element_text(size=18,face="bold")) + 
  theme(panel.margin = unit(1, "lines")) # just to space 'em out and avoid the 25 obscuring the 0
#png(file = 'NdistDecline_Plot_HS_Concurrent.png', width = 1200, height = 800)

print(z)


# calculate median annual % change in Nnotdist and CLs
Nnotdist_dec <- matrix(0, nrow=years, ncol=nboot)
Nnotdist_dec <-dat.out[start_year:24,5,1:nboot]
Nnotdist_dec_med<-rowMedians(Nnotdist_dec)

# Nnotdist_dec_mean <- rowMeans(Nnotdist_dec)
# df=data.frame(x=time,y=Nnotdist_dec_med)

for (t in 1:years){limits[t,] <- quantile(Nnotdist_dec[t,], c(0.10,0.90))}
lcl <- limits[,1]
ucl <- limits[,2]

df = data.frame(x = time, y = Nnotdist_dec_med, lcl = limits[,1], ucl = limits[,2])
# calculate range for y-axis values
lower <- min(limits[,1]) - 1
upper <- max(limits[,2]) + 1

# plot % annual change in Nnotdist
z <- ggplot(aes(x, y), data=df) + 
  geom_line(method="loess") + 
  geom_line(aes(x, lcl), lty = 'dotted') +
  geom_line(aes(x, ucl), lty = 'dotted') +
  scale_y_continuous(limits = c(lower,upper)) +
  labs(x = "Year", y = "Annual % Change in Population Size") + 
  theme(axis.text=element_text(size=12),axis.title=element_text(size=18,face="bold")) + 
  theme(panel.margin = unit(1, "lines")) # just to space 'em out and avoid the 25 obscuring the 0
#png(file = 'NdistDecline_Plot_HS_Concurrent.png', width = 1200, height = 800)

print(z)

