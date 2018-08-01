library(splines)
library(R2jags)
###############################################################################
###############################################################################
###  Ivo M. Foppa, Aug 2018---National RSV excess mort. analysis          #####
###############################################################################
###############################################################################
bfolder <- "C:/Users/vor1/Dropbox/Misc work/RSV/RSV git project/"

setwd(paste0(bfolder,'RSVData'))
###############################################################################
fromyr <- 2010
toyr <- 2016
fromwk <- 40
toseaswk <- 17
towk <- 52
nseas <- toyr - fromyr
### Define 'global N'

mmwrdat<-read.table("mmwrweeks.txt",header=T)
mmwrdat<-mmwrdat[which((mmwrdat$wyear==fromyr&mmwrdat$week>=fromwk)|(mmwrdat$wyear>fromyr&mmwrdat$wyear<toyr)|(mmwrdat$wyear==toyr&mmwrdat$week<=towk)),]

DVDun <- mmwrdat$dvdweek

N <- length(mmwrdat$dvdweek)
seasbeg <- which(mmwrdat$week==fromwk)
seasend <- which(mmwrdat$week==toseaswk)
time <- (1:N)/N
#######################################################################################
for (k in 1:6){
  assign(paste("seas",k,sep=""),c(max(seasbeg[k],3),min(seasend[k],N)))
}
#######################################################################
niter <- 5000
nadapt <- 1000
###################################################################################################
###################################################################################################
###################################################################################################
filename <- paste0('data_National_2010_16_5ag.dat')
setwd(paste0(bfolder,'RSVData'))
datarr <- data.frame(read.table(file = filename,header = T))
###################################################################################################
###################################################################################################
###  NS, flu mort indicator #######################################################################
###################################################################################################
qual1 <- 'flumort'
model1.file <- paste0('model NS ',qual1,'.txt')
ag <- 5
###################################################################################################
for (ag in rev(1:5)){
  agdata <- datarr[which(datarr$age==ag),]
  RSV <- agdata$RSV
  flumort <- agdata$flu
  mort <- agdata$rcu - flumort
  N <- length(mort)
  pop <- agdata$pop
  
  N <- length(mort)
  time <- (1:N)/N
  ###################################################################################################
  variables7 <- c('EMRSV1tot','EMRSV2tot','EMRSV3tot','EMRSV4tot','EMRSV5tot','EMRSV6tot','EMRSVtot',
                  'EMflu1tot','EMflu2tot','EMflu3tot','EMflu4tot','EMflu5tot','EMflu6tot','EMflutot')
  
  nknots <- 4
  #######################################################################################
  ###################################################################################################
  ndf <- round((nknots + 1) * nseas/(nseas*52.5)*N)
  nsarr <- ns(time,df = ndf)[,]  ## defining basis matrix
  #######################################################################################
  mod <- lm(mort ~ ns(time, df = ndf))
  smod <- summary(mod)
  coeffls <- as.numeric(smod$coefficients[,1])/mean(pop)
  nsinit <- coeffls[2:(ndf + 1)]
  b0init <- coeffls[1]
  #######################################################################################
  data <- list('N'=N,'ndf'=ndf,'ns'=nsarr,'y'=mort,'RSV'=RSV,'pop'=pop, 'flumort'=flumort, 
               'seas1'=seas1,'seas2'=seas2,'seas3'=seas3,'seas4'=seas4,'seas5'=seas5,'seas6'=seas6)
  #######################################################################################
  rinit <- (mort/(b0init*pop))[3:N]
  sigmainit <- sapply(1:N,function(x) ifelse(flumort[x] > 0, flumort[x]/pop[x],.1/pop[x]))
  inits <- function() {
    list(
      b0=b0init,
      
      b10=0,
      b11=0,
      b12=0,
      
      b20=0,
      b21=0,
      b22=0,
      
      sigma=sigmainit,
      
      b=nsinit,
      
      logalpha=0.0,  
      r=rinit
    )}
  
  setwd(paste0(bfolder,'RSVmodels'))
  
  j.model <- jags.model(file=model1.file,data=data, inits=inits, n.adapt=nadapt, n.chains=3)
  j.samples<-coda.samples(j.model, variable.names=variables7, n.iter=niter, thin = 5) 
  
  codaarr <- rbind(j.samples[[1]],j.samples[[2]],j.samples[[3]],deparse.level=0)
  assign(paste0("codaarr",ag),codaarr)
  
  setwd(paste0(bfolder,'RSVMCMCoutput'))
  fname <- paste0('codaarr',ag,' RSV ',qual1,' ',nknots,' knots ps ',nadapt,' ',round(niter/5*3),'.RData')
  setwd(paste0(bfolder,'RSVMCMCoutput'))
  save(codaarr,file = fname)
  
  cat(paste0('\nAge group ',ag,': done\n'))
}
  
###################################################################################################
###################################################################################################
###################################################################################################
###  NS, vir indicator  ###########################################################################
###################################################################################################
qual2 <- 'vir'
model2.file <- paste0('model NS ',qual2,'.txt')
#######################################################################################
for (ag in rev(1:5)){
  AH1P <- agdata$AH1P
  AH3 <- agdata$AH3
  B <- agdata$B
  ###################################################################################################
  data <- list('N'=N,'ndf'=ndf,'ns'=nsarr,'y'=mort,'RSV'=RSV,'pop'=pop, 
               'AH1P'=AH1P,'AH3'=AH3,'B'=B, 
               'seas1'=seas1,'seas2'=seas2,'seas3'=seas3,'seas4'=seas4,'seas5'=seas5,'seas6'=seas6)
  #######################################################################################
  inits <- function() {
    list(
      b0=b0init,
      
      b30=0,
      b31=0,
      b32=0,
      
      b40=0,
      b41=0,
      b42=0,
      
      b50=0,
      b51=0,
      b52=0,
      
      b70=0,
      b71=0,
      b72=0,
      
      b=nsinit,
      
      logalpha=0.0,  
      r=rinit
    )}
  
  setwd(paste0(bfolder,'RSVmodels'))
  
  j.model <- jags.model(file=model2.file,data=data, inits=inits, n.adapt=nadapt, n.chains=3)
  j.samples<-coda.samples(j.model, variable.names=variables7, n.iter=niter, thin = 5) 
  
  codals <- rbind(j.samples[[1]],j.samples[[2]],j.samples[[3]],deparse.level=0)
  codaarr <- rbind(j.samples[[1]],j.samples[[2]],j.samples[[3]],deparse.level=0)
  assign(paste0("codaarr",ag),codaarr)
  
  setwd(paste0(bfolder,'RSVMCMCoutput'))
  fname <- paste0('codaarr',ag,' RSV ',qual2,' ',nknots,' knots ps ',nadapt,' ',round(niter/5*3),'.RData')
  setwd(paste0(bfolder,'RSVMCMCoutput'))
  save(codaarr,file = fname)
  
  cat(paste0('\nAge group ',ag,': done\n'))
}


###################################################################################################
###################################################################################################
setwd(paste0(bfolder1,'RSVresults'))

for (k in 1:7){
  assign(paste('el',k,sep = ''), round(as.vector(quantile(codaarr[,k],prob=c(.5,.05,0.975)))))
}
newline <- paste0('Influenza viral indicator\t',el1[1],' (',el1[2],',',el1[3],')\t',
                  el2[1],' (',el2[2],',',el2[3],')\t',
                  el3[1],' (',el3[2],',',el3[3],')\t',
                  el4[1],' (',el4[2],',',el4[3],')\t',
                  el5[1],' (',el5[2],',',el5[3],')\t',
                  el6[1],' (',el6[2],',',el6[3],')\t',
                  el7[1],' (',el7[2],',',el7[3],')')

write.table(newline,outfile,append=T,row.names=FALSE,col.names=FALSE, quote=FALSE)
###################################################################################################
plot(RSV,type = 'l',lwd = 2,xlab = 'Index Week')
lines(AH3/max(AH3)*max(RSV)*0.8,col = 'red',lwd = 2)
lines(flumort/max(flumort)*max(RSV)*0.8,col = 'blue',lwd = 2)
par(mar=c(0,0,0,0))

legend('topleft',c('RSV','Influenza A(H3N2'), col = c('black','red'), lwd = 2,y.intersp = .8,bty = 'n')
###################################################################################################plot(RSV,type = 'l',lwd = 2,xlab = 'Index Week')
plot(mort,type = 'l',lwd = 2,xlab = 'Index Week',ylim = c(0,max(mort)))
lines(RSV/max(RSV)*max(mort)*0.8,col = 'red',lwd = 2)
lines(flumort/max(flumort)*max(mort)*0.8,col = 'blue',lwd = 2)
par(mar=c(0,0,0,0))

legend('topleft',c('R&C mortality','RSV','Influenza mort.'), col = c('black','red','blue'), lwd = 2,y.intersp = .8,bty = 'n')
###################################################################################################