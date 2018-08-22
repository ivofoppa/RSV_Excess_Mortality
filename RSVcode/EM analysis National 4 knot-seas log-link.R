library(splines)
library(R2jags)
###############################################################################
###############################################################################
###  Ivo M. Foppa, Aug 2018---National RSV excess mort. analysis          #####
###############################################################################
###############################################################################
setwd(paste0(bfolder,'RSVdata'))
datarr <- data.frame(read.table(file = filename,header = T))
###################################################################################################
###################################################################################################
###  NS, flu mort indicator #######################################################################
###################################################################################################
###################################################################################################
ag <- 5
for (ag in rev(4:5)){
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
  variables7 <- c('b0','b[1:31]','b10','b11','b12','b20','b21','b22')
  #######################################################################################
  ###################################################################################################
  ndf <- round((nknots + 1) * nseas/(nseas*52.5)*N)
  nsarr <- ns(time,df = ndf)[,]  ## defining basis matrix
  #######################################################################################
  mod <- lm(log(mort/pop) ~ ns(time, df = ndf))
  smod <- summary(mod)
  coeffls <- as.numeric(smod$coefficients[,1])
  nsinit <- coeffls[2:(ndf + 1)]
  b0init <- coeffls[1]
  #######################################################################################
  data <- list('N'=N,'ndf'=ndf,'ns'=nsarr,'y'=mort,'RSV'=RSV,'pop'=pop, 'flumort'=flumort, 
               'seas1'=seas1,'seas2'=seas2,'seas3'=seas3,'seas4'=seas4,'seas5'=seas5,'seas6'=seas6)
  #######################################################################################
  # rinit <- (mort/(exp(b0init)*pop))[3:N]
  
  inits <- function() {
    list(
      b0=b0init,
      
      b10=0,
      b11=0,
      b12=0,
      
      b20=0,
      b21=0,
      b22=0,
      
      b=nsinit
    )}
  
  setwd(paste0(bfolder,'RSVmodels'))
  
  j.model <- jags.model(file=model3.file,data=data, inits=inits, n.adapt=nadapt, n.chains=3)
  j.samples<-coda.samples(j.model, variable.names=variables7, n.iter=niter, thin = 5) 
  
  codaarr <- rbind(j.samples[[1]],j.samples[[2]],j.samples[[3]],deparse.level=0)
  assign(paste0("codaarr",ag),codaarr)
  
  setwd(paste0(bfolder,'RSVMCMCoutput'))
  fname <- paste0('codaarr',ag,' RSV ',qual3,' ',nknots,' knots ps ',nadapt,' ',round(niter/5*3),'.RData')
  setwd(paste0(bfolder,'RSVMCMCoutput'))
  save(codaarr,file = fname)
  
  cat(paste0('\nAge group ',ag,': done\n'))
}
  
###################################################################################################