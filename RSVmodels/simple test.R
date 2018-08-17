library(splines)
library(R2jags)
###############################################################################
###############################################################################
###  Ivo M. Foppa, Aug 2018---Trouble shooting Bayesian-frequentist comparison 
###############################################################################
###############################################################################
setwd(paste0(bfolder,'RSVData'))
datarr <- data.frame(read.table(file = filename,header = T))
###################################################################################################
###################################################################################################
###  NS, flu mort indicator #######################################################################
###################################################################################################
###################################################################################################
ag <- 5

agdata <- datarr[which(datarr$age==ag),]
RSV <- agdata$RSV
flumort <- agdata$flu
mort <- agdata$rcu - flumort
N <- length(mort)
pop <- agdata$pop

N <- length(mort)
time <- (1:N)/N
###################################################################################################
variables7 <- c('b0','b10','b11','b12','b20','b21','b22')
#######################################################################################
###################################################################################################
ndf <- round((nknots + 1) * nseas/(nseas*52.5)*N)
nsarr <- ns(time,df = ndf)[,]  ## defining basis matrix
#######################################################################################
mod0 <- lm(log(mort/pop) ~ ns(time, df = ndf))
smod <- summary(mod0)
coeffls <- as.numeric(smod$coefficients[,1])
nsinit <- coeffls[2:(ndf + 1)]
b0init <- coeffls[1]
#######################################################################################
data <- list('N'=N,'y'=mort,'RSV'=RSV,'pop'=pop, 'flumort'=flumort, 
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
    b22=0
  )}

setwd(paste0(bfolder,'RSVmodels'))

j.model <- jags.model(file='test model.txt',data=data, inits=inits, n.adapt=nadapt, n.chains=3)
j.samples<-coda.samples(j.model, variable.names=variables7, n.iter=niter, thin = 5) 

#### Frequentist equivalent
agdata <- datarr[which(datarr$age==ag),]
RSV <- agdata$RSV
RSV0 <- RSV[-c(1:2)]
RSV1 <- RSV[-c(1,length(RSV))]
RSV2 <- RSV[-c((length(RSV)-1):length(RSV))]

flumort <- agdata$flu
flumort0 <- flumort[-c(1:2)]
flumort1 <- flumort[-c(1,length(flumort))]
flumort2 <- flumort[-c((length(flumort)-1):length(flumort))]

mort <- (agdata$rcu - flumort)[-c(1:2)]
pop <- agdata$pop[-c(1:2)]

N <- length(mort)
time <- (1:N)/N
#######################################################################################
###################################################################################################
ndf <- round((nknots + 1) * nseas/(nseas*52.5)*N)
### using basis matrix from Bayesian analysis
nsarr1 <- nsarr[1:length(mort),]
#######################################################################################
data <- data.frame(mort,time=time,flumort0,flumort1,flumort2,RSV0,RSV1,RSV2,pop)

mod <- glm(mort ~ flumort0 + flumort1 + flumort2 + RSV0 + RSV1 + RSV2 + offset(log(pop)),
           data = data,family = poisson(link = 'log'))
