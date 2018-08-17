library(splines)
###############################################################################
###############################################################################
###  Ivo M. Foppa, Aug 2018---National RSV excess mort. analysis          #####
###############################################################################
###############################################################################
setwd(paste0(bfolder,'RSVData'))
###############################################################################
###################################################################################################
###################################################################################################
setwd(paste0(bfolder,'RSVData'))
datarr <- data.frame(read.table(file = filename,header = T))
###################################################################################################
###################################################################################################
###  NS, flu mort indicator #######################################################################
###################################################################################################
ag <- 5
agls <- c('<5','5-17','18-49','50-64','65+','All')
###################################################################################################
for (ag in rev(1:5)){
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
  data <- data.frame(mort,time=time,flumort0,flumort1,flumort2,RSV0,RSV1,RSV2,pop,nsarr1)

  mod <- glm(mort ~ ns(time,df=ndf) + flumort0 + flumort1 + flumort2 + RSV0 + RSV1 + RSV2 + offset(log(pop)),
             data = data,family = poisson(link = 'log'))
 
  dataRSV0 <- data.frame(mort,time,flumort0,flumort1,flumort2,
                         RSV0=RSV0*0,RSV1=RSV1*0,RSV2=RSV2*0,pop,nsarr1)
  dataflu0 <- data.frame(mort,time,flumort0=flumort0*0,flumort1=flumort1*0,flumort2=flumort2*0,
                         RSV0,RSV1,RSV2,pop,nsarr1)
  
  mort1 <- predict(mod, type = 'response')
  mortRSV0 <- predict(mod, newdata = dataRSV0, type = 'response')
  mortflu0 <- predict(mod, newdata = dataflu0, type = 'response')
  cat(paste0(agls[ag],':\tRSV: ',round(sum(mort-mortRSV0)),'\t Flu: ', round(sum(mort + flumort[-c(1:2)]-mortflu0)),'\n'))
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
  agdata <- datarr[which(datarr$age==ag),]

  AH1P <- agdata$AH1P
  AH3 <- agdata$AH3
  B <- agdata$B
  ###################################################################################################
  #######################################################################################
  RSV <- agdata$RSV
  RSV0 <- RSV[-c(1:2)]
  RSV1 <- RSV[-c(1,length(RSV))]
  RSV2 <- RSV[-c((length(RSV)-1):length(RSV))]
  
  AH1P <- agdata$AH1P
  AH1P0 <- AH1P[-c(1:2)]
  AH1P1 <- AH1P[-c(1,length(AH1P))]
  AH1P2 <- AH1P[-c((length(AH1P)-1):length(AH1P))]
  
  AH3 <- agdata$AH3
  AH30 <- AH3[-c(1:2)]
  AH31 <- AH3[-c(1,length(AH3))]
  AH32 <- AH3[-c((length(AH3)-1):length(AH3))]
  
  B <- agdata$B
  B0 <- B[-c(1:2)]
  B1 <- B[-c(1,length(B))]
  B2 <- B[-c((length(B)-1):length(B))]
  
  mort <- agdata$rcu[-c(1:2)]
  pop <- agdata$pop[-c(1:2)]
  
  N <- length(mort)
  time <- ((1:N)/N)
  ###################################################################################################
  nknots <- 4
  #######################################################################################
  ###################################################################################################
  ndf <- round((nknots + 1) * nseas/(nseas*52.5)*N)
  #######################################################################################
  data <- data.frame(mort,time,
                     AH1P0,AH1P1,AH1P2,AH30,AH31,AH32,B0,B1,B2,
                     RSV0,RSV1,RSV2,pop)
  
mod <- glm(mort ~ ns(time,df = ndf) + AH1P0 + AH1P1 + AH1P2 + AH30 + 
             AH31 + AH32 + B0 + B1 + B2 + RSV0 + RSV1 + RSV2 + offset(log(pop)),
             data = data,family = poisson(link = 'log'))
  
data0 <- data.frame(mort,time,
                    AH1P0=AH1P0*0,AH1P1=AH1P1*0,AH1P2=AH1P2*0,
                    AH30=AH30*0,AH31=AH31*0,AH32=AH32*0,
                    B0=B0*0,B1=B1*0,B2=B2*0,
                    RSV0,RSV1,RSV2,
                    pop)

mort1 <- predict(mod, type = 'response')
mort0 <- predict(mod, newdata = data0, type = 'response')
sum(mort-mort0)
}
###################################################################################################
###################################################################################################