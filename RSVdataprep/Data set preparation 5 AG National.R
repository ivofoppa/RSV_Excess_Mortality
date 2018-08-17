#### Creating National data set 2010-16; for five age groups
### Ivo M. Foppa, 8/2018
#########################################################################################
#########################################################################################
### The directory containing all datasets has to be defined here:
setwd(paste0(bfolder,'RSVdata'))
#########################################################################################
###  VERY IMPORTANT: Global parameters have to be set here:   ###########################
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
library(stringr)
library(mgcv)

### Relevant time markers
fromyr <- 2010
toyr <- 2016
fromwk <- 40
towk <- 52

vPHdat <- data.frame(head(read.csv('PHdata.csv',header=T),-1)) ### virol. data from Public Health labs
DVD <- as.matrix(vPHdat$DVD)
selind <- which(DVD >= as.numeric(paste0(fromyr,fromwk)) & DVD <= as.numeric(paste(toyr,towk,sep="")))
vPHdat <- vPHdat[selind,]
#########################################################################################
DVD <- as.matrix(vPHdat$DVD)
DVDun <- as.numeric(unique(DVD))
#########################################################################################
datarr <- c()
missvals <- c()
for (dvd in DVDun){
  selind <- which(vPHdat$DVD==dvd)
  datblock <- vPHdat[selind,]
  datblock <- datblock[order(datblock$HHS),]
  missind <- which(!c(1:10) %in% datblock$HHS)
  if (length(missind)>0){
    datblock1 <- rbind(datblock[1:(missind-1),],c(dvd,c(1:10)[missind],0,0,0,0,0,0,0),
                       datblock[missind:9,],deparse.level = 0)
    datblock <- datblock1
    missvals <- rbind(missvals,c(dvd,missind))
  }
  datarr <- rbind(datarr,colSums(datblock),deparse.level = 0)
}
vPHdat <- data.frame(datarr[,-2])
vPHdat$DVD <- DVDun
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
vnonPHdat <- data.frame(head(read.csv('non-PHdata.csv',header=T),-1)) ### virologic data from clinical labs
DVD <- as.matrix(vnonPHdat$DVD)
selind <- which(DVD >= as.numeric(paste(fromyr,fromwk,sep="")) & DVD <= as.numeric(paste(toyr,towk,sep="")))
vnonPHdat <- vnonPHdat[selind,]
#########################################################################################
#########################################################################################
datarr <- c()
missvals <- c()
for (dvd in DVDun){
  selind <- which(vnonPHdat$DVD==dvd)
  datblock <- vnonPHdat[selind,]
  datblock <- datblock[order(datblock$HHS),]
  missind <- which(!c(1:10) %in% datblock$HHS)
  if (length(missind)>0){
    datblock1 <- rbind(datblock[1:(missind-1),],c(dvd,c(1:10)[missind],0,0,0,0,0,0,0),
                       datblock[missind:9,],deparse.level = 0)
    datblock <- datblock1
    missvals <- rbind(missvals,c(dvd,missind))
  }
  datarr <- rbind(datarr,colSums(datblock),deparse.level = 0)
}
vnonPHdat <- data.frame(datarr[,-2])
vnonPHdat$DVD <- DVDun
#########################################################################################
#########################################################################################
#########################################################################################
###  Read in ILInet data  ###############################################################
#########################################################################################
ILIdat <- data.frame(read.csv('ILINet.csv',header=T)) ### ILInet data
begtmselind <- head(which(ILIdat$YEAR==fromyr & ILIdat$WEEK==fromwk),1)
endtmselind <- tail(which(ILIdat$YEAR==toyr & ILIdat$WEEK==towk),1)
ILIdat <- ILIdat[begtmselind:endtmselind,]

datarr <- c()
for (dvd in as.numeric(DVDun)){
  datblock <- NULL
  for (reg in 1:10){
    selind <- which(ILIdat$YEAR==as.numeric(substr(dvd,1,4)) & ILIdat$WEEK==as.numeric(substr(dvd,5,6)) &
                      ILIdat$HHS==reg)
    datrow <- ILIdat[selind,]
    datrow <- c(DVD=dvd,ilitot=as.numeric(datrow$ILITOTAL),
                  numprov=as.numeric(datrow$NUM..OF.PROVIDERS),
                  total=as.numeric(datrow$TOTAL.PATIENTS)) ## converting categorical to numeric
    datblock <- rbind(datblock,datrow,deparse.level = 0)
  }
  datarr <- rbind(datarr,colSums(datblock),deparse.level = 0)
}
ILIdat <- data.frame(datarr)
ILIdat$DVD <- DVDun
#########################################################################################
#########################################################################################
mmwrdat0<-read.table("mmwrweeks.txt",header=T)
mmwrdat<-mmwrdat0[which((mmwrdat0$wyear==fromyr&mmwrdat0$week>=fromwk)|(mmwrdat0$wyear>fromyr&mmwrdat0$wyear<toyr)|(mmwrdat0$wyear==toyr&mmwrdat0$week<=towk)),]
#########################################################################################
#########################################################################################
begweek <- head(mmwrdat$studywk,1)
endweek <- tail(mmwrdat$studywk,1)
#########################################################################################
#########################################################################################
strtwk<-as.Date(mmwrdat$mmwrstrt,format="%m/%d/%Y")
endwk<-as.Date(mmwrdat$mmwrend,format="%m/%d/%Y")

censusdates<-as.Date(c("2011-07-01","2012-07-01","2013-07-01","2014-07-01","2015-07-01","2016-07-01"))

# function for getting weeknumbers for census estimates
weekenum <- function(censusdate){
  as.numeric(mmwrdat$studywk[which(strtwk<=censusdate&endwk>=censusdate)])
}

nodes <- as.numeric(sapply(censusdates,weekenum))-head(mmwrdat$studywk,1)+1
#as.numeric(nodes)

### function to create weekly pop numbers
spf <- function(k,x,agepop){
  fn <- splinefun(nodes,as.numeric(agepop[k,]),method="natural")
  fn(x)
}

mmwrdvd <- mmwrdat$dvdweek
weeks <- mmwrdat$studywk
weeks <- weeks-weeks[1]+1
#########################################################################################
###  creating virol surveillance data   #################################################
#########################################################################################
infilersv <- "RSV.csv"
rsvdat <- read.csv(infilersv)
rsvdat$DVD <- sapply(rsvdat$RepWeekCode,function(x) as.numeric(paste0('20',x)))
rsvdat <- rsvdat[which(rsvdat$DVD >= head(DVDun,1) & rsvdat$DVD <= tail(DVDun,1)),]
rsvdat <- rsvdat[,c(7,6,4,5)]

selind <- which(rsvdat$DVD >= as.numeric(DVDun[1]) & rsvdat$DVD <= as.numeric(tail(DVDun,1)))
rsvdat <- rsvdat[selind,]

DVD <- HHS <- RSVtest <- RSVpos <- NULL

for (wk in DVDun){
  selindwk <- which(rsvdat$DVD==wk)
  if (length(selindwk) > 1){
    RSVtest <- c(RSVtest,sum(rsvdat$RSVtest[selindwk]))
    RSVpos <- c(RSVpos,sum(rsvdat$RSVpos[selindwk]))
  } else {
    RSVtest <- c(RSVtest,0)
    RSVpos <- c(RSVpos,0)
  }
}

rsvdat <- data.frame(DVD=DVDun,RSVtest,RSVpos)
#########################################################################################
#########################################################################################

#
datarr <- NULL
#
DVD <- vnonPHdat$DVD
nonPHA <-
  vnonPHdat$totalA
spec <-
  sapply(vnonPHdat$total, function(x)
    max(x,1)) # to avoid division by 0

propAH3 <-
  vPHdat$H3N2 / sapply(vPHdat$H3N2 + vPHdat$H1N1pdm, function(x)
    max(x,1))

propAH1pdm <-
  vPHdat$H1N1pdm / sapply(vPHdat$H3N2 + vPHdat$H1N1pdm, function(x)
    max(x,1))

AH3num <- nonPHA * propAH3
AH1Pnum <- nonPHA * propAH1pdm
Bnum <- vnonPHdat$totalB
### Creating percent-positive indicator
AH3pp <- AH3num / spec
AH1Ppp <- AH1Pnum / spec
Bpp <- Bnum / spec
#########################################################################################
###  RSV data prep  #####################################################################
#########################################################################################
rsvspec <- rsvdat$RSVtest
rsvspec <- sapply(rsvspec, function(x) ifelse(x==0,1,x))
rsvpos <- rsvdat$RSVpos

#########################################################################################
###  ILInet data  #######################################################################
#########################################################################################
total <- as.numeric(as.vector(ILIdat$total))
ILI <- as.numeric(as.vector(ILIdat$ilitot))
ILIperc <- ILI / total
#########################################################################################
#########################################################################################
###  Creating flu proxy: multiply pp by ILI%    #########################################
#########################################################################################
#########################################################################################
AH3 <- AH3pp * ILIperc ; AH3 <- ifelse(is.na(AH3),0,AH3)
AH1P <- AH1Ppp * ILIperc ; AH1P <- ifelse(is.na(AH1P),0,AH1P)
B <- Bpp * ILIperc ; B <- ifelse(is.na(B),0,B)

AH3pp <- ifelse(is.na(AH3pp),0,AH3pp)
AH1Ppp <- ifelse(is.na(AH1Ppp),0,AH1Ppp)
Bpp <- ifelse(is.na(Bpp),0,Bpp)

rsvpp <- rsvpos / rsvspec
rsv <- rsvpp * ILIperc
#########################################################################################
#########################################################################################
#########################################################################################

virdatarr <-  data.frame(DVD=DVDun,rsv,AH1P,AH3,B)

colnames(virdatarr) <- c('DVD','RSV','AH1P','AH3','B')
#########################################################################################
#########################################################################################
###  Reading in mortality data                 ##########################################
#########################################################################################
#########################################################################################
infile <- paste0('mort2010_16_HHS_5ag.csv')

mortdat <- read.csv(infile, header=T)
mortdat <- data.frame(mortdat)
###########################################################################################################
### converting the week variable to DVD for compatibility
###########################################################################################################
mortdvd <- sapply(1:length(mortdat$studywk), function(x) {
  selind <- which(mmwrdat0$studywk==mortdat$studywk[x])
  mmwrdat0$dvdweek[selind]})

colnames(mortdat) <- c('DVD','agecat','HHS','rcu','pi','flu','allcause')
mortdat$DVD <- mortdvd
###########################################################################################################
###########################################################################################################
### making sure only relevant time period included
selind <- which(mortdat$DVD >= as.numeric(paste0(fromyr,fromwk)) & mortdat$DVD <= as.numeric(paste0(toyr,towk)))
mortdat <- mortdat[selind,]
mortdat <- data.frame(mortdat)
###################################################################################################
###################################################################################################
nag <- 5 ### # of age groups

poptotdat <- array(0,dim = c(nag,6))

for (HHS in 1:10){
  infile <- paste0('pop_reg',HHSreg,'.dat')
  popdat <- read.table(infile,header = F)
  popdat <- popdat[,-1]
  poptotdat <- poptotdat + popdat
}  

popag <- c()
for (k in 1:nag) {
  poprow <- c()
  for (x in weeks) {
    el <- round(spf(k,x,poptotdat))
    poprow <- c(poprow,el)
  }
  popag <- cbind(popag,poprow,deparse.level = 0)
}
#############################################################################
#############################################################################
DVD <-  rcu <- pi <- flu <- allcause <- AH1P <- AH3 <- B <- RSV <- pop <- age <- NULL
for (wk in DVDun){

  selindwk2 <- which(virdatarr$DVD==wk)
  virdatwk <- virdatarr[selindwk2,]
  
  selindwk3 <- which(mmwrdvd==wk)
  popwk <- popag[selindwk3,]
  
  for (ag in 1:nag){
    selindwk1 <- which(mortdat$DVD==wk&mortdat$agecat==ag)
    mortdatwk <- mortdat[selindwk1,]
    
    rcuel <- sum(mortdatwk$rcu)
    piel <- sum(mortdatwk$pi)
    allcauseel <- sum(mortdatwk$allcause)
    fluel <- sum(mortdatwk$flu)
    
    rcu <- c(rcu,rcuel); pi <- c(pi,piel); 
    flu <- c(flu,fluel); allcause <- c(allcause,allcauseel); 
    age <- c(age,ag)
  }
  DVD <- c(DVD,rep(wk,nag)); 
  pop <- c(pop,popwk); 
  AH1P <- c(AH1P,rep(virdatwk$AH1P,nag)); 
  AH3 <- c(AH3,rep(virdatwk$AH3,nag)); 
  B <- c(B,rep(virdatwk$B,nag)); 
  RSV <- c(RSV,rep(virdatwk$RSV,nag))
}


datcoll <- data.frame(DVD,age,rcu,pi,flu,allcause,pop,AH1P,AH3,B,RSV)
###########################################################################################################
###########################################################################################################
filename <- paste("data_National_2010_16_5ag",".dat",sep = "")

setwd(paste0(bfolder,"RSVdata"))
write.table(datcoll,file = filename,row.names = F,col.names = T,quote = F,sep = "\t")

