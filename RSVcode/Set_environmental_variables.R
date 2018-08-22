###################################################################################################
### Set directory path here ... remember to use forward slashes ... ###############################
###################################################################################################
bfolder <- "C:/Users/VOR1/Documents/Git projects/RSV_Excess_Mortality/"  #this is an example ...

#bfolder <- "C:/Users/.../.../" ### put "/" at the end
###################################################################################################
###  Define censdatafile (census data file) and filename (name of analysis file to be created)#####
###################################################################################################
filename <- paste0('data_National_2010_16_5ag.dat')
###################################################################################################
###  Now we are defining all time variables; make sure mmrweeks is present ########################
###################################################################################################
fromyr <- 2010
toyr <- 2016
fromwk <- 40
toseaswk <- 17
towk <- 52

nseas <- toyr - fromyr
### Define 'global N'

setwd(paste0(bfolder,'RSVData'))

mmwrdat<-read.table("mmwrweeks.txt",header=T)
mmwrdat<-mmwrdat[which((mmwrdat$wyear==fromyr&mmwrdat$week>=fromwk)|(mmwrdat$wyear>fromyr&mmwrdat$wyear<toyr)|(mmwrdat$wyear==toyr&mmwrdat$week<=towk)),]

DVDun <- mmwrdat$dvdweek

N <- length(mmwrdat$dvdweek)
seasbeg <- which(mmwrdat$week==fromwk)
seasend <- which(mmwrdat$week==toseaswk)
time <- (1:N)/N

for (k in 1:6){
  assign(paste("seas",k,sep=""),c(max(seasbeg[k],3),min(seasend[k],N)))
}
###################################################################################################
###  Length of burn-in period nadapt and the number of iterations used niter ######################
###################################################################################################
nadapt <- 1000
niter <- 5000
###################################################################################################
### Defining other variables relevant the MCMC; make sure that model files exist ##################
###################################################################################################
ag <- 5
nknots <- 4

qual1 <- 'flumort'
model1.file <- paste0('model NS ',qual1,'.txt')

qual2 <- 'vir'
model2.file <- paste0('model NS ',qual2,'.txt')

qual3 <- 'flumort log-link'
model3.file <- paste0('model NS ',qual3,'.txt')

