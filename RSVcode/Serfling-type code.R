setwd('C:/Users/VOR1/Documents/Git projects/RSV_Excess_Mortality')

filename <- paste0('data_National_2010_16_5ag.dat')
setwd('../') ## change to root directory
setwd('./RSVdata')
datarr <- data.frame(read.table(file = filename,header = T))
###################################################################################################
###  Serfling-type model, virologic indicators (Thompson et al., 2009)  ###########################
###################################################################################################
agdata <- datarr[which(datarr$age==ag),]

AH1P <- agdata$AH1P
AH3 <- agdata$AH3
B <- agdata$B

mort <- agdata$rcu
pop <- agdata$pop

N <- length(mort)
time <- (1:N)/N
t <- time * N * 2* pi/52.25
###################################################################################################
data <- data.frame(mort,time,t,AH1P,AH3,B,pop)
###################################################################################################
Poimod <- glm(mort ~ time + time^2 + time^3 + sin(t) + cos(t) + AH1P + AH3 + B + offset(log(pop)), data = data,family = poisson(link = 'log'))

dataflu0 <- data.frame(mort,time,t,AH1P=AH1P*0,AH3=AH3*0,B=B*0,pop)
mortflu0 <- predict(Poimod, newdata = dataflu0, type = 'response')

EM <- sum(mort - mortflu0)
