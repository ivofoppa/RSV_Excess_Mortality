###################################################################################################
###  Post-processing of coda files                     ############################################
###  Ivo M. Foppa, Sept. 2018                          ############################################
###################################################################################################
setwd('..')
setwd('./RSVMCMCoutput')
ag <- 5
# nknots <- x ## number of knots per year
# nadapt <- y ## duration of burn-in period
# niter <- z ## number of iterations
fname <- paste0('codaarr',ag,' RSV ',qual1,' ',nknots,' knots ps ',nadapt,' ',round(niter/5*3),'.RData')

load(fname)

codaarr <- data.frame(codaarr)
RSVEM <- codaarr$EMRSVtot
fluEM <- codaarr$EMflutot

RSVci <- round(as.numeric(quantile(RSVEM, probs = c(.5,.025,.975))))
fluci <- round(as.numeric(quantile(fluEM, probs = c(.5,.025,.975))))

cat(paste0(fluci[1],' (',fluci[2],',',fluci[3],')'))
