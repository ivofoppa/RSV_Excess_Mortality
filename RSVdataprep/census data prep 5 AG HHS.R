setwd(paste0(bfolder,'RSVdata'))

censdat <- data.frame(read.csv(censdatafile, header = T))

st_list <- list(c('Connecticut','Maine','Massachusetts','New Hampshire','Rhode Island','Vermont'),
                c('New Jersey','New York','Puerto Rico','The Virgin Islands'), 
                c('Delaware','District of Columbia','Maryland','Pennsylvania','Virginia','West Virginia'),
                c('Alabama','Florida','Georgia','Kentucky','Mississippi','North Carolina','South Carolina','Tennessee'),
                c('Illinois','Indiana','Michigan','Minnesota','Ohio','Wisconsin'),
                c('Arkansas','Louisiana','New Mexico','Oklahoma','Texas'),
                c('Iowa','Kansas','Missouri','Nebraska'),
                c('Colorado','Montana','North Dakota','South Dakota','Utah','Wyoming'),
                c('Arizona','California','Hawaii','Nevada','American Samoa','Commonwealth of the Northern Mariana Islands','Federated States of Micronesia','Guam','Marshall Islands','Republic of Palau'),
                c('Alaska','Idaho','Oregon','Washington'))


age_list <- list(c(0:4),c(5:17),c(18:49),c(50:64),c(65:85))

for (HHSreg in 1:10){
  agepop <- c()
  regselind <- which(censdat$NAME %in% st_list[[HHSreg]] & censdat$SEX == 0)
  tempdat <- censdat[regselind,]
  for (ag in 1:5){
    
    agselind <- which(tempdat$AGE %in% age_list[[ag]]&tempdat$SEX==0)
    agepop <- rbind(agepop,c(sum(tempdat$POPESTIMATE2010[agselind]),
                             sum(tempdat$POPESTIMATE2011[agselind]),
                             sum(tempdat$POPESTIMATE2012[agselind]),
                             sum(tempdat$POPESTIMATE2013[agselind]),
                             sum(tempdat$POPESTIMATE2014[agselind]),
                             sum(tempdat$POPESTIMATE2015[agselind]),
                             sum(tempdat$POPESTIMATE2016[agselind])),deparse.level = 0)
  }
  setwd(paste0(bfolder,'RSVdata'))
  outfilename <- paste0('pop_reg',HHSreg,'.dat')
  write.table(agepop,outfilename,quote = F,col.names = F,row.names = F)
}
