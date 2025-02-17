## Analysis is day-by-day

## Define baseline by:
##    Integrate NAIS negative ions 1.5-2nm
##    Calculate noise level for measurement 0:00-4:00
##    noise=mean(abs(diff(conc)))

## Threshold for signal TH = noise*20 (Factor 20 is just to make sure I’m selecting only signal)

## Starting time: Last data point below threshold that has 10 consecutively following measurements above threshold
## Stop time: First data point above threshold that has 10 consecutively following  measurements bellow threshold

m1 <- nais.n
log.dp <- log10(nais.diam[-1])-log10(nais.diam[-28])
log.dp <- c(log.dp,log.dp[27])

a <- data.frame(mapply(`*`,m1[,-c(1,2)],log.dp))
m1 <- cbind.data.frame(m1[,1:2],a)
rm(a,log.dp)

m1$date.start <- m1[,1]-120

#a <- split(m1[,2],list(m1[,ncol(m1)]-as.numeric(m1[,ncol(m1)])%%1800)) #CHECK SPLITTING
m1.ave <- aggregate(m1[,-c(1,ncol(m1))],by=list(m1[,ncol(m1)]-as.numeric(m1[,ncol(m1)])%%360),function(x) ifelse(length(x[!is.na(x)])==0,NA,mean(x,na.rm=T)) )

names(m1.ave)[1] <- "date"
m1.ave[,1] <- m1.ave[,1]+360

if(any(apply(m1.ave[,-1],2,is.nan))){ #if any nan, replace
    a <- apply(m1.ave[,-1],2,function(x) replace(x,is.nan(x),NA))
    a <- cbind.data.frame(m1.ave[,1],data.frame(matrix(unlist(a),nrow=nrow(m1.ave),byrow=F),stringsAsFactors=F))
    names(a) <- names(m1.ave)
    m1.ave <- a
    rm(a)
}


m1.int <- cbind.data.frame(m1.ave[,1],rowSums(m1.ave[,c(7:8)]))
names(m1.int)[2] <- "int.1.44.1.66"

a <- split(m1.int,format(m1.int[,1],"%j"))
npf.timing.ave <- data.frame(matrix(NA,0,2))
names(npf.timing.ave) <- c("b[start.npf, 1]","b[stop.npf, 1]")

for (l in 1:length(a)){
    b <- a[[l]]
    noise <- mean(abs(diff(b[as.numeric(format(b[,1],"%H"))<4 & as.numeric(format(b[,1],"%H"))>=0,2])),na.rm=T)
    if (l==6) noise <- mean(abs(diff(b[as.numeric(format(b[,1],"%H"))>21,2])),na.rm=T)  #Dec 3rd there is NPF but no morning data. Noise estimated on nighttime data.

    start.npf <- NA
    stop.npf <- NA
    if(is.nan(noise)) next ### checks for NPF
    
    th.ratio <- 20 
    if (any(l==15,l==19)) th.ratio <- 10
    th <- th.ratio*noise
    count.fact <- 5
    above.th <- 0 #counts the number of data above th
    b <- b[!is.na(b[,2]),] #remove NA
    
### forward
    for (k in 1:nrow(b)){
        if(b[k,2]>th) {
            above.th <- above.th+1 #count one
###        print(paste("count:", above.th,"row:",k))
            if(above.th > count.fact){ #if it's above count-threshold
                start.npf <- k-count.fact-1 #set the start 
                break} #and exit
        } else {
            above.th <- 0 # if we found a value below th, restart the counter
        }
    }
    
### backward
    th.ratio <- 15
    if (l==15) th.ratio <- 10
    above.th <- 0 #control the sequence of NA
    b <- b[!is.na(b[,2]),]
    for (k in nrow(b):1){
        if(b[k,2]>th) {
            above.th <- above.th+1
###        print(paste("count:", above.th,"row:",k))
            if(above.th > count.fact){
                stop.npf <- k+count.fact+1
                break}
        } else {
            above.th <- 0      
        }
    }
    t <- cbind.data.frame(b[start.npf,1],b[stop.npf,1])
    npf.timing.ave <- rbind.data.frame(npf.timing.ave,t)
    
    png(paste("npf_start_stop_nais_neg",format(b[1,1],"%Y%m%d"),".png",sep=""),res=300,width=8*300,height=8*300)
    plot(b[,2]~b[,1],cex=.3,type="o",ylab="Total count for 1.66+1.44 nm negative ions",xlab="",main=paste("NAIS neg ions ",format(b[1,1],"%Y%m%d"),sep=""))
    points(b[b[,2]>th,2]~b[b[,2]>th,1],cex=.8,pch=20,type="p",col=2)
    if(!is.na(start.npf)) abline(v=b[start.npf,1],lty=2,col=2)
    if(!is.na(stop.npf)) abline(v=b[stop.npf,1],lty=2,col=3)
    dev.off()
    
    rm(count.fact,th.ratio,above.th,b,k,th,noise,stop.npf,start.npf,t)
}

rm(l)
npf.timing.ave <- npf.timing.ave[rowSums(is.na(npf.timing.ave))!=2, ]
rownames(npf.timing.ave) <- 1:nrow(npf.timing.ave)
npf.timing.ave <- npf.timing.ave[-match(c("20141128","20141223","20141212"),format(npf.timing.ave[,1],"%Y%m%d")),] #remove non nucleation days

### IMPORTANT NOTE ###
## The algorithm skipped 4 NPF days: Dec 3, Dec 8, Dec 12 and Dec 16.
## Dec 3 had a missing morning and the algorithm failed (fixed)
## Dec 8 no clue
## Dec 12 and Dec 16 had a very low growth (or nothing) (manual timing)
## all statistics (polarplots, mean event etc...) do not take into account these events. Event of 12 and 16 is no problem (too little). Event of Dec 3 no problem (half day was missing). Event 8 Dec, too bad...

## NPF times for 2 and 16 dec
a <- split(psm.6min,format(psm.6min[,1],"%j"))
npf.timing.psm <- data.frame(matrix(NA,0,2))
names(npf.timing.psm) <- c("b[start.npf, 1]","b[stop.npf, 1]")

for (i in length(a)){
#b <- psm.6min[2471:2614,] #take only 16 dic data
#b <- psm.6min[455:598,] #take only 2 dic data

    b <- a[[i]]
    noise <- mean(abs(diff(b[as.numeric(format(b[,1],"%H"))<4 & as.numeric(format(b[,1],"%H"))>=0,7])),na.rm=T)

    start.npf <- NA
    stop.npf <- NA
    if(is.nan(noise)) next ### checks for NPF

    th.ratio <- 40 
    th <- th.ratio*noise
    count.fact <- 5 
    above.th <- 0 #counts the number of data above th
    b <- b[!is.na(b[,7]),] #remove NA
    
### forward
    for (k in 1:nrow(b)){
        if(b[k,7]>th) {
            above.th <- above.th+1 #count one
###        print(paste("count:", above.th,"row:",k))
            if(above.th > count.fact){ #if it's above count-threshold
                start.npf <- k-count.fact-1 #set the start 
                break} #and exit
        } else {
            above.th <- 0 # if we found a value below th, restart the counter
        }
    }
    
### backward
    th.ratio <- 40
    th <- th.ratio*noise
    above.th <- 0 #control the sequence of NA
    b <- b[!is.na(b[,7]),]
    for (k in nrow(b):1){
    if(b[k,7]>th) {
        above.th <- above.th+1
###        print(paste("count:", above.th,"row:",k))
        if(above.th > count.fact){
            stop.npf <- k+count.fact+1
            break}
    } else {
        above.th <- 0      
    }
    }

    if (i==2) {
        start.npf <- 60
        stop.npf <- 93}

    if (i==3) {
        start.npf <- 51
        stop.npf <- 92}

    if (i==4) {
        start.npf <- 67
        stop.npf <- 100}
     
    t <- cbind.data.frame(b[start.npf,1],b[stop.npf,1])
    npf.timing.psm <- rbind.data.frame(npf.timing.psm,t)

#    png(paste("npf_start_stop_PSM",format(b[1,1],"%Y%m%d"),".png",sep=""),res=300,width=8*300,height=8*300)
    plot(b[,7]~b[,1],cex=.3,type="o",ylab="Total count for PSM 1.1 nm particles",xlab="",main=paste("PSM 1.1 nm ",format(b[1,1],"%Y%m%d"),sep=""))
    identify(b[,7]~b[,1])

    points(b[b[,7]>th,7]~b[b[,7]>th,1],cex=.8,pch=20,type="p",col=2)
    if(!is.na(start.npf)) abline(v=b[start.npf,1],lty=2,col=2)
    if(!is.na(stop.npf)) abline(v=b[stop.npf,1],lty=2,col=3)
#    dev.off()  

}


rm(l)
npf.timing.psm <- npf.timing.psm[rowSums(is.na(npf.timing.psm))!=2, ]
rownames(npf.timing.psm) <- 1:nrow(npf.timing.psm)


## create continuous list of nucleation time
m1 <- data.frame(seq.POSIXt(as.POSIXct("2014-11-29 00:00:00", tz="Asia/Katmandu"),as.POSIXct("2014-12-26 00:00:00", tz="Asia/Katmandu"),by="min"))

i <- 1
nucl.time <- data.frame(m1[grep(format(npf.timing.ave[i,1],"%Y%m%d"),format(m1[,1],"%Y%m%d")),]) #extract day of nucleation
nucl.time$nucl.time <- nucl.time[,1]-npf.timing.ave[i,1] #add column with nucleation time

for (i in 2:nrow(npf.timing.ave)){
    dum <- data.frame(m1[grep(format(npf.timing.ave[i,1],"%Y%m%d"),format(m1[,1],"%Y%m%d")),]) #extract day of nucleation
#    dum$real.time <- dum[,1] # copy real time of measurements
    dum$nucl.time <- dum[,1]-npf.timing.ave[i,1] #add column with nucleation time
    nucl.time <- rbind.data.frame(nucl.time,dum) #cbind real time, nucleation time and NAIS data
    rm(dum)
}

names(nucl.time)[1] <- "date"


## create continuous list of nucleation time
## NEGATIVE IONS
m1 <- data.frame(seq.POSIXt(as.POSIXct("2014-11-29 00:00:00", tz="Asia/Katmandu"),as.POSIXct("2014-12-26 00:00:00", tz="Asia/Katmandu"),by="min"))

i <- 4
nucl.time.neg <- data.frame(m1[grep(format(npf.timing.ave[i,1],"%Y%m%d"),format(m1[,1],"%Y%m%d")),]) #extract day of nucleation
nucl.time.neg$nucl.time <- nucl.time.neg[,1]-npf.timing.ave[i,1] #add column with nucleation time

for (i in c(5,6,7,17,18,19)){
    dum <- data.frame(m1[grep(format(npf.timing.ave[i,1],"%Y%m%d"),format(m1[,1],"%Y%m%d")),]) #extract day of nucleation
#    dum$real.time <- dum[,1] # copy real time of measurements
    dum$nucl.time <- dum[,1]-npf.timing.ave[i,1] #add column with nucleation time
    nucl.time.neg <- rbind.data.frame(nucl.time.neg,dum) #cbind real time, nucleation time and NAIS data
    rm(dum)
}

names(nucl.time.neg)[1] <- "date"


## create continuous list of nucleation time
## POSITIVE IONS
m1 <- data.frame(seq.POSIXt(as.POSIXct("2014-11-29 00:00:00", tz="Asia/Katmandu"),as.POSIXct("2014-12-26 00:00:00", tz="Asia/Katmandu"),by="min"))

i <- 1
nucl.time.pos <- data.frame(m1[grep(format(npf.timing.ave[i,1],"%Y%m%d"),format(m1[,1],"%Y%m%d")),]) #extract day of nucleation
nucl.time.pos$nucl.time <- nucl.time.pos[,1]-npf.timing.ave[i,1] #add column with nucleation time

for (i in c(2,3,8:16)){
    dum <- data.frame(m1[grep(format(npf.timing.ave[i,1],"%Y%m%d"),format(m1[,1],"%Y%m%d")),]) #extract day of nucleation
#    dum$real.time <- dum[,1] # copy real time of measurements
    dum$nucl.time <- dum[,1]-npf.timing.ave[i,1] #add column with nucleation time
    nucl.time.pos <- rbind.data.frame(nucl.time.pos,dum) #cbind real time, nucleation time and NAIS data
    rm(dum)
}

names(nucl.time.pos)[1] <- "date"

## m1 <- lapply(split(m1,format(m1$date,"%Y-%m-%d")),function(x){
##     f <- data.frame(seq.POSIXt(round.POSIXt(min(x$date),"day"),round.POSIXt(max(x$date),"day")-60,by="min"))
##     f$nucl.time.approx <- f[,1]-x$date[which(x$nucl.time==0)]
##     f
## }
## )


## m1.ave <- nais.n.6min

## m1.ave.tnorm <- m1.ave[,0]
## m1.ave.tnorm$real.time <- m1.ave[,1]

## for (i in 1:nrow(npf.timing.ave)){
##     dum <- m1.ave[grep(format(npf.timing.ave[i,1],"%Y%m%d"),format(m1.ave[,1],"%Y%m%d")),] #extract day of nucleation
##     dum$real.time <- dum[,1] # copy real time of measurements
##     dum[,1] <- dum[,1]-npf.timing.ave[i,1] #add column with nucleation time
##     m1.ave.tnorm <- rbind.data.frame(m1.ave.tnorm,dum) #cbind real time, nucleation time and NAIS data
## }

## #720 * round(as.numeric(m1.ave.tnorm[1:10,1])/720)

## ## b <- lapply(split(m1.ave.tnorm,m1.ave.tnorm[,1]),function(x){#group NAIS data according to the diff time (i.e. seconds before/after nucleation start) and average
## ##     colMeans(x[,-c(1,ncol(x))],na.rm=T)
## ##     }
## ##     )

## ## b1 <- as.data.frame(do.call("rbind",b))
## ## b1$nuc.time <- as.numeric(names(split(m1.ave.tnorm,m1.ave.tnorm[,1]))) #add nucleation time
## ## nais.n.mean.tnorm <- b1

## CONTROL NPF TIMING 26 Feb 2019
npf.timing.fixed.2019 <- npf.timing.fixed

plot(nucl.time[,2]~nucl.time[,1],type="o",cex=.1,ylim=c(-20000,20000))
abline(h=0)
abline(v=npf.timing.fixed.2019[,1],col=2)

plot(nucl.time.neg[,2]~nucl.time.neg[,1],type="o",cex=.1,ylim=c(-20000,20000))
abline(h=0)
abline(v=npf.timing.fixed.2019[,1],col=2)

plot(nucl.time.pos[,2]~nucl.time.pos[,1],type="o",cex=.1,ylim=c(-20000,20000))
abline(h=0)
abline(v=npf.timing.fixed.2019[,1],col=2)

## CHANGE SOME NPF TIMING on 26 Feb 2019
## manually changed NPF start time on 9/12 to 9:33. It was 10:25 in nucl.time and 11:40 on npf.time
## manually changed NPF start time on 4/12 to 8:15. It was 9:21 in nucl.time and 10:50 on npf.time
##
png(paste("npf_start_stop_nais_neg20141209_bis.png",sep=""),res=300,width=8*300,height=8*300)
plot(m1.int[2370:2611,2]~m1.int[2370:2611,1],type="o",cex=.2)
identify(m1.int[2370:2611,2]~m1.int[2370:2611,1],type="o",cex=.2)
abline(v=m1.int[2370+97,1],col=2)
dev.off()
npf.timing.fixed.2019[11,1] <- m1.int[2370+97,1]

png(paste("npf_start_stop_nais_neg20141204_bis.png",sep=""),res=300,width=8*300,height=8*300)
plot(m1.int[1175:1410,2]~m1.int[1175:1410,1],type="o",cex=.2)
identify(m1.int[1175:1410,2]~m1.int[1175:1410,1])
abline(v=m1.int[1175+79,1],col=2)
dev.off()
npf.timing.fixed.2019[6,1] <- m1.int[1175+79,1]

## create a new continuous list of nucleation time
## NEGATIVE IONS
m1 <- data.frame(seq.POSIXt(as.POSIXct("2014-11-29 00:00:00", tz="Asia/Katmandu"),as.POSIXct("2014-12-26 00:00:00", tz="Asia/Katmandu"),by="sec"))
nucl.time.neg.2019 <- data.frame(matrix(ncol=2,nrow=0)) ## empty data frame

for (i in which(npf.timing.fixed.2019[,3]=="neg")){
    dum <- data.frame(m1[grep(format(npf.timing.fixed.2019[i,1],"%Y%m%d"),format(m1[,1],"%Y%m%d")),]) #extract day of nucleation
#    dum$real.time <- dum[,1] # copy real time of measurements
    dum$nucl.time <- dum[,1]-npf.timing.fixed.2019[i,1] #add column with nucleation time
    nucl.time.neg.2019 <- rbind.data.frame(nucl.time.neg.2019,dum) #cbind real time, nucleation time and NAIS data
    rm(dum)
}

names(nucl.time.neg.2019)[1] <- "date"

## create a new continuous list of nucleation time
## POSITIVE IONS
m1 <- data.frame(seq.POSIXt(as.POSIXct("2014-11-29 00:00:00", tz="Asia/Katmandu"),as.POSIXct("2014-12-26 00:00:00", tz="Asia/Katmandu"),by="sec"))
nucl.time.pos.2019 <- data.frame(matrix(ncol=2,nrow=0)) ## empty data frame

for (i in which(npf.timing.fixed.2019[,3]=="pos")){
    dum <- data.frame(m1[grep(format(npf.timing.fixed.2019[i,1],"%Y%m%d"),format(m1[,1],"%Y%m%d")),]) #extract day of nucleation
#    dum$real.time <- dum[,1] # copy real time of measurements
    dum$nucl.time <- dum[,1]-npf.timing.fixed.2019[i,1] #add column with nucleation time
    nucl.time.pos.2019 <- rbind.data.frame(nucl.time.pos.2019,dum) #cbind real time, nucleation time and NAIS data
    rm(dum)
}

names(nucl.time.pos.2019)[1] <- "date"

