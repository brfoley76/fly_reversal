###################
## Read in the raw data
## Process it
## Reshape it for analysis

setwd("/home/brad/Documents/reversal")
library("lme4")
library("arm")
library("nplr")
library("nls2")
library("plyr")
library("car")

rev=read.csv("fly_reversal.csv")
aversive=rep(0, nrow(rev))
aversive[which(rev[,"treat"]=="gc" & rev[,"juice"]=="g")]=1
aversive[which(rev[,"treat"]=="pc" & rev[,"juice"]=="p")]=1
rev=cbind(rev, aversive)
revID=paste(rev[,"geno"], rev[,"exp"], sep="")
revDay=paste(rev[,"mon"], rev[,"day"], sep="_")
pop=rep(0, nrow(rev))
pop[which(rev[,"geno"]=="X")]=1
pop[which(rev[,"geno"]=="Y")]=1
pop[which(rev[,"geno"]=="Z")]=1
rev=cbind(rev, revID, revDay, pop)

#Analysis 1
#Take the raw data, and aggregate it by time, bitterness, and period.

aggDat=ddply(rev, .(exp, geno, revID, revDay, pop, period, min, aversive), summarise, sum=sum(n_flies))

wideAgg=reshape(aggDat, v.names="sum", 
                idvar=c("exp", "geno", "period", "min"), 
                direction="wide", timevar="aversive")

totAgg=c(wideAgg[,"sum.0"]+wideAgg[,"sum.1"])
wideAgg=cbind(wideAgg, totAgg)
index=c(1:nrow(wideAgg))
wideAgg=cbind(index, wideAgg)

maxAgg=ddply(wideAgg, .(exp, geno, revID), summarise, max=max(totAgg, na.rm=T))
wideAgg=merge(wideAgg, maxAgg, by=c("revID", "exp", "geno"), all=T)
wideAgg=wideAgg[order(wideAgg$exp, wideAgg$geno,wideAgg$period,wideAgg$min),]
head(wideAgg)

#In a few cases, almost no flies went on food at all.
# Remove those trials (this really makes little difference in the end)
wideAgg=wideAgg[which(wideAgg[,"revID"]!="Y11"),]
wideAgg=wideAgg[which(wideAgg[,"revID"]!="Z24"),]
#write.csv(wideAgg, "wag.csv", row.names=F)

################################
## Create preference scores

## get the maximum on-food for each period ("scale")
# this is the value preference will be scaled by
# accounts for dead flies, sick flies, escaped flies, or mis-counts
uniExps=unique(wideAgg[,"revID"])
uniExps
scale=rep(0, nrow(wideAgg))
for(i in uniExps){
    max1=max(wideAgg[which(wideAgg[,"revID"]==i & wideAgg[,"period"]==1), "totAgg"], na.rm=T)
    max2=max(wideAgg[which(wideAgg[,"revID"]==i & wideAgg[,"period"]==2), "totAgg"], na.rm=T)
    scale[which(wideAgg[,"revID"]==i & wideAgg[,"period"]==1)]=max1
    scale[which(wideAgg[,"revID"]==i & wideAgg[,"period"]==2)]=max2
}

wideAgg=cbind(wideAgg, scale)
wideAgg[,"oddsYum"]=log((wideAgg[,"sum.0"]+0.5)/(wideAgg[,"scale"]-wideAgg[,"sum.0"]+0.5))
wideAgg[,"oddsPref"]=log((wideAgg[,"sum.0"]+0.5)/(wideAgg[,"sum.1"]+0.5))
wideAgg[,"oddsAvers"]=log((wideAgg[,"sum.1"]+0.5)/(wideAgg[,"scale"]-wideAgg[,"sum.1"]+0.5))
wideAgg=wideAgg[which(wideAgg[,"scale"]>2),]

# get a vector to do time-point comparisons
write.csv(wideAgg, "wide_pref_data.csv", row.names=F)

uniVec=unique(wideAgg[,"min"])
uniVec=uniVec[order(uniVec)]

########################################
#####################################
###### These are the time-point comparisons
######  in Figure 1, aversive
######

#Tests the probability of being on the aversive food
#the p-values of interest are in column 3, the ests are in column 2
contrasts=c("contr.sum", "contr.poly")
uniVec=unique(wideAgg[,"min"])
uniVec=uniVec[order(uniVec)]
out1=c()
for(i in uniVec){
    #print(i)
    data=wideAgg[which(wideAgg[,"min"]==i),]
    myMod=glm(sum.1~period+geno+revDay, data=data, family="poisson")
    est=myMod$coefficients[2]
    sig=Anova(myMod, type=2)
    #pval=sig[1,"Pr(>Chisq)"]
    line=c(i, est, Anova(myMod)$"Pr(>Chisq)"[1],Anova(myMod)$"Pr(>Chisq)"[2],Anova(myMod)$"Pr(>Chisq)"[3])
    out1=rbind(out1, line)
}
out1


#####################################
#####################################
###### These are the time-point comparisons
######  in Figure 1, non-aversive
######

#Tests the probability of being on the non-aversive food
#Tests the probability of being on the aversive food
#the p-values of interest are in column 3, the ests are in column 2

uniVec=unique(wideAgg[,"min"])
uniVec=uniVec[order(uniVec)]
out2=data.frame(min=NULL, est=NULL, p=NULL)
for(i in uniVec){
    data=wideAgg[which(wideAgg[,"min"]==i),]
    myMod=glm(sum.0~period+geno+revDay, data=data, family="poisson")
    est=myMod$coefficients[2]
    sig=Anova(myMod, type=2)
    line=c(i, est, Anova(myMod)$"Pr(>Chisq)"[1],Anova(myMod)$"Pr(>Chisq)"[2],Anova(myMod)$"Pr(>Chisq)"[3])
    out2=rbind(out2, line)
}
out2


#####################################
#####################################
###### This generates figure 1
######

pdf(file="boxed_fly.pdf", width=6.2, height=6.1)
par(xpd=T)
symbs=c(21, 22, 23)
par(mfrow=c(2,2))
par(mar=c(4,4.5,3,1.5))
#TableLabels=c(rep("0.0008", 3), rep("0.013", 3))
boxplot(sum.1~min, data=wideAgg[which(wideAgg[,"period"]==1),], ylim=c(0,6), ylab="n flies",main="training, aversive", cex.lab=1.5)
points(c(1,2,3,4,5,6, 22,23, 24, 25, 26, 27, 28), rep(-1.8, 13), pch="*")
points(c(26, 27, 28), rep(-2, 3), pch="*")
points(c(27, 28), rep(-2.2, 2), pch="*")

boxplot(sum.0~min, data=wideAgg[which(wideAgg[,"period"]==1),], ylim=c(0,6), main="training, non-aversive", cex.lab=1.5)
points(c(2,5, 25, 26, 27, 28), rep(-1.8, 6), pch="*")
points(c(26,27, 28), rep(-2, 3), pch="*")
points(c(27, 28), rep(-2.2, 2), pch="*")

boxplot(sum.1~min, data=wideAgg[which(wideAgg[,"period"]==2),], ylim=c(0,6),  ylab="n flies", xlab="minutes", main="reversal, aversive", cex.lab=1.5)
boxplot(sum.0~min, data=wideAgg[which(wideAgg[,"period"]==2),], ylim=c(0,6), main="reversal, non-aversive", xlab="minutes", cex.lab=1.5)
dev.off()

##########################################################
##########################################################
##########################################################
#####                    #################################
##### Do analysis using  #################################
#####  estimate of       #################################
#####  midpoint from     #################################
#####  loess             #################################
#####                    #################################
##########################################################
##########################################################
uniExps=unique(wideAgg[,"revID"])
uniExps

options(contrasts=c("contr.sum", "contr.poly")) #for proper type III SS 
library(car)

# an identifier for doing per-rep, per period analysis
repIt=paste(wideAgg[,"revID"], wideAgg[,"period"], sep="_")
wideAgg=cbind(wideAgg, repIt)


getMid=function(data, metric, span=NULL){
    #this function takes in a chunk of timecourse data=data
    # the name of a preference measure=metric
    # and an optional span/smoothing parameter = span.
    # note that the R default span is 0.75
    if(is.null(span)){
        span=0.75
    }
    #fit a loess curve
    #write.csv(data, "what_is_it.csv", row.names=F)
    fitCurve=loess(data[,metric]~log(data[,"min"]+1), span=span)
    curvAture=fitCurve$fitted
    len=length(curvAture)
    # Find the first local minimum, and take that as "baseline" == "minnie"
    # The ultimate learned preference will be the maximum estimated preference == "maxie"
    # The first point at which estimated preference is halfway between "minnie" and "maxie" is the learning rate proxy
    
    localMinVec=c(curvAture[2:(len-1)]<curvAture[1:(len-2)] & curvAture[2:(len-1)]<curvAture[3:(len)])
    localMinVec=c(curvAture[1]<curvAture[2], localMinVec, curvAture[len]<curvAture[(len-1)])
    trunc=c(1:len)
    trunc=trunc[localMinVec]
    trunc=trunc[1]
    curvRemain=curvAture[trunc:len]
    #minnie=min(curvRemain)
    minnie=curvAture[trunc]
    maxie=max(curvRemain)
    middie=(minnie+maxie)/2
    #nadir=which(curvAture==minnie)    
    # create 2 vectors that we can use to grab the flanking times/prefValues
    # and get an estimate of the midpoint
    preMidScreenRight=c(curvAture>middie)
    #I want to omit anything before the nadir from consideration
    preMidScreenRight[1:trunc]=FALSE
    preMidScreenLeft=preMidScreenRight[2:length(preMidScreenRight)]
    preMidScreenLeft=c(preMidScreenLeft, TRUE)
    leftPref=curvAture[preMidScreenLeft]
    leftPref=leftPref[1]
    rightPref=curvAture[preMidScreenRight]
    rightPref=rightPref[1]
    leftTime=data[preMidScreenLeft,"min"]
    leftTime=log(leftTime[1]+1)
    rightTime=data[preMidScreenRight,"min"]
    rightTime=log(rightTime[1]+1)
    intervalScale=(middie-leftPref)/(rightPref-leftPref)
    timeEst=intervalScale*(rightTime-leftTime)+leftTime
    
    #print(data[1,"repIt"])
    
    return(data.frame(pref=middie, time=timeEst))
}

#############################################################################
#############################################################################
### The following is a god-awful crude way to identify the optimal smoothing
### Basically, I just took a vector of arbitrary length, of span parameters
### did loess regression using those span parameters 
### and looked to see which ones gave the most explanatory midpoint-estimates
### I used "pseudorsq", or proportion deviance explained, as the choice metric

par(mfrow=c(1,1))
int1=0.4 #left limit of the span values
int2=2.2 #right limit of the span values
testVec=c(1:100) #a vector of 25 test-vlues
testVec=(testVec/100)*(int2-int1) #scale to the limits
testVec=testVec+int1
res=data.frame("span"=numeric(), "rsq"=numeric()) 
for(i in 1:100){
    span=testVec[i]
    # run the curve-getting function with a given span
    realRes=ddply(wideAgg[which(wideAgg[,"max"]>2),], .(geno, revID, pop,revDay,period,repIt), function(data) getMid(data, "oddsYum", span=span))
    # do the full model using the span metric
    modFit=glm(time~period*geno+revDay, data=realRes)
    signess=summary(modFit)        
    # get an appropriate fit statistic
    pseudorsq=1-signess$deviance/signess$null.deviance
    res[i,"span"]=span    
    res[i,"rsq"]=pseudorsq

}
res
# Check: how does it look?
plot(res[,"rsq"]~res[,"span"])
# Answer: pretty much flat, from about span=1.2 on.

#######################################
#######################################
###
### Preference analysis, empirical data
###
#######################################
#######################################
span=1.45

realRes=ddply(wideAgg[which(wideAgg[,"scale"]>2),], .(geno, revID, pop,revDay,period,repIt), function(data) getMid(data, "oddsYum", span=span))

modFit=glm(time~as.factor(period)+as.factor(revDay)+as.factor(geno), data=realRes)
summary(modFit)    
logLik(modFit)
signess=Anova(modFit, type=2) #no interaction terms. Also, note we set the contrasts properly at the top: contrasts=c("contr.sum", "contr.poly")
signess

modFit=glm(time~period*geno+revDay, data=realRes) #interaction term
summary(modFit)    
logLik(modFit)
signess=Anova(modFit, type=3)
signess
# the period*geno effect is not even close to significance
# nothing I could do got it anywhere near approaching significance



#########################################################################
#########################################################################
### SanityCheck: look through the plots 
###   and see if the position of the "midpoint" makes sense
###
###   This is the code that generates the figure in the appendix
###
#########################################################################

uniVec=unique(wideAgg[,"repIt"])
i=1

span=1.45
uniVec=uniVec[sample(uniVec, length(uniVec))]
#pdf(file="sanity.pdf", width=12, height=8)
par(mfrow=c(2,4))
for(j in 1:8){
    data=wideAgg[which(wideAgg[,"repIt"]==uniVec[j]),]
    plot(data[,"oddsYum"]~log(data[,"min"]+1),  ylab="preference", xlab="log mins", cex=1.5, pch=19,cex.lab=1.5)
    fitCurve=loess(data[,"oddsYum"]~log(data[,"min"]+1), span=span)
    lines(log(data[,"min"]+1), fitCurve$fitted,col=rgb(0,0,200, max=255), lty=2, cex=2)
    points(realRes[which(realRes[,"repIt"]==uniVec[j]),"time"], realRes[which(realRes[,"repIt"]==uniVec[j]),"pref"], pch=19, col="red", cex=1.8)
}
#dev.off()



########################################
########################################
########################################
#####
##### Permutation analysis to evaluate 
#####   evaluate model fit
#####
########################################
########################################

### Need to get rid of all observations which are missing a morning-night pair

head(realRes)
uniRev1=unique(realRes[which(realRes[,"period"]==1),"revID"])
length(uniRev1)
uniRev2=unique(realRes[which(realRes[,"period"]==2),"revID"])
length(uniRev2)
uniRev=intersect(uniRev1, uniRev2)
length(uniRev)

realResPruned=realRes[which(is.element(realRes[,"revID"], uniRev)),]

#######################################################
#######################################################
###   permutation keeping morning-evening pairs intact
#######################################################
realResPruned=realResPruned[order(realResPruned$period, realResPruned$revID),]
#head(realResPruned)
rownames(realResPruned)=c(1:226) ### R sorts by row name, not index, I learned to my chagrin
realResPermed=realResPruned
modFit=glm(time~geno+period+revDay, data=realResPruned)

signess=Anova(modFit, type=2)
#signess

indexVec=c(1:113) 
permUteVec=indexVec
outLike=data.frame("it"=integer(), "logLik"=numeric(), "rsq"=numeric(), "AIC"=numeric())
outLike
## For validation sake, just 100 times
#### for the paper, 100000
for(i in 1:100){
    sampVec=c(permUteVec, permUteVec+113) ## apply the same permutation to the first period, and 2nd period
    realResPermed[,"time"]=realResPruned[sampVec,"time"]
    modFit=glm(time~geno+period+revDay, data=realResPermed)
    lik=logLik(modFit)
    
    signess=summary(modFit)  
    pseudorsq=1-signess$deviance/signess$null.deviance
    #print(pseudorsq)
    aic=AIC(modFit)
    outLike[i,"it"]=i
    outLike[i,"logLik"]=lik[1]
    outLike[i,"rsq"]=pseudorsq
    outLike[i,"AIC"]=aic
    
    permUteVec=sample(indexVec, length(indexVec), replace=F)
}

### Examine the distribution of (pseudo) r-squareds
par(mfrow=c(1,1))
modFit=glm(time~geno+period+revDay, data=realResPruned)
signess=summary(modFit)
pseudorsq=1-signess$deviance/signess$null.deviance
pseudorsq
nrow(outLike)
nrow(outLike[which(outLike[,"rsq"]>pseudorsq),])
hist(outLike[,"rsq"], breaks=100, xlim=c(0,0.25))
abline(v=pseudorsq, col="red", cex=2 )


#############################################################
#############################################################
###   permutation, breaking learning-reversal pair associations
###   but preserving learning-reversal identity
#############################################################
modFit=glm(time~geno+period+revDay, data=realResPruned)
logLik(modFit)
signess=Anova(modFit, type=2)
signess

indexVec1=c(1:113)
indexVec2=c(114:226)

permUteVec1=indexVec1
permUteVec2=indexVec2
outLike2=data.frame("it"=integer(), "logLik"=numeric(), "rsq"=numeric(), "AIC"=numeric())
## For validation sake, just 100 times
#### for the paper, 100000
for(i in 1:100){
    sampVec=c(permUteVec1, permUteVec2)
    realResPermed[,"time"]=realResPruned[sampVec,"time"]
    modFit=glm(time~geno+period+revDay, data=realResPermed)
    lik=logLik(modFit)
    signess=summary(modFit)  
    pseudorsq=1-signess$deviance/signess$null.deviance
    #print(pseudorsq)
    aic=AIC(modFit)
    outLike2[i,"it"]=i
    outLike2[i,"logLik"]=lik[1]
    outLike2[i,"rsq"]=pseudorsq
    outLike2[i,"AIC"]=aic
    
    permUteVec1=sample(indexVec1, length(indexVec1), replace=F)
    permUteVec2=sample(indexVec2, length(indexVec2), replace=F)
}

### examine distribution of (pseudo) r-squareds
modFit=glm(time~geno+period+revDay, data=realResPruned)
signess=summary(modFit)
pseudorsq=1-signess$deviance/signess$null.deviance
pseudorsq
nrow(outLike2)
nrow(outLike2[which(outLike2[,"rsq"]>pseudorsq),])
hist(outLike2[,"rsq"], breaks=100, xlim=c(0,0.25))
abline(v=pseudorsq, col="red", cex=2 )

#########################################################
#########################################################
##
## permutation analysis destroying every association
##
#########################################################
modFit=glm(time~geno+period+revDay, data=realResPruned)
logLik(modFit)
signess=Anova(modFit, type=2)
signess

indexVec=c(1:226)

permUteVec=indexVec
outLike3=data.frame("it"=integer(), "logLik"=numeric(), "rsq"=numeric(), "AIC"=numeric())

for(i in 1:100){
    sampVec=c(permUteVec)
    realResPermed[,"time"]=realResPruned[sampVec,"time"]
    modFit=glm(time~geno+period+revDay, data=realResPermed)
    lik=logLik(modFit)
    signess=summary(modFit)  
    pseudorsq=1-signess$deviance/signess$null.deviance
    #print(pseudorsq)
    aic=AIC(modFit)
    outLike3[i,"it"]=i
    outLike3[i,"logLik"]=lik[1]
    outLike3[i,"rsq"]=pseudorsq
    outLike3[i,"AIC"]=aic
    
    permUteVec=sample(indexVec, length(indexVec), replace=F)
}

modFit=glm(time~geno+period+revDay, data=realResPruned)
signess=summary(modFit)
pseudorsq=1-signess$deviance/signess$null.deviance
pseudorsq
nrow(outLike3)
nrow(outLike3[which(outLike3[,"rsq"]>pseudorsq),])
hist(outLike3[,"rsq"], breaks=100, xlim=c(0,0.25))
abline(v=pseudorsq, col="red", cex=2 )
