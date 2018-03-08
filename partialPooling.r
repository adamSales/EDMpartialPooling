library(ggplot2)
library(dplyr)
library(rstanarm)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

#### benchmarks
dat <- read.csv('ASSISTments_dataset_22_experiments.csv')

dat$sampleSize <- as.vector(table(dat$problem_set)[as.character(dat$problem_set)])

sdat <- split(dat,as.factor(dat$User.ID))

studTable <- table(dat$User.ID)

for(stud in names(studTable)[studTable>1])
    sdat[[stud]] <- sdat[[stud]][which.min(sdat[[stud]]$sampleSize),]

dat <- do.call('rbind',sdat)

dat$problem_set <- as.factor(dat$problem_set)


partialPooling <- function(Yname,covs=FALSE){
    fam <- ifelse(length(unique(na.omit(dat[[Yname]])))==2,binomial,gaussian)
    form <- paste(Yname,'~condition+(condition|problem_set)')
    if(covs) form <- paste(form,'+Prior.Percent.Correct+Prior.Percent.Completion')

    ppMod <- stan_glmer(as.formula(form),data=dat,family=fam)

    draws <- as.matrix(ppMod)

    ddd <- draws[,paste0('b[conditionE problem_set:',levels(dat$problem_set),']')]

    effs <- sweep(ddd,1,draws[,'conditionE'],'+')

    results <- cbind(est=colMeans(effs),se=apply(effs,2,sd))
    results <- cbind(results, pval=2*pnorm(-abs(results[,'est']/results[,'se'])))

    rownames(results) <- gsub('b\\[conditionE problem_set\\:|\\]','',rownames(results))

    list(results=results,mod=ppMod)
}


noPooling <- function(Yname,covs=FALSE){
    fam <- ifelse(length(unique(na.omit(dat[[Yname]])))==2,binomial,gaussian)
    form <- paste(Yname,'~condition')
    if(covs) form <- paste(form,'+Prior.Percent.Correct+Prior.Percent.Completion')

    out <- do.call('rbind',lapply(levels(dat$problem_set),function(ps)
        summary(glm(as.formula(form),family=fam,data=dat,subset=problem_set==ps))$coef['conditionE',]))
    rownames(out) <- levels(dat$problem_set)
    out
}

plotit <- function(np,pp){
    pp <- pp$results
    plotDat <- as.data.frame(rbind(pp[order(np[,1]),],np[order(np[,1]),-3]))
    names(plotDat) <- c('est','se','p')
    plotDat$ps <- c(1:22,1:22+0.2)
    plotDat$type=rep(c('Partial Pooling','No Pooling'),each=22)

    p1 <- ggplot(plotDat,aes(ps,est,fill=type,color=type))+geom_point()+
        geom_errorbar(aes(ymin=est-2*se,ymax=est+2*se))+
        geom_hline(yintercept=0,linetype='dotted')+
        xlab('Experiment')+
        labs(color=NULL,fill=NULL)+
        theme(legend.position='bottom',text=element_text(size=7.5),
              axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank())

    plotDat2 <- data.frame(pp=pp[,'se'],np=np[,'Std. Error'],n=as.vector(table(dat$problem_set)))

    lims <- range(c(pp[,2],np[,2]))

    p2 <- qplot(np, pp,data=plotDat2,ylim=lims,xlim=lims,xlab='SEs: No Pooling', ylab='SEs: Partial Pooling',size=n)+theme(legend.pos='none',aspect.ratio=1)+
        geom_abline(intercept=0,slope=1)+
        scale_size_continuous(range=c(0.3,3))

    list(p1,p2)
}

ppComplete <- partialPooling('complete')
npComplete <- noPooling('complete')
completePlots <- plotit(npComplete,ppComplete)

ggsave(completePlots[[1]]+labs(y='Estimated Effect\n(Log Odds Ratio)',title='Outcome: Complete'),file='completeEst1.jpg',width=3,height=3)
ggsave(completePlots[[2]]+labs(title='Outcome: Complete'),file='completeSE1.jpg',width=3,height=3)

npCompleteCovs <- noPooling('complete',TRUE)
ppCompleteCovs <- partialPooling('complete',TRUE)

dat$speed <- ifelse(dat$complete==1,1/(dat$problem_count),0)

ppSpeed <- partialPooling('speed')
npSpeed <- noPooling('speed')

speedPlots <- plotit(npSpeed,ppSpeed)

ggsave(speedPlots[[1]]+labs(y='Estimated Effect\n(1/# Problems)',title='Outcome: Mastery Speed'),file='speedEst1.jpg',width=3,height=3)
ggsave(speedPlots[[2]]+labs(title='Outcome: Mastery Speed'),file='speedSE1.jpg',width=3,height=3)
