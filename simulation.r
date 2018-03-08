library(rstan)
library(rstanarm)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

datSimp <- function(sigTau){
    ngroup <- 20
    ns <- seq(20,by=5,length=ngroup)
    tau <- rnorm(ngroup,0,sigTau)
    group <- rep(1:ngroup,ns)
    groupEff <- rnorm(ngroup)
    Z <- rep(c(0,1),sum(ns)/2)
    Y <- rnorm(sum(ns))+groupEff[group]+Z*tau[group]
    list(data=data.frame(group=as.factor(group),Z=Z,Y=Y),tau=tau)
}

noPooling <- function(data)
    sapply(levels(data[[1]]$group),
           function(grp) summary(lm(Y~Z,data=data[[1]],subset=group==grp))$coef['Z',])

partialPooling <- function(data){
    mod <- stan_lmer(Y~Z+(Z|group),data=data[[1]],iter=1000,warmup=300,chains=6,
                     refresh=-1,open_progress=FALSE)
    mod <- as.matrix(mod)
    ddd <- mod[,paste0('b[Z group:',1:20,']')]
    effs <- sweep(ddd,1,mod[,'Z'],'+')
    ps <- apply(effs,2,function(x) mean(x>0))
    ps <- ifelse(ps>.5,1-ps,ps)
    rbind(est=colMeans(effs),se=apply(effs,2,sd),pval=2*ps)
}

sim <- function(){
    res <- list()
    for(sigTau in c(0,0.1,0.2,0.5,1)){
        print(sigTau)
        datasets <- replicate(500,datSimp(sigTau),simplify=FALSE)
        np <- lapply(datasets,noPooling)
        pp <- lapply(datasets,partialPooling)
        res[[as.character(sigTau)]] <- list(datasets=datasets,np=np,pp=pp)
        save(res,file='simResults.RData')
    }
    res
}

reForm <- function(sig){
    site <- res[[sig]]
    tau <- lapply(site$datasets,function(x) x$tau)
    tau <- do.call('c',tau)
    np <- lapply(site$np,t)
    np <- do.call('rbind',np)[,-3]
    pp <- lapply(site$pp,t)
    pp <- do.call('rbind',pp)
    run <- rep(1:length(site$pp),each=20)
    out <- data.frame(tau,np,pp,run,sig)
    names(out) <- c('tau','NPest','NPse','NPp','PPest','PPse','PPp','run','sig')
    out
}

#### summary stats
## level
lev <- function(run){
    byRun <- sapply(unique(run$run),function(x) with(subset(run,run==x),c(pp=sum(PPp<0.05),np=sum(NPp<0.05))))
    apply(byRun,1,function(x) mean(x>0))
}

se <- function(run){
    c(pp=mean(run$PPse),np=mean(run$NPse))
}

bias <- function(run){
    with(run,c(pp=mean((PPest-tau)/sign(tau)),np=mean((NPest-tau)/sign(tau))))
}

mse <- function(run){
    c(pp=sqrt(mean((run$PPest-run$tau)^2)),np=sqrt(mean((run$NPest-run$tau)^2)))
}

coverage <- function(run){
    with(run, c(pp=mean(abs(PPest-tau)<2*PPse),np=mean(abs(NPest-tau)<2*NPse)))
}

typeSerror <- function(run){
    c(pp=with(subset(run,PPp<0.05), mean((tau>0 & PPest<0)|(tau<0 & PPest>0))),
      np=with(subset(run,NPp<0.05), mean((tau>0 & NPest<0)|(tau<0 & NPest>0))))
}

typeMerror <- function(run){
    c(pp=with(subset(run,PPp<0.05), mean(abs(PPest/tau))),
      np=with(subset(run,NPp<0.05), mean(abs(NPest/tau))))
}


tab1 <- rbind(se=sapply(res,se),bias=sapply(res,bias),mse=sapply(res,mse),coverage=sapply(res,coverage))

sink('table1.tex')
for(bbb in c('se','bias','mse','coverage')){
    ttt <- sapply(res,get(bbb))
    cat('\\multirow{2}{*}{',bbb,'}& Partial Poooling &',paste(round(ttt['pp',],2),collapse='&'),'\\\\ \n')
    cat('& No Pooling &',paste(round(ttt['np',],2),collapse='&'),'\\\\ \n')
    cat('\\hline')
}


plotDat <- lapply(res, function(rr) with(rr,
                                         data.frame(eff=c(tau,tau),
                                                    est=c(NPest,PPest),
                                                    se=c(NPse,PPse),
                                                    run=c(run,run),
                                                    tau=sig,
                                                    type=rep(c('No Pooling','Partial Pooling'),
                                                             each=10000))))

plotDat <- do.call('rbind',plotDat)


ggplot(plotDat,aes(tau,est-eff,fill=type))+geom_boxplot()+xlab(expression(tau))+ylab('Estimation Error')+theme(legend.pos='top')+labs(fill=NULL)
ggsave('simEstError.jpg',height=3,width=3)
