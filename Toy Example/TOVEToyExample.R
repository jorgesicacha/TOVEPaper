### Very simple example for TOVE model ##
library(sp)
library(rgeos)
library(INLA)
library(raster)
library(tiff)
library(spatstat)
library(maptools)

## Region of study ##

studyregion <- shapefile("studyregion1.shp")

## Covariate data ##

cov <- raster("covariate.tif")

## Observed counts ##

## y1 and y2: observed counts in Norway
## y3 and y4: observed counts in Sweden

data <- shapefile("obscounts.shp")

## Fitting the three models proposed 

fit <- function(data,cov){
  
  ## Making the mesh and A-matrix
  datalocs <- data@coords
  mesh = inla.mesh.2d(loc = datalocs, max.edge = c(10000,50000),cutoff = 10000, offset=c(5000, -0.2)) 
  A=inla.spde.make.A(mesh,loc=as.matrix(datalocs))
  
  ## Prior on spatial hyperparameters
  spde1=inla.spde2.pcmatern(mesh=mesh,alpha = 2, prior.range = c(5000,0.05), prior.sigma = c(3,0.01)) 
  spde1_sept=spde1
  spde1_noli=spde1
  spde1_seli=spde1
  ## w_2 in the manuscript
  spde1_noli.alt=inla.spde2.pcmatern(mesh=mesh,alpha = 2, prior.range = c(1000,0.01), prior.sigma = c(5,0.01))  
  
  index.spde1=inla.spde.make.index(name="field.nopt",n.spde=spde1$n.spde)
  index.spde1_sept=inla.spde.make.index(name="field.sept",n.spde=spde1_sept$n.spde)
  index.spde1_noli=inla.spde.make.index(name="field.noli",n.spde=spde1_noli$n.spde)
  index.spde1_seli=inla.spde.make.index(name="field.seli",n.spde=spde1_seli$n.spde)
  index.spde1_noli.alt=inla.spde.make.index(name="field.noli.alt",n.spde=spde1_noli.alt$n.spde)
  
  cov_locs <- extract(cov,datalocs)
  m <- length(data)
  
  ## Making the stacks for each data type ##
  stack_nopt=inla.stack(data=list(y=cbind(data$y1_obs,NA,NA,NA)), 
                        effects=list(index.spde1,data.frame(intercept_nopt=rep(1,m),
                                                            cov_nopt = cov_locs,
                                                            beta1_nopt = rep(1,m))),
                        A=list(A,1),tag="nopt")
  
  stack_sept=inla.stack(data=list(y=cbind(NA,data$y2_obs,NA,NA)), 
                        effects=list(index.spde1_sept,data.frame(intercept_nopt=rep(1,m),
                                                                 zeta1 = rep(1,m),
                                                                 cov_nopt = cov_locs,
                                                                 beta1_sept = rep(1,m))),
                        A=list(A,1),tag="sept")
  
  stack_noli=inla.stack(data=list(y=cbind(NA,NA,data$y3_obs,NA)), 
                        effects=list(index.spde1_noli,index.spde1_noli.alt,data.frame(intercept_nopt=rep(1,m),
                                                                                      zeta2 = rep(1,m),
                                                                                      cov_nopt = cov_locs,
                                                                                      beta1_noli = rep(1,m))),
                        A=list(A,A,1),tag="noli")
  
  stack_seli=inla.stack(data=list(y=cbind(NA,NA,NA,data$y4_obs)), 
                        effects=list(index.spde1_seli,data.frame(intercept_nopt=rep(1,m),
                                                                 zeta3 = rep(1,m),
                                                                 cov_nopt = cov_locs,
                                                                 beta1_seli = rep(1,m))),
                        A=list(A,1),tag="seli")
  
  #Join stacks.
  stack.obs=inla.stack(stack_nopt,stack_sept,stack_noli,stack_seli) 
  
  ## prior distribution for random intercepts
  prec.prior <- list(prec = list(prior = "loggamma", param = c(1, 0.00005)))

  ## Model 1 specification 
  formula1 <- y ~ -1 + intercept_nopt + 
    cov_nopt +
    f(zeta1, model="iid",hyper=prec.prior) +
    f(zeta2, model="iid",hyper=prec.prior) +
    f(zeta3, model="iid",hyper=prec.prior) +
    f(field.nopt,model=spde1)+
    f(field.sept, copy="field.nopt",fixed=TRUE) +
    f(field.noli, copy="field.nopt",fixed=TRUE) +
    f(field.seli, copy="field.nopt",fixed=TRUE)
  
  ## Model 2 specification
  formula2 <- y ~ -1 + intercept_nopt + 
    cov_nopt +
    f(zeta1, model="iid",hyper=prec.prior) +
    f(zeta2, model="iid",hyper=prec.prior) +
    f(zeta3, model="iid",hyper=prec.prior) +
    f(field.nopt,model=spde1)+
    f(field.sept, copy="field.nopt",fixed=FALSE) +
    f(field.noli, copy="field.nopt",fixed=FALSE) +
    f(field.seli, copy="field.nopt",fixed=FALSE)
  
  ## Model 3 specification
  formula3 <- y ~ -1 + intercept_nopt + 
    cov_nopt +
    f(zeta1, model="iid",hyper=prec.prior) +
    f(zeta2, model="iid",hyper=prec.prior) +
    f(zeta3, model="iid",hyper=prec.prior)+
    f(field.nopt,model=spde1)+
    f(field.sept, copy="field.nopt",fixed=FALSE) +
    f(field.noli, copy="field.nopt",fixed=FALSE) +
    f(field.seli, copy="field.nopt",fixed=FALSE)+
    f(field.noli.alt,model=spde1_noli.alt)
  
  
  ## Our models
  mod1=inla(formula1,
            family=c("Poisson","Poisson","Poisson","Poisson"),
            data=inla.stack.data(stack.obs),control.predictor=list(A=inla.stack.A(stack.obs),
                                                                   compute=TRUE),
            control.compute = list(cpo = TRUE, dic = TRUE,waic=TRUE,config=TRUE),
            verbose = F)
  mod2=inla(formula2,
            family=c("Poisson","Poisson","Poisson","Poisson"),
            data=inla.stack.data(stack.obs),control.predictor=list(A=inla.stack.A(stack.obs),
                                                                   compute=TRUE),
            control.compute = list(cpo = TRUE, dic = TRUE,waic=TRUE,config=TRUE),
            verbose = F)
  mod3=inla(formula3,
            family=c("Poisson","Poisson","Poisson","Poisson"),
            data=inla.stack.data(stack.obs),control.predictor=list(A=inla.stack.A(stack.obs),
                                                                   compute=TRUE),
            control.compute = list(cpo = TRUE, dic = TRUE,waic=TRUE,config=TRUE),
            verbose = F)
  
  return(list(mod1=mod1,
              mod2=mod2,
              mod3=mod3))
}

results <- fit(data,cov)





