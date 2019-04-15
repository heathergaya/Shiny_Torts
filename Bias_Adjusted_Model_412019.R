
##check to make sure all the needed packages are installed and install them if need be ####
list.of.packages <- c('rjags', 'runjags', 'mixAK', 'miscF', 'Rdistance', "lme4", "Rcpp")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages,function(x){library(x,character.only=TRUE)})


#### Models ##########
modelstring.IC = "
model
{
  for (i in 1:(nind +nz)) {
  w[i] ~ dbern(psi)					#augmentation 
  x[i] ~ dunif(0,Bx)					#distance from line; Bx = max(distances) 
  
  z[i] ~ dgamma(a[i],c.beta[i])T(4,65)		#covariate                
  a[i] <- shape[clust[i]]		## a depends on which cluster the size comes from
  c.beta[i] <- betaClust[clust[i]]	#precision is allowed to vary between clusters
  clust[i] ~ dcat( pClust[1:Nclust] )	## clusters are defined as categories, assuming Nclust clusters
  
  sigma[i] <- exp(sigma.int+sigma.beta*z[i])  #sigma of detection
  logp[i] <- -((x[i]*x[i])/(2*sigma[i]*sigma[i]))	#This is the normal distribution with an adjustment for covs
  p[i] <- exp(logp[i])*xi[i]
  xi[i] <- ifelse(z[i] < b.point, m*z[i]+intercept, 1) #if less than b.point, probability is linear; if larger than b.point, perfect detection on line
  
  mu[i] <- w[i]*p[i] 					# probabilty of seeing it for all the ones we DID see (where w[i] = 1)
  
  y[i] ~ dbern(mu[i])         #found vs. missed
  
  o[i]~ dbin(o2[i], 1)        #occupancy
  logit(o2[i]) <- o.int + z.beta*z[i]						
  }
  
  
  p.online ~ dunif(pmin, pmax)		#estimated detection on line for 4 cm burrows
  b.point ~ dunif(10,25)
  m <- (1-p.online)/(b.point-4) 	#slope for detection on the line for smaller burrows		
  intercept <- p.online-(4*m)	## finding intercept via the detection of the 4 cm burrow 
  
  
  sigma.int~ dnorm(0,s.int.tau)T(0,)
  s.int.tau <- 1/(s.int.sd*s.int.sd)
  s.int.sd ~ dunif(.00001,5)	
  sigma.beta~ dnorm(0,s.beta.tau)T(0,)
  s.beta.tau <- 1/(s.beta.sd*s.beta.sd)
  s.beta.sd ~ dunif(.00001,5)
  o.int~ dnorm(0,o.int.tau)
  o.int.tau <- 1/(o.int.sd*o.int.sd)
  o.int.sd ~ dunif(.00001,5)		
  z.beta~ dnorm(0,z.beta.tau)
  z.beta.tau <- 1/(z.beta.sd*z.beta.sd)
  z.beta.sd ~ dunif(.00001,5)	
  
  for (clustIdx in 1: Nclust) {
  shape[clustIdx] ~ dunif(1,80)
  betaClust[clustIdx] ~ dunif(0.2,2)
  }
  
  pClust[1:Nclust] ~ ddirch(onesRepNclust) # probability of each cluster = the probability of that category in the ddirch distribution
  
  
  psi~ dunif(0,1)			#exists or not		
  
  
  Occ <- sum(o)/(nind+nz)
  N <- sum(w)	
  D <- N/(2*L*Bx)   #all indiviudals density
  Nt <- N*Occ
  Dt <- Nt/(2*L*Bx)	#occupied density
  
  juvi1 <- sum(z < 13)/(nind+nz)
  juvi2 <- (sum(z < 21)- sum(z < 13))/(nind+nz)
  juvi3 <- sum(z > 21)/(nind+nz)
  
  
}
"

modelstring.PC = "
model
{
  for (i in 1:(nind +nz)) {
  w[i] ~ dbern(psi)					#augmentation 
  x[i] ~ dunif(0,Bx)					#distance from line; Bx = max(distances) 
  
  z[i] ~ dgamma(a[i],c.beta[i])T(4,65)		#covariate                
  a[i] <- shape[clust[i]]		## a depends on which cluster the size comes from
  c.beta[i] <- betaClust[clust[i]]	#precision is allowed to vary between clusters
  clust[i] ~ dcat( pClust[1:Nclust] )	## clusters are defined as categories, assuming Nclust clusters
  
  sigma[i] <- exp(sigma.int+sigma.beta*z[i])  #sigma of detection
  logp[i] <- -((x[i]*x[i])/(2*sigma[i]*sigma[i]))	#This is the normal distribution with an adjustment for covs
  p[i] <- exp(logp[i])

  mu[i] <- w[i]*p[i] 					# probabilty of seeing it for all the ones we DID see (where w[i] = 1)
  
  y[i] ~ dbern(mu[i])         #found vs. missed
  
  o[i]~ dbin(o2[i], 1)        #occupancy
  logit(o2[i]) <- o.int + z.beta*z[i]						
  }
  
  
  p.online <-1		#estimated detection on line for 4 cm burrows
  b.point <- 4
  
  sigma.int~ dnorm(0,s.int.tau)T(0,)
  s.int.tau <- 1/(s.int.sd*s.int.sd)
  s.int.sd ~ dunif(.00001,5)	
  sigma.beta~ dnorm(0,s.beta.tau)T(0,)
  s.beta.tau <- 1/(s.beta.sd*s.beta.sd)
  s.beta.sd ~ dunif(.00001,5)
  o.int~ dnorm(0,o.int.tau)
  o.int.tau <- 1/(o.int.sd*o.int.sd)
  o.int.sd ~ dunif(.00001,5)		
  z.beta~ dnorm(0,z.beta.tau)
  z.beta.tau <- 1/(z.beta.sd*z.beta.sd)
  z.beta.sd ~ dunif(.00001,5)	
  
  for (clustIdx in 1: Nclust) {
  shape[clustIdx] ~ dunif(1,80)
  betaClust[clustIdx] ~ dunif(0.2,2)
  }
  
  pClust[1:Nclust] ~ ddirch(onesRepNclust) # probability of each cluster = the probability of that category in the ddirch distribution
  
  psi~ dunif(0,1)			#exists or not		
  
  
  Occ <- sum(o)/(nind+nz)
  N <- sum(w)	
  D <- N/(2*L*Bx)   #all indiviudals density
  Nt <- N*Occ
  Dt <- Nt/(2*L*Bx)	#occupied density
  
  juvi1 <- sum(z < 13)/(nind+nz)  
  juvi2 <- (sum(z < 21)- sum(z <13))/(nind+nz)
  juvi3 <- sum(z > 21)/(nind+nz)
}
"

modelstring.ICV = "
model
{
  for (i in 1:(nind +nz)) {
  w[i] ~ dbern(psi)					#augmentation 
  x[i] ~ dunif(0,Bx)					#distance from line; Bx = max(distances) 
  
  z[i] ~ dgamma(a[i],c.beta[i])T(4,65)		#covariate                
  a[i] <- shape[clust[i]]		## a depends on which cluster the size comes from
  c.beta[i] <- betaClust[clust[i]]	#precision is allowed to vary between clusters
  clust[i] ~ dcat( pClust[1:Nclust] )	## clusters are defined as categories, assuming Nclust clusters
  
  v[i] ~ dbeta(d,e)T(.05,.95)				## vegetation is unknown so beta prior
  sigma[i] <- exp(sigma.int+sigma.beta*z[i]+sigma.gamma*v[i])  #sigma of detection
  logp[i] <- -((x[i]*x[i])/(2*sigma[i]*sigma[i]))	#This is the normal distribution with an adjustment for covs
  p[i] <- exp(logp[i])*xi[i]
  xi[i] <- ifelse(z[i] < b.point, m*z[i]+intercept, 1) #if less than b.point, probability is linear; if larger than b.point, perfect detection on line
  
  mu[i] <- w[i]*p[i] 					# probabilty of seeing it for all the ones we DID see (where w[i] = 1)
  
  y[i] ~ dbern(mu[i])         #found vs. missed
  
  o[i]~ dbin(o2[i], 1)        #occupancy
  logit(o2[i]) <- o.int + z.beta*z[i]						
  }
  
  for (q in 1:q) {
    v2[q] ~ dbeta (d,e)T(.05,.95)		# all veg measurements at the site
  }
  
  p.online ~ dunif(pmin, pmax)		#estimated detection on line for 4 cm burrows
  b.point ~ dunif(10,25)
  m <- (1-p.online)/(b.point-4) 	#slope for detection on the line for smaller burrows		
  intercept <- p.online-(4*m)	## finding intercept via the detection of the 4 cm burrow 
  
  
  sigma.int~ dnorm(0,s.int.tau)T(0,)
  s.int.tau <- 1/(s.int.sd*s.int.sd)
  s.int.sd ~ dunif(.00001,5)	
  sigma.beta~ dnorm(0,s.beta.tau)T(0,)
  s.beta.tau <- 1/(s.beta.sd*s.beta.sd)
  s.beta.sd ~ dunif(.00001,5)
  o.int~ dnorm(0,o.int.tau)
  o.int.tau <- 1/(o.int.sd*o.int.sd)
  o.int.sd ~ dunif(.00001,5)		
  z.beta~ dnorm(0,z.beta.tau)
  z.beta.tau <- 1/(z.beta.sd*z.beta.sd)
  z.beta.sd ~ dunif(.00001,5)	
  sigma.gamma~ dnorm(0,s.gamma.tau)T(,0)
  s.gamma.tau <- 1/(s.gamma.sd*s.gamma.sd)
  s.gamma.sd ~ dunif(.00001,5)
  
  d~dunif(.1,40)
  e~dunif(.1,40)
  
  for (clustIdx in 1: Nclust) {
  shape[clustIdx] ~ dunif(1,80)
  betaClust[clustIdx] ~ dunif(0.2,2)
  }
  
  pClust[1:Nclust] ~ ddirch(onesRepNclust) # probability of each cluster = the probability of that category in the ddirch distribution
  
  
  psi~ dunif(0,1)			#exists or not		
  
  
  Occ <- sum(o)/(nind+nz)
  N <- sum(w)	
  D <- N/(2*L*Bx)   #all indiviudals density
  Nt <- N*Occ
  Dt <- Nt/(2*L*Bx)	#occupied density
  
  juvi1 <- sum(z < 13)/(nind+nz)  
  juvi2 <-  (sum(z < 21)- sum(z <13))/(nind+nz)
  juvi3 <- sum(z > 21)/(nind+nz)
}
"

modelstring.PCV = "
model
{
  for (i in 1:(nind +nz)) {
  w[i] ~ dbern(psi)					#augmentation 
  x[i] ~ dunif(0,Bx)					#distance from line; Bx = max(distances) 
  
  z[i] ~ dgamma(a[i],c.beta[i])T(4,65)		#covariate                
  a[i] <- shape[clust[i]]		## a depends on which cluster the size comes from
  c.beta[i] <- betaClust[clust[i]]	#precision is allowed to vary between clusters
  clust[i] ~ dcat( pClust[1:Nclust] )	## clusters are defined as categories, assuming Nclust clusters
  
  v[i] ~ dbeta(d,e)T(.05,.95)				## vegetation is unknown so beta prior
  sigma[i] <- exp(sigma.int+sigma.beta*z[i]+sigma.gamma*v[i])  #sigma of detection
  logp[i] <- -((x[i]*x[i])/(2*sigma[i]*sigma[i]))	#This is the normal distribution with an adjustment for covs
  p[i] <- exp(logp[i])
  
 
  mu[i] <- w[i]*p[i] 					# probabilty of seeing it for all the ones we DID see (where w[i] = 1)
  
  y[i] ~ dbern(mu[i])         #found vs. missed
  
  o[i]~ dbin(o2[i], 1)        #occupancy
  logit(o2[i]) <- o.int + z.beta*z[i]						
  }
  
  for (q in 1:q) {
  v2[q] ~ dbeta (d,e)T(.05,.95)		# all veg measurements at the site
  }
  
  p.online <- 1		#estimated detection on line for 4 cm burrows
  b.point <- 4
  
  sigma.int~ dnorm(0,s.int.tau)T(0,)
  s.int.tau <- 1/(s.int.sd*s.int.sd)
  s.int.sd ~ dunif(.00001,5)	
  sigma.beta~ dnorm(0,s.beta.tau)T(0,)
  s.beta.tau <- 1/(s.beta.sd*s.beta.sd)
  s.beta.sd ~ dunif(.00001,5)
  o.int~ dnorm(0,o.int.tau)
  o.int.tau <- 1/(o.int.sd*o.int.sd)
  o.int.sd ~ dunif(.00001,5)		
  z.beta~ dnorm(0,z.beta.tau)
  z.beta.tau <- 1/(z.beta.sd*z.beta.sd)
  z.beta.sd ~ dunif(.00001,5)	
  sigma.gamma~ dnorm(0,s.gamma.tau)T(,0)
  s.gamma.tau <- 1/(s.gamma.sd*s.gamma.sd)
  s.gamma.sd ~ dunif(.00001,5)
  
  d~dunif(.1,40)
  e~dunif(.1,40)
  
  for (clustIdx in 1: Nclust) {
  shape[clustIdx] ~ dunif(1,80)
  betaClust[clustIdx] ~ dunif(0.2,2)
  }
  
  pClust[1:Nclust] ~ ddirch(onesRepNclust) # probability of each cluster = the probability of that category in the ddirch distribution
  
  
  psi~ dunif(0,1)			#exists or not		
  
  
  Occ <- sum(o)/(nind+nz)
  N <- sum(w)	
  D <- N/(2*L*Bx)   #all indiviudals density
  Nt <- N*Occ
  Dt <- Nt/(2*L*Bx)	#occupied density
  
  juvi1 <- sum(z < 13)/(nind+nz)  
  juvi2 <-  (sum(z < 21)- sum(z <13))/(nind+nz)
  juvi3 <- sum(z > 21)/(nind+nz)
}
"

modelstring.IEV = "
  model
{
  ## Loop over burrows detected in LTDS (excluding ISU boxes) plus augmented burrows
  for (i in (nind.IS+1):(nind.IS+nind.LT + nz)) {
  w[i] ~ dbern(psi)					##augmentation 
  x[i] ~ dunif(0,Bx)					#distance from line for the missed ones; Bx = max(distances) 
  
  z[i] ~ dgamma(a[i],c.beta[i])T(4,65)		#covariate                
  a[i] <- shape[clust[i]]		## a depends on which cluster the size comes from
  c.beta[i] <- betaClust[clust[i]]	#precision is allowed to vary between clusters
  clust[i] ~ dcat( pClust[1:Nclust] )	## clusters are defined as categories, assuming Nclust clusters
  
  v[i] ~ dbeta(d,e)T(.05,.95)				## vegetation is unknown so beta prior
  sigma[i] <- exp(sigma.int+sigma.beta*z[i]+sigma.gamma*v[i])  #sigma of detection
  logp[i] <- -((x[i]*x[i])/(2*sigma[i]*sigma[i]))	#This is the normal distribution with an adjustment for covs
  p[i] <- exp(logp[i])*xi[i]
  xi[i] <- ifelse(z[i] < b.point, m*z[i]+intercept, 1) #if less than b.point, probability is linear; if larger than b.point, perfect detection on line
  
  mu[i] <- w[i]*p[i] 					## probabilty of seeing it for all the ones we DID see (where w[i] = 1)
  y[i] ~ dbern(mu[i])         ## found or missed
  
  o[i]~ dbin(o2[i], 1)        #occupancy
  logit(o2[i]) <- o.int + z.beta*z[i]	
  o.real[i] <- o[i]*w[i]
  }
  
  ## Loop over burrows detected in ISUs
  for (i in 1:nind.IS) {  
  w[i] ~ dbern(1)       ## forcing w[i] to take value 1
  x[i] ~ dunif(0,Bx)					#distance from line for the missed ones; Bx = max(distances) 
  
  z[i] ~ dgamma(a[i],c.beta[i])T(4,65)		                
  a[i] <- shape[clust[i]]		## a depends on which cluster the size comes from
  c.beta[i] <- betaClust[clust[i]]	#precision is allowed to vary between clusters
  clust[i] ~ dcat( pClust[1:Nclust] )	## clusters are defined as categories, assuming Nclust clusters
  
  v[i] ~ dbeta(d,e)T(.05,.95)				## vegetation is unknown so beta prior
  sigma[i] <- exp(sigma.int+sigma.beta*z[i]+sigma.gamma*v[i])	#log link for size
  logp[i] <- -((x[i]*x[i])/(2*sigma[i]*sigma[i]))	#This is the normal distribution with an adjustment for covs
  p[i] <- exp(logp[i])*xi[i]
  xi[i] <- ifelse(z[i] < b.point, m*z[i]+intercept, 1) #if less than b.point, probability is linear; if larger than b.point, perfect detection on line
  mu[i] <- p[i] 					## probabilty of seeing it for all the ones we DID see (where w[i] = 1)
  
  y[i] ~ dbern(mu[i])
  
  o[i]~ dbin(o2[i], 1)
  logit(o2[i]) <- o.int + z.beta*z[i]
  
  o.real[i] <- o[i]*w[i]
  }
  
  p.online ~ dunif(pmin, pmax)		#estimated detection on line for 4 cm burrows
  b.point ~ dunif(10,25)
  m <- (1-p.online)/(b.point-4) 	#slope for detection on the line for smaller burrows		
  intercept <- p.online-(4*m)	## finding intercept via the detection of the 4 cm burrow 
  
  for (q in 1:q) {
  v2[q] ~ dbeta (d,e)T(.05,.95)		# all veg measurements at the site
  }
  
  sigma.int~ dnorm(0,s.int.tau)T(,)
  s.int.tau <- 1/(s.int.sd*s.int.sd)
  s.int.sd ~ dunif(.00001,5)
  sigma.beta~ dnorm(0,s.beta.tau)T(0,)
  s.beta.tau <- 1/(s.beta.sd*s.beta.sd)
  s.beta.sd ~ dunif(.00001,5)
  o.int~ dnorm(0,o.int.tau)
  o.int.tau <- 1/(o.int.sd*o.int.sd)
  o.int.sd ~ dunif(.00001,5)		
  z.beta~ dnorm(0,z.beta.tau)
  z.beta.tau <- 1/(z.beta.sd*z.beta.sd)
  z.beta.sd ~ dunif(.00001,5)	
  sigma.gamma~ dnorm(0,s.gamma.tau)T(,0)
  s.gamma.tau <- 1/(s.gamma.sd*s.gamma.sd)
  s.gamma.sd ~ dunif(.00001,5)
  
  d~dunif(.1,40)
  e~dunif(.1,40)
  
  

  for (clustIdx in 1: Nclust) {
  shape[clustIdx] ~ dunif(1,80)
  betaClust[clustIdx] ~ dunif(0.2,2)
  }
  
  pClust[1:Nclust] ~ ddirch(onesRepNclust) ## probability of each cluster = the probability of that category in the ddirch distribution
  psi~ dunif(0,1)			#exists or not		
  
  N <- sum(w)      ## w vector has (nind.IS + nind.LT) ones, plus psi fraction of nz's
  Occ <- sum(o.real)/N
  D <- N/(2*L*Bx)
  Nt <- N*Occ
  Dt <- D*Occ		#tort density
  
  juvi1 <- sum(z < 13)/(nind.IS+nind.LT + nz)
  juvi2 <-  (sum(z < 21)- sum(z <13))/(nind.IS+nind.LT + nz)
  juvi3 <- sum(z > 22)/(nind.IS+nind.LT + nz)
  
}
"

modelstring.IE = "
  model
{
  ## Loop over burrows detected in LTDS (excluding ISU boxes) plus augmented burrows
  for (i in (nind.IS+1):(nind.IS+nind.LT + nz)) {
  w[i] ~ dbern(psi)					##augmentation 
  x[i] ~ dunif(0,Bx)					#distance from line for the missed ones; Bx = max(distances) 
  
  z[i] ~ dgamma(a[i],c.beta[i])T(4,65)		#covariate                
  a[i] <- shape[clust[i]]		## a depends on which cluster the size comes from
  c.beta[i] <- betaClust[clust[i]]	#precision is allowed to vary between clusters
  clust[i] ~ dcat( pClust[1:Nclust] )	## clusters are defined as categories, assuming Nclust clusters
  
  sigma[i] <- exp(sigma.int+sigma.beta*z[i])  #sigma of detection
  logp[i] <- -((x[i]*x[i])/(2*sigma[i]*sigma[i]))	#This is the normal distribution with an adjustment for covs
  p[i] <- exp(logp[i])*xi[i]
  xi[i] <- ifelse(z[i] < b.point, m*z[i]+intercept, 1) #if less than b.point, probability is linear; if larger than b.point, perfect detection on line
  
  mu[i] <- w[i]*p[i] 					## probabilty of seeing it for all the ones we DID see (where w[i] = 1)
  y[i] ~ dbern(mu[i])         ## found or missed
  
  o[i]~ dbin(o2[i], 1)        #occupancy
  logit(o2[i]) <- o.int + z.beta*z[i]	
  o.real[i] <- o[i]*w[i]
  }
  
  ## Loop over burrows detected in ISUs
  for (i in 1:nind.IS) {  
  w[i] ~ dbern(1)       ## forcing w[i] to take value 1
  x[i] ~ dunif(0,Bx)					#distance from line for the missed ones; Bx = max(distances) 
  
  z[i] ~ dgamma(a[i],c.beta[i])T(4,65)		                
  a[i] <- shape[clust[i]]		## a depends on which cluster the size comes from
  c.beta[i] <- betaClust[clust[i]]	#precision is allowed to vary between clusters
  clust[i] ~ dcat( pClust[1:Nclust] )	## clusters are defined as categories, assuming Nclust clusters
  
  sigma[i] <- exp(sigma.int+sigma.beta*z[i])	#log link for size
  logp[i] <- -((x[i]*x[i])/(2*sigma[i]*sigma[i]))	#This is the normal distribution with an adjustment for covs
  p[i] <- exp(logp[i])*xi[i]
  xi[i] <- ifelse(z[i] < b.point, m*z[i]+intercept, 1) #if less than b.point, probability is linear; if larger than b.point, perfect detection on line
  mu[i] <- p[i] 					## probabilty of seeing it for all the ones we DID see (where w[i] = 1)
  
  y[i] ~ dbern(mu[i])
  
  o[i]~ dbin(o2[i], 1)
  logit(o2[i]) <- o.int + z.beta*z[i]
  
  o.real[i] <- o[i]*w[i]
  }
  
  p.online ~ dunif(pmin, pmax)		#estimated detection on line for 4 cm burrows
  b.point ~ dunif(10,25)
  m <- (1-p.online)/(b.point-4) 	#slope for detection on the line for smaller burrows		
  intercept <- p.online-(4*m)	## finding intercept via the detection of the 4 cm burrow 
  
  
  sigma.int~ dnorm(0,s.int.tau)T(,)
  s.int.tau <- 1/(s.int.sd*s.int.sd)
  s.int.sd ~ dunif(.00001,5)
  sigma.beta~ dnorm(0,s.beta.tau)T(0,)
  s.beta.tau <- 1/(s.beta.sd*s.beta.sd)
  s.beta.sd ~ dunif(.00001,5)
  o.int~ dnorm(0,o.int.tau)
  o.int.tau <- 1/(o.int.sd*o.int.sd)
  o.int.sd ~ dunif(.00001,5)		
  z.beta~ dnorm(0,z.beta.tau)
  z.beta.tau <- 1/(z.beta.sd*z.beta.sd)
  z.beta.sd ~ dunif(.00001,5)	
  
  
  
  for (clustIdx in 1: Nclust) {
  shape[clustIdx] ~ dunif(1,80)
  betaClust[clustIdx] ~ dunif(0.2,2)
  }
  
  pClust[1:Nclust] ~ ddirch(onesRepNclust) ## probability of each cluster = the probability of that category in the ddirch distribution
  psi~ dunif(0,1)			#exists or not		
  
  N <- sum(w)      ## w vector has (nind.IS + nind.LT) ones, plus psi fraction of nz's
  Occ <- sum(o.real)/N
  D <- N/(2*L*Bx)
  Nt <- N*Occ
  Dt <- D*Occ		#tort density
  
  juvi1 <- sum(z < 13)/(nind.IS+nind.LT + nz)
  juvi2 <-  (sum(z < 21)- sum(z <13))/(nind.IS+nind.LT + nz)
  juvi3 <- sum(z > 22)/(nind.IS+nind.LT + nz)
  
}
"
modelstring.PEV = "
  model
{
  ## Loop over burrows detected in LTDS (excluding ISU boxes) plus augmented burrows
  for (i in (nind.IS+1):(nind.IS+nind.LT + nz)) {
  w[i] ~ dbern(psi)					##augmentation 
  x[i] ~ dunif(0,Bx)					#distance from line for the missed ones; Bx = max(distances) 
  
  z[i] ~ dgamma(a[i],c.beta[i])T(4,65)		#covariate                
  a[i] <- shape[clust[i]]		## a depends on which cluster the size comes from
  c.beta[i] <- betaClust[clust[i]]	#precision is allowed to vary between clusters
  clust[i] ~ dcat( pClust[1:Nclust] )	## clusters are defined as categories, assuming Nclust clusters
  
  v[i] ~ dbeta(d,e)T(.05,.95)				## vegetation is unknown so beta prior
  sigma[i] <- exp(sigma.int+sigma.beta*z[i]+sigma.gamma*v[i])  #sigma of detection
  logp[i] <- -((x[i]*x[i])/(2*sigma[i]*sigma[i]))	#This is the normal distribution with an adjustment for covs
  p[i] <- exp(logp[i])
  
  mu[i] <- w[i]*p[i] 					## probabilty of seeing it for all the ones we DID see (where w[i] = 1)
  y[i] ~ dbern(mu[i])         ## found or missed
  
  o[i]~ dbin(o2[i], 1)        #occupancy
  logit(o2[i]) <- o.int + z.beta*z[i]	
  o.real[i] <- o[i]*w[i]
  }
  
  ## Loop over burrows detected in ISUs
  for (i in 1:nind.IS) {  
  w[i] ~ dbern(1)       ## forcing w[i] to take value 1
  x[i] ~ dunif(0,Bx)					#distance from line for the missed ones; Bx = max(distances) 
  
  z[i] ~ dgamma(a[i],c.beta[i])T(4,65)		                
  a[i] <- shape[clust[i]]		## a depends on which cluster the size comes from
  c.beta[i] <- betaClust[clust[i]]	#precision is allowed to vary between clusters
  clust[i] ~ dcat( pClust[1:Nclust] )	## clusters are defined as categories, assuming Nclust clusters
  
  v[i] ~ dbeta(d,e)T(.05,.95)				## vegetation is unknown so beta prior
  sigma[i] <- exp(sigma.int+sigma.beta*z[i]+sigma.gamma*v[i])	#log link for size
  logp[i] <- -((x[i]*x[i])/(2*sigma[i]*sigma[i]))	#This is the normal distribution with an adjustment for covs
  p[i] <- exp(logp[i])
  mu[i] <- p[i] 					## probabilty of seeing it for all the ones we DID see (where w[i] = 1)
  
  y[i] ~ dbern(mu[i])
  
  o[i]~ dbin(o2[i], 1)
  logit(o2[i]) <- o.int + z.beta*z[i]
  
  o.real[i] <- o[i]*w[i]
  }
  
  p.online <- 1		#estimated detection on line for 4 cm burrows
  b.point <- 4

  for (q in 1:q) {
  v2[q] ~ dbeta (d,e)T(.05,.95)		# all veg measurements at the site
  }
  
  sigma.int~ dnorm(0,s.int.tau)T(,)
  s.int.tau <- 1/(s.int.sd*s.int.sd)
  s.int.sd ~ dunif(.00001,5)
  sigma.beta~ dnorm(0,s.beta.tau)T(0,)
  s.beta.tau <- 1/(s.beta.sd*s.beta.sd)
  s.beta.sd ~ dunif(.00001,5)
  o.int~ dnorm(0,o.int.tau)
  o.int.tau <- 1/(o.int.sd*o.int.sd)
  o.int.sd ~ dunif(.00001,5)		
  z.beta~ dnorm(0,z.beta.tau)
  z.beta.tau <- 1/(z.beta.sd*z.beta.sd)
  z.beta.sd ~ dunif(.00001,5)	
  sigma.gamma~ dnorm(0,s.gamma.tau)T(,0)
  s.gamma.tau <- 1/(s.gamma.sd*s.gamma.sd)
  s.gamma.sd ~ dunif(.00001,5)
  
  d~dunif(.1,40)
  e~dunif(.1,40)
  
  
  
  for (clustIdx in 1: Nclust) {
  shape[clustIdx] ~ dunif(1,80)
  betaClust[clustIdx] ~ dunif(0.2,2)
  }
  
  pClust[1:Nclust] ~ ddirch(onesRepNclust) ## probability of each cluster = the probability of that category in the ddirch distribution
  psi~ dunif(0,1)			#exists or not		
  
  N <- sum(w)      ## w vector has (nind.IS + nind.LT) ones, plus psi fraction of nz's
  Occ <- sum(o.real)/N
  D <- N/(2*L*Bx)
  Nt <- N*Occ
  Dt <- D*Occ		#tort density
  
  juvi1 <- sum(z < 13)/(nind.IS+nind.LT + nz)
  juvi2 <- (sum(z < 21)- sum(z <13))/(nind.IS+nind.LT + nz)
  juvi3 <- sum(z > 22)/(nind.IS+nind.LT + nz)
  
}
"

modelstring.PE = "
  model
{
  ## Loop over burrows detected in LTDS (excluding ISU boxes) plus augmented burrows
  for (i in (nind.IS+1):(nind.IS+nind.LT + nz)) {
  w[i] ~ dbern(psi)					##augmentation 
  x[i] ~ dunif(0,Bx)					#distance from line for the missed ones; Bx = max(distances) 
  
  z[i] ~ dgamma(a[i],c.beta[i])T(4,65)		#covariate                
  a[i] <- shape[clust[i]]		## a depends on which cluster the size comes from
  c.beta[i] <- betaClust[clust[i]]	#precision is allowed to vary between clusters
  clust[i] ~ dcat( pClust[1:Nclust] )	## clusters are defined as categories, assuming Nclust clusters
  
  sigma[i] <- exp(sigma.int+sigma.beta*z[i])  #sigma of detection
  logp[i] <- -((x[i]*x[i])/(2*sigma[i]*sigma[i]))	#This is the normal distribution with an adjustment for covs
  p[i] <- exp(logp[i])
  
  mu[i] <- w[i]*p[i] 					## probabilty of seeing it for all the ones we DID see (where w[i] = 1)
  y[i] ~ dbern(mu[i])         ## found or missed
  
  o[i]~ dbin(o2[i], 1)        #occupancy
  logit(o2[i]) <- o.int + z.beta*z[i]	
  o.real[i] <- o[i]*w[i]
  }
  
  ## Loop over burrows detected in ISUs
  for (i in 1:nind.IS) {  
  w[i] ~ dbern(1)       ## forcing w[i] to take value 1
  x[i] ~ dunif(0,Bx)					#distance from line for the missed ones; Bx = max(distances) 
  
  z[i] ~ dgamma(a[i],c.beta[i])T(4,65)		                
  a[i] <- shape[clust[i]]		## a depends on which cluster the size comes from
  c.beta[i] <- betaClust[clust[i]]	#precision is allowed to vary between clusters
  clust[i] ~ dcat( pClust[1:Nclust] )	## clusters are defined as categories, assuming Nclust clusters
  
  sigma[i] <- exp(sigma.int+sigma.beta*z[i])	#log link for size
  logp[i] <- -((x[i]*x[i])/(2*sigma[i]*sigma[i]))	#This is the normal distribution with an adjustment for covs
  p[i] <- exp(logp[i])
  mu[i] <- p[i] 					## probabilty of seeing it for all the ones we DID see (where w[i] = 1)
  
  y[i] ~ dbern(mu[i])
  
  o[i]~ dbin(o2[i], 1)
  logit(o2[i]) <- o.int + z.beta*z[i]
  
  o.real[i] <- o[i]*w[i]
  }
  
  p.online <- 1		#estimated detection on line for 4 cm burrows
  b.point <- 4
  
  sigma.int~ dnorm(0,s.int.tau)T(,)
  s.int.tau <- 1/(s.int.sd*s.int.sd)
  s.int.sd ~ dunif(.00001,5)
  sigma.beta~ dnorm(0,s.beta.tau)T(0,)
  s.beta.tau <- 1/(s.beta.sd*s.beta.sd)
  s.beta.sd ~ dunif(.00001,5)
  o.int~ dnorm(0,o.int.tau)
  o.int.tau <- 1/(o.int.sd*o.int.sd)
  o.int.sd ~ dunif(.00001,5)		
  z.beta~ dnorm(0,z.beta.tau)
  z.beta.tau <- 1/(z.beta.sd*z.beta.sd)
  z.beta.sd ~ dunif(.00001,5)	
  
  for (clustIdx in 1: Nclust) {
  shape[clustIdx] ~ dunif(1,80)
  betaClust[clustIdx] ~ dunif(0.2,2)
  }
  
  pClust[1:Nclust] ~ ddirch(onesRepNclust) ## probability of each cluster = the probability of that category in the ddirch distribution
  psi~ dunif(0,1)			#exists or not		
  
  N <- sum(w)      ## w vector has (nind.IS + nind.LT) ones, plus psi fraction of nz's
  Occ <- sum(o.real)/N
  D <- N/(2*L*Bx)
  Nt <- N*Occ
  Dt <- D*Occ		#tort density
  
  juvi1 <- sum(z < 13)/(nind.IS+nind.LT + nz)
  juvi2 <- (sum(z < 21)- sum(z <13))/(nind.IS+nind.LT + nz)
  juvi3 <- sum(z > 22)/(nind.IS+nind.LT + nz)
  
}
"

############################# Create data, initials, etc.########################################
#source('RJMCMC.R')
Burrs <- data.frame(dist = Burrows$Dist, size = Burrows$Diameter,
                     Occ = Burrows$Occ, Found = Burrows$Found)
ifelse((chosenmodel == "111" | chosenmodel == "211" | chosenmodel == "121" |chosenmodel == "221"),
 Burrs$Veg <- Burrows$Veg,
 Burrs$Veg <- rep(.1,nrow(Burrs)))
Burrs <- subset(Burrs, Burrs$size >0)

Run.me <- function(q1){
  
  starttime <- Sys.time()
  Transects <- Trans
  nz <- nrow(Burrs)*3
  
  if(Freq == "Enhanced LTDS") {Found1 <- subset(Burrs, Burrs$Found != 1)  ## ISUS
  Found2 <- subset(Burrs, Burrs$Found == 1)}
  
  if((chosenmodel == "111" | chosenmodel == "211" | chosenmodel == "121" |chosenmodel == "221")){
    v2 <- Veg$Veg
    v2 <- ifelse(v2 < 0.05, .05, v2)
    v2 <- ifelse(v2 > 0.95, .95, v2)
    q = length(v2)}
  
  if(Freq == "Conventional LTDS"){    #Conventional LTDS
    x <- c(Burrs$dist, rep(NA,nz))
    y <- c(rep(1,nrow(Burrs)), rep(0,nz))
    z <- c(Burrs$size, rep(NA,nz))
    z <- ifelse(z > 65, 65, z)
    z <- ifelse(z < 4, 4, z)
    o <- c(Burrs$Occ, rep (NA, nz))
    L <- sum(Transects$Length)*10^-4
    v <- c(Burrs$Veg, rep(NA,nz))
    v <- ifelse(v < .05, .05 , ifelse(v > .95, .95, v))
    
    nind = nrow(Burrs)
    Bx = max(Burrs$dist, na.rm = TRUE)
    clust <- rep(NA,(nind +nz)) 	# no idea which cluster anything is in, so unknown	
  }
  
  if(Freq == "Enhanced LTDS"){  #ELTDS
    x <- c(Found1$dist, Found2$dist, rep(NA, nz))
    y <- c(Found1$Found, Found2$Found,rep(0, nz))
    z <- c(Found1$size, Found2$size, rep(NA, nz))
    z <- ifelse(z > 65, 65, z)
    z <- ifelse(z < 4, 4, z)
    o <- c(Found1$Occ,Found2$Occ, rep(NA, nz))
    v <- c(Found1$Veg, Found2$Veg, rep(NA,nz))
    v <- ifelse(v < .05, .05 , ifelse(v > .95, .95, v))
    nind.LT <- nrow(Found2)   ## calling these counts nind.LT and nind.IS to keep things straight
    nind.IS <- nrow(Found1) ## Things in ISU boxes
    Bx <- max(Found2$dist)
    clust = rep(NA, (nind.IS+nind.LT+nz)) 
    L <- sum(Transects$Length)*10^-4
  }
  
  # w.clust <- c(0.50, 0.1, 0.33, .03, .04)
  # shape.clust <- c(10, 20, 30, 40, 66)
  # rate.clust <- c(1,.7,1, 2, 2)
  # y.clust <- Burrs$size
  # Z <- do.call(cbind, lapply(1:5, function(j)
  #  w.clust[j]*dgamma(y.clust, shape.clust[j], rate = rate.clust[j])))
  # Z <- apply(Z, 1, function(x) which(x==max(x))[1])
  # print("Estimating clusters for size mixture model")
  # res <- mixgam.rjmcmc(y = y.clust, nsweep = 80000, kmax = 10, k = 5, w = w.clust,
  #                     shape = shape.clust, rate = rate.clust, Z,
  #                     delta=1.75, xi=NULL, kappa=NULL, alpha=NULL,
  #                     beta=NULL, g=NULL, h=NULL, verbose=TRUE)
  # 
  # ksave <- res$k.save
  # groups <- round(table(ksave[-(1:40000)])/40000,4)
  # Nclust <<- as.numeric(names(groups)[which(groups == max(groups))])+1
  # if(length(Nclust) != 1){Nclust <<- sample(Nclust,1)}
  Nclust <- 4

  clust[which.min(z)]=1 # smallest value assigned to cluster 1; generally represents juvis
  clust[which.max(z)]=Nclust # highest value assigned to largest cluster value; generally represents large adults


  ### Run the model you want to run - initial values and inputs will change, as will desired outputs

  if(selected.model == "modelstring.PCV" | selected.model == "modelstring.ICV") {#conventional with veg
    jd.test <<- list(nind= nind, nz = nz, L = L, Bx = Bx, x=x,y=y,z=z, Nclust = Nclust, q=q, o=o, v=v,
                  clust= clust, onesRepNclust = rep(1,Nclust),  pmin = p.min, pmax = p.max)
    ji.test1 <<-list(pClust = runif(Nclust, 0.1, 0.9), shape = runif(Nclust, 10, 80),
                     betaClust = runif(Nclust, 0.2, 1.9), sigma.int = .01, sigma.beta = .05,
                     w = c(rep(1, nind), rbinom(nz, 1,.5)))
    ji.test2 <<-list(pClust = runif(Nclust, 0.1, 0.9), shape = runif(Nclust, 10, 80),
                     betaClust = runif(Nclust, 0.2, 1.9), sigma.int = .01, sigma.beta = .05,
                     w = c(rep(1, nind), rbinom(nz, 1,.5)))
    ji.test3 <<-list(pClust = runif(Nclust, 0.1, 0.9), shape = runif(Nclust, 10, 80),
                     betaClust = runif(Nclust, 0.2, 1.9), sigma.int = .02, sigma.beta = .05,
                     w = c(rep(1, nind), rbinom(nz, 1,.5)))}
  if(selected.model == "modelstring.PC" | selected.model == "modelstring.IC") {#conventional no veg
    jd.test <<- list(nind= nind, nz = nz, L = L, Bx = Bx, x=x,y=y,z=z, Nclust = Nclust, o=o,
                    clust= clust, onesRepNclust = rep(1,Nclust),  pmin = p.min, pmax = p.max)
    ji.test1 <<-list(pClust = runif(Nclust, 0.1, 0.9), shape = runif(Nclust, 10, 80),
                     betaClust = runif(Nclust, 0.2, 1.9), sigma.int = .01, sigma.beta = .05,
                     w = c(rep(1, nind), rbinom(nz, 1,.5)))
    ji.test2 <<-list(pClust = runif(Nclust, 0.1, 0.9), shape = runif(Nclust, 10, 80),
                     betaClust = runif(Nclust, 0.2, 1.9), sigma.int = .01, sigma.beta = .05,
                     w = c(rep(1, nind), rbinom(nz, 1,.5)))
    ji.test3 <<-list(pClust = runif(Nclust, 0.1, 0.9), shape = runif(Nclust, 10, 80),
                     betaClust = runif(Nclust, 0.2, 1.9), sigma.int = .02, sigma.beta = .05,
                     w = c(rep(1, nind), rbinom(nz, 1,.5)))}
  if(selected.model == "modelstring.PEV" | selected.model == "modelstring.IEV") {#enhanced with veg
    jd.test <<- list(nind.LT= nind.LT, nz=nz, L = L, x=x,y=y,z=z, o=o, Bx = Bx, Nclust = Nclust,
                 clust= clust, onesRepNclust = rep(1,Nclust), nind.IS= nind.IS, q= q, v= v,  pmin = p.min, pmax = p.max)
    ji.test1 <<-list(pClust = runif(Nclust, 0.1, 0.9), shape = runif(Nclust, 10, 80),
                     betaClust = runif(Nclust, 0.2, 1.9), sigma.int = .01, sigma.beta = .05,
                     w = c(rep(1, nind.IS+nind.LT), rbinom(nz, 1,.5)))
    ji.test2 <<-list(pClust = runif(Nclust, 0.1, 0.9), shape = runif(Nclust, 10, 80),
                     betaClust = runif(Nclust, 0.2, 1.9), sigma.int = .01, sigma.beta = .05,
                     w = c(rep(1, nind.IS+nind.LT), rbinom(nz, 1,.5)))
    ji.test3 <<-list(pClust = runif(Nclust, 0.1, 0.9), shape = runif(Nclust, 10, 80),
                     betaClust = runif(Nclust, 0.2, 1.9), sigma.int = .02, sigma.beta = .05,
                     w = c(rep(1, nind.IS+nind.LT), rbinom(nz, 1,.5)))}
  if(selected.model == "modelstring.PE" | selected.model == "modelstring.IE") {#enhanced no veg
    jd.test <<- list(nind.LT= nind.LT, nz = nz, L = L, Bx = Bx, x=x,y=y,z=z, Nclust = Nclust, o=o,
                 clust= clust, onesRepNclust = rep(1,Nclust), nind.IS= nind.IS, pmin = p.min, pmax = p.max)
    ji.test1 <<-list(pClust = runif(Nclust, 0.1, 0.9), shape = runif(Nclust, 10, 80),
                     betaClust = runif(Nclust, 0.2, 1.9), sigma.int = .01, sigma.beta = .05,
                     w = c(rep(1, nind.IS+nind.LT), rbinom(nz, 1,.5)))
    ji.test2 <<-list(pClust = runif(Nclust, 0.1, 0.9), shape = runif(Nclust, 10, 80),
                     betaClust = runif(Nclust, 0.2, 1.9), sigma.int = .01, sigma.beta = .05,
                     w = c(rep(1, nind.IS+nind.LT), rbinom(nz, 1,.5)))
    ji.test3 <<-list(pClust = runif(Nclust, 0.1, 0.9), shape = runif(Nclust, 10, 80),
                     betaClust = runif(Nclust, 0.2, 1.9), sigma.int = .02, sigma.beta = .05,
                     w = c(rep(1, nind.IS+nind.LT), rbinom(nz, 1,.5)))
    }


  

  jp.test <<- c("D","Dt", "N","Nt", "betaClust", "pClust","shape",
                "sigma.int", "sigma.beta", "p.online", "b.point", "Occ",
                "o.int", "z.beta", "juvi1", "juvi2", "juvi3")

  if(selected.model == "modelstring.PEV" | selected.model == "modelstring.IEV" |
     selected.model == "modelstring.PCV" | selected.model == "modelstring.ICV" ){
    jp.test <<- c("D","Dt", "N","Nt", "betaClust", "pClust","shape",
               "sigma.int", "sigma.beta", "p.online", "b.point", "Occ",
               "o.int", "z.beta", "sigma.gamma", "d", "e", "juvi1", "juvi2", "juvi3")}

  print("Time to actually run the model!")
  Foo <<- suppressWarnings(autorun.jags(eval(parse(text = selected.model)), data = jd.test, inits = list(ji.test1, ji.test2, ji.test3), monitor= jp.test,
                                       startburnin = 4000, startsample = 5000, adapt = 2000, method = 'parallel', n.chain = 3, max.time= "1s",
                                       psrf.target = 1.1,silent.jags = FALSE, summarise = TRUE))
  print("Results coming soon!")
  Covmod <- summary(Foo)
  Ext <<- Foo

  results = list()
  results$Site <- paste(selected.model)
  results$LTDS <- Covmod
  results$Time <- paste("Above Results Created:", Sys.time())

  lapply(results, function(x) write.table(data.frame(x), 'results_Bias_Adjusted_LTDS.csv'  ,
                                          append= T, sep=',' ,row.names = T, col.names = NA ))

  return(Foo)
}


Extend.me <- function(q2){
  try2 <- extract.runjags(Ext, 'end.state')
  Foo <<- suppressWarnings(autorun.jags(eval(parse(text = selected.model)), data = jd.test, inits = try2, monitor= jp.test,
                                       startburnin = 10000, startsample = 8000, adapt = 200, method = 'parallel', n.chain = 3, max.time=maxusrtime,
                                       psrf.target = 1.1,silent.jags = FALSE, summarise = TRUE))
  print("New and improved results coming soon!")
  Covmod <- summary(Foo)
  Ext <<- Foo

  results = list()
  results$Site <- paste(selected.model)
  results$LTDS <- Covmod
  results$Time <- paste("Above Results Created:", Sys.time())


  lapply(results, function(x) write.table(data.frame(x), 'results_Bias_Adjusted_LTDS.csv'  ,
                                          append= T, sep=',' ,row.names = T, col.names = NA ))
}

#Plot.me <- function(q3){
#     Covmod <- summary(Foo)
#     ## occupancy plot
#     Occ.overall <- Covmod["Occ", 4]
#     z.beta <- Covmod["z.beta", 4]
#     o.int <- Covmod["o.int",4]
#     xs <- seq(4, 65, by = .01)
#     par(mar = c(5, 5, 5, 1))
#     plot(xs, exp(o.int + z.beta*xs)/(1+exp(o.int + z.beta*xs)), type = "l", 
#          main = paste("Probability of Occupancy by Size", "\nAverage Occupancy = ", floor(Occ.overall*100), "%"), 
#                       xlab = "Burrow Width in cm", cex.main = 1.5,
#          ylab = "Probability of Occupancy", lwd = 2, ylim = c(0,1), col = "blue",
#          cex.lab = 1.5, cex.axis = 1.5
#          )
#     points(Burrows$size, Burrows$Occ, pch = 19)
#     text(60, .1, paste("y = ", round(o.int, 2), round(z.beta, 2), "*x", sep = ""), cex = 1.25)
#   
#     ##detection curve plot 
#     sigma.int <- Covmod["sigma.int", 4]
#     sigma.beta <- Covmod["sigma.beta", 4]
#     p.online <- Covmod["p.online", 4]
#     b.point <- Covmod["b.point", 4]
#     m <- (1-p.online)/(b.point-4)	
#     intercept <- p.online-(4*m)
#     sigma <- exp(sigma.int+sigma.beta*c(4, 15, 25, 35))  #sigma of detection
#     xi <- ifelse(c(4,15,25,35) < b.point, m*c(4,15,25,35)+intercept, 1)
#     xs <- seq(0,60, by = .01)
#     p.4 <- exp(-((xs*xs)/(2*sigma[1]*sigma[1])))*xi[1]	#This is the normal distribution with an adjustment for covs
#     p.15 <- exp(-((xs*xs)/(2*sigma[2]*sigma[2])))*xi[2]
#     p.25 <- exp(-((xs*xs)/(2*sigma[3]*sigma[3])))*xi[3]
#     p.35 <- exp(-((xs*xs)/(2*sigma[4]*sigma[4])))*xi[4]
#     plot(xs, p.35, xlim = c(0, 40), ylim = c(0,1), main = "Detection Curves By Size", 
#          xlab = "Distance From Transect in Meters", ylab = "Detection Probability", 
#          type = "l", col = "black", lwd =2,  cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
#     lines(xs, p.25, lwd = 2, col = "blue")
#     lines(xs, p.15, lwd = 2, col = "purple")
#     lines(xs, p.4, lwd = 2, col = "darkgreen")
#     legend("topright", c("4 cm", "15 cm", "25 cm", "35 cm"), col = c("darkgreen", "purple", "blue", "black"), cex = 1.5, lty = c(1,1,1, 1), lwd = 2)
#     
#     ### Size curves 
#     xs <<- seq(0,65, by = .01)
#           ifelse(nrow(Covmod) == 32, #6 clusters
#                          {
#                            k1 <- dgamma(xs,Covmod[17,4], rate = Covmod[5,4])*Covmod[11,4]
#                            k2 <- dgamma(xs,Covmod[18,4], rate = Covmod[6,4])*Covmod[12,4]
#                            k3 <- dgamma(xs,Covmod[19,4], rate = Covmod[7,4])*Covmod[13,4]
#                            k4 <- dgamma(xs,Covmod[20,4], rate = Covmod[8,4])*Covmod[14,4]
#                            k5 <- dgamma(xs,Covmod[21,4], rate = Covmod[9,4])*Covmod[15,4]
#                            k6 <- dgamma(xs,Covmod[22,4], rate = Covmod[10,4])*Covmod[16,4]
#                            Kall <<- k1+k2+k3+k4+k5+k6
#                          },
#                          ifelse(nrow(Covmod) == 29, #5 clusters
#                                 {
#                                   k1 <- dgamma(xs,Covmod[15,4], rate = Covmod[5,4])*Covmod[10,4]
#                                   k2 <- dgamma(xs,Covmod[16,4], rate = Covmod[6,4])*Covmod[11,4]
#                                   k3 <- dgamma(xs,Covmod[17,4], rate = Covmod[7,4])*Covmod[12,4]
#                                   k4 <- dgamma(xs,Covmod[18,4], rate = Covmod[8,4])*Covmod[13,4]
#                                   k5 <- dgamma(xs,Covmod[19,4], rate = Covmod[9,4])*Covmod[14,4]
#                                   Kall <<- k1+k2+k3+k4+k5
#                                 },
#                                 ifelse(nrow(Covmod) == 26, #4 clusters
#                                        {
#                                          k1 <- dgamma(xs,Covmod[13,4], rate = Covmod[5,4])*Covmod[9,4]
#                                          k2 <- dgamma(xs,Covmod[14,4], rate = Covmod[6,4])*Covmod[10,4]
#                                          k3 <- dgamma(xs,Covmod[15,4], rate = Covmod[7,4])*Covmod[11,4]
#                                          k4 <- dgamma(xs,Covmod[16,4], rate = Covmod[8,4])*Covmod[12,4]
#                                          Kall <<- k1+k2+k3+k4
#                                        },
#                                        { k1 <- dgamma(xs,Covmod[11,4], rate = Covmod[5,4])*Covmod[8,4]
#                                        k2 <- dgamma(xs,Covmod[12,4], rate = Covmod[6,4])*Covmod[9,4]
#                                        k3 <- dgamma(xs,Covmod[13,4], rate = Covmod[7,4])*Covmod[10,4]
#                                        Kall <<- k1+k2+k3
#                                        }
#                                 )))
#     Occ <<- exp(o.int + z.beta*xs)/(1+exp(o.int + z.beta*xs)) 
#     par(mfrow = c(1,1), mar = c(5,2,3,2))
#     plot(xs, Kall*Occ, type = "l", main = "Estimated Size Curve for Population",
#          col = "black", lwd =2,  cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, xlab = "Burrow Width in cm",
#          ylab = "", yaxt = "n", ylim = c(-.0002, max(Kall*Occ)+.0005))
#     percents <- Covmod[c("juvi1", "juvi2", "juvi3"), 4]
#     legend("topright", c(paste(round(percents[1],2)*100, "%", "  ", "< 13 cm"), 
#                          paste(round(percents[2],2)*100, "%", "  ", "13-21 cm"), 
#                          paste(round(percents[3],2)*100, "%", "  ", "> 21 cm")), 
#                          text.col = c("darkgreen", "purple", "blue"), cex = 1.35)
#     arrows(4, -.0002, 12.9, -.0002, length = 0.075, col = "darkgreen", code = 3, angle = 30, lwd = 2)
#     arrows(13, -.0002, 21.9, -.0002, length = 0.075, col = "purple", code = 3, angle = 30, lwd = 2)
#     arrows(22, -.0002, 65, -.0002, length = 0.075, col = "blue", code = 3, angle = 30, lwd = 2)
#     #juvi 1, 2, and 3 are sizes < 13, 13-21, >21
# }

Plot.me.shiny <- function(q3){
  Covmod <<- summary(Foo)
  ## occupancy plot
  Occ.overall <- Covmod["Occ", 4]
  z.beta <- Covmod["z.beta", 4]
  o.int <- Covmod["o.int",4]
  xs <<- seq(0,65, by = .01)
  ##detection curve plot 
  sigma.int <- Covmod["sigma.int", 4]
  sigma.beta <- Covmod["sigma.beta", 4]
  p.online <- Covmod["p.online", 4]
  b.point <- Covmod["b.point", 4]
  m <- (1-p.online)/(b.point-4)	
  intercept <- p.online-(4*m)
  sigma <- exp(sigma.int+sigma.beta*c(4, 15, 25, 35))  #sigma of detection
  xi <- ifelse(c(4,15,25,35) < b.point, m*c(4,15,25,35)+intercept, 1)
  p.4 <- exp(-((xs*xs)/(2*sigma[1]*sigma[1])))*xi[1]	#This is the normal distribution with an adjustment for covs
  p.15 <- exp(-((xs*xs)/(2*sigma[2]*sigma[2])))*xi[2]
  p.25 <- exp(-((xs*xs)/(2*sigma[3]*sigma[3])))*xi[3]
  p.35 <- exp(-((xs*xs)/(2*sigma[4]*sigma[4])))*xi[4]
  ### Size curves 
  ifelse(nrow(Covmod) == 32, #6 clusters
         {
           k1 <- dgamma(xs,Covmod[17,4], rate = Covmod[5,4])*Covmod[11,4]
           k2 <- dgamma(xs,Covmod[18,4], rate = Covmod[6,4])*Covmod[12,4]
           k3 <- dgamma(xs,Covmod[19,4], rate = Covmod[7,4])*Covmod[13,4]
           k4 <- dgamma(xs,Covmod[20,4], rate = Covmod[8,4])*Covmod[14,4]
           k5 <- dgamma(xs,Covmod[21,4], rate = Covmod[9,4])*Covmod[15,4]
           k6 <- dgamma(xs,Covmod[22,4], rate = Covmod[10,4])*Covmod[16,4]
           Kall <<- k1+k2+k3+k4+k5+k6
         },
         ifelse(nrow(Covmod) == 29, #5 clusters
                {
                  k1 <- dgamma(xs,Covmod[15,4], rate = Covmod[5,4])*Covmod[10,4]
                  k2 <- dgamma(xs,Covmod[16,4], rate = Covmod[6,4])*Covmod[11,4]
                  k3 <- dgamma(xs,Covmod[17,4], rate = Covmod[7,4])*Covmod[12,4]
                  k4 <- dgamma(xs,Covmod[18,4], rate = Covmod[8,4])*Covmod[13,4]
                  k5 <- dgamma(xs,Covmod[19,4], rate = Covmod[9,4])*Covmod[14,4]
                  Kall <<- k1+k2+k3+k4+k5
                },
                ifelse(nrow(Covmod) == 26, #4 clusters
                       {
                         k1 <- dgamma(xs,Covmod[13,4], rate = Covmod[5,4])*Covmod[9,4]
                         k2 <- dgamma(xs,Covmod[14,4], rate = Covmod[6,4])*Covmod[10,4]
                         k3 <- dgamma(xs,Covmod[15,4], rate = Covmod[7,4])*Covmod[11,4]
                         k4 <- dgamma(xs,Covmod[16,4], rate = Covmod[8,4])*Covmod[12,4]
                         Kall <<- k1+k2+k3+k4
                       },
                       { k1 <- dgamma(xs,Covmod[11,4], rate = Covmod[5,4])*Covmod[8,4]
                       k2 <- dgamma(xs,Covmod[12,4], rate = Covmod[6,4])*Covmod[9,4]
                       k3 <- dgamma(xs,Covmod[13,4], rate = Covmod[7,4])*Covmod[10,4]
                       Kall <<- k1+k2+k3
                       }
                )))
  Occ <<- exp(o.int + z.beta*xs)/(1+exp(o.int + z.beta*xs)) 
  plot(xs, Kall*Occ, type = "l", main = "Estimated Size Curve for Population",
       col = "black", lwd =2,  cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, xlab = "Burrow Width in cm",
       ylab = "", yaxt = "n", ylim = c(-.0002, max(Kall*Occ)+.0008))
  percents <- Covmod[c("juvi1", "juvi2", "juvi3"), 4]
  legend("topright", c(paste(round(percents[1],2)*100, "%", "  ", "< 13 cm"), 
                       paste(round(percents[2],2)*100, "%", "  ", "13-21 cm"), 
                       paste(round(percents[3],2)*100, "%", "  ", "> 21 cm")), 
         text.col = c("darkgreen", "purple", "blue"), cex = 1.35)
  arrows(4, -.0002, 12.9, -.0002, length = 0.075, col = "darkgreen", code = 3, angle = 30, lwd = 2)
  arrows(13, -.0002, 21.9, -.0002, length = 0.075, col = "purple", code = 3, angle = 30, lwd = 2)
  arrows(22, -.0002, 65, -.0002, length = 0.075, col = "blue", code = 3, angle = 30, lwd = 2)
  #juvi 1, 2, and 3 are sizes < 13, 13-21, >21
}
  
