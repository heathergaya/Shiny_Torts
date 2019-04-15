####################################################################
### Reversible Jump MCMC Code for Mixtures of Gammas             ###
### Code modified from the uvnm.rjmcmc function in package MiscF ###
### Last Edited 5/20/2018                                        ###
####################################################################


#decide whether to split-combine and birth-death
decideScBd <- function(k, kmax){
  if(k<= 2)
    action <- 1
  else{
    if(k==kmax)
      action <- 0
    else
      action <- ifelse(runif(1) <0.5, 1, 0)
  }
}


#get the probabilities for z in the split and combine cases

#getP <- function(yjstar, muj1, muj2, sigma2j1, sigma2j2, wj1, wj2){
#p <- cbind(dnorm(yjstar, muj1, sqrt(sigma2j1)),
#           dnorm(yjstar, muj2, sqrt(sigma2j2)))

getP <- function(yjstar, shapej1, shapej2, ratej1, ratej2, wj1, wj2){  
p <- cbind(dgamma(yjstar, shape = shapej1, rate = ratej1), 
           dgamma(yjstar, shape = shapej2, rate = ratej2))
  p[,1] <- p[,1] * wj1
  p[,2] <- p[,2] * wj2                 
  s <- rowSums(p)
  p[which(s == 0 | !is.finite(s)),] <- 1
  p / rowSums(p)
  
}

#get acceptance probabilities for split-combine
#getAsc <- function(yjstar, j1.pos, j2.pos, k,
#                   wj1, wj2, wjstar,
#                   muj1, muj2, mujstar,
#                   sigma2j1, sigma2j2, sigma2jstar,
#                   delta, xi, kappa, alpha, beta,
#                   p, b, d, u1, u2, u3
#){
getAsc <- function(yjstar, j1.pos, j2.pos, k,
                     wj1, wj2, wjstar,
                     shapej1, shapej2, shapejstar,
                     ratej1, ratej2, ratejstar,
                     delta, xi, kappa, alpha, beta,
                     p, b, d, u1, u2, u3){
  
  #browser()
  #calculate acceptance probability
  ## take log first then transform to exponential b/c there are
  ## multiplication of larger number of densities and so on.
  
#likelihood ratio

#log.lr <- sum(dnorm(yjstar[j1.pos], mean=muj1, sd=sqrt(sigma2j1), log=TRUE)) +
#    sum(dnorm(yjstar[j2.pos], mean=muj2, sd=sqrt(sigma2j2), log=TRUE)) -
#    sum(dnorm(yjstar, mean=mujstar, sd=sqrt(sigma2jstar), log=TRUE))
log.lr <- sum(dgamma(yjstar[j1.pos], shape = shapej1, rate = ratej1, log=TRUE)) +
          sum(dgamma(yjstar[j2.pos], shape = shapej2, rate = ratej2, log=TRUE)) -
          sum(dgamma(yjstar, shape=shapejstar, rate = ratejstar, log=TRUE))
l1 <- length(j1.pos)
l2 <- length(j2.pos)
  
#term1-3: ratio between two states
  #term1 includes various ratios
  log.term1 <- log.lr + #likelihood ratio
    log(1) + #ratio of p(k+1)/p(k)
    (delta-1+l1)*log(wj1) + (delta-1+l2)*log(wj2) - (delta-1+l1+l2)*log(wjstar) - #ratio of z
    lbeta(delta, k*delta) #ratio of w
  #  (lgamma(k*delta+n) + lgamma(delta+l1) +lgamma(delta+l2) - lgamma((k+1)*delta + n) - lgamma(delta+l1+l2))   
  
  #term2 is the ratio of mu
  ## mu = shape/rate = alpha/beta for gamma
  ## sd = sqrt(alpha/beta^2) for gamma

  log.term2 <- log(k+1) +
    0.5*log(kappa/(2*pi)) -
    0.5*kappa*((shapej1-xi)*2 + (shapej2-xi)*2 - (shapejstar-xi)*2) 
  
  
    #term3 is the ratio of sigma
  #log.term3 <- alpha*log(beta) - lgamma(alpha) +
  #  (-alpha-1)*(log(sigma2j1) + log(sigma2j2) - log(sigma2jstar)) -
  #  beta*(1/sigma2j1 + 1/sigma2j2 - 1/sigma2jstar)
  
  log.term3 <- alpha*log(beta) - lgamma(alpha) +
    (-alpha-1)*(log(ratej1) + log(ratej2) - log(ratejstar)) -
    beta*(ratej1+ratej2- ratejstar)
  
  #term4 is the transform ratio between two states
  log.Palloc <- sum(log(p[,1][j1.pos])) + sum(log(p[,2][j2.pos]))
  log.term4 <- log(d[k+1]) - log(b[k]) - log.Palloc -
    dbeta(u1, 2, 2, log=TRUE) - dbeta(u2, 2, 2, log=TRUE) - dbeta(u3, 1, 1, log=TRUE)
  
  #term5 is the Jacobian
  #log.term5 <- log(wjstar) +log(abs(muj1-muj2)) + log(sigma2j1) + log(sigma2j2) -
  #  log(u2) -log(1-u2^2) - log(u3) - log(1-u3) - log(sigma2jstar)
  ##################################################################################
  ### TERM 5 is our current problem areas ###############################
  
  log.term5 <- log(wjstar) +log(abs(shapej1-shapej2)) + log(ratej1) + log(ratej2) -
    log(u2) -log(1-u2^2) - log(u3) - log(1-u3) - log(shapejstar)
  
  terms <- exp(log.term1 + log.term2 + log.term3 + log.term4 + log.term5)
  ifelse(terms == "Inf", 10^100, terms)
}


split <- function(y, k, w, shape, rate, Z, b, d, delta, xi, kappa, alpha, beta){
  #choose a component to split
  jstar <- sample(1:k, 1)
  
  #generate intermidiate parameters
  u1 <- rbeta(1, 2, 2)
  u2 <- rbeta(1, 2, 2)
  u3 <- rbeta(1, 1, 1)
  
  #generate two new ws
  wjstar <- w[jstar]
  wj1 <- wjstar * u1
  wj2 <- wjstar * (1-u1)
  
  
  #generate two new mus
  shapejstar <- shape[jstar]
  ratejstar <- rate[jstar]

  shapej1 <- shapejstar - u2*shapejstar*ratejstar*sqrt(wj2/wj1)
  shapej1 <- ifelse(shapej1 < .5, runif(1,.5,110), shapej1)
  shapej2 <- shapejstar + u2*shapejstar*ratejstar*sqrt(wj1/wj2)
  shapej2 <- ifelse(shapej2 > 110, runif(1,.5,110), shapej2)
  
    
  #generate two new rates
    ratej2 <-  ratejstar+sqrt(u3 * (1-u2^2) * sqrt(ratejstar)* wjstar / wj1)
    ratej1 <- sqrt((1-u3) * (1-u2^2) * sqrt(ratejstar) * wjstar / wj2)
    ratej2 <- ifelse(ratej2 >4, 4, ratej2)
    ratej1 <- ifelse(ratej1 <.5, .5, ratej1)
  
  #check order of means
    newshape.part <- c((shape[jstar-1])/rate[jstar-1], shapej1/ratej1, shapej2/ratej2, shape[jstar+1]/rate[jstar+1])
    if(all(order(newshape.part) == 1:length(newshape.part))) {
    #allocate z_i=jstar
    jstar.pos <- which(Z==jstar)
    if(length(jstar.pos) > 0){
      yjstar <- y[jstar.pos]
      p <- getP(yjstar, shapej1, shapej2, ratej1, ratej2, wj1, wj2)
      zj12 <- rMultinom(p)
      j1.pos <- which(zj12==1)
      j2.pos <- which(zj12==2)

      A <- getAsc(yjstar, j1.pos, j2.pos, k,
                  wj1, wj2, wjstar,
                  shapej1, shapej2, shapejstar,
                  ratej1, ratej2, ratejstar,
                  delta, xi, kappa, alpha, beta,
                  p, b, d, u1, u2, u3)
      
      if(runif(1) < min(1, A)){
        indicator <- rep(0,k)
        indicator[jstar] <- 1
        indicator[which((1:k) > jstar)] <- 2
        ind0 <- which(indicator==0)
        ind2 <- which(indicator==2)
        #generate new w vector
        w <- c(w[ind0], wj1, wj2, w[ind2])
        #generate new mu vector
        shape <- c(shape[ind0],shapej1, shapej2, shape[ind2])
        #generate new sigma2 vector
        rate <- c(rate[ind0], ratej1, ratej2, rate[ind2])
        #gernerate new Z matrix
        larger <- which(Z > jstar)
        Z[larger] <- Z[larger] + 1
        Z[jstar.pos[j1.pos]] <- jstar
        Z[jstar.pos[j2.pos]] <- jstar + 1
        
        }
      }
    }
  list(w=w, shape=shape, rate=rate, Z=Z)
}

combine <- function(y, k, w, shape, rate, Z, b, d, delta, xi, kappa, alpha, beta){
  #choose a pair of components to combine
  j1 <- sample(1:(k-1), 1)
  j2 <- j1 + 1
  
  #generate new parameters
  wj1 <- w[j1]
  wj2 <- w[j2]
  shapej1 <- shape[j1]
  shapej2 <- shape[j2]
  ratej1 <- rate[j1]
  ratej2 <- rate[j2]
  wjstar <- w[j1] + w[j2]
  shapejstar <- (wj1*shapej1 + wj2*shapej2) / wjstar
  shapejstar <- ifelse(shapejstar > 110 | shapejstar < .5, runif(1,shapej1,shapej2), shapejstar)
  ratejstar <- (wj1*ratej1 + wj2*ratej2)/wjstar
  ratejstar <- ifelse(ratejstar >4 |ratejstar<.5, runif(1,min(c(ratej1, ratej2)),max(c(ratej1,ratej2))), ratejstar)
  
  #calculate acceptance probability
  #likelihood ratio
  jstar.pos <- which(Z==j1 | Z==j2)
  j1.pos <- match(which(Z==j1), jstar.pos)
  j2.pos <- match(which(Z==j2), jstar.pos)
  yjstar <- y[jstar.pos]
  p <- getP(yjstar, shapej1, shapej2, ratej1, ratej2, wj1, wj2)
  
  #generate intermidiate parameters
  #u1 u2 u3 are derived from split move (equations below eqn (10) on page 739
  #u1 <- rbeta(1, 2, 2)
  #u2 <- rbeta(1, 2, 2)
  #u3 <- rbeta(1, 1, 1)
  u1 <- wj1 / wjstar
  
  #sqrt(sigma2jstar) <- (shapejstar/(ratejstar)^2)
  #muj2 <- shapej2/ratej2
  u2 = abs(shapej2-shapej1)* (sqrt(wj1*wj2)/ratejstar)/ wjstar/100
  u2 = max(u2,.02)
  u2 = min(u2,.98)
  u3 = rbeta(1,1,1)
  #u2 <- u3 <- .4
  A <- getAsc(yjstar, j1.pos, j2.pos, k-1,
              wj1, wj2, wjstar,
              shapej1, shapej2, shapejstar,
              ratej1, ratej2, ratejstar,
              delta, xi, kappa, alpha, beta,
              p, b, d, u1, u2, u3)
  A.c <- ifelse(is.na(1/A), 0, 1/A)
  
  if(runif(1) < min(1, A.c)){
    #browser()
    indicator <- rep(0,k)
    indicator[c(j1,j2)] <- 1
    indicator[which((1:k) > j2)] <- 2
    ind0 <- which(indicator==0)
    ind2 <- which(indicator==2)
    #generate new w vector
    w <- c(w[ind0], wjstar, w[ind2])
    #generate new shape vector
    shape <- c(shape[ind0], shapejstar, shape[ind2])
    #generate new sigma2 vector
    rate <- c(rate[ind0], ratejstar, rate[ind2])
    #gernerate new Z matrix
    Z[which(Z==j2)] <- j1
    large <- which(Z > j2)
    Z[large] <- Z[large] - 1
  }
  list(w=w, shape=shape, rate=rate, Z=Z)
}

getAbd <- function(n, k, k0, delta, wjstar, b, d){
  log.term1 <- log(1) - lbeta(k*delta, delta) + (delta-1)*log(wjstar) +
    (n+k*delta-k)*log(1-wjstar) + log(k+1)
  #note that there is an error in the original paper:
  # (1-wjstar)^(k-1) instead of (1-wjstar)^k
  log.term2 <- log(d[k+1]) - log(k0+1) - log(b[k]) - dbeta(wjstar,1,k, log=TRUE) + (k-1)*log(1-wjstar)
  exp(log.term1 + log.term2)
}


birth <- function(n, k, w, shape, rate, Z, delta, xi, kappa, alpha, beta, b, d){
  
  wjstar <- rbeta(1, 1, k)
  k0 <- sum(unlist(lapply(1:k, function(i) sum(Z==i)))==0)
  A <- getAbd(n, k, k0, delta, wjstar, b, d)
  if(runif(1) < min(1, A)){
    
    shapejstar <- rgamma(1, shape=xi, rate = kappa)
    shapejstar <- ifelse(shapejstar >110, 110, ifelse(shapejstar < .5, .5, shapejstar))
    #shapejstar <- runif(1, 1, 110)
    ratejstar <- rgamma(1, shape=alpha, rate=beta)
    ratejstar <- ifelse(ratejstar >4, 4, ifelse(ratejstar<.5, .5, ratejstar))
    #ratejstar <- runif(1,.1,4)
    
    w <- w*(1-wjstar)
    jstar.pos <- which(shape > shapejstar)[1]
    if(is.na(jstar.pos)){
      w <- c(w, wjstar)
      shape <- c(shape, shapejstar)
      rate <- c(rate, ratejstar)
    }
    else{
      indicator <- rep(0,k)
      indicator[which((1:k) >= jstar.pos)] <- 1
      ind0 <- which(indicator==0)
      ind1 <- which(indicator==1)
      w <- c(w[ind0], wjstar, w[ind1])
      shape <- c(shape[ind0], shapejstar, shape[ind1])
      rate <- c(rate[ind0], ratejstar, rate[ind1])
      larger <- which(Z >= jstar.pos)
      Z[larger] <- Z[larger] + 1
    }
  }
  list(w=w, shape=shape, rate=rate, Z=Z)
}


death <- function(n, k, w, shape, rate, Z, delta, xi, kappa, alpha, beta, b, d){
  d.candidate <- which(unlist(lapply(1:k, function(i) sum(Z==i)==0)))
  if(length(d.candidate) > 0){
    d.pos <- sample(1:length(d.candidate),1)
    d.pos <- d.candidate[d.pos]
    wjstar <- w[d.pos]
    k0 <- length(d.candidate) - 1
    A <- getAbd(n, k-1, k0, delta, wjstar, b, d)
    if(runif(1) < min(1, 1/A)){
      w <- w[-d.pos]
      w <- w / sum(w)
      shape <- shape[-d.pos]
      rate <- rate[-d.pos]
      larger <- which(Z > d.pos)
      Z[larger] <- Z[larger] - 1
    }
  }
  list(w=w, shape=shape, rate=rate, Z=Z)
}

#uvnm.rjmcmc <- function(y, nsweep, kmax, k, w, mu, sigma2, Z,
#                        delta=1, xi=NULL, kappa=NULL, alpha=2,
#                        beta=NULL, g=0.2, h=NULL, verbose=TRUE){

mixgam.rjmcmc <- function(y, nsweep, kmax, k, w, shape, rate, Z,
                          delta=1.5, xi=NULL, kappa=NULL, alpha=2,
                          beta=NULL, g=0.2, h=NULL, verbose=TRUE){
  #Error checking
  if(nsweep <= 0)
    stop("The number of sweeps has to be positive.")
  if(kmax < 0 || k < 0)
    stop("The number of components have to be positive.")
  if(kmax < k)
    stop("The maximum number of components allowed is larger than
         the intitial value of components.")
  if(length(w) != k){
    w <- rep(w, len=k)
    warning("The length of 'w' was not equal to k and
            was forced to be k by being cut off or recycled.")
  }
  if(any(w <=0))
    stop()
  if(sum(w) != 1){
    
  }
  
  #########################
  if(length(shape) != k){
    shape <- rep(shape, len=k)
    warning("The length of 'shape' was not equal to k and
            was forced to be k by being cut off or recycled.")
  }
  if(length(rate) != k){
    rate <- rep(rate, len=k)
    warning("The length of 'rate' was not equal to k and
            was forced to be k by being cut off or recycled.")
  }
  if(any(rate <= 0))
    stop
  if(any(shape <= 0))
    stop
  if(length(Z) != length(y)){
    Z <- rep(Z, len=length(y))
    warning("The length of 'Z' was not equal to the length of 'y' and
            was forced to be equal by being cut off or recycled.")
  }
  
  R <- diff(range(y))
  if(is.null(xi))
    xi <- 5
  if(is.null(kappa))
    kappa <- .1
  if(is.null(delta))
    delta <- runif(1, 2, 4)
  if(is.null(alpha))
    alpha <- 2
  if(is.null(g))
    g <- 6
  if(is.null(h))
    h <- .03
  if(is.null(beta))
    beta <- .1
  
  n=length(y)
  
  #split probabilities
  b <- rep(0.5, kmax)
  b[kmax] <- 0
  #combine probabilities
  d <- rep(0.5, kmax)
  d[c(1:2)] <- 0  #can't combine when there are only 1 or 2 clusts
  
  k.save <- rep(0, nsweep)
  w.save <- shape.save <- rate.save <-  vector("list", nsweep)
  Z.save <- matrix(0, nrow=n, ncol=nsweep)
  
  for(i in 1:nsweep){
    k <- length(shape)
    
    #update w|...
    Z.expand <- do.call(cbind, lapply(1:k, function(i) ifelse(Z==i, 1, 0)))
    Nj <- colSums(Z.expand)
    w <- rDirichlet(1, delta + Nj)
  
    #update shape|...
    sumYj <- colSums(y * Z.expand)
    shape.temp <- shape +Nj/2
    shape.temp <- ifelse(shape.temp > 110, 110, ifelse(shape.temp < .5, .5, shape.temp))
    #rate.temp <- rate + sqrt(1/(kappa + sumYj/2))
    #rate.temp <- ifelse(rate.temp >4, 4, ifelse(rate.temp <.1, .1, rate.temp))
    shape.new <- rgamma(k, shape=shape.temp, rate=rate)
    shape.new <- ifelse(shape.new > 110 | shape.new < .5,runif(1,.5,110), shape.new)

    

    #update rate|...
    Diff2j <- outer(y, shape/rate, `-`) ^ 2
    sumDiff2j <- colSums(Diff2j * Z.expand)
    rate.new <- rgamma(k, shape=alpha + Nj/2, rate=beta + sqrt(sumDiff2j/2))
    rate.new <- ifelse(rate.new >4, 4, ifelse(rate.new <.5, .5, rate.new))
    if(all(order(shape.new/rate.new) == 1:k)){
      shape <- shape.new
      rate <- rate.new
    }
    #update Z
    #p <- exp(- Diff2j %*% diag(1/(2*sigma2), nrow=k, ncol=k))
    #p <- p %*% diag(w/sqrt(sigma2), nrow=k, ncol=k)
    #p <- do.call(cbind, lapply(1:k, function(i)
    #  w[i] * dnorm(y, mu[i], sqrt(sigma2[i]))))
    p <- do.call(cbind, lapply(1:k, function(i)
      w[i]*dgamma(y, shape = shape[i], rate = rate[i])))
    s <- rowSums(p)
    p[which(s == 0 | !is.finite(s)),] <- 1
    p <- p / rowSums(p)
    if(ncol(p) > 1){
      Z <- rMultinom(p)
    }
    else{
      Z <- rep(1, n)
    }
    
    #update beta
    beta <- rgamma(1, shape=g + k*alpha, rate=h + rate)

    #combine or split
    action <- decideScBd(k, kmax)
    if(action==1){
      split.results <- split(y, k, w, shape, rate, Z, b, d, delta, xi, kappa, alpha, beta)
      w <- split.results$w
      shape <- split.results$shape
      rate <- split.results$rate
      Z <- split.results$Z
      k <- length(w)
    }
    else{
      combine.results <- combine(y, k, w, shape, rate, Z, b, d, delta, xi, kappa, alpha, beta)
      w <- combine.results$w
      shape <- combine.results$shape
      rate <- combine.results$rate
      Z <- combine.results$Z
      k <- length(w)
    }
    
    #birth-death
    action <- decideScBd(k, kmax)
    if(action==1){
      birth.results <- birth(n, k, w, shape, rate, Z, delta,
                             xi, kappa, alpha, beta, b, d)
      w <- birth.results$w
      shape <- birth.results$shape
      rate <- birth.results$rate
      Z <- birth.results$Z
      k <- length(w)
    }
    else{
      
      death.results <- death(n, k, w, shape, rate, Z, delta,
                             xi, kappa, alpha, beta, b, d)
      w <- death.results$w
      shape <- death.results$shape
      rate <- death.results$rate
      Z <- death.results$Z
      k <- length(w)
    }
    
    k.save[i] <- k
    w.save[[i]] <- w
    shape.save[[i]] <- shape
    rate.save[[i]]<- rate
    Z.save[,i] <- Z
    if (verbose && i %% 1000 == 0){
      cat(paste(i, " sweeps", " have finished.\n", sep=""))
    }
    
  }
  list(k.save=k.save, w.save=w.save, shape.save= shape.save,
       rate.save=rate.save, Z.save=Z.save)
}
###########################################################################
######################### Time to test out the code... ####################
############################################################################
library(mixAK)
library(miscF)

#y.clust <- c(rgamma(30, shape = 80, rate =.8),
#             rgamma(30, shape = 20, rate = 2),
#             rgamma(35, shape = 60, rate = 1.5), 
#             rgamma(30, shape = 100, rate = 1.5))
#hist(y.clust, breaks = seq(0,150, by = 2.5))
#k = 9
#w.clust <- rep(1,k)
#shape.clust <- c(5,15, 40, 60, 10, 20, 20, 80, 10)
#rate.clust <- rep(1,k)
##################
#Z <- do.call(cbind, lapply(1:k, function(j)
#  w.clust[j]*dgamma(y.clust, shape.clust[j], rate = rate.clust[j])))
#Z <- apply(Z, 1, function(x) which(x==max(x))[1])
#nsweep = 50000
#result <- mixgam.rjmcmc(y = y.clust, nsweep = nsweep, kmax = 10, k = k, w = w.clust, 
#                      shape = shape.clust, rate = rate.clust, Z,
#                      delta=1, xi=50, kappa=.5, alpha=1,
#                      beta=NULL, g=20, h=NULL, verbose=TRUE)
#
#ksave <- result$k.save
#groups <- round(table(ksave)/nsweep,4)
#groups <- round(table(ksave[-(1:10000)])/(nsweep-10000),4)
#groups
#Nclust <- as.numeric(names(groups)[which(groups == max(groups))])
#if(length(Nclust) != 1){Nclust <- max(Nclust)}
#Nclust
