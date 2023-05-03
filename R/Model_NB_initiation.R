#' Probability to acquire 2 mutations during exponential expansion (compare Haeno et al., 2007)
#'
#' @param lambda1 proliferation rate during expansion; set to 1
#' @param delta1 loss rate during expansion
#' @param N maximal precursor count
#' @param t time
#' @param lambda2 proliferation rate during contraction; set to 1
#' @param mu1 rate of acquiring 1st oncogenic event
#' @param mu2 rate of acquiring 2nd oncogenic event
#' @param p.surv survival probability of 2nd mutant
#' P.both.hits.during.expansion()
#' @export

P.both.hits.during.expansion <- function(lambda1, delta1, N, t, lambda2, mu1, mu2, p.surv){

  ## stem cell count at time t
  M <- exp((lambda1-delta1)*t)
  ## compute a2, the proliferation rate of the 2nd mutant (lambda1*s) assuming that selection acts on proliferation. Then p.surv = 1-delta1/(s*lambda1) --> s*lambda1 = delta1/(1-p.surv)
  a2 <- delta1/(1-p.surv)

  ## require a2 >= 1 (else it would be a selective disadvantage
  if(a2 < 1){
    a2 <- 1
  }

  beta <- (1-delta1/lambda1)*mu1/(1-delta1/lambda1)

  func <- function(a2, delta1, lambda1){
    integrand <- function(lambda1, delta1, a2, z){
      (1-delta1/a2)/(1-(delta1/a2)*z^((a2-delta1)/(lambda1-delta1)))
    }
    integrate(f = integrand, lower = 0, upper = 1, lambda1=lambda1, delta1=delta1, a2=a2)$value
  }

  integrand <- function(M, lambda1, delta1, beta, mu2, x){
    taux <- log(M)/(lambda1-delta1) - log(x)/(lambda1-delta1)
    y <- exp((lambda1-delta1)*taux)
    res <- exp(-beta*(x-1))*(1-exp(-beta))*(1-exp(-mu2*y*func(a2, delta1, lambda1)/(1-delta1/lambda1)))
    res
  }
  if(is.infinite(integrand(M, lambda1, delta1, beta, mu2, 1))){
    return(10^11)
  }
  res <- sum(sapply(seq(1, log10(M)), function(bins){
    integrate(integrand, lower = 10^(bins-1), upper = 10^bins, M=M, lambda1=lambda1, delta1=delta1, beta = beta, mu2=mu2)$value
  }))
  res
}



#' Probability to acquire the first hit during exponential expansion and the second hit during contraction
#'
#' @param lambda1 proliferation rate during expansion; set to 1
#' @param delta1 loss rate during expansion
#' @param N maximal precursor count
#' @param t time
#' @param lambda2 proliferation rate during contraction; set to 1
#' @param delta2 loss rate during contraction, set to arbitrary value if modeling homeostasis rather than contraction
#' @param p.surv survival probability of 2nd mutant
#' @param mu1 rate of acquiring 1st oncogenic event
#' @param mu2 rate of acquiring 2nd oncogenic event
#' @param r selective advantage associated with 1st oncogenic event
#' @param model inidicates whether expansion and homeostasis or expansion and contraction is modeled
#' P.2nd.driver.1st.hit.during.expansion()
#' @export


P.2nd.driver.1st.hit.during.expansion <- function(lambda1, delta1, N, t, lambda2, delta2=NA, p.surv, mu1, mu2, r, model="contraction"){

  ## time point at which the population peaks
  t.max <- log(N)/(lambda1-delta1)
  if(model=="contraction"){
    P <-   1 - exp(-mu1*mu2*p.surv*lambda1*lambda2/(lambda2-delta2*r)*N*t.max*(exp((lambda2-delta2*r)*t) - 1))
  }else if(model=="homeostasis"){
    P <- 1 - exp(-mu1*mu2*p.surv*lambda1/(1-r)*N*t.max*(exp(lambda2*(1-r)*t) - 1))
  }
  P
}


#' Probability to acquire both hits during contraction
#'
#' @param lambda1 proliferation rate during expansion; set to 1
#' @param delta1 loss rate during expansion
#' @param N maximal precursor count
#' @param t time
#' @param lambda2 proliferation rate during contraction; set to 1
#' @param delta2 loss rate during contraction, set to arbitrary value if modeling homeostasis rather than contraction
#' @param p.surv survival probability of 2nd mutant
#' @param mu1 rate of acquiring 1st oncogenic event
#' @param mu2 rate of acquiring 2nd oncogenic event
#' @param r selective advantage associated with 1st oncogenic event
#' @param model inidicates whether expansion and homeostasis or expansion and contraction is modeled
#' P.2nd.driver.1st.hit.during.2nd.phase()
#' @export

P.2nd.driver.1st.hit.during.2nd.phase <- function(lambda1, delta1, N, t, lambda2, delta2=NA, p.surv, mu1, mu2, r, model="contraction"){

  if(model=="contraction"){
    P <- 1-exp(-mu1*mu2*lambda2^2*p.surv*N/(delta2*(r-1))*(
      1/(lambda2-delta2)*(exp((lambda2-delta2)*t) - 1) - 1/(lambda2 - r*delta2)*(exp((lambda2 - r*delta2)*t)-1)
    ))
  }else if(model=="homeostasis"){
    P <- 1-exp(mu1*mu2*lambda2*p.surv*N/(1-r)*(1/(lambda2*(1-r))*(1-exp(lambda2*(1-r)*t))+t))
  }


  return(P)
}


#' Probability to acquire two hits at time t
#'
#' @param parms a list containing delta1 (loss rate during expansion), delta2 (loss rate during contraction, set to arbitrary value if modeling homeostasis), N (log10 of maximal precursor count), p.surv (survival probability associated with 2nd event), mu1 (rate of acquiring the first oncogenic event, log10), mu2 (rate of acquiring the second oncogenic event, log10), mu (neutral mutation rate), r (selective advantage associated with first oncogenic event)
#' @param t time
#' @param model inidicates whether expansion and homeostasis or expansion and contraction is modeled
#' P.2nd.driver()
#' @export


P.2nd.driver <- function(parms, t, model="contraction"){

  delta1 <- parms$delta1
  N <- 10^parms$N
  if(model=="homeostasis"){
    delta2 <- NA
  }else{
    delta2 <- parms$delta2
  }
  p.surv <- parms$psurv
  mu1 <- 10^parms$muD1
  mu2 <- 10^parms$muD2
  r <- parms$r
  ## set lambda1, lambda2 to 1
  lambda1 <- 1
  lambda2 <- 1

  ## time point at which the population peaks
  t.turning.point <- log(N)/(lambda1-delta1)

  ## return Inf if the maximum population size is not reached after 2000 steps
  if(log10(N)/(lambda1-delta1) > 2000){
    return(10^14)
  }

  if(t < t.turning.point){
    res <-  P.both.hits.during.expansion(lambda1, delta1, N, t, lambda2, mu1, mu2, p.surv)
  }else{
    res <- 1 -(1-P.both.hits.during.expansion(lambda1, delta1, N, t.turning.point, lambda2, mu1, mu2, p.surv))*
      (1- P.2nd.driver.1st.hit.during.expansion(lambda1, delta1, N, t - t.turning.point, lambda2, delta2, p.surv, mu1, mu2, r, model=model))*
      (1-P.2nd.driver.1st.hit.during.2nd.phase(lambda1, delta1, N, t - t.turning.point, lambda2, delta2, p.surv, mu1, mu2, r, model=model))
  }

  res
}


#' Probability to acquire the first hit at time t. Only written for the contraction case
#'
#' @param parms a list containing delta1 (loss rate during expansion), delta2 (loss rate during contraction, set to arbitrary value if modeling contraction), N (log10 of maximal precursor count), mu1 (rate of acquiring the first oncogenic event, log10), mu (neutral mutation rate), r (selective advantage associated with first oncogenic event)
#' @param t time
#' P.1st.driver()
#' @export


P.1st.driver <- function(parms, t){

  ## shut off the generation of first-mutants after 2000 mutations (approximately the maximum measured mutation count at MRCA)
  if(t > 2500/parms$mu){
    return(0)
  }

  delta1 <- parms$delta1
  N <- 10^parms$N
  delta2 <- parms$delta2
  mu1 <- 10^parms$muD1
  r <- parms$r
  lambda1 <- 1
  lambda2 <- 1


  ## in case that r*delta2 > lambda2: the first hit cannot expand during involution. Assuming that a sizeable tumor must make up
  ## at least 1% of the tissue, it must've started its growth before N reached 1% of its maximal size.
  if(lambda2 < delta2*r){
    t.turning.point <- log(N/100)/(lambda1-delta1)
  }else{
    t.turning.point <- log(N)/(lambda1-delta1)
  }

  ## return Inf if the population size is not met
  if(log10(N)/(lambda1-delta1) > 2000){
    return(10^14)
  }

  ## survival probability of the first mutant during contraction
  p.surv <- 1-delta2*r/lambda2
  if(p.surv < 0){p.surv <- 0}

  if(t < t.turning.point){
    ## survival probability of the first mutant during expansion (same as for normal cells)
    p.surv <- 1-delta1/lambda1
    res <- 1-exp(-mu1*p.surv/(lambda1-delta1)*(exp((lambda1-delta1)*t)-1))
  }else{
    ## survival probability of the first mutant during expansion (same as for normal cells)
    p.surv.1 <- 1-delta1/lambda1

    res <- 1 - exp(-mu1*p.surv.1/(lambda1-delta1)*(exp((lambda1-delta1)*t.turning.point)-1))*
       exp(-mu1*p.surv*N/(lambda2-delta2)*(exp((lambda2-delta2)*(t-t.turning.point)) -1))

  }

  res

}


#' Probability of the first event happening at time t1 given the second event happened at t2
#'
#' @param parms a list containing delta1 (loss rate during expansion), delta2 (loss rate during contraction), N (log10 of maximal precursor count), p.surv (survival probability associated with 2nd event), mu1 (rate of acquiring the first oncogenic event, log10), mu2 (rate of acquiring the second oncogenic event, log10), mu (neutral mutation rate), r (selective advantage associated with first oncogenic event)
#' @param t1 time point of first mutation
#' @param t2 time point of second mutation
#' @param model inidicates whether expansion and homeostasis or expansion and contraction is modeled
#' P.eca.at.t()
#' @export

P.eca.at.t <- function(parms, t1, t2, model="contraction"){

  delta1 <- parms$delta1
  N <- 10^parms$N
  if(model=="homeostasis"){
    delta2 <- NA
  }else{
    delta2 <- parms$delta2
  }
  p.surv <- parms$psurv
  mu1 <- 10^parms$muD1
  mu2 <- 10^parms$muD2
  r <- parms$r
  lambda1 <- 1
  lambda2 <- 1

  ## time point at which the probability peaks
  t.turning.point <- log(N)/(lambda1-delta1)

  if(model=="contraction"){
    ## for cases where t2 > t.turning.point: the number of expected type-1-cells generated during expansion and contraction
    total.count.of.first.hit.during.expansion <- mu1*N*exp((lambda2 - delta2*r)*(t2 - t.turning.point))*t.turning.point
    total.count.of.first.hit.during.decay <-mu1*N/(delta2*(r-1))*(exp((lambda2-delta2)*(t2-t.turning.point)) - exp((lambda2 - delta2*r)*(t2-t.turning.point)))

    if(total.count.of.first.hit.during.decay + total.count.of.first.hit.during.expansion ==0){
      return(0)
    }
    if(is.infinite(total.count.of.first.hit.during.decay + total.count.of.first.hit.during.expansion)){
      return(0)
    }

    if(t2 < t.turning.point){
      res <- t1/t2
    }else if(t1 <= t.turning.point){
      res <- mu1*N*exp((lambda2-delta2*r)*(t2-t.turning.point))*t1/
        (total.count.of.first.hit.during.expansion + total.count.of.first.hit.during.decay)
    }else if(t1 <= t2){
      res <- mu1*N*exp((lambda2-delta2*r)*(t2-t.turning.point))*t.turning.point/
        (total.count.of.first.hit.during.expansion + total.count.of.first.hit.during.decay) +
        mu1*N/(delta2*(r-1))*(-exp((lambda2 - delta2*r)*(t2-t.turning.point)) + exp(lambda2*(t2-t.turning.point) +
                                                                                      delta2*(t.turning.point + (r-1)*t1 - r*t2)))/
        (total.count.of.first.hit.during.expansion + total.count.of.first.hit.during.decay)

    }else{
      res <-0
    }
  }else if(model=="homeostasis"){
    ## Time point at which the population size peaks
    t.turning.point <- log(N)/(lambda1-delta1)

    total.count.of.first.hit.during.expansion <- mu1*N*exp(lambda2*(1-r)*(t2 - t.turning.point))*t.turning.point
    total.count.of.first.hit.during.homeostasis <-mu1*N/(r-1)*(1-exp(lambda2*(1-r)*(t2-t.turning.point)))

    if(total.count.of.first.hit.during.expansion + total.count.of.first.hit.during.homeostasis ==0){
      return(0)
    }
    if(is.infinite(total.count.of.first.hit.during.homeostasis + total.count.of.first.hit.during.expansion)){
      return(0)
    }

    if(t2 < t.turning.point){
      res <- t1/t2
    }else if(t1 <= t.turning.point){
      res <- mu1*N*exp(lambda2*(1-r)*(t2-t.turning.point))*t1/
        (total.count.of.first.hit.during.expansion + total.count.of.first.hit.during.homeostasis)
    }else if(t1 <= t2){
      res <- mu1*N*exp(lambda2*(1-r)*(t2-t.turning.point))*t.turning.point/
        (total.count.of.first.hit.during.expansion + total.count.of.first.hit.during.homeostasis) +
        mu1*N/(r-1)*(exp(lambda2*(1-r)*(t2-t1))-exp(lambda2*(1-r)*(t2-t.turning.point)))/
        (total.count.of.first.hit.during.expansion + total.count.of.first.hit.during.homeostasis)

    }else{
      res <-0
    }
  }

  res
}

#' Sample the time point of tumor initiation
#'
#' @param parms a list containing delta1 (loss rate during expansion), delta2 (loss rate during contraction), N (log10 of maximal precursor count), p.surv (survival probability associated with 2nd event), mu1 (rate of acquiring the first oncogenic event, log10), mu2 (rate of acquiring the second oncogenic event, log10), mu (neutral mutation rate), r (selective advantage associated with first oncogenic event)
#' @param model inidicates whether expansion and homeostasis or expansion and contraction is modeled
#' Sample.timepoint()
#' @export

Sample.timepoint <- function(parms, model="contraction"){

  # First, sample a uniformly distributed number; assume total incidence of 10^-5
  x <- runif(1, 0, 1)*10^-5
  ## Then, choose the associated probability
  tmp <- function(x, parms, t, model){
    P.2nd.driver(parms, t, model=model) - x
  }

  if((P.2nd.driver(parms, 10000, model=model) >x & P.2nd.driver(parms, 0, model=model)  > x)|
     (P.2nd.driver(parms, 10000, model=model) <x & P.2nd.driver(parms, 0, model=model)  < x)){
    return(100000000)
  }

  t <- uniroot(tmp, c(0, 10000), x=x, parms=parms, model=model)$root
  t
}


#' Sample the number of mutations present in the founer cell at time t
#'
#' @param parms a list containing delta1 (loss rate during expansion), delta2 (loss rate during contraction), N (log10 of maximal precursor count), p.surv (survival probability associated with 2nd event), mu1 (rate of acquiring the first oncogenic event, log10), mu2 (rate of acquiring the second oncogenic event, log10), mu (neutral mutation rate), r (selective advantage associated with first oncogenic event)
#' @param t time
#' Sample.mutations()
#' @export

Sample.mutations <- function(parms, t){

  n.mutations <- rpois(n=1, lambda=parms$mu*t)
  n.mutations
}


#' Sample the time point of the ECA
#'
#' @param parms a list containing delta1 (loss rate during expansion), delta2 (loss rate during contraction), N (log10 of maximal precursor count), p.surv (survival probability associated with 2nd event), mu1 (rate of acquiring the first oncogenic event, log10), mu2 (rate of acquiring the second oncogenic event, log10), mu (neutral mutation rate), r (selective advantage associated with first oncogenic event)
#' @param t2 time point of the MRCA
#' @param model inidicates whether expansion and homeostasis or expansion and contraction is modeled
#' Sample.timepoint.eca()
#' @export

Sample.timepoint.eca <- function(parms, t2, model="contraction"){
  # First, sample a uniformly distributed number
  x <- runif(1, 0, 1)

  ## Then, choose the associated probability
  tmp <- function(x, parms, t1, t2, model){
    P.eca.at.t(parms, t1, t2, model) - x
  }
  if(P.eca.at.t(parms, t2, t2, model)==0){
    return(10^8)
  }

  t <- uniroot(tmp, c(0, t2), x=x, parms=parms, t2 = t2, model=model)$root
  t
}


#' Simulate mutation densities and their 95% CI using population genetic model of neuroblastoma initiation.
#'
#' This method takes as input the parameter estimates obtained by fitting a population genetics model to NB initiation. It then simulates the model using the parameter estimates, which can be useful for plotting the fitting result.
#' @param parameter.samples a data frame with the parameter samples
#' @param P.MRCA a data frame with the measured densities at MRCA among high-risk tumors (TMM). Rownames must identify the individual tumors. Must contain 1 column named density in units #SNVs/haploid genome
#' @param P.MRCA.lr a data frame with the measured densities at MRCA among low-risk tumors (no TMM). Rownames must identify the individual tumors. Must contain 1 column named density in units #SNVs/haploid genome
#' @param P.ECA  a data frame with the measured densities at ECA among high-risk tumors (TMM). Rownames must identify the individual tumors. Must contain 1 column named density in units #SNVs/haploid genome
#' @param mode Whether tumor initiation should be modeled in a homeostatic tissue ("homeostasis") or in a transient tissue ("decay")
#' @return A list containing lower and upper bounds of the 95% CI of the simulated probabilities at the measured incidences of ECA and MRCA
#' simulateCI()
#' @export


simulateCI <- function(parameter.samples, measured.mutation.times, measured.mutation.times.lr=NA,
                        measured.mutation.times.eca, mode="expansion_only"){

  if(mode=="decay"){
    model <- "contraction"
  }else{
    model <- "homeostasis"
  }
  ## HR tumors MRCA
  simulated.incidence <- matrix(0, nrow=nrow(parameter.samples), ncol=nrow(measured.mutation.times))
  for(i in 1:nrow(parameter.samples)){
    set.seed(Sys.time())
    parms <- as.list(parameter.samples[i,])

    ## rescale measured mutations into time (gamma-distributed)
    time.points <- sapply(sort(measured.mutation.times$Density), function(x){rgamma(n = 1, shape = x, rate = parms$mu)})

    ## then at each time point, take the probability of having acquired the 2 hits
    P <- sapply(time.points, function(t){
      P.2nd.driver(parms, t, model=model)
    })

    simulated.incidence[i,] <- P

  }

  min.probabilities <- apply(simulated.incidence, 2, quantile, p=0.025)
  max.probabilities <- apply(simulated.incidence, 2, quantile, p=0.975)

  ## LR tumors MRCA; only done for the decay scenario
  if(mode=="decay" & !is.na(measured.mutation.times.lr)){
    simulated.incidence <- matrix(0, nrow=nrow(parameter.samples), ncol=length(c(measured.mutation.times.lr$Density,
                                                                                 measured.mutation.times$Density)))
    for(i in 1:nrow(parameter.samples)){
      set.seed(Sys.time())
      parms <- as.list(parameter.samples[i,])

      ## rescale measured mutations into time
      time.points <- sapply(sort(c(measured.mutation.times.lr$Density, measured.mutation.times$Density)), function(x){rgamma(n = 1, shape = x, rate = parms$mu)})
      time.points <- sort(time.points)
      ## then at each time point, take the probability of having acquired the 1st hit

      P <- sapply(time.points, function(t){
        P.1st.driver(parms, t)
      })

      simulated.incidence[i,] <- P/max(P)

    }

    min.probabilities.lr <- apply(simulated.incidence, 2, quantile, p=0.025)
    max.probabilities.lr <- apply(simulated.incidence, 2, quantile, p=0.975)

  }else{
    min.probabilities.lr <- NA
    max.probabilities.lr <- NA
  }

  ## ECA
  simulated.incidence <- matrix(0, nrow=nrow(parameter.samples)*1, ncol=length(measured.mutation.times.eca$Density))
  j <- 1
  for(i in rep(1:nrow(parameter.samples), each=1)){
    set.seed(Sys.time())
    parms <- as.list(parameter.samples[i,])

    ## rescale measured mutations into time; MRCA (gamma distributed)
    time.points <- sapply(sort(measured.mutation.times[rownames(measured.mutation.times.eca),]$Density), function(x){rgamma(n = 1, shape = x, rate = parms$mu)})
    ## and ECA
    time.points.eca <- sapply(sort(measured.mutation.times.eca$Density), function(x){rgamma(n = 1, shape = x, rate = parms$mu)})

    ## then sample a time point for the ECA
    x <- runif(length(time.points), 0, 1)

    ## Then, choose the associated probability
    tmp <- function(x, parms, t1, t2){
      P.eca.at.t(parms, t1, t2, model=model) - x
    }

    t <- apply(rbind(time.points, x), 2, function(x){
      uniroot(tmp, c(0, x[1]), x=x[2], parms=parms, t2 = x[1])$root
    })

    incidence <- sapply(sort(time.points.eca), function(x){
      sum(t <= x)
    })/length(time.points.eca)


    simulated.incidence[j,] <- incidence
    j <- j + 1
  }

  min.probabilities.eca <- apply(simulated.incidence, 2, quantile, p=0.025)
  max.probabilities.eca <- apply(simulated.incidence, 2, quantile, p=0.975)

  return(list(min.mrca = min.probabilities,
              max.mrca = max.probabilities,
              min.lr = min.probabilities.lr,
              max.lr = max.probabilities.lr,
              min.eca = min.probabilities.eca,
              max.eca = max.probabilities.eca))

}

