# Random Walk of a Poisson Distribution. 

#install.packages('tseries')
#install.packages('circular')

library(tseries)
library(circular)

####### my function from question 1 #######

gen_pois <- function(lambda, n) {
  # This loop will create n iid poison distributed values with mean lambda
  result <- rep(NA,n)
  for (i in 1:n) {
    x <- 0
    prod <- 1
    while (TRUE) {
      unif1 <- runif(1, min=0, max=1)
      prod <- prod*unif1
      if (prod>exp(-lambda)) {
        x <- x+1
      }
      else {
        result[i] <- x
        break
      }
    }
  }
  return(result)
}

# Created a loop for increasing sample sizes of the random walk
ns <- c(10, 100, 1000, 5000)
jb_test <- c() # save p-values from JB test 
sw_test <- c() # save p-values from SW test 
ks_test <-c() # save p-values from KS test 
index <- 1 # loop index
for (N in ns) {
     # mean of each variable 
    lambda <- 10
    sum_lambda <- N*lambda # The mean of the sums
    
    # Find the theoretical CDF values 
    k <- seq(0, 2*sum_lambda, 1) # k values to put into the pdf
    pdf1 <- dpois(k, sum_lambda) # calculate PDF values
    cdf1 <- rep(0,2*sum_lambda) # calculate CDF values
    for (i in 1:(2*sum_lambda)) {
      cdf1[i] <- sum(pdf1[1:i])
    }
    
    # Monte Carlo simulation to find Empiracal CDF values
    Niter <- 10000
    ZN <- rep(NA,Niter)
    for (i in 1:Niter) {
      ZN[i] <- sum(gen_pois(lambda, N))
    } 
    
    # Plot the Empiracal pdf histogram
    hist(ZN, breaks='FD', probability='TRUE', xlab='ZN values', col='purple', main=paste('Histogram of ZN (lambda = ', lambda, ', N = ', N, ')'))
    
    # plot the Theoretical and Empiracal CDFs together
    plot(ecdf(ZN), main=paste('Theoretical and Empirical CDF (lambda = ', lambda, ', N = ', N, ')'), xlab='K', col='blue')
    lines(cdf1, col='red')
    legend("bottomright", c("Empirical CDF","Theoretical CDF"), fill=c("blue","red"))
    
    ########## normalisation #############
    
    # For poisson distributions, E[X]=Var[X]=lambda.
    UN = (ZN - N*lambda)/sqrt(N*lambda)
    
    # Compute the theoretical pdf of a standard normal distribution - N(0,1) 
    sigma <- 1
    mu <- 0
    x = seq(-10, 10, 0.1)
    norm_pdf <- dnorm(x)
    
    hist(UN, breaks='FD', probability='TRUE', main=paste('UN converging to a standard normal distribution (lambda = ', lambda, ', N = ', N, ')'), col='red')
    # lines(x, norm_pdf, col='blue', type='h')
    lines(x, norm_pdf, col='blue', type='l')
    legend("bottomright", c("Empirical Values (UN)","Theoretical pdf"), fill=c("red", "blue"))
    
    # Kolmogorov-Smirnov test
    KS <- ks.test(UN, "pnorm")$p.value
    ks_test[index] <- KS
    
    # Jarque Bera Test of normality - p-value close to 0 means it suggests normality
    JB <- jarque.bera.test(UN)$p.value
    jb_test[index] <- JB
    
    # PP plot of UN
    pp.plot(UN, ref.line=TRUE)
    text(0.5, 0.95, labels=paste('PP plot for N =', N), cex=2)
    
    # Shapiro wilk test for normality (can only use a sample size up to 5000 for this test )
    SW <- shapiro.test(UN[1:5000])$p.value
    sw_test[index] <- SW
    
    index <- index+1
}

