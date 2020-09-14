# Creating a function that will generate the Poisson Distribution from the uniform distribution

gen_pois <- function(lambda, n, p_value) {
  
  # This loop will create n iid poison distributed values with mean lambda
  result <- c()  # the values will be saved into this vector
  
  # algorithm from Devroye
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
  
  # Caculate the expected theoretical probabilities for this distribution 
  x <- seq(0, max(result), length=(max(result)+1))
  y <- dpois(x, lambda)        # pdf values by hand (exp(-lambda)*(lambda**x))/(factorial(x)) DELETE IF NOT NEEDED

  # plot the histogram of the observed values with the line of the expected values 
  hist(result, probability = TRUE, breaks=max(result)-min(result), xlab='Observed Values', col='orange', main = paste('Histogram of Generating', n, 'Poisson Variables with Rate', lambda))
  lines(x,y)
  
  # Chi-Square goodness-of-fit test: calculate the observed frequencies and the expected frequencies.
  
  expected_vals <- y*n # expected frequencies up to the max observed value
  expected_vals[max(result)+2] <- (1-sum(y)) * n # last bin of expected values in the expected frequency of more than the max observed value 
  
  observed_vals <- c()
  for (i in 0:(max(result))) {observed_vals[i+1] <- length(which(result==i))} # find the frequencies of each number up to the max result
  observed_vals[max(result)+2] <- 0 # 0 observed over the max value of the result

  
  # Calculate the Chi Suqared value
  chi_sq <- sum(((observed_vals-expected_vals)^2)/expected_vals)
  print(chi_sq)
  
  # Table value
  chi_q <- qchisq(p=p_value, df=(length(observed_vals)-2))
  print(chi_q)
  
  # hypothesis test, if 
  test_result <- chi_sq - chi_q
  print(test_result)
  
  if (test_result > 0) {
    print(paste('For sample size', n, 'and lambda =', lambda, 'We reject the null hypothesis and conclude there is a significant difference between these 2 distributionsat a significance level of p =', p_value))
    }
  else {
    print(paste('For sample size', n, 'and lambda =', lambda,'We accept the null hypothesis and conclude there is no significant difference between these 2 distributions at a significance level of p =', p_value))
    }
  
  return(test_result)
}


# test multiple lambdas with differing sample sizes to see 
lambdas <- c(10, 25, 50, 100)
ns <- c(10000, 100000, 1000000)
test <- rep(0, 12)
i <- 1
for (lam in lambdas) {
  for (n in ns) {
    test[i] <- gen_pois(lam, n, 0.1)
    i <- i+1
  }
}
test

