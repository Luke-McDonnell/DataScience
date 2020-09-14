# 2 Users, 2 Websites. 
# 1 step includes:
# - 1 marble randomly selected (prob of 1/2)
# - Then decide if the selected marble moves to other pot (1/3 prob it does move)

# state space (n=2 and g=2) is S = {(0,2), (1,1), (2,0)}

# transition probability matrix.
library(expm)

probs <- matrix(rep(0,9),nrow=3,ncol=3)
probs[1,1] <- 2/3
probs[1,2] <- 1/3
probs[2,1] <- 1/6
probs[2,2] <- 2/3
probs[2,3] <- 1/6
probs[3,2] <- 1/3
probs[3,3] <- 2/3
probs

Equil_dist <- matrix(rep(0,3),nrow=3,ncol=1)
Equil_dist[1] <- 1/4 
Equil_dist[2] <- 1/2
Equil_dist[3] <- 1/4
Equil_dist


X <- rep(1,2) # individual descriptions: initial state with both users in website 1
Niter <- 100000 # number of iterations
K <- rep(0,Niter) # vector of amount of people on site 1

for (i in 1:Niter) {
  k <- sample(2, size = 1) # selects a person at random 
  if (runif(1)<1/3){ # probability of 1/3 moving
    if (X[k] == 1) {  
      X[k] <- 0  # person moves from site 1 to site 2
    }
    else { 
      X[k] <- 1 # person moves from site 2 to site 1
    }
  }
  K[i] <- sum(X) # Save total number of users in site 1
}

freq <- rep(0,3)
for (j in 1:3) {
  freq[j] <- length(which(K==j-1))/Niter
}
freq


########## Plots #############

# Difference in Empirical and Invariant Distr 

state = c(0:2)
plot(state,freq,xlab="State",ylab="Probability of state", col="blue", main="Difference in Empirical and Invariant Distr")
lines(state,Equil_dist,type="h", col="red")
legend("topright", c("Empirical distribution","Theoretical invariant distribution."), fill=c("blue","red"))


######### Running Mean (Converges to theoretical mean as number of steps increase)

m <- ((Equil_dist[1]*0)+(Equil_dist[2]*1)+(Equil_dist[3]*2)) # expected value of k (sum of the states*probs)
n <- 10
Mean = seq(0,Niter/n-1)
for (i in 1:Niter/n) {Mean[i]=mean(K[1:i])}
TMean = rep(m,Niter/n -1)


plot (1:10000, Mean, xlab="Number of Steps", ylab="Mean State", col="blue", main="Mean state with increasing number of steps")
lines(TMean, col="red")
legend("topright", c("Mean States","Expected Value"), fill=c("blue","red"))

# Running variance (sampled every n Monte Carlo steps)

Var = seq(0,Niter/n-1)
for (i in 1:Niter/n) {Var[i]=var(K[1:i])}
v <- (((0-m)^2)*Equil_dist[1]) + (((1-m)^2)*Equil_dist[2]) + (((2-m)^2)*Equil_dist[3]) # theoretical variance of k 
TVar = rep(v,Niter/n-1)
plot(1:10000, Var, xlab="Number of Steps", ylab="Variance", col="red", main="Variance with increasing number of steps")
lines(TVar, col='blue')
legend("bottomright", c("Empirical Variance","Expected Variance"), fill=c("red","blue"))

