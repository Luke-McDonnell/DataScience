
kTJ = 1.9 # temp
betaJ = 1/kTJ # inverse temp
Niter = 100000 # number of iterations 

s <- rep(1,5) # Vector of states with initial condition s[i]=1

H = -betaJ*(s[1]*s[2]+s[2]*s[3]+s[3]*s[4]+s[4]*s[5]) # initial energy

E <- rep(0,Niter+1) # Energy Vector
M <- rep(0,Niter+1) # Magnetisation vector
E[1] <- H #initial energy
M[1] <- sum(s) #initial magnetization

for (i in 1:Niter) {
  k <- sample(5, size = 1) # selects a spin
  Hs = -(s[1]*s[2]+s[2]*s[3]+s[3]*s[4]+s[4]*s[5]) #old hamiltonian
  s[k]<- -s[k] # the spin is flipped
  Hy = -(s[1]*s[2]+s[2]*s[3]+s[3]*s[4]+s[4]*s[5]) #new hamiltonian
  DeltaH <- -betaJ*(Hy - Hs) 
  alpha <- min(1,exp(DeltaH)) # Generate a random uniform number in [0,1]. If less than alpha then accept; otherwise reject.
  U <- runif(1)
  if (U<=alpha) # the move is accepted, energy and magnetization are updated
  {E[i+1]<-E[i]+DeltaH
  M[i+1]<-sum(s)}
  else # The move is rejected and the new energies and magnetisations stay the same
  {s[k] <- -s[k]
  E[i+1]<- E[i]
  M[i+1]<- M[i]}
}


# Average absolute magnetization per spin
memp <- mean(abs(M[5000:10000]))/5


# Running average
n <- 10
Mean = seq(0,Niter/n-1)
for (i in 1:Niter/n) {Mean[i]=mean(abs(M[1:i]))/5} #Calculating the mean
Mean

m <- 0.57 # expected value of m
TMean = rep(m,Niter/n -1) #Theoretical mean
plot (Mean,xlab="MC steps per spin",ylab="|m| (kT/J = 1.9)", main = "Running Average vs Expected Value", col='red')
lines(TMean, col = "blue")
legend("topright", c("Expected value","Running mean."), fill=c("blue", "red"))


# Plot of results
kTJ <- c(0.1,0.3,0.5,0.7,0.9,1,1.1,1.3,1.5,1.7,1.9)
m <- (2*exp(4/kTJ) + (24/5) + (8/5)*exp(-2/kTJ) + (16/5)*exp(2/kTJ) + (2/5)*exp(-4/kTJ)) / (8*exp(2/kTJ) + 8*exp(-2/kTJ) + 12 + 2*exp(4/kTJ)  + 2*exp(-4/kTJ) )
memp <- round(m, digits=2)
plot(kTJ,memp,xlab="kT/J",ylab="|m|",col="red", main = "Change in temperature against magnetisation")
lines(kTJ,m, col='blue')
legend("topright", c("Expected value","Empirical mean"), fill=c("blue", "red"))

