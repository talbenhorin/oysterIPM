rm(list=ls(all=TRUE))

library(ggplot2)
require(deSolve)

##2.2 Numerical Integration of SIR Model
sirmod = function(t, y, parms) {
  # Pull state variables from y vector
  S = y[1]
  I = y[2]
  R = y[3]
  # Pull parameter values from parms vector
  beta = parms["beta"]
  mu = parms["mu"]
  gamma = parms["gamma"]
  N = parms["N"]
  # Define equations
  dS = mu * (N - S) - beta * S * I/N
  dI = beta * S * I/N - (mu + gamma) * I
  dR = gamma * I - mu * R
  res = c(dS, dI, dR)
  # Return list of gradients
  list(res)
}

times = seq(0, 26, by = 1/10)
parms = c(mu = 0, N = 1, beta = 2, gamma = 1/2)
start = c(S = 0.999, I = 0.001, R = 0)
##start is proportion of population in each group

out=ode(y=start, times=times, func=sirmod, parms=
          parms)
out=as.data.frame(out)

p <- ggplot(out, aes(x = time)) +
  geom_line(aes(y = S), color = "darkgoldenrod") +
  geom_line(aes(y = I), color = "darkorchid") +
  geom_line(aes(y = R), color = "firebrick") +
  theme_classic()


#Calculate R0
R0=parms["beta"]/(parms["gamma"]+parms["mu"])
#Adjust margins to accommodate a second right axis
par(mar = c(5,5,2,5))
#Plot state variables
plot(x=out$time, y=out$S, ylab="Fraction", xlab="Time",
     type="l")
lines(x=out$time, y=out$I, col="red")
lines(x=out$time, y=out$R, col="green")
#Add vertical line at turnover point
xx=out$time[which.max(out$I)]
lines(c(xx,xx), c(1/R0,max(out$I)), lty=3)
#prepare to superimpose 2nd plot
par(new=TRUE)
#plot effective reproductive ratio (w/o axes)
plot(x=out$time, y=R0*out$S, type="l", lty=2, lwd=2,
     col="black", axes=FALSE, xlab=NA, ylab=NA,
     ylim=c(-.5, 4.5))
lines(c(xx, 26), c(1,1), lty=3)
#Add right-hand axis for RE
axis(side = 4)
mtext(side = 4, line = 4, expression(R[E]))
#Add legend
legend("right", legend=c("S", "I", "R",
                         expression(R[E])), lty=c(1,1,1, 2),
       col=c("black", "red", "green", "black"))


##2.3 Final Epidemic Size

require(rootSolve)
equil=runsteady(y=c(S=1-1E-5, I=1E-5, R=0),
                times=c(0,1E5), func=sirmod, parms=parms)
round(equil$y, 3)
##2 percent of S escape infection and 98 percent become infected by the end of the epidemic

#Candidate values for R0 and beta
R0 = seq(0.1, 5, length=50)
betas= R0 * 1/2
#Vector of NAs to be filled with numbers
f = rep(NA, 50)
#Loop over i from 1, 2, ..., 50
for(i in seq(from=1, to=50, by=1)){
  equil=runsteady(y=c(S=1-1E-5, I=1E-5,
                      R=0), times=c(0,1E5), func=sirmod,
                  parms=c(mu=0, N=1, beta=betas[i], gamma=1/2))
  f[i]=equil$y["R"]
}
plot(R0, f, type="l", xlab=expression(R[0]))
curve(1-exp(-x), from=1, to=5, add=TRUE, col="red")

##I don't understand what "f" represents in the above equation/graph - figure 2.3 in book
##final epidemic size is a function of R0


#Define function - fraction of susceptibles that escape infection
fn=function(x, R0){
  exp(-(R0*(1-x))) - x
}

1-uniroot(fn, lower = 0, upper = 1-1E-9,
          tol = 1e-9, R0=2)$root
## [1] 0.7968121
#For R0 = 2, final epidemic size is 79.6%
#check accuracy of approximation:
exp(-2)-uniroot(fn, lower = 0, upper = 1-1E-9,
                tol = 1e-9, R0=2)$root
## [1] -0.06785259
#approximation from 1-uniroot function. Approximation is offf by 6.7% points


##2.4 Open epidemic
##mu > 0 - recruitment of new susceptibles

times = seq(0, 52*50, by=.1)
#life expectancy of 50 years and stable population size gives mu a weekly birth rate of 1/(50*52)
#integrate model for 50 years
parms = c(mu = 1/(50*52), N = 1, beta = 2,
          gamma = 1/2)
start = c(S=0.19, I=0.01, R = 0.8)
#19% of initial pop is susceptible, 1% infected
out = as.data.frame(ode(y=start, times=times,
                        func=sirmod, parms=parms))
par(mfrow=c(1,2)) #Make room for side-by-side plots
plot(times, out$I, ylab="Fraction", xlab="Time",
     type="l")
plot(out$S, out$I, type="l", xlab="Susceptible",
     ylab="Infected")



##2.5 Phase Analysis

simod = function(t, y, parameters) {
  S = y[1]
  I = y[2]
  beta = parameters["beta"]
  mu = parameters["mu"]
  gamma = parameters["gamma"]
  N = parameters["N"]
  dS = mu * (N - S) - beta * S * I/N
  dI = beta * S * I/N - (mu + gamma) * I
  res = c(dS, dI)
  list(res)
}


require(phaseR)
#Plot vector field
fld=flowField(simod, x.lim=c(0.15,0.35), y.lim=c(0,.01),
              parameters=parms, system="two.dim", add=FALSE,
              ylab="I", xlab="S")
#Add trajectory
out = as.data.frame(ode(y = c(S=0.19, I=0.01), times=
                          seq(0, 52*100, by=.1), func=simod, parms=parms))
lines(out$S, out$I, col="red")

#Add S-isocline
curve(parms["mu"]*(1/x-1)/parms["beta"], 0.15, 0.35,
      xlab="S", ylab="I", add=TRUE)
#Add I-isocline
shat=(parms["gamma"]+parms["mu"])/parms["beta"]
lines(rep(shat, 2), c(0,0.01))
legend("topright", legend=c("Transient", "Isoclines"),
       lty=c(1, 1), col=c("red", "black"))

##I don't understand isoclines
##flowfield arrows are not showing up on my plot


#2.7

chainSIR=function(t, logx, params){
  x=exp(logx)
  u=params["u"]
  S=x[1]
  I=x[2:(u+1)]
  R=x[u+2]
  with(as.list(params),{
    dS = mu * (N - S) - sum(beta * S * I) / N
    dI = rep(0, u)
    dI[1] = sum(beta * S * I) / N - (mu + u*gamma) * I[1]
    if(u>1){
      for(i in 2:u){
        dI[i]= u*gamma * I[i-1] - (mu+u*gamma)* I[i]
      }
    }
    dR = u*gamma * I[u] - mu * R
    res=c(dS/S, dI/I, dR/R)
    list(res)
  })
}

times = seq(0, 10, by=1/52)
paras2 = c(mu = 1/75, N = 1, beta = 625,
           gamma = 365/14, u=1)
xstart2 = log(c(S=.06, I=c(0.001, rep(0.0001,
                                      paras2["u"]-1)), R = 0.0001))
out = as.data.frame(ode(xstart2, times, chainSIR,
                        paras2))
plot(times, exp(out[,3]), ylab="Infected", xlab=
       "Time", ylim=c(0, 0.01), type=’l’)
paras2["u"] =2
xstart2 = log(c(S=.06, I=c(0.001, rep(0.0001/
                                        paras2["u"], paras2["u"]-1)), R = 0.0001))
out2 = as.data.frame(ode(xstart2, times, chainSIR,
                         paras2))
lines(times, apply(exp(out2[,-c(1:2,length(out2))]),
                   1 ,sum), col=’blue’)
paras2["u"] =73
xstart2 = log(c(S=.06, I=c(0.001, rep(0.0001/
                                        paras2["u"], paras2["u"]-1)), R = 0.0001))
out3 = as.data.frame(ode(xstart2, times, chainSIR,
                         paras2))
lines(times, apply(exp(out3[,-c(1:2,length(out3))]),
                   1, sum), col=’red’, lwd=2, lty=2)
paras2["u"] =500
xstart2 = log(c(S=.06, I=c(0.001, rep(0.0001/
                                        paras2["u"], paras2["u"]-1)), R = 0.0001))
out4 = as.data.frame(ode(xstart2, times, chainSIR,
                         paras2))
lines(times, apply(exp(out4[,-c(1:2,length(out4))]),
                   1,sum, na.rm=TRUE), col=’green’)
legend("topright", legend=c("SIR", "u=2", "u=500",
                            "u=73 (H-S)"), lty=c(1,1,1,2), lwd=c(1,1,1, 2),
       col=c("black", "blue", "green", "red"))