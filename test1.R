rm(list=ls(all=TRUE))

library(ggplot2)
require(deSolve)

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

out=ode(y=start, times=times, func=sirmod, parms=
          parms)
out=as.data.frame(out)

p <- ggplot(out, aes(x = time)) +
  geom_line(aes(y = S), color = "darkgoldenrod") +
  geom_line(aes(y = I), color = "darkorchid") +
  geom_line(aes(y = R), color = "firebrick") +
  theme_classic()


