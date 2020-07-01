library(statnet)
library(bbmle)
library(epimdr)

data(niamey)
head(niamey[, 1:5])

par(mar = c(5,5,2,5))
plot(niamey$absweek, niamey$tot_cases, type="b",
     xlab="Week", ylab="Incidence")
par(new=T)
plot(niamey$absweek, niamey$cum_cases, type="l",
     col="red", axes=FALSE, xlab=NA, ylab=NA, log="y")
axis(side = 4)
mtext(side = 4, line = 4, "Cumulative incidence")
legend("topleft", legend=c("Cases", "Cumulative"),
       lty=c(1,1), pch=c(1,NA), col=c("black", "red"))

#serial interval for measels is 10-12 days
#data is weekly
#V = 1.5 - 1.8 weeks 
#V = serial interval, average time between infection and reinfection

##Calculate R0 - errors in code?

fit=lm(log(cum_cases)~absweek, subset=absweek<7,data=niamey)
r=fit$coef["absweek"]
V=c(1.5, 1.8)
V*r+1


##Calculate R0
V = c(1.5, 1.8)
f = (5/7)/V
V * r + 1 + f * (1 - f) * (V * r)^2

##Maximum likelihood chain binomial model
llik.cb = function(S0, beta, I) {
  n = length(I)
  S = floor(S0 - cumsum(I[-n]))
  p = 1 - exp(-beta * (I[-n])/S0)
  L = -sum(dbinom(I[-1], S, p, log = TRUE))
  return(L)}

twoweek = rep(1:15, each = 2)
y = sapply(split(niamey$cases_1[1:30], twoweek), sum)
sum(y)

#result is total number of cases

S0cand=6500
betacand=seq(0,10, by=.1)
ll=rep(NA, length(betacand))
for(i in 1:length(betacand)){
  ll[i]=llik.cb(S0=S0cand, beta=betacand[i], I=y)
}
plot(ll~betacand, ylab="Neg log-lik", xlab=
       expression(beta))
betacand[which.min(ll)]

#If SO estimation is correct, beta should be 2.3

betacand=2.3
S0cand=seq(5920,8000, length=101)
ll=rep(NA, length=101)
for(i in 1:101){
  ll[i]=llik.cb(S0=S0cand[i], beta=betacand, I=y)
}
plot(ll~S0cand, ylab="Neg log-lik", xlab=
       expression(S[0]))

##
require(bbmle)
fit=mle2(llik.cb, start=list(S0=7085, beta=2.3),
         method="Nelder-Mead",data = list(I = y))
summary(fit)
## Maximum likelihood estimation
##
## Call:
## mle2(minuslogl = llik.cb, start = list(S0 = 7085,
## beta = 2.3), beta = 2), data = list(I = y))
##
## Coefficients:
## Estimate Std. Error z value Pr(z)
## S0 7.8158e+03 1.3022e+02 60.019 < 2.2e-16 ***
## beta 1.8931e+00 3.6968e-02 51.209 < 2.2e-16 ***
## ---
## Signif. codes:
## 0 ’***’ 0.001 ’**’ 0.01 ’*’ 0.05 ’.’ 0.1 ’ ’ 1
##
## -2 log L: 841.831

confint(fit)
#gives confidence interval

#correlation matrix
cov2cor(vcov(fit))