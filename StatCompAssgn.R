#Calling Necessary Libraries
library(alr4)
library(mosaic)
library(DT)
library(pander)
#Loading the data
datatable(Challeng)
df<-Challeng
#Logistic Regression (Frequentist)
chall.glm <- glm(fail>0 ~ temp+pres , data=Challeng, family=binomial)
summary(chall.glm) %>% pander()
#Plotting the fitted curve
plot( fail>0 ~ temp, data=Challeng, xlab="Outside Temperature at Time of Launch (Fahrenheit)", ylab='Probability of At least One O-ring Failing', pch=16, main="NASA Shuttle Launch Data from 1981 to 1985", xlim=c(30,85))
curve(exp(15.043-0.232*x)/(1+exp(15.043-0.232*x)), add=TRUE)
abline(v=31, lty=2, col="lightgray")
text(31,0.3,"Outside Temp on Jan. 28, 1986 was 31", pos=4, cex=0.7, col="lightgray")
abline(h=c(0,1), lty=1, col=rgb(.2,.2,.2,.2))
legend("right", bty="n", legend=c("Previous Launches"), pch=16, col="black", cex=0.7)
#Prediction of paati logit
library(ResourceSelection)
hoslem.test(chall.glm$y, chall.glm$fitted, g=6) %>% pander()
pred <- predict(chall.glm, data.frame(temp=c(31,45,55,65,75,85,95),pres=50), type='response')
pred
plot(c(31,45,55,65,75,85,95),pred,type="l",xlab="Temperature",
     ylab="Probability of failure")
#Eibar Bayesian Logistic
x1<-Challeng$temp
x2<-Challeng$pres
y<-Challeng$fail
for(i in 1:23)
{
  if(y[i]>0)
    y[i]<-1
  else
    y[i]<-0
}
y
library(MCMCpack)
#Model 1 (two explanatory variables)
out = MCMClogit(y~x1+x2, burnin=1000, mcmc=21000, b0=0, B0=.001)
densityplot(as.numeric(out[,1]),xlab="Beta_0")
densityplot(as.numeric(out[,2]),xlab="Beta_1")
densityplot(as.numeric(out[,3]),xlab="Beta_2")
hist(as.numeric(out[,2]),probability=T)
hist(as.numeric(out[,3]),probability=T)
#Plotting posterior of p for given x1,x2(intuitive)
library(boot)
temp=60
pres=15
fail.prob=inv.logit(out[,1]+temp*out[,2]+pres*out[,3])
p<-as.numeric(fail.prob)
p
par(mfrow=c(3,2))
densityplot(p)
hist(p,probability=T)
#Model 2 (one explanatory variable)
out = MCMClogit(y~x1, burnin=1000, mcmc=21000, b0=0, B0=.001)
plot(out)
hist(as.numeric(out[,1]),probability=T,col="cyan",xlab="Beta_0",main=paste("Histogram of" , "Intercept"))
hist(as.numeric(out[,2]),probability=T,col="cyan",xlab="Beta_1",main=paste("Histogram of" , "Slope"))
densityplot(as.numeric(out[,1]),xlab="Beta_0",col="red")
densityplot(as.numeric(out[,2]),xlab="Beta_1")


library(boot)
temp=87
fail.prob=inv.logit(out[,1]+temp*out[,2])
p<-as.numeric(fail.prob)
p
densityplot(p)
hist(p,probability = T,
     col="cyan",xlab="p",
     main=paste("Histogram for" , temp,"degree Farenheit"))
#Model 3: Probit Model
out = MCMCprobit(y~x1, burnin=1000, mcmc=21000, b0=0, B0=.001)
plot(out)
densityplot(as.numeric(out[,1]),xlab="Beta_0")
densityplot(as.numeric(out[,2]),xlab="Beta_1")
hist(as.numeric(out[,1]),probability=T,col="cyan",xlab="Beta_0",main=paste("Histogram of" , "Intercept"))
hist(as.numeric(out[,2]),probability=T,col="cyan",xlab="Beta_1",main=paste("Histogram of" , "Slope"))
#Plotting posterior of p for given x1(intuitive)
library(boot)
temp=87
fail.prob=pnorm(out[,1]+temp*out[,2])
p<-as.numeric(fail.prob)
hist(p,probability = T,
     col="cyan",xlab="p",
     main=paste("Histogram for" , temp,"degree Farenheit"))

#The above models were fitted using some R packages, only comparing 
#the two models are left.Still I have written some codes from scratch
#to do the whole work.
rm(list = ls())
library(alr4)

# Data loading
Datas <- Challeng
Datas <- Datas[,1:3]
rownames(Datas) <- NULL

normalize <- function(x)
{
  return((x-mean(x))/sd(x))
}


Datas <- data.frame(lapply(Datas[,c(1,2)], normalize),Datas$fail)
Datas$Datas.fail[which(Datas$Datas.fail == 2)] <- 1

## MCMC Algorithm
## Model 1
likelihood1 <- function(X,y,beta,lambda = 0.0001,M = 100)
{
  beta = matrix(beta,nrow = 1)
  a = exp(-(lambda/2)*beta%*%t(beta))
  b = exp(y*(beta%*%t(X)))
  c = exp(beta%*%t(X))
  return(M*a*prod(b/(1+c)))
}

MCMC.Sampler1 <- function(X,y,beta0,B,sg = c(1,1,1))
{
  X = cbind(rep(1,nrow(X)),X)
  beta0 = matrix(beta0,nrow = 1)
  post.sample = c(0,0,0)
  beta1 = beta0
  beta2 = matrix(c(0,0,0),nrow = 1)
  for(i in 1:B)
  {
    beta2[1] = beta1[1] + rnorm(1,0,sg[1])
    beta2[2] = beta1[2] + rnorm(1,0,sg[2])
    beta2[3] = beta1[3] + rnorm(1,0,sg[3])
    ratio = likelihood1(X,y,beta = beta2)/likelihood1(X,y,beta1)
    unif = runif(1)
    if(unif <= min(1,ratio)) beta1=beta2
    post.sample = rbind(post.sample,beta1)
  }
  return(post.sample)
}

# MCMC parameters
B = 5*10^4
n.thin = 3

# First we fit a logit model to get an initial guess
library(glmnet)
probit.mod <- glm(formula = Datas.fail ~ temp+pres,
                 data = Datas,family = binomial(link = "probit"))
summary(probit.mod)
summary(logit.mod)
# Running the MCMC sampler
Post.Sample1 = MCMC.Sampler1(X = Datas[,1:2],y = Datas$Datas.fail,beta0 = c(probit.mod$coefficients[1],probit.mod$coefficients[2],probit.mod$coefficients[3]),B,sg = c(3,3,3))
Post.Sample1 = (Post.Sample1)[-(1:(B/10)),]
n.length = nrow(Post.Sample1)
batch.size = floor(n.length/n.thin)
Post.Sample1 = Post.Sample1[n.thin*(1:batch.size),]
Post.Samp1 = data.frame(Post.Sample1)
names(Post.Samp1) <- c('b0','b1','b2')

# Posterior distributions of beta0,beta1
hist(Post.Samp1$b0,probability = TRUE)
hist(Post.Samp1$b1,probability = TRUE)
hist(Post.Samp1$b2,probability = TRUE)
densityplot(Post.Samp1$b0,xlab="Beta_0")
densityplot(Post.Samp1$b1,xlab="Beta_1")
densityplot(Post.Samp1$b2,xlab="Beta_2")

# means
mean(Post.Samp1$b0)
mean(Post.Samp1$b1)
mean(Post.Samp1$b2)

# plotting the mean cumulatively w.r.t sample size
b0.mean.cum <- cumsum(Post.Samp1$b0)/(1:nrow(Post.Samp1))
b1.mean.cum <- cumsum(Post.Samp1$b1)/(1:nrow(Post.Samp1))
b2.mean.cum <- cumsum(Post.Samp1$b2)/(1:nrow(Post.Samp1))
plot(b0.mean.cum,type = "l")
plot(b1.mean.cum,type = "l")
plot(b2.mean.cum,type = "l")

# Joint posterior density
library(ggplot2)

# b0,b1
ggplot(Post.Samp1, aes(x = b0, y = b1, fill = ..level..)) +
  stat_density_2d(geom = "polygon")
# b0,b2
ggplot(Post.Samp1, aes(x = b0, y = b2, fill = ..level..)) +
  stat_density_2d(geom = "polygon")
# b1,b2
ggplot(Post.Samp1, aes(x = b1, y = b2, fill = ..level..)) +
  stat_density_2d(geom = "polygon")

# Plotting the probability
Post.Prob1 <- function(x.point)
{
  x.norm = NULL
  x.norm[1] = (x.point[1] - mean(Challeng$temp))/sd(Challeng$temp)
  x.norm[2] = (x.point[2] - mean(Challeng$pres))/sd(Challeng$pres)
  x_val = matrix(c(1,x.norm),nrow = 1)
  y_reg = x_val%*%t(Post.Samp1)
  y_reg = as.vector(y_reg)
  Pi.Posterior <- exp(y_reg)/(1+exp(y_reg))
  return(list("samples" = Pi.Posterior,"post.mean" = mean(Pi.Posterior)))
}

Samples.PRob <- Post.Prob1(x.point = c(59,5))$samples
hist(Samples.PRob,probability = TRUE)
densityplot(Samples.PRob,cex.main="Posterior Density of p")

## Model 2
likelihood2 <- function(X,y,beta,lambda = 0.0001,M = 100)
{
  beta = matrix(beta,nrow = 1)
  a = exp(-(lambda/2)*beta%*%t(beta))
  b = exp(y*(beta%*%t(X)))
  c = exp(beta%*%t(X))
  return(M*a*prod(b/(1+c)))
}

MCMC.Sampler2 <- function(X,y,beta0,B,sg = c(1,1,1))
{
  X = cbind(rep(1,nrow(X)),X)
  beta0 = matrix(beta0,nrow = 1)
  post.sample = c(0,0,0)
  beta1 = beta0
  beta2 = matrix(c(0,0,0),nrow = 1)
  for(i in 1:B)
  {
    beta2[1] = beta1[1] + rnorm(1,0,sg[2])
    beta2[2] = beta1[2] + rnorm(1,0,sg[2])
    beta2[3] = 0
    ratio = likelihood1(X,y,beta = beta2)/likelihood1(X,y,beta1)
    unif = runif(1)
    if(unif <= min(1,ratio)) beta1=beta2
    post.sample = rbind(post.sample,beta1)
  }
  return(post.sample)
}

# MCMC parameters
B = 5*10^4
n.thin = 2

# First we fit a logit model to get an initial guess
library(glmnet)
logit.mod <- glm(formula = Datas.fail ~ temp ,
                 data = Datas,family = "binomial")
# Running the MCMC sampler
Post.Sample2 = MCMC.Sampler2(X = Datas[,1:2],y = Datas$Datas.fail,beta0 = c(logit.mod$coefficients[1],logit.mod$coefficients[2],0),B,sg = c(3,3,3))
Post.Sample2 = (Post.Sample2)[-(1:(B/10)),]
n.length = nrow(Post.Sample2)
batch.size = floor(n.length/n.thin)
Post.Sample2 = Post.Sample2[n.thin*(1:batch.size),]
Post.Samp2 = data.frame(Post.Sample2)
names(Post.Samp2) <- c('b0','b1','b2')

# Posterior distributions of beta0,beta1
hist(Post.Samp2$b0,probability = TRUE)
hist(Post.Samp2$b1,probability = TRUE)


# means
mean(Post.Samp2$b0)
mean(Post.Samp2$b1)


# plotting the mean cumulatively w.r.t sample size
b0.mean.cum <- cumsum(Post.Samp2$b0)/(1:nrow(Post.Samp2))
b1.mean.cum <- cumsum(Post.Samp2$b1)/(1:nrow(Post.Samp2))
b2.mean.cum <- cumsum(Post.Samp2$b2)/(1:nrow(Post.Samp2))
plot(b0.mean.cum,type = "l")
plot(b1.mean.cum,type = "l")
plot(b2.mean.cum,type = "l")

# Joint posterior density
library(ggplot2)

# b0,b1
ggplot(Post.Samp2, aes(x = b0, y = b1, fill = ..level..)) +
  stat_density_2d(geom = "polygon")

# Plotting the probability
Post.Prob <- function(x.point)
{
  x.norm = (x.point - mean(Datas$temp))/sd(Datas$temp)
  x_val = matrix(c(1,x.norm),nrow = 1)
  y_reg = x_val%*%t(Post.Samp2)
  y_reg = as.vector(y_reg)
  Pi.Posterior <- exp(y_reg)/(1+exp(y_reg))
  return(list("samples" = Pi.Posterior,"post.mean" = mean(Pi.Posterior)))
}

# Plotting
Post.Prob2 <- function(x.point)
{
  x.norm = NULL
  x.norm[1] = (x.point[1] - mean(Challeng$temp))/sd(Challeng$temp)
  x_val = matrix(c(1,x.norm),nrow = 1)
  y_reg = x_val%*%t(Post.Samp2)
  y_reg = as.vector(y_reg)
  Pi.Posterior <- exp(y_reg)/(1+exp(y_reg))
  return(list("samples" = Pi.Posterior,"post.mean" = mean(Pi.Posterior)))
}

Samples.PRob <- Post.Prob1(x.point = c(60,40))$samples
hist(Samples.PRob,probability = TRUE)
densityplot(Samples.PRob,cex.main="Posterior Density of p")



## Model Selection :Bayes Factor
m_logit <- function(X,y,lambda = 0.01,N = 10^3,null = TRUE)
{
  beta = c()
  Total = 0
  X = cbind(rep(1,nrow(X)),X)
  for(i in 1:N)
  {
    beta[1] = rnorm(1,0,sd = 1/lambda)
    beta[2] = rnorm(1,0,sd = 1/lambda)
    beta[3] = rnorm(1,0,sd = 1/lambda)
    
    if(null){
      beta[3] = 0
    }
    
    beta = matrix(beta,byrow = TRUE,nrow = 1)
    b = exp(y*(beta%*%t(X)))
    c = exp(beta%*%t(X))
    M = prod(b/(1+c))
    Total = M + Total
  }
  return(Total/N)
}

m_probit <- function(X,y,lambda = 0.01,N = 10^3,null = TRUE)
{
  beta = c()
  Total = 0
  X = cbind(rep(1,nrow(X)),X)
  for(i in 1:N)
  {
    beta[1] = rnorm(1,0,sd = 1/lambda)
    beta[2] = rnorm(1,0,sd = 1/lambda)
    beta[3] = rnorm(1,0,sd = 1/lambda)
    
    if(null){
      beta[3] = 0
    }
    
    beta = matrix(beta,byrow = TRUE,nrow = 1)
    b = pnorm(beta%*%t(X))^y
    c = (1-pnorm(beta%*%t(X)))^(1-y)
    M = prod(b*c)
    Total = M + Total
  }
  return(Total/N)
}

set.seed(481)
a = m_logit(X = Datas[,1:2],y = Datas$Datas.fail,N = 10^4,lambda = 1000,null = TRUE)
set.seed(481)
b = m_probit(X = Datas[,1:2],y = Datas$Datas.fail,N = 10^4,lambda = 1000,null = FALSE)
(BF=a/b)
densityplot(Post.Samp1$b1)
(Savage_Dickey_BF=sqrt(2*pi)*0.02)
