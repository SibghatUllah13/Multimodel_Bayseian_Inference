

# Initialize --------------------------------------------------------------
set.seed(1)
dugong = read.table('dugong-data.txt',header = TRUE)
Xi = dugong$Age
Yi = dugong$Length
number_of_sample = length(Yi)
nsims = 10000
n_size=1
# Generate Alpha ----------------------------------------------------------
alpha_sd=15
pa=pnorm(1,0,alpha_sd)
pb=pnorm(Inf,0,alpha_sd)
alpha=alpha_sd*qnorm(runif(n_size,pa,pb))
# Generate Beta -----------------------------------------------------------
beta_sd=10
pa_b=pnorm(1,0,beta_sd)
pb_b=pnorm(Inf,0,beta_sd)
beta=beta_sd*qnorm(runif(n_size,pa_b,pb_b))
beta

# Generate Gamma ----------------------------------------------------------
gam=runif(n_size,0,1)
# Generate Tau ------------------------------------------------------------
tau=1/rgamma(n_size,2,3)

# Generate mu -------------------------------------------------------------
gen_mu = function(alpha,beta,gam,Xi) 
{  return(alpha-(beta*((gam)^Xi)))
  
}
mu_samples = gen_mu(alpha,beta,gam,Xi)
mu_samples

# Set the Parameters ------------------------------------------------------
alpha_pos = array(alpha,dim = nsims)
beta_pos = array(beta, dim = nsims)
gamma_pos = array(gam,dim = nsims)
tau_pos = array(tau,dim = nsims)


# Likelihood --------------------------------------------------------------
log_likelihood = function(y_vec,x_vec,alp,bet,gam,tau)
{
  total = 0
  for (i in 1:length(y_vec))
  {
    temp = gam ^ x_vec[i]
    temp = bet * temp
    temp = alp - temp
    temp = temp **2
    temp = y_vec[i] - temp
    total= total+temp
  }
  total / -(2 * (tau^2))
}



# Run the Simulation of Metropolis to Get Gamma Posterior ------------------------------------------------------
n = 10000
current_gam = 0.634 
samps = rep(NA,n)
for (i in 1:n)
{proposed = rnorm(1, current_gam, 0.4) 
  while(proposed<0 | proposed>1){
    proposed = rnorm(1, current_gam, 0.4) 
  }
  logr = log_likelihood(Yi,Xi,alpha,beta,proposed,tau)-log_likelihood(Yi,Xi,alpha,beta,current_gam,tau)
  if (log(runif(1))< logr) {current_gam = proposed}
  samps[i] = current_gam
}
#gamma
acf(samps,col='blue',lwd=3)
gamma_pos=samps


# Gibbs Sampling to Simulate alpha,beta and Tau ---------------------------
require(truncnorm)
require(truncnorm)
for (i in 2:nsims){
  alpha_var=alpha_sd^2
  alpha_val=sum(Yi+beta_pos[i-1]*(gamma_pos[i-1]^Xi))
  alpha_pos_mean=(alpha_var*alpha_val)/(alpha_var*number_of_sample+tau_pos[i-1])
  alpha_pos_std=tau_pos[i-1]*alpha_var/(alpha_var*number_of_sample+tau_pos[i-1])
  k=rnorm(1,alpha_pos_mean,alpha_pos_std)
  alpha_pos[i] = k 
  
  beta_var = beta_sd^2
  beta_val = sum((alpha_pos[i]-Yi)*(gamma_pos[i-1]^Xi))
  beta_pos_mean = beta_var*beta_val/(beta_var+tau_pos[i-1])
  beta_pos_std = tau_pos[i-1]*beta_var/(beta_var+tau_pos[i-1])
  k_b=rtruncnorm(1,1,Inf,beta_pos_mean,beta_pos_std)
  beta_pos[i] = k_b 
  
  tau_pos_shape = (2+number_of_sample)
  tau_pos_rate = (6+sum((Yi-(alpha_pos[i]-beta_pos[i]*(gamma_pos[i-1]^Xi)))^2))/2
  tau_pos[i] = 1/rgamma(1,tau_pos_shape,tau_pos_rate)}



# Part d, Trace Plots -----------------------------------------------------
# For Alpha ---------------------------------------------------------------
par(mfrow=c(4,1))
plot(alpha_pos[1:100],type="l",main=paste("1:100 Iterations"),col='blue')
title(sub=paste("Mean estimate =",round(mean(alpha_pos),4),"Variance estimate =",round(var(alpha_pos),4)))
plot(alpha_pos[500:length(alpha_pos)],type="l",main="Burnt the first 500 simulations",col='red')
title(sub=paste("Mean estimate =",round(mean(alpha_pos[501:length(alpha_pos)]),4),"Variance estimate =",round(var(alpha_pos[501:length(alpha_pos)]),4)))
hist(alpha_pos[500:length(alpha_pos)],freq=F,col='orchid')
lines(density(alpha_pos[500:length(alpha_pos)]),col='blue')
acf(alpha_pos,col='blue')
# For Beta ----------------------------------------------------------------
par(mfrow=c(4,1))
plot(beta_pos[1:100],type="l",main=paste("1:100 Iterations"),col='blue')
title(sub=paste("Mean estimate =",round(mean(beta_pos),4),"Variance estimate =",round(var(beta_pos),4)))
plot(beta_pos[500:length(beta_pos)],type="l",main="Burnt the first 500 simulations",col='red')
title(sub=paste("Mean estimate =",round(mean(beta_pos[501:length(beta_pos)]),4),"Variance estimate =",round(var(beta_pos[501:length(alpha_pos)]),4)))
hist(beta_pos[500:length(beta_pos)],freq=F,col='orchid')
lines(density(beta_pos[500:length(beta_pos)]),col='blue')
acf(beta_pos,col='blue')
# For Gamma ---------------------------------------------------------------
par(mfrow=c(4,1))
plot(gamma_pos[1:100],type="l",main=paste("1:100 Iterations"),col='blue')
title(sub=paste("Mean estimate =",round(mean(gamma_pos),4),"Variance estimate =",round(var(gamma_pos),4)))
plot(gamma_pos[500:length(gamma_pos)],type="l",main="Burnt the first 500 simulations",col='red')
title(sub=paste("Mean estimate =",round(mean(gamma_pos[501:length(gamma_pos)]),4),"Variance estimate =",round(var(gamma_pos[501:length(alpha_pos)]),4)))
hist(gamma_pos[500:length(gamma_pos)],freq=F,col='orchid')
lines(density(gamma_pos[500:length(gamma_pos)]),col='blue')
acf(gamma_pos,col='blue')
# For Tau -----------------------------------------------------------------
par(mfrow=c(4,1))
plot(tau_pos[1:100],type="l",main=paste("1:100 Iterations"),col='blue')
title(sub=paste("Mean estimate =",round(mean(tau_pos),4),"Variance estimate =",round(var(tau_pos),4)))
plot(tau_pos[500:length(beta_pos)],type="l",main="Burnt the first 500 simulations",col='red')
title(sub=paste("Mean estimate =",round(mean(tau_pos[501:length(tau_pos)]),4),"Variance estimate =",round(var(tau_pos[501:length(alpha_pos)]),4)))
hist(tau_pos[500:length(tau_pos)],freq=F,col='orchid')
lines(density(tau_pos[500:length(tau_pos)]),col='blue')
acf(tau_pos,col='blue')







# Part e ------------------------------------------------------------------
par(mfrow=c(1,1))
# Empirical Behaviour of Alpha with Growing T -----------------------------
sum_vector=cumsum(alpha_pos)
v = rep(NA,length(alpha_pos))
for (i in 1:length(v))
{
  v[i] = sum_vector[i] / i
}
length_vector = seq(1,nsims)
plot(length_vector,v,type = 'p',xlab = 't = 1....T',ylab='empirical_mean of alpha',col='orchid')



# Empricial Bheaviour of Beta with Growing T ------------------------------
sum_vector=cumsum(beta_pos)
v = rep(NA,length(beta_pos))
for (i in 1:length(v))
{
  v[i] = sum_vector[i] / i
}
length_vector = seq(1,nsims)
plot(length_vector,v,type = 'p',xlab = 't = 1....T',ylab='empirical_mean for beta',col='orchid')


# Empirical Behaviour of Gamma with Growing T -----------------------------
sum_vector=cumsum(gamma_pos)
v = rep(NA,length(gamma_pos))
for (i in 1:length(v))
{
  v[i] = sum_vector[i] / i
}
length_vector = seq(1,nsims)
plot(length_vector,v,type = 'p',xlab = 't = 1....T',ylab='empirical_mean for gamma',col='orchid')


# Empirical Behaviour of Tau with Growing Tau -----------------------------
sum_vector=cumsum(tau_pos)
v = rep(NA,length(tau_pos))
for (i in 1:length(v))
{
  v[i] = sum_vector[i] / i
}
length_vector = seq(1,nsims)
plot(length_vector,v,type = 'p',xlab = 't = 1....T',ylab='empirical_mean for tau',col='orchid')


# Part f ------------------------------------------------------------------
# Estimate for alpha and its error ----------------------------------------
maximum_uncertainty = rep(NA,4)
batch_function = function(chain,bins){
  overall_mean = mean(chain)
  result = split(chain, cut(seq(chain), bins, labels = FALSE))
  
}
batch_size = nsims / 10
mean_array = rep(NA,10)
alpha_error = 0
index = 1
for (i in 1:10){
  current_batch = alpha_pos[index:(i * batch_size) ]
  mean_array[i] = mean(current_batch)
  mu_i = mean(current_batch) - mean(alpha_pos)
  mu_i = mu_i **2
  alpha_error=alpha_error+mu_i
  index = index+batch_size
}
mean(mean_array)
alpha_error
upper_val = mean(alpha_pos)+qnorm(0.95)*sqrt(alpha_error)/sqrt(nsims)
lower_val = mean(alpha_pos)-qnorm(0.95)*sqrt(alpha_error)/sqrt(nsims)
cat(lower_val,upper_val)
maximum_uncertainty[1]=upper_val-lower_val
# Estimate for Beta and Its Error -----------------------------------------
batch_size = nsims / 10
mean_array = rep(NA,10)
beta_error = 0
index = 1
for (i in 1:10){
  current_batch = beta_pos[index:(i * batch_size) ]
  mean_array[i] = mean(current_batch)
  mu_i = mean(current_batch) - mean(beta_pos)
  mu_i = mu_i **2
  beta_error=beta_error+mu_i
  index = index+batch_size
}
mean(mean_array)
beta_error
upper_val = mean(beta_pos)+qnorm(0.95)*sqrt(beta_error)/sqrt(nsims)
lower_val = mean(beta_pos)-qnorm(0.95)*sqrt(beta_error)/sqrt(nsims)
cat(lower_val,upper_val)
maximum_uncertainty[2]=upper_val-lower_val


# Estimate for Gamma and Its Error ----------------------------------------
batch_size = nsims / 10
mean_array = rep(NA,10)
gamma_error = 0
index = 1
for (i in 1:10){
  current_batch = gamma_pos[index:(i * batch_size) ]
  mean_array[i] = mean(current_batch)
  mu_i = mean(current_batch) - mean(gamma_pos)
  mu_i = mu_i **2
  gamma_error=gamma_error+mu_i
  index = index+batch_size
}
mean(mean_array)
gamma_error 
upper_val = mean(gamma_pos)+qnorm(0.95)*sqrt(gamma_error)/sqrt(nsims)
lower_val = mean(gamma_pos)-qnorm(0.95)*sqrt(gamma_error)/sqrt(nsims)
cat(lower_val,upper_val)
maximum_uncertainty[3]=upper_val-lower_val
# Estimate for Tau and Its Error ------------------------------------------
batch_size = nsims / 10
mean_array = rep(NA,10)
tau_error = 0
index = 1
for (i in 1:10){
  current_batch = tau_pos[index:(i * batch_size) ]
  mean_array[i] = mean(current_batch)
  mu_i = mean(current_batch) - mean(tau_pos)
  mu_i = mu_i **2
  tau_error=tau_error+mu_i
  index = index+batch_size
}
mean(mean_array)
tau_error 
upper_val = mean(tau_pos)+qnorm(0.95)*sqrt(tau_error)/sqrt(nsims)
lower_val = mean(tau_pos)-qnorm(0.95)*sqrt(tau_error)/sqrt(nsims)
cat(lower_val,upper_val)
maximum_uncertainty[4]=upper_val-lower_val
# Part g ------------------------------------------------------------------
which.max(maximum_uncertainty)
cat("The Maximum Uncertainty is for Alpha")

# Part h ------------------------------------------------------------------
Corr = rep(NA,6)
Corr[1]=cor(alpha_pos[500:10000],beta_pos[500:10000]) 
Corr[2]=cor(alpha_pos[500:10000],gamma_pos[500:10000])
Corr[3]=cor(alpha_pos[500:10000],tau_pos[500:10000]) 
Corr[4]=cor(beta_pos[500:10000],gamma_pos[500:10000]) 
Corr[5]=cor(beta_pos[500:10000],tau_pos[500:10000]) 
Corr[6]=cor(gamma_pos[500:10000],tau_pos[500:10000]) 
which.max(Corr)
cat("Alpha and Gamma have highest Correlation")

# Part i ------------------------------------------------------------------
Precision = rep(NA,2)
par(mfrow = c(2,1))
x_new = 20
updated_mean = alpha_pos-beta_pos*(gamma_pos^x_new)
updated_tau = tau_pos
draw_sample = 0
for (i in 1:length(updated_mean)){
  draw_sample = draw_sample+ rnorm(100,updated_mean[i],updated_tau[i])
}
post_pred_samp = draw_sample/nsims
hist(post_pred_samp,col='orchid')
Precision[1]=var(post_pred_samp)
# Part J ------------------------------------------------------------------
x_30 = 30
updated_mean_30 = alpha_pos-beta_pos*(gamma_pos^x_30)
updated_tau_30 = tau_pos
num = 1000
draw_sample = 0
for (i in 1:length(updated_mean_30)){
  draw_sample = draw_sample+ rnorm(100,updated_mean_30[i],updated_tau_30[i])}
post_pred_samp = draw_sample/nsims
hist(post_pred_samp,col='orchid')
Precision[2]=var(post_pred_samp)
# Part k -------------------------------------------------------------------
which.max(Precision)
cat("Distribution for Dugong with 30 Age is less Precise")

