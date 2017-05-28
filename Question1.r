


# Some Suitable Options for the Choice of Hyper Parameters Alongside their  --------
#Possible Distribution for alpha, Note that We can only Use the 
# Part of Distribution from 1 on Wards
curve(dnorm(x,0,10),from = -10,to=30, col='orchid',lwd=4)
#Possible Distribution for beta, Note that We can only Use the 
# Part of Distribution from 1 on Wards
curve(dnorm(x,0,10),from = -10,to=30, col='orchid',lwd=4)
#tau2 
curve(1/dgamma(x,5,5), col='orchid',lwd=4,from = 25, to = 30)
# Part d ------------------------------------------------------------------
set.seed(1)
n_size=1
dugong = read.table('dugong-data.txt',header = TRUE)
Xi = dugong$Age
Yi = dugong$Length

# Generate alpha ----------------------------------------------------------
alpha_sd=10
a=1
b = Inf
pa=pnorm(a,0,alpha_sd)
pb=pnorm(b,0,alpha_sd)
alpha=alpha_sd*qnorm(runif(n_size,pa,pb))


# Generate Beta -----------------------------------------------------------
beta_sd=10
pa_b=pnorm(a,0,beta_sd)
pb_b=pnorm(b,0,beta_sd)
beta=beta_sd*qnorm(runif(n_size,pa_b,pb_b))


# Generate Gamma ----------------------------------------------------------
gam=runif(n_size,0,1)
# Generate Tau Square -----------------------------------------------------
tau=1/rgamma(n_size,5,5)



# Generate mu -------------------------------------------------------------
gen_mu_samples = function(alpha,beta,gam,Xi) 
{  return(alpha-(beta*((gam)^Xi)))
  
}
mu_samples = gen_mu_samples(alpha,beta,gam,Xi)


# Likelihood --------------------------------------------------------------
likelihood = function(vec,y_vec,x_i){
  alpha_vec=vec[1]
  beta_vec = vec[2]
  gam_vec = vec[3]
  tau=vec[4]
  mu = alpha_vec-(beta_vec*((gam_vec)^x_i))
  term_1 = sum((y_vec-mu)^2)/(tau^2)
  term_2 = -term_1 * (0.5)
  const = -0.5*length(y_vec)*log(2*pi*(tau^2))
  expres = -(const + term_2)
  return(expres)
}


# MLE Estimates -----------------------------------------------------------
mle=optim(c(1.5,1.5,0.5,1.5),likelihood,gr = NULL,Yi,Xi,hessian = TRUE)
alpha_estimated = mle$par[1]
beta_estimated=mle$par[2]
gam_estimatede=mle$par[3]
tau_estimated=mle$par[4]
estimates = c(alpha_estimated,beta_estimated,gam_estimatede,tau_estimated)
estimates
# Maximum a Posteriori  ---------------------------------------------------
log_prior = function(vec,y_vec,x_i){
  alpha_vec=vec[1]
  beta_vec = vec[2]
  gam_vec = vec[3]
  tau=vec[4]
  Xi = x_i
  mu = alpha_vec-(beta_vec*((gam_vec)^Xi))
  term_1 = sum((y_vec-mu)^2)/(tau^2)
  term_2 = -term_1 * (0.5)
  const = -0.5*length(y_vec)*log(2*pi*(tau^2))
  expres = -(const + term_2)
  val = -log( dnorm(alpha_vec,0,5)*dnorm(beta_vec,0,2)*dunif(gam_vec,0,1)*(1/dgamma(tau,2,3))) 
  return(expres+val)
}

# Use Maximum_Posteriori --------------------------------------------------
Map=optim(c(1.5,1.5,0.5,1.5),log_prior,gr = NULL,Yi,Xi,hessian = TRUE)
alpha_estimated = Map$par[1]
beta_estimated=Map$par[2]
gam_estimated=Map$par[3]
tau_estimated=Map$par[4]
Map_estimates = c(alpha_estimated,beta_estimated,gam_estimated,tau_estimated)

# Compare Both ------------------------------------------------------------
cat ("Maximum a Posteriori Estimates",Map_estimates)
cat ("Maximum Likelihood Estimates",estimates)



