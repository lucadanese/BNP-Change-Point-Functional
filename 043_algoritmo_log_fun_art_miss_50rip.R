#
rm(list = ls())
#

# LIBRARIES  --------

library("Rcpp")
library("RcppArmadillo")
library("BNPmix")
library("gmp")
library("coda")
library("LaplacesDemon")
library("MASS")
library("FAdist")
library("combinat")
library("mvtnorm")

# SEED --------

set.seed(123)

timestamp()

list_of_partitions <- list()
list_of_gamma <- list()
list_of_sigma <- list()
list_of_theta <- list()
VI_vec <- as.numeric()

for(sim in 1:1){
  
# DATA SIMULATION --------

gamma_sim_1 = 0.5
gamma_sim_2 = 0.5
gamma_sim_3 = 0.5

mu_1 <- c(-2,-2,-2)
mu_2 <- c(3,3,3)
mu_3 <- c(-1,-1,-1)


sigma_1 = sigma_2 = sigma_3 = matrix(0, nrow = 3, ncol = 3)


diag(sigma_1) = c(0.5,0.5,0.5)

#sigma_1[1,2]  = 0.2
#sigma_1[1,3]  = 0.3
#sigma_1[2,1]  = 0.2
#sigma_1[3,1]  = 0.3


diag(sigma_2) <- c(0.5,0.5,0.5)

diag(sigma_3) <- c(0.5,0.5,0.5)
#sigma_3[2,3]  = 0.1
#sigma_3[3,2]  = 0.1

data_scenario_1 <- as.data.frame(matrix(nrow = 300, ncol = 3))

data_scenario_1[1,] = mu_1

for(i in 2:100){ 
  
  data_scenario_1[i,] = gamma_sim_1*data_scenario_1[i-1,] + (1-gamma_sim_1)*mu_1 + mvrnorm(n = 1,mu = rep(0,3), Sigma =  sigma_1)
  
}


data_scenario_1[101,] = mu_2

for(i in 102:200){ 
  
  data_scenario_1[i,] = gamma_sim_2*data_scenario_1[i-1,] + (1-gamma_sim_2)*mu_2 + mvrnorm(n = 1,  mu = rep(0,3), Sigma =  sigma_2)
  
}

data_scenario_1[201,] = mu_3

for(i in 202:300){ 
  
  data_scenario_1[i,] = gamma_sim_3*data_scenario_1[i-1,] + (1-gamma_sim_3)*mu_3 + mvrnorm(n = 1, mu = rep(0,3), Sigma =  sigma_3)
  
}

data <- data_scenario_1


# PARAMETERS --------

#
m_0 = rep(0,3)

k_0 = 0.25

nu_0 = 4 

phi_0 = diag(1,nrow = 3, ncol = 3)
#

a = b = 1
c = 0.1

y <- data

n <- nrow(y)

rho_n_0 = c(150,150) # initial partition

k = length(rho_n_0)

frequenze = rep(0, n) #vettore frequenze punti di cambio


## IMPORT FUNZIONI C++

sourceCpp("wade.cpp")
sourceCpp("cpp_funz.cpp")

##

prior_c = 1 

prior_d = 1 


gamma_prova <- as.numeric()
accettato_gamma <- as.numeric()
media_gamma <- as.numeric()

sigma_prova <- as.numeric()
accettato_sigma <- as.numeric()
media_sigma <- as.numeric()

theta_prova <- as.numeric()
accettato <- as.numeric()
media_theta <- as.numeric()

# Missing Imputation 

miss_index = c()


z_missing = rep(0,n)
z_missing[miss_index] = 1

# SOURCE --------

source("01_funzioni_log_fun_art_miss.R")
source("02_verosomiglianza_log_fun_art_miss.R")
source("03_alpha_log_fun_art_miss.R")

# MCMC ALGORITHM --------

q <- runif(1)

Nsim = 6*10^3

rho_n <- rho_n_0 #all'inzio il vettore n coincide con la partitizione specificata in partenza


# inizializzazione gamma 
gamma_prova[1] = runif(1) 
accettato_gamma[1] <- 1
media_gamma[1] <- gamma_prova[1]
#


# inizializzazione sigma 
sigma_prova[1] = rbeta(1,1,1) 
accettato_sigma[1] <- 1
media_sigma[1] <- sigma_prova[1]
#

# inizializzazione theta
theta_prova[1] <- r_shifted_gamma(prior_c,prior_d,sigma_prova[1])
accettato[1] <- 1
media_theta[1] <- theta_prova[1]
#


partitions <- list()

effective_sample <- as.numeric()

for(step in 1:Nsim){
  
  gamma <- gamma_prova[step]
  sigma <- sigma_prova[step]
  theta <- theta_prova[step]
  
  u <- runif(1)
  
  if (u <= q*indicator_1(k,n) + indicator_2(k)){          
    
    
    if (length(rho_n) != n){      #SPLIT FUNCTION
      
      output_list <- split(k,rho_n) # split proposal
      
      rho_n_proposal <- output_list[[1]]
      
      j_proposal <- output_list[[2]]
      
      u <- log(runif(1))
      
      alpha = alpha_split(y,j_proposal,rho_n_proposal,rho_n)
      
      if (u <= alpha){rho_n <- rho_n_proposal}
      
      
      k <- length(rho_n)
      
      
    }
    
    
    
  }
  
  if (u > q*indicator_1(k,n) + indicator_2(k)){     #MERGE FUNCTION
    
    output_list <- merge(k,rho_n) # merge proposal
    
    rho_n_proposal <- output_list[[1]]
    
    j_proposal <- output_list[[2]]
    
    u <- log(runif(1))
    
    alpha = alpha_merge(y,j_proposal,rho_n_proposal,rho_n)
    
    if (u <= alpha){rho_n <- rho_n_proposal}
    
    k <- length(rho_n)
    
  }
  
  
  if (k > 1){     #SHUFFLE FUNCTION
    
    rho_n_proposal <- shuffle(k,rho_n) #  shuffle proposal
    
    u <- log(runif(1))
    
    alpha = alpha_shuffle(y,rho_n_proposal,rho_n)
    
    if (u <= alpha){rho_n <- rho_n_proposal}
    
    k <- length(rho_n)
    
  }
  
  partitions[[step]] = rho_n
  
  
  # INFERENZA SU GAMMA 
  
  gamma_k <- gamma_splitting_MULTIVARIATE(y,rho_n)
  
  tau = log(gamma/(1-gamma))
  
  sigma_tau = 0.45 # CAPIRE COME MODIFICARE 
  
  tau_star = tau + rnorm(1,0,sigma_tau) 
  
  gamma_star = exp(tau_star)/(1+exp(tau_star))
  
  deriv_inv_star = abs(exp(tau_star) / (1+exp(tau_star))^2)
  
  deriv_inv = abs(exp(tau) / (1+exp(tau))^2)
  
  alpha_MH <- full_conditional_gamma(n, k, rho_n, sigma, theta, gamma_star, gamma_k) + log(deriv_inv_star) - full_conditional_gamma(n, k, rho_n, sigma, theta, gamma, gamma_k) - log(deriv_inv) 
  
  
  if(log(runif(1))<= min(alpha_MH,log(1))){
    gamma_prova[step+1] = gamma_star
    accettato_gamma[step+1] = 1
  }else{
    gamma_prova[step+1] = gamma_prova[step]
    accettato_gamma[step+1] = 0  
  }
  
  gamma <-gamma_prova[step+1]
  
  # INFERENZA SU SIGMA 
  
  candidato <- rbeta(1,1,1)
  
  alpha_MH <- full_conditional_sigma(candidato,theta,k,rho_n,1,1,1,1) - full_conditional_sigma(sigma_prova[step],theta,k,rho_n,1,1,1,1)
  
  if(log(runif(1)) <= min(alpha_MH,log(1))){
    sigma_prova[step+1] = candidato
    accettato_sigma[step+1] = 1
  }else{
    sigma_prova[step+1] = sigma_prova[step]
    accettato_sigma[step+1] = 0  
  }
  
  
  sigma <- sigma_prova[step+1]
  
  
  # INFERENZA SU THETA 
  
  theta_prova[step+1] = full_conditional_theta(prior_c,prior_d,theta_prova[step],k)
  theta <- theta_prova[step+1]
  
  # MISSING IMPUTATION
  
  miss_index = which(z_missing == 1)
  
  # assign the correspondent group to each observation
  
  y_grouped <- y
  
  for(i in 1:length(rho_n)){
    
    if(i == 1){
      y_grouped[1:cumsum(rho_n)[i],"group"] = 1
    }
    
    if(i != 1){
      y_grouped[(cumsum(rho_n)[i-1] + 1):(cumsum(rho_n)[i]),"group"] = i
    }
    
    
  }
  
  # imputation cycle 
  
  for(miss in miss_index){
    
    
    
    j = y_grouped[miss,"group"] 
    
    
    # extract parameters for the Normal Inverse Wishart 
    
    
    gamma_k = gamma_splitting_MULTIVARIATE(y,rho_n)
    
    d <- ncol(as.data.frame(gamma_k[j]))
    
    n_j = rho_n[j]
    
    k_nj = compute_k_n(k_0, gamma, n_j)
    
    gamma_k_j = as.matrix(as.data.frame(gamma_k[j]))
    
    m_nj = compute_m_n(d,m_0,k_0,gamma,n_j,gamma_k_j,nbasis)
    
    V_nj = nu_0 + n_j
    
    S_nj = compute_s_n(d, m_0, k_0, nu_0, phi_0, gamma, n_j, gamma_k_j, nbasis)
    
    
    mu_j = rnorminvwishart(n=1, t(m_nj), k_nj, S_nj, V_nj)$mu
    
    lambda_j = rnorminvwishart(n=1, t(m_nj), k_nj, S_nj, V_nj)$Sigma
    
    # lambda and phi 
    
    m_j = ifelse(j!=1,cumsum(rho_n)[j-1],0)
    
    if(miss == m_j + 1){
      
      lambda_miss = gamma*y[miss + 1,] + (1-gamma)*mu_j
      
      phi_miss = lambda_j * (1- gamma^2)
      
    }
    
    
    if(miss < m_j + n_j & miss > m_j + 1){
      
      lambda_miss = (gamma*(y[miss + 1,] + y[miss - 1,]) / (1 + gamma^2)) + (1-gamma^2)*mu_j/(1+gamma^2)
      
      phi_miss = lambda_j * ((1- gamma^2)/(1+gamma^2))
      
    }
    
    if(miss == m_j + n_j & m_j != 0){
      
      lambda_miss = gamma*(y[n_j + m_j - 1,]) + (1-gamma)*mu_j
      
      phi_miss = lambda_j * (1- gamma^2)
      
    }
    
    y[miss,] = mvrnorm(1,t(lambda_miss), phi_miss)
    
    
  }
  
  # CREAZIONE VETTORE FREQUENZA PUNTI DI CAMBIO 
  
  if(step > 20000){
    
    frequenze[cumsum(rho_n)] = frequenze[cumsum(rho_n)] + 1
    
  }
  
  
  
  # EFFECTIVE SAMPLE SIZE 
  
  effective_sample[step] = length(rho_n) - 1 
  
  
}


# VI LOSS --------

partitions_no_burn_in <- partitions[(1*10^3 + 1):length(partitions)] #remove burn-in period
#partitions_no_burn_in <- partitions 

X <- matrix(data = NA, nrow = length(partitions) - 1*10^3, ncol = 300)

for (i in 1:length(partitions_no_burn_in)){
  for (j in 1:length(partitions_no_burn_in[[i]])){
    if (j==1){
      X[i,c(1:partitions_no_burn_in[[i]][j])] = j
    }
    if (j!=1){
      X[i,c((sum(partitions_no_burn_in[[i]][1:j-1])+1):sum(partitions_no_burn_in[[i]][1:j]))] = j
    }
  }
}

M <- psm(X)

dists <- VI_LB(clean_partition(X), psm_mat = M)
a <- X[which.min(dists),]

final_partition_VI <- as.numeric(table(a)); final_partition_VI
#


# GAMMA ESTIMATE --------

mean(gamma_prova) # estimated value 

# SIGMA ESTIMATE --------

mean(sigma_prova) # estimated value 

# THETA ESTIMATE --------

mean(theta_prova) # estimated value 



list_of_partitions[[sim]] <- final_partition_VI
list_of_gamma[[sim]] <- gamma_prova
list_of_sigma[[sim]] <- sigma_prova
list_of_theta[[sim]] <- theta_prova

timestamp()
print(sim)

}


vec_true <- as.character()

vec_true[1:100] <- "1"
vec_true[101:200] <- "2"
vec_true[201:300] <- "3"

for(j in 1:length(list_of_partitions)){
  partition <- cumsum(list_of_partitions[[j]])
  vec <- rep("0", n)
  if(length(partition) != 1){
    vec[1:partition[1]] <- "1"
    
    counter = 2
    for(i in 1:(length(partition)-1)){
      vec[(partition[i]+1):partition[i+1]] = as.character(counter)
      counter = counter + 1
    }
  }
  VI_vec[j] = VI(vec,vec_true)
} 



mean(unlist(lapply(list_of_partitions, FUN = length))) #AVERAGE NUMBER OF BLOCKS
mean(VI_vec) #AVERAGE VALUE OF VI 

print(mean(unlist(list_of_gamma))) #ESTIMATED GAMMA
print(mean(unlist(list_of_sigma))) #ESTIMATED SIGMA
print(mean(unlist(list_of_theta))) #ESTIMATED THETA





