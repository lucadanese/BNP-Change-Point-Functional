n_s <- function(j,rho_n){ # number of observations inside group j
  
  n_s <- rho_n[j]
  
  return(n_s)
  
}

n_s_1 <- function(j,rho_n){  # number of observations inside group j+1
  
  n_s_1 <- rho_n[j + 1]
  
  return(n_s_1)
}

n_g_k <- function(rho_n){ 
  
  n_g_k <- length(which(rho_n > 1))
  
  return(n_g_k)
  
}

gamma_splitting_MULTIVARIATE <- function(y,rho_n){ # it takes Y and rho_n and split the observations in groups
  
  gamma_k <- list()
  
  for (i in 1:length(rho_n)){
    ifelse(i != 1, {gamma_k[[i]] <- y[(sum(rho_n[1:i-1]) + 1) : (rho_n[i] + sum(rho_n[1:i-1])),]},
           {gamma_k[[i]] <- y[1:rho_n[i],]})
  }
  
  return(gamma_k)
}


# ALPHA SPLIT -------------------------------------------------------------

alpha_split <- function(y,j,rho_n_proposal,rho_n,m_0){  
  
  
  gamma_k_proposal <- gamma_splitting_MULTIVARIATE(y,rho_n_proposal)
  k_proposal <- length(rho_n_proposal)
    
  gamma_k <- gamma_splitting_MULTIVARIATE(y,rho_n)
  k <- length(rho_n)
  
  if (k < n & k > 1){
    
    a1 = log((1-q)/q) + (posterior(k_proposal, gamma_k_proposal,rho_n_proposal) - posterior(k,gamma_k, rho_n)) + log(((n_g_k(rho_n) * (n_s(j,rho_n) - 1) )/(k)))
      
    a2 = log(1)
    
    alpha = min(a1,a2)
    
  } 
  
  if (k == 1){  #Potrebbe essere inutile questa parte perché è inglobata sopra 
    
    a1 = log((1-q) * (n-1)) + (posterior(2,gamma_k_proposal,rho_n_proposal) - posterior(1,gamma_k,rho_n))
      
    a2 = log(1)
    
    alpha = min(a1,a2)
    
  }
  
  return(alpha)
  
}


# ALPHA MERGE -------------------------------------------------------------

alpha_merge <- function(y,j,rho_n_proposal,rho_n,m_0){  #INSERIRE INPUT
  
  gamma_k_proposal <- gamma_splitting_MULTIVARIATE(y,rho_n_proposal)
  k_proposal <- length(rho_n_proposal)
  
  gamma_k <- gamma_splitting_MULTIVARIATE(y,rho_n)
  k <- length(rho_n)

  if (k < n & k > 1){
    a1 = log((q/(1-q))) + (posterior(k_proposal, gamma_k_proposal,rho_n_proposal)-posterior(k,gamma_k,rho_n)) + log(((k-1)/((n_g_k(rho_n) + 1) * (n_s(j,rho_n) + n_s_1(j,rho_n) - 1))))
    
    a2 = log(1)
    
    alpha = min(a1,a2)
  } 
  
  if (k == n){
    a1 = log(q * (n-1)) + (posterior(k_proposal, gamma_k_proposal,rho_n_proposal) - posterior(k,gamma_k,rho_n))
    
    a2 = log(1)
    
    alpha = min(a1,a2)
  }  
  
  return(alpha)
}


# ALPHA SHUFFLE -------------------------------------------------------------

alpha_shuffle <- function(y,rho_n_proposal,rho_n){    
  
  gamma_k_proposal <- gamma_splitting_MULTIVARIATE(y,rho_n_proposal)
  k_proposal <- length(rho_n_proposal)
  
  gamma_k <- gamma_splitting_MULTIVARIATE(y,rho_n)
  k <- length(rho_n)
  
  a1 = {posterior(k_proposal,gamma_k_proposal,rho_n_proposal) - posterior(k,gamma_k,rho_n)}
  
  a2 = log(1)
  
  return(min(a1,a2))
}










