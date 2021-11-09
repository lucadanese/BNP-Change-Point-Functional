#include "RcppArmadillo.h"
#include <math.h> 
// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
float lgamma_multi(int d, int a){
  
  arma::vec vettore(d, arma::fill::zeros);
  
  //float ret;
  
  for(int j = 0; j < d ; j++){
    
    float j_float = j;
    
    vettore[j] = std::lgamma(a + (-j_float/2)); // da 1-j a -j perchè sarebbe (1-(j+1))
    
  }
  
  return log(std::pow(M_PI,(d*(d-1)/4))) + arma::accu(vettore);
  
}

// [[Rcpp::export]]
float inner_likelihood(int d, arma::vec m_0, float k_0, int V_0, arma::mat phi_0, float gamma, int n_gamma, arma::mat gamma_k, int nbasis){
  
  float V_n;
  float k_n;
  float value;
  
  arma::mat vettore_m_n(nbasis, n_gamma - 1, arma::fill::zeros); // inizializziamo la matrice m_n

  arma::mat vettore_phi_n(nbasis, nbasis, arma::fill::zeros);     // inizializziamo la matrice Phi_n
  
  arma::vec m_n(d);
  
  arma::mat phi_n(d,d, arma::fill::zeros);
  
  V_n = V_0 + (n_gamma);
  
  k_n = k_0 + (std::pow((1-gamma),2)/(1 - std::pow(gamma,2)))*(n_gamma-1) + 1 ;
  
  // m_n
  
  if(n_gamma != 1){
    
    for(int i = 1; i < n_gamma; i++){

      vettore_m_n.col(i-1) = gamma_k.row(i).t() - (gamma * gamma_k.row(i-1).t());
      
     
    }

    m_n = (1/k_n)*( (m_0*k_0) + gamma_k.row(0).t() + ((1-gamma)/(1 - std::pow(gamma,2))) * sum(vettore_m_n,1));
    
    
    
  }
  else{
    
    m_n = (1/k_n)*((m_0*k_0) + gamma_k.row(0).t());
    
  }
  
  // Phi_n
  
  if(n_gamma != 1){
    
    for(int i = 1; i < n_gamma; i++){
     
     
        vettore_phi_n = vettore_phi_n + ((gamma_k.row(i).t() - gamma * gamma_k.row(i-1).t()) * (gamma_k.row(i).t() - gamma * gamma_k.row(i-1).t()).t())/(1-std::pow(gamma,2));
      
      }

        phi_n = phi_0 + (gamma_k.row(0).t() * gamma_k.row(0)) + vettore_phi_n + ((m_0*m_0.t())*k_0) - ((m_n * m_n.t())*k_n) ;
    
      
    }
    else{

       phi_n = phi_0 + (gamma_k.row(0).t() * gamma_k.row(0)) + ((m_0 * m_0.t())*k_0) - ((m_n * m_n.t())*k_n) ;
      
    }
    
    
    value = log(pow(k_n,(-d/2))) - log(std::pow(k_0,(-d/2))) + (V_0/2) * log(arma::det(phi_0)) + lgamma_multi(d, (V_n/2)) - (V_n/2)*log(arma::det(phi_n)) - lgamma_multi(d,(V_0/2)) + log(std::pow(M_PI,(-n_gamma*d/2))) - log(pow((1-pow(gamma,2)),((n_gamma-1)/2)));
  
  return value;
}

// /*** R
// inner_likelihood(d,m_0,k_0,V_0,phi_0,gamma,n_gamma)
// */


// [[Rcpp::export]]
arma::mat prodottino(int d, arma::mat phi_0){
  
  arma::mat final; 
  
  final = phi_0.row(1) - d * phi_0.row(1);
  
  return final;
  
}


// [[Rcpp::export]]
float compute_k_n(float k_0, float gamma, int n_gamma){
 
  float k_n;
  float value;

  k_n = k_0 + (std::pow((1-gamma),2)/(1 - std::pow(gamma,2)))*(n_gamma-1) + 1 ;
    
  value = k_n;
  
  return value;
}


// [[Rcpp::export]]
arma::mat compute_s_n(int d, arma::vec m_0, float k_0, int V_0, arma::mat phi_0, float gamma, int n_gamma, arma::mat gamma_k, int nbasis){
  
  float k_n;
  
  arma::mat vettore_m_n(nbasis, n_gamma - 1, arma::fill::zeros); // inizializziamo la matrice m_n

  arma::mat vettore_phi_n(nbasis, nbasis, arma::fill::zeros);     // inizializziamo la matrice Phi_n
  
  arma::vec m_n(d);
  
  arma::mat phi_n(d,d, arma::fill::zeros);
  
  k_n = k_0 + (std::pow((1-gamma),2)/(1 - std::pow(gamma,2)))*(n_gamma-1) + 1 ;
  
  // m_n
  
  if(n_gamma != 1){
    
    for(int i = 1; i < n_gamma; i++){

      vettore_m_n.col(i-1) = gamma_k.row(i).t() - (gamma * gamma_k.row(i-1).t());
      
     
    }

    m_n = (1/k_n)*( (m_0*k_0) + gamma_k.row(0).t() + ((1-gamma)/(1 - std::pow(gamma,2))) * sum(vettore_m_n,1));
    
    
    
  }
  else{
    
    m_n = (1/k_n)*((m_0*k_0) + gamma_k.row(0).t());
    
  }
  
  // Phi_n
  
  if(n_gamma != 1){
    
    for(int i = 1; i < n_gamma; i++){
     
     
        vettore_phi_n = vettore_phi_n + ((gamma_k.row(i).t() - gamma * gamma_k.row(i-1).t()) * (gamma_k.row(i).t() - gamma * gamma_k.row(i-1).t()).t())/(1-std::pow(gamma,2));
      
      }

        phi_n = phi_0 + (gamma_k.row(0).t() * gamma_k.row(0)) + vettore_phi_n + ((m_0*m_0.t())*k_0) - ((m_n * m_n.t())*k_n) ;
    
      
    }
    else{

       phi_n = phi_0 + (gamma_k.row(0).t() * gamma_k.row(0)) + ((m_0 * m_0.t())*k_0) - ((m_n * m_n.t())*k_n) ;
      
    }
    

  
  return phi_n;
}



// [[Rcpp::export]]
arma::vec compute_m_n(int d, arma::vec m_0, float k_0, float gamma, int n_gamma, arma::mat gamma_k, int nbasis){
  
  float k_n;
  
  arma::mat vettore_m_n(nbasis, n_gamma - 1, arma::fill::zeros); // inizializziamo la matrice m_n
  
  arma::vec m_n(d);
  
  arma::mat phi_n(d,d, arma::fill::zeros);
  
  k_n = k_0 + (std::pow((1-gamma),2)/(1 - std::pow(gamma,2)))*(n_gamma-1) + 1 ;
  
  // m_n
  
  if(n_gamma != 1){
    
    for(int i = 1; i < n_gamma; i++){

      vettore_m_n.col(i-1) = gamma_k.row(i).t() - (gamma * gamma_k.row(i-1).t());
      
     
    }

    m_n = (1/k_n)*( (m_0*k_0) + gamma_k.row(0).t() + ((1-gamma)/(1 - std::pow(gamma,2))) * sum(vettore_m_n,1));
    
    
    
  }
  else{
    
    m_n = (1/k_n)*((m_0*k_0) + gamma_k.row(0).t());
    
  }
  
  return m_n;
}


