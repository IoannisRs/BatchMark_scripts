  #include <RcppArmadilloExtensions/sample.h>
  #include <random>
  #include <iostream>
  #include <cmath>
  #include <vector>
  #include <numeric>
  #include <algorithm>
  #include <stdio.h>    
  #include <math.h>
  #include<Rmath.h>
  using namespace Rcpp;

const double log2pi = std::log(2.0 * M_PI);


// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::vec Mahalanobis(arma::mat x, arma::rowvec center, arma::mat cov){
  int n = x.n_rows;
  arma::mat x_cen;
  x_cen.copy_size(x);
  for (int i=0; i < n; i++) {
    x_cen.row(i) = x.row(i) - center;
  }
  return sum((x_cen * cov.i()) % x_cen, 1);    
}

// [[Rcpp::export]]
arma::vec dmvnorm_arma(arma::mat x,  arma::rowvec mean,  arma::mat sigma, bool log = false) { 
  arma::vec distval = Mahalanobis(x,  mean, sigma);
  double logdet = sum(arma::log(arma::eig_sym(sigma)));
  double log2pi = std::log(2.0 * M_PI);
  arma::vec logretval = -( (x.n_cols * log2pi + logdet + distval)/2  ) ;
  
  if (log){ 
    return(logretval);
  }else { 
    return(exp(logretval));
  }
}

// [[Rcpp::export]]

double log_dbinomial(double x, double N, double p){
  
  return R::lchoose(N,x) + x*std::log(p) + (N-x)*std::log(p);
  
}

// [[Rcpp::export]]

double p_lapl(double x, double m, double s){
  
  if(x <= m){
    return 0.5*exp((x-m)/s);
  }else{
    return 1-0.5*exp(-(x-m)/s);
  }
}

// [[Rcpp::export]]
arma::vec IndAvail(arma::mat n, int K){
  arma::vec N(K);
  N.zeros();
  for(int k=0; k<K; ++k){
    for(int j=0; j<(k+1); ++j){
      for(int i=0; i<(K-k); ++i){
        N[k] += n(i,j);
      }
    }
  }
  return N;
}

// [[Rcpp::export]]
arma::vec IndAvel(arma::mat n, int K){
  
  arma::vec N(K);
  N.zeros();
  
  for(int k=0; k<K; ++k){
    
    for(int j=0; j<(k+1); ++j){
      for(int i=0; i<(K-k); ++i){
        N[k] += n(i,j);
      }
    }
    
  }
  
  return N;
  
}

// [[Rcpp::export]]

double lgamma_2(double x){
  if(x != 0) {
    return lgamma(x);
  } else {
    return 1;
  }
}

// [[Rcpp::export]]
arma::mat Weights_n_0(int K, double p){
  
  arma::mat Omega(K+1,K+1);
  Omega.ones();
  
  for(int idx=0; idx<K; ++idx){
    
    for(int c=0; c<(idx+1); ++c){
      for(int r=0; r<(K-idx); ++r){
        Omega(r,c) *= (1-p);
      }
    }
    
  }
  return Omega; 
}


// [[Rcpp::export]]


arma::vec Pois_Lik_unif(List m_cup, List loss, arma::mat n_cand, arma::mat m_old, arma::mat O, int w, int K, double p){
  
  double log_cand=0;
  double log_old=0;
  arma::vec res(2);
  
  arma::mat n(K+1,K+1);
  n.zeros();
  n = n_cand;
  arma::mat m(K+1,K+1);
  m.zeros();
  m = m_old;
  
  arma::mat Weights(K+1,K+1);
  Weights.zeros();
  Weights = Weights_n_0(K,  p);
  
  // for(int idx=0; idx<K; ++idx){
    //  n += as<arma::mat>(m_cup[idx]) + as<arma::mat>(loss[idx]); 
    //    m += as<arma::mat>(m_cup[idx]) + as<arma::mat>(loss[idx]);
    //  }
  
  int k = 0; 
  for(int i=0; i<K; ++i){
    for(int j=0; j<(K-k); ++j){
      log_cand += R::dpois(n(i,j),w*O(i,j)*Weights(i,j),1);
      log_old +=  R::dpois(m(i,j),w*O(i,j)*Weights(i,j),1);
    }
    k = k + 1;
  }
  
  res[0] = log_cand;
  res[1] = log_old;
  
  return res;
}

// [[Rcpp::export]]
int PlusOrMinusOne() {
  return (rand() % 2) * 2 - 1;
}

// [[Rcpp::export]]
arma::mat Prop_Mat(arma::mat n, int m, int K){
  
  int r;
  int c;
  Rcpp::Environment pkg = Rcpp::Environment::namespace_env("extraDistr");
  Rcpp::Function f = pkg["rdunif"];
  SEXP val = f(1,-m,m);
  int l;
  l = Rcpp::as<int>(val);
  r = (rand()%(K));    
  c = (rand()%(K-r));
  n(r,c) = n(r,c) + (l*PlusOrMinusOne());
  return n;
}


// [[Rcpp::export]]
arma::mat Prop_Mat_cup(arma::mat n, arma::vec m, int K, int k){
  
  Rcpp::Environment pkg = Rcpp::Environment::namespace_env("extraDistr");
  Rcpp::Function f = pkg["rdunif"];
  SEXP val = f(1,-m[k],m[k]);
  int l;
  l = Rcpp::as<int>(val);
  int r1;
  int c1;
  int r2;
  int c2;
  r1 = (rand()%(K-k));    
  c1 = (rand()%(k+1));
  r2 = (rand()%(K-k));    
  c2 = (rand()%(k+1));
  n(r1,c1) = n(r1,c1) + l;
  n(r2,c2) = n(r2,c2) - l;
  return n;
}

// [[Rcpp::export]]
arma::mat Prop_Mat_cup_l(arma::mat n, arma::vec m, int K, int k){
  
  Rcpp::Environment pkg = Rcpp::Environment::namespace_env("extraDistr");
  Rcpp::Function f = pkg["rdunif"];
  SEXP val = f(1,-m[k],m[k]);
  int l;
  l = Rcpp::as<int>(val);
  int r1;
  int c1;
  int r2;
  int c2;
  r1 = K-1-k;
  r2 = K-1-k;
  c1 = (rand()%(k+1));
  c2 = (rand()%(k+1));
  n(r1,c1) = n(r1,c1) + l;
  n(r2,c2) = n(r2,c2) - l;
  
  return n;
}

// [[Rcpp::export]]  
double chooseC(double n, double k) {
  return Rf_choose(n, k);
}

// [[Rcpp::export]]
double Captured_Likelihood_unif(List m_cup, arma::mat m_obs, int k, int K, double p){
  
  arma::vec Like(K);
  Like.zeros();
  arma::vec M;
  M.zeros();
  M = IndAvel(m_cup[k],K);
  
  double s;
  s=0;
  
  for(int i=k+1; i<K; ++i){
    Like[i] = R::dbinom(m_obs(k,i),M[i],p,1);
  }
  s = accu(Like);
  
  
  return s;
}

// [[Rcpp::export]]

arma::vec Unmarked_Available_unif(List m_cup, List loss, arma::mat m_uncup, int K){
  
  arma::vec N_0(K);
  N_0.zeros();
  N_0 = IndAvel(m_uncup, K);
  
  arma::vec U(K);
  U.zeros();
  
  arma::mat N(K,K);
  N.zeros();
  
  for(int k=0; k<(K); ++k){
    N.submat(0,k,K-1,k) = (IndAvel(as<arma::mat>(m_cup[k]),K))  +  (IndAvel(as<arma::mat>(loss[k]),K));
    
  }
  
  for(int i=0; i<(K); ++i){
    U[i] = N_0[i] + accu(N.submat(i,i,i,K-1));
  }
  
  
  
  return U;
}




// [[Rcpp::export]]
arma::vec Avail_p(List m_cup, List loss, arma::mat m_uncup, int K){
  
  arma::vec N_0(K);
  N_0.zeros();
  N_0 = IndAvel(m_uncup, K);
  
  for(int i=0; i<K; ++i){
    N_0 += IndAvel(as<arma::mat>(m_cup[i]),K);
  }
  for(int i=(K-1); i>=0; --i){
    N_0.subvec(0,i) += IndAvel(as<arma::mat>(loss[i]),K).subvec(0,i);
  }
  
  return N_0;
}




// [[Rcpp::export]]

double Unmarked_Likelihood_unif(List m_cup, List loss, arma::mat m_uncup, arma::vec unmarked_obs, int K, double p){
  
  arma::vec U(K);
  U.zeros();
  U = Unmarked_Available_unif(m_cup, loss, m_uncup, K);
  
  
  double Like=0;
  
  for(int k=0; k<K; ++k){
    Like += R::dbinom(unmarked_obs[k],U[k],p,1);
  }
  
  
  return Like;
}










// [[Rcpp::export]]


int IndAvelSingle(arma::mat n, int K){
  
  int N=0;
  
  
  N = accu(IndAvel(n,K));
  
  
  return N;
}






// [[Rcpp::export]]
arma::mat matrix2vector(arma::mat m, const bool byrow=false){
  if (byrow) {
    return m.as_row();
  } else {
    return m.as_col();
  }
}


// [[Rcpp::export]]


arma::mat Normalization(arma::mat W){
  
  arma::mat Y;
  
  Y = W / accu(W);  
  return Y;
}


// [[Rcpp::export]]

arma::mat List_sum(List n, int K){
  
  arma::mat m(K+1,K+1);
  m.zeros();
  
  for(int k=0; k<K; ++k){
    m += as<arma::mat>(n[k]);
  }
  
  return m;
}


// [[Rcpp::export]]
List Weights(int K, double p){
  
  List weights(K);
  
  arma::mat W(K+1,K+1);
  
  for(int idx=0; idx<K; ++idx){
    W.ones();
    for(int c=0; c<(idx+1); ++c){
      for(int r=0; r<(K-idx); ++r){
        if(idx==c){
          W(r,c) *= p;
        }else{
          W(r,c) *= pow(1-p,idx-c)*p;
        }
      }
    }
    
    weights[idx] = W;
  }
  
  
  return weights;
}

// [[Rcpp::export]]
double Available_p(List M, List L, arma::mat n_0, int K){
  double res=0;
  for(int r=0; r<(K-1); r++){
    res += accu(IndAvail(as<arma::mat>(M[r]),K).subvec(r+1,K-1));
    
  }
  res += accu(Unmarked_Available_unif(M,L,n_0,K));
  return res;
}

// [[Rcpp::export]]
List grid_prob_fixed(arma::vec V, arma::vec W, int K){
  arma::mat Omega(K+1,K+1);
  Omega.ones();
  arma::mat A(K+1 , K+1);
  A.ones();
  for(int idx=0; idx<(K+1); ++idx){
    for(int r=0; r<(K+1-idx); ++r){
      Omega(r,idx) *= V[idx];
    }
  }
  int b;
  for(int idx=0; idx<K; ++idx){
    A.ones();
    int k_winver = 0;
    for(int i=idx; i<(K+1); ++i){
      if(i==idx){
        b = 0; 
      }else{
        b = 1;
      }
      for(int j=0; j<(K-k_winver+b-idx); ++j){
        if(i==idx){
          A(i,j) = W[idx];
        }else{
          A(i,j) = 1-W[idx];
        }
      } 
      k_winver = k_winver + 1;
    }
    Omega = Omega % A;  
  }
  return List::create(_["Omega"]=Omega,_["Vj"]=V,_["Vij"]=W);
}

// [[Rcpp::export]]
double Captured_Likelihood_unif_p(List m_cup, arma::mat m_obs, int k, int K, arma::vec p){
  arma::vec Like(K);
  Like.zeros();
  arma::vec M;
  M.zeros();
  M = IndAvail(as<arma::mat>(m_cup[k]),K);
  double s;
  s=0;
  for(int i=k+1; i<K; ++i){
    Like[i] = R::dbinom(m_obs(k,i),M[i],p[i],1);
  }
  s = accu(Like);
  return s;
}

// [[Rcpp::export]]
double Unmarked_Likelihood_unif_p(List m_cup, List loss, arma::mat m_uncup, arma::vec unmarked_obs, int K, arma::vec p){
  arma::vec U(K);
  U.zeros();
  U = Unmarked_Available_unif(m_cup, loss, m_uncup, K);
  double Like=0;
  for(int k=0; k<K; ++k){
    Like += R::dbinom(unmarked_obs[k],U[k],p[k],1);
  }
  return Like;
}


// [[Rcpp::export]]
List Weights_p(int K, arma::vec p, double lambda){
  List weights(K);
  arma::vec On(K);
  On.ones();
  arma::mat W(K+1,K+1);
  for(int idx=0; idx<K; ++idx){
    W.ones();
    for(int c=0; c<(idx+1); ++c){
      for(int r=0; r<(K-idx); ++r){
        if(idx==c){
          W(r,c) *= p[idx]*(1-lambda);
        }else{
          W(r,c) *= prod(On.subvec(c,idx-1)-p.subvec(c,idx-1))*p[idx]*(1-lambda);
        }
      }
    }
    weights[idx] = W;
  }
  return weights;
}

// [[Rcpp::export]]
List Weights_loss_p(int K, arma::vec p, double lambda){
  List weights(K);
  arma::vec On(K);
  On.ones();
  arma::mat W(K+1,K+1);
  for(int idx=0; idx<K; ++idx){
    W.ones();
    for(int c=0; c<(idx+1); ++c){
      for(int r=0; r<(K-idx); ++r){
        if(idx==c){
          W(r,c) *= p[idx]*lambda;
        }else{
          W(r,c) *= prod(On.subvec(c,idx-1)-p.subvec(c,idx-1))*p[idx]*lambda;
        }
      }
    }
    weights[idx] = W;
  }
  return weights;
}

// [[Rcpp::export]]
List Weights_loss_p_mult_l(int K, arma::vec p, arma::vec lambda){
  List weights(K);
  arma::vec On(K);
  On.ones();
  arma::mat W(K+1,K+1);
  for(int idx=0; idx<K; ++idx){
    W.ones();
    for(int c=0; c<(idx+1); ++c){
      for(int r=0; r<(K-idx); ++r){
        if(idx==c){
          W(r,c) *= p[idx]*lambda[idx];
        }else{
          W(r,c) *= prod(On.subvec(c,idx-1)-p.subvec(c,idx-1))*p[idx]*lambda[idx];
        }
      }
    }
    weights[idx] = W;
  }
  return weights;
}


// [[Rcpp::export]]
List Weights_p_mult_l(int K, arma::vec p, arma::vec lambda){
  List weights(K);
  
  arma::vec On(K);
  On.ones();
  
  arma::mat W(K+1,K+1);
  
  for(int idx=0; idx<K; ++idx){
    W.ones();
    for(int c=0; c<(idx+1); ++c){
      for(int r=0; r<(K-idx); ++r){
        if(idx==c){
          W(r,c) *= p[idx]*(1-lambda[idx]);
        }else{
          W(r,c) *= prod(On.subvec(c,idx-1)-p.subvec(c,idx-1))*p[idx]*(1-lambda[idx]);
        }
      }
    }
    
    weights[idx] = W;
  }
  return weights;
}


// [[Rcpp::export]]
double dmultinom_rcpp(int N, arma::vec x, arma::vec p){
  double den = 0;
  den = lgamma(N+1);
  int n = x.size();
  for(int i=0; i<n; ++i){
    den += x[i]*log(p[i]) - lgamma(x[i]+1);
  }
  return den;
}

// [[Rcpp::export]]
List MH_caps_p(List n_cand, List m_old, List loss, arma::mat m_0, 
               arma::mat O,  arma::vec unmarked_obs, arma::mat marked, int K,  arma::vec p, double lambda){
  
  List n_output(K);
  n_output = clone(m_old);
  List n_like(K);
  n_like = clone(m_old);
  
  List W = Weights_p(K,p,lambda);
  
  arma::mat x;
  arma::mat y;
  
  arma::vec t;
  arma::vec z;
  arma::vec w;
  arma::vec w1;
  
  arma::vec accept(K);
  accept.zeros();
  
  arma::mat weights(K+1,K+1);
  weights.zeros();
  
  
  for(int k=0; k<K; ++k){
    weights = as<arma::mat>(W[k]);
    x = as<arma::mat>(n_cand[k]);
    y =  as<arma::mat>(m_old[k]) ;
    double cand=0;
    double old=0;
    
    if( accu(as<arma::mat>(n_cand[k]) < 0) > 0  ){
      n_output[k] = as<arma::mat>(m_old[k]);
      n_like[k] = as<arma::mat>(m_old[k]);
    }else{
      n_like[k] = as<arma::mat>(n_cand[k]);
      t = matrix2vector(x.submat(0,0,K-k-1,k));
      z = matrix2vector(y.submat(0,0,K-k-1,k));
      w = matrix2vector(O.submat(0,0,K-k-1,k) / accu(O.submat(0,0,K-k-1,k)));
      w1 = matrix2vector(weights.submat(0,0,K-k-1,k));
      cand =  dmultinom_rcpp(accu(t), t,  w%w1)  + Captured_Likelihood_unif_p(n_like, marked, k, K, p) +
        Unmarked_Likelihood_unif_p(n_like, loss, m_0, unmarked_obs, K, p);
      
      old =  dmultinom_rcpp(accu(z), z,  w%w1)  +  Captured_Likelihood_unif_p(n_output, marked, k, K, p) +
        Unmarked_Likelihood_unif_p(n_output, loss, m_0, unmarked_obs, K, p);
      
      if( (R::runif(0, 1) < exp(cand - old)) & (std::isnan(cand)!=TRUE) ){
        
        n_output[k] = as<arma::mat>(n_cand[k]);
        accept[k] = 1;
        n_like[k] = as<arma::mat>(n_cand[k]);
        
      }else{
        
        n_output[k] = as<arma::mat>(m_old[k]);
        n_like[k] = as<arma::mat>(m_old[k]);
        
      }
      
    }
    
  }
  
  
  
  return List::create(_["n_output"] = n_output, _["accept"] = accept);
  
}

// [[Rcpp::export]]
List MH_loss_p(List m_cup, List loss_cand, List loss_old, arma::mat m_0, arma::mat O,
               arma::vec unmarked_obs, int K,  arma::vec p, double lambda){
  List n_output(K);
  n_output = clone(loss_old);
  List W = Weights_loss_p(K,p,lambda);//Weights_p(K,p);
  arma::mat x;
  arma::mat y;
  arma::vec t;
  arma::vec z;
  arma::vec w;
  arma::vec w1;
  arma::mat n_trans(K+1,K+1);
  n_trans.zeros();
  n_trans = n_trans + m_0;
  for(int idx=0; idx<K; ++idx){
    n_trans += as<arma::mat>(loss_old[idx]) + as<arma::mat>(m_cup[idx]);
  }
  arma::mat n_new(K+1,K+1); 
  n_new.zeros();
  List n_like(K);
  arma::vec accept(K);
  accept.zeros();
  n_like = clone(loss_old);
  arma::mat weights(K+1,K+1);
  weights.zeros();
  for(int k=0; k<K; ++k){
    weights = as<arma::mat>(W[k]);
    n_like[k] = as<arma::mat>(loss_cand[k]);
    n_new = n_trans - as<arma::mat>(loss_old[k]) + as<arma::mat>(n_like[k]);
    x = as<arma::mat>(loss_cand[k]);
    y = as<arma::mat>(loss_old[k]);
    double cand=0;
    double old=0;
    
    if( accu(as<arma::mat>(loss_cand[k]) < 0) > 0  ){
      
      n_output[k] = as<arma::mat>(loss_old[k]);
      n_like[k] = as<arma::mat>(loss_old[k]);
    }else{
      
      t = matrix2vector(x.submat(0,0,K-k-1,k)); 
      z = matrix2vector(y.submat(0,0,K-k-1,k));
      w = matrix2vector(O.submat(0,0,K-k-1,k) / accu(O.submat(0,0,K-k-1,k)));
      w1 = matrix2vector(weights.submat(0,0,K-k-1,k));
      
      cand =  dmultinom_rcpp(accu(t),t, w%w1)  + Unmarked_Likelihood_unif_p(m_cup, n_like, m_0, unmarked_obs, K, p);  
      
      old =  dmultinom_rcpp(accu(z),z, w%w1)  +  Unmarked_Likelihood_unif_p(m_cup, n_output, m_0, unmarked_obs, K, p);  
      
      if( (R::runif(0, 1) < exp(cand - old)) & (std::isnan(cand)!=TRUE) ){
        
        n_output[k] = as<arma::mat>(loss_cand[k]);
        accept[k] = 1;
        n_like[k] = as<arma::mat>(loss_cand[k]);
        
      }else{
        
        n_output[k] = as<arma::mat>(loss_old[k]);
        n_like[k] = as<arma::mat>(loss_old[k]);
        
      }
      
    }
    
  }
  
  
  
  return  List::create(_["n_output"] = n_output, _["accept"] = accept)  ;
  
}


// [[Rcpp::export]]
arma::mat Weights_n_0_p(int K, arma::vec p){
  arma::mat Omega(K+1,K+1);
  Omega.ones();
  for(int idx=0; idx<K; ++idx){
    for(int c=0; c<(idx+1); ++c){
      for(int r=0; r<(K-idx); ++r){
        Omega(r,c) *= (1-p[idx]);
      }
    }
  }
  return Omega; 
}


// [[Rcpp::export]]
arma::vec Pois_prior(List m_cup, List loss, arma::mat n_cand, arma::mat m_old, arma::mat O, double w, int K, arma::vec p){
  
  double log_cand=0;
  double log_old=0;
  arma::vec res(2);
  res.zeros();
  arma::mat Weights(K+1,K+1);
  Weights.zeros();
  Weights = Weights_n_0_p(K,  p);
  for(int i=0; i<(K); ++i){
    for(int j=0; j<(K-i); ++j){
      log_cand += R::dpois(n_cand(i,j),w*O(i,j)*Weights(i,j),1);
      log_old +=  R::dpois(m_old(i,j),w*O(i,j)*Weights(i,j),1);
    }
  }
  
  res[0] = log_cand;
  res[1] = log_old;
  
  return res;
}




// [[Rcpp::export]]

arma::mat MH_mis_p(List m_cup, List loss, arma::vec unmarked_obs, arma::mat n_cand, arma::mat n_old,
                   arma::mat O, int K, double w, arma::vec p){
  
  arma::vec log_Pois;
  log_Pois.zeros();
  
  double r=0;
  
  arma::mat n_output(K+1,K+1);
  n_output.zeros();
  
  if( accu(n_cand < 0) > 0  ){
    n_output = n_old;
  }else{
    log_Pois = Pois_prior(m_cup, loss, n_cand,n_old,O,w,K, p);
    r = (log_Pois[0] +  Unmarked_Likelihood_unif_p(m_cup, loss, n_cand, unmarked_obs, K, p) ) -
      ( log_Pois[1] + Unmarked_Likelihood_unif_p(m_cup, loss, n_old, unmarked_obs, K, p)  );
    if( (R::runif(0, 1) < exp(r)) & (std::isnan(r)!=TRUE) ){
      n_output = n_cand;
    }else{
      n_output = n_old;
    }
  }
  for(int i=0; i<(K+1); ++i){
    for(int j=0; j<(K+1); ++j){
      if( (i+j)==K ){
        n_output(i,j) = R::rpois(w*O(i,j));
      }
    }
  }
  return n_output;
}    


// [[Rcpp::export]]

arma::vec Available_mult_p(List M, List L, arma::mat n_0, int K){
  
  arma::vec res(K);
  res.zeros();
  res = Unmarked_Available_unif(M,L,n_0,K);
  for(int idx=1; idx<K; ++idx){
    for(int i=idx; i<K; ++i){
      res[i] += IndAvail(as<arma::mat>(M[idx-1]),K)[i]; 
    }
  }
  
  return res;
}


// [[Rcpp::export]]
List Generate_mult_p(arma::mat nij, int K, arma::vec l,arma::vec p){
  arma::mat u_trans(K+1,K+1);
  arma::mat r_trans(K+1,K+1);
  arma::mat n_0(K+1,K+1);
  n_0.zeros();
  arma::mat n_trans(K+1,K+1);
  n_trans = nij;
  List M(K);
  List R(K);
  arma::vec Unmarked(K);
  Unmarked.zeros();
  arma::vec lost(K);
  lost.zeros();
  for(int k=0; k<K; ++k){
    u_trans.zeros();
    r_trans.zeros();  
    for(int j=0; j<(k+1); ++j){
      for(int i=0; i<(K-k); ++i){
        u_trans(i,j) = R::rbinom(n_trans(i,j),p[k]);
        r_trans(i,j) = R::rbinom(u_trans(i,j),l[k]);
      }
    } 
    M[k] = u_trans-r_trans;
    Unmarked[k] = accu(u_trans);
    R[k] = r_trans;
    lost[k] = accu(r_trans);
    n_trans = n_trans - u_trans;
  }
  n_0 = n_trans;
  arma::mat Data(K,K);
  Data.zeros();
  for(int k=0; k<K; ++k){
    for(int i=k; i<K; ++i){
      if( i!=k ){
        Data(k,i) = R::rbinom(IndAvel(as<arma::mat>(M[k]),K)[i], p[i]);
      }else{
        Data(k,k) = accu(as<arma::mat>(M[k]));
      }
    }
  }
  return List::create(_["M"] = M, _["R"] = R, _["n_0"] = n_0, _["Data"] = Data, _["Unmarked"] = Unmarked, _["r"] = lost);
}


// [[Rcpp::export]]
arma::mat grid_cell_dir(arma::mat n, int K){
  int N = ( (K+1)*(K+2) ) / 2;
  arma::vec n_vec(N);
  int idx = 0 ;
  for(int c=0; c<(K+1); ++c){
    for(int r=0; r<(K+1-c); ++r){
      n_vec[idx] = n(r,c);
      idx += 1;
    }
  } 
  arma::vec wij(N);
  for(int i=0; i<N; ++i){
    wij[i] = R::rgamma(n_vec[i]+1,1);
  }
  wij = wij/accu(wij);
  arma::mat Omega(K+1,K+1);
  idx = 0;
  for(int c=0; c<(K+1); ++c){
    for(int r=0; r<(K+1-c); ++r){
      Omega(r,c) = wij[idx];
      idx += 1;
    }
  }   
  return Omega;
}

