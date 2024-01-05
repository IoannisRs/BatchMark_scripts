//#include <Rcpp.h>
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
arma::vec Rep_Avail(arma::vec N, arma::vec TT){
  arma::vec Rep(sum(TT));
  int n_2 = size(N)[0];
  int u = 0;
  for(int i=0; i<n_2; ++i){
    for(int j=0; j<TT[i]; ++j){
      Rep[u] = N[i];
      u += 1;
    }
  }
  return Rep;
  
}



// [[Rcpp::export]]
List Generate_Unmarked(arma::mat nij, arma::vec T_K, arma::vec T_K_0,
                       arma::vec p, int K){
  arma::vec U(sum(T_K));
  int idx = 0;
  List Cup_List(K);
  arma::mat n_trans = nij;
  arma::mat u_trans(K+1,K+1);
  for(int k = 0; k<K; ++k){
    arma::cube m_mat(K+1,K+1,T_K[k]);
    Cup_List[k] = m_mat;
    for(int t=0; t<T_K[k]; ++t){
      u_trans.zeros();
      for(int j=0; j<(k+1); ++j){
        for(int i=0; i<(K-k); ++i){
          u_trans(i,j) = R::rbinom(n_trans(i,j),p[idx]);//
        }
      }
      m_mat.slice(t) = u_trans;
      U[idx] = accu(u_trans);
      n_trans -= u_trans;
      idx += 1;
    }
    Cup_List[k] = m_mat;
  }
  arma::mat n_0(K+1,K+1);
  n_0 = n_trans;
  return List::create(_["n_0"] = n_0, _["Unmarked"] = U,
                       _["Cup_List"]=Cup_List);
  
}



// [[Rcpp::export]]
List grid_prob_unif_unif_univariate_t(arma::mat n, int K){
  arma::mat Omega(K+1,K+1);
  Omega.ones();
  arma::vec L(K); 
  arma::vec V(K);
  arma::vec n_V_1(K);
  arma::vec n_V_2(K);
  arma::mat alpha_arr(K,2);
  alpha_arr(0,0) = 1;
  alpha_arr(0,1) = 1;
  alpha_arr(1,0) = 1;
  alpha_arr(1,1) = 1;
  alpha_arr(2,0) = 1;
  alpha_arr(2,1) = 1;
  alpha_arr(3,0) = 1;
  alpha_arr(3,1) = 1;
  alpha_arr(4,0) = 1;
  alpha_arr(4,1) = 1;
  alpha_arr(5,0) = 1;
  alpha_arr(5,1) = 1;
  
  for(int idx=0; idx<K; ++idx){
    for(int c=idx; c<(K+1); ++c){
      for(int r=0; r<(K+1-c); ++r){
        if(c==idx){
          n_V_1[idx] += n(r,c);
        }else{
          n_V_2[idx] += n(r,c);
        }
      }
    }
    V[idx] =  R::rbeta(n_V_1[idx]+alpha_arr(idx,0),n_V_2[idx]+alpha_arr(idx,1));
  }
  for(int idx=0; idx<K; ++idx){
    for(int c=idx; c<(K+1); ++c){
      for(int r=0; r<(K+1-c); ++r){
        if(c==idx){
          Omega(r,c) *= V[idx];
        }else{
          Omega(r,c) *= (1-V[idx]);
        }
      }
    }
  }
  arma::vec n_L_1(K);
  arma::vec n_L_2(K);
  arma::mat alpha_dep(K,2);
  alpha_dep(0,0) = 1;
  alpha_dep(0,1) = 1;
  alpha_dep(1,0) = 1;
  alpha_dep(1,1) = 1;
  alpha_dep(2,0) = 1;
  alpha_dep(2,1) = 1;
  alpha_dep(3,0) = 1;
  alpha_dep(3,1) = 1;
  alpha_dep(4,0) = 1;
  alpha_dep(4,1) = 1;
  alpha_dep(5,0) = 1;//alpha_dep(5,0) = 10;
  alpha_dep(5,1) = 1;//alpha_dep(5,1) = 10;
  for(int idx=0; idx<K; ++idx){
    for(int r=(K-idx); r>-1; --r){
      for(int c=0; c<(idx+1); ++c){
        if(r==(K-idx)){
          n_L_1[idx] += n(r,c);
        }else{
          n_L_2[idx] += n(r,c);
        }
      }
    }
    L[idx] = R::rbeta(n_L_1[idx]+alpha_dep(idx,0),n_L_2[idx]+alpha_dep(idx,1));
  }  
  for(int idx=0; idx<K; ++idx){
    for(int r=(K-idx); r>-1; --r){
      for(int c=0; c<(idx+1); ++c){
        if(r==(K-idx)){
          Omega(r,c) *= L[idx];
        }else{
          Omega(r,c) *= (1-L[idx]);
        }
      }
    }
  }
  return List::create(_["Omega"]=Omega,_["V"]=V,_["L"]=L);
}



// [[Rcpp::export]]
List grid_prob_unif_unif_univariate_t_fixed(arma::vec V, arma::vec L, int K){
  arma::mat Omega(K+1,K+1);
  Omega.ones();
  
  
  for(int idx=0; idx<K; ++idx){
    for(int c=idx; c<(K+1); ++c){
      for(int r=0; r<(K+1-c); ++r){
        if(c==idx){
          Omega(r,c) *= V[idx];
        }else{
          Omega(r,c) *= (1-V[idx]);
        }
      }
    }
  }
  
  for(int idx=0; idx<K; ++idx){
    for(int r=(K-idx); r>-1; --r){
      for(int c=0; c<(idx+1); ++c){
        if(r==(K-idx)){
          Omega(r,c) *= L[idx];
        }else{
          Omega(r,c) *= (1-L[idx]);
        }
      }
    }
  }    
  return List::create(_["Omega"]=Omega,_["V"]=V,_["L"]=L);
}





// [[Rcpp::export]]
arma::vec Unmarked_Avail_T_K(arma::vec T_K, List m, arma::mat n_0, int k, int K){
  int n = T_K[k];
  arma::vec N(n);
  for(int t=0; t<n; ++t){
    for(int tt=t; tt<n; ++tt){
      N[t] += accu((as<arma::cube>(m[k])).slice(tt));
    }
  }
  for(int kk=(k+1); kk<K; ++kk){
    n = T_K[kk];
    for(int t=0; t<n; ++t){
      N += IndAvail(as<arma::cube>(m[kk]).slice(t),K)[k];
    }
  }
  N += IndAvail(n_0,K)[k];
  return N;
}


// [[Rcpp::export]]
arma::vec Marked_Avail_T_K(arma::vec T_K, arma::vec T_K_0, int k,
                           arma::cube m, int K){
  arma::vec N(sum(T_K));
  arma::mat m_tot(K+1,K+1);
  arma::vec N_avail(sum(T_K));
  for(int t=(sum(T_K_0.subvec(0,k))); t<sum(T_K_0); ++t){
    if(t==(sum(T_K_0.subvec(0,k)))){
      N[t] = 0;
      m_tot += m.slice(t-(sum(T_K_0.subvec(0,k))));
    }else{
      N_avail = Rep_Avail(IndAvail(m_tot,K),T_K);
      N[t] = N_avail[t];
      if(t<(sum(T_K_0.subvec(0,k+1)))){
        m_tot += m.slice(t-(sum(T_K_0.subvec(0,k))));
      }
    }
  }
  return N;
}






// [[Rcpp::export]]
double Likelihood(List M, List U, List Cup_list, arma::mat n_0, arma::vec p, 
                  arma::vec T_K, arma::vec T_K_0, int K){
  double Like = 0;
  arma::vec N_marked(sum(T_K));
  arma::vec M_marked(sum(T_K));
  for(int k=0; k<K; ++k){
    arma::vec U_unmarked(T_K[k]);
    arma::vec N_unmarked(T_K[k]);
    N_marked = Marked_Avail_T_K(T_K,T_K_0,k,as<arma::cube>(Cup_list[k]),K);
    N_unmarked = Unmarked_Avail_T_K(T_K, Cup_list,  n_0,  k,  K);
    U_unmarked = as<arma::vec>(U[k]);
    M_marked = as<arma::vec>(M[k]);
    for(int t=(sum(T_K_0.subvec(0,k))); t<sum(T_K_0); ++t){
      if(t!=(sum(T_K_0.subvec(0,k)))){
        if(t<(sum(T_K_0.subvec(0,k+1)))){
          Like += R::dbinom(M_marked[t],N_marked[t],p[t],1);
          Like += R::dbinom(U_unmarked[t-(sum(T_K_0.subvec(0,k)))],N_unmarked[t-(sum(T_K_0.subvec(0,k)))],p[t],1);
        }else{
          Like += R::dbinom(M_marked[t],N_marked[t],p[t],1);
        }
      }else{
        Like += R::dbinom(U_unmarked[t-(sum(T_K_0.subvec(0,k)))],N_unmarked[t-(sum(T_K_0.subvec(0,k)))],p[t],1);
      }
      
    }
  }
  return Like;
}

// [[Rcpp::export]]
arma::vec Unmarked_Available(List m, arma::mat n_0, 
                             arma::vec T_K, int K){
  arma::vec Unmarked_Avail(sum(T_K));
  Unmarked_Avail = Rep_Avail(IndAvail(n_0,K),T_K);
  int u = 0;
  for(int k=0; k<K; ++k){
    for(int t=0; t<T_K[k]; ++t){
      Unmarked_Avail.subvec(0,u) += 
        Rep_Avail(IndAvail((as<arma::cube>(m[k])).slice(t),K),T_K).subvec(0,u);
      u += 1;
    }
  }
  return Unmarked_Avail;
}


// [[Rcpp::export]]
double Unmarked_Likelihood(arma::vec U, List m, arma::mat n_0, 
                           arma::vec T_K,arma::vec p, int K){
  arma::vec Unmarked_Avail(sum(T_K));
  Unmarked_Avail = Rep_Avail(IndAvail(n_0,K),T_K);
  int u = 0;
  for(int k=0; k<K; ++k){
    for(int t=0; t<T_K[k]; ++t){
      Unmarked_Avail.subvec(0,u) += 
        Rep_Avail(IndAvail((as<arma::cube>(m[k])).slice(t),K),T_K).subvec(0,u);
      u += 1;
    }
  }
  double Like = 0;
  for(int t=0; t<sum(T_K); ++t){
    Like += R::dbinom(U[t],Unmarked_Avail[t],p[t],1);
  }
  return Like;
}


// [[Rcpp::export]]
double Unmarked_Likelihood_p(arma::vec U, List m, arma::mat n_0, 
                             arma::vec T_K,double p, int K){
  arma::vec Unmarked_Avail(sum(T_K));
  Unmarked_Avail = Rep_Avail(IndAvail(n_0,K),T_K);
  int u = 0;
  for(int k=0; k<K; ++k){
    for(int t=0; t<T_K[k]; ++t){
      Unmarked_Avail.subvec(0,u) += 
        Rep_Avail(IndAvail((as<arma::cube>(m[k])).slice(t),K),T_K).subvec(0,u);
      u += 1;
    }
  }
  double Like = 0;
  for(int t=0; t<sum(T_K); ++t){
    Like += R::dbinom(U[t],Unmarked_Avail[t],p,1);
  }
  return Like;
}

// [[Rcpp::export]]
arma::vec Unmarked_Avail_p(arma::vec T_K, int K, arma::mat n_0, List m){
  arma::vec Unmarked_Avail(sum(T_K));
  Unmarked_Avail = Rep_Avail(IndAvail(n_0,K),T_K);
  int u = 0;
  for(int k=0; k<K; ++k){
    for(int t=0; t<T_K[k]; ++t){
      Unmarked_Avail.subvec(0,u) += 
        Rep_Avail(IndAvail((as<arma::cube>(m[k])).slice(t),K),T_K).subvec(0,u);
      u += 1;
    }
  }
  return Unmarked_Avail;
}



// [[Rcpp::export]]
double Marked_Likelihood(List M, List m, arma::vec T_K_0,  
                         arma::vec T_K,arma::vec p, int K){
  arma::vec Marked_Avail(sum(T_K));
  double Like = 0;
  int u = 1;
  for(int k=0; k<K; ++k){
    Marked_Avail =   Marked_Avail_T_K( T_K, T_K_0, k, as<arma::cube>(m[k]),  K);
    for(int t=u; t<(sum(T_K)); ++t){
      Like += R::dbinom(as<arma::vec>(M[k])[t],Marked_Avail[t],p[t],1);
    }
    u += 2;
  }
  return Like;
}

// [[Rcpp::export]]
double Marked_Likelihood_p(List M, List m, arma::vec T_K_0,  
                           arma::vec T_K,double p, int K){
  arma::vec Marked_Avail(sum(T_K));
  double Like = 0;
  int u = 1;
  for(int k=0; k<K; ++k){
    Marked_Avail =   Marked_Avail_T_K( T_K, T_K_0, k, as<arma::cube>(m[k]),  K);
    for(int t=u; t<(sum(T_K)); ++t){
      Like += R::dbinom(as<arma::vec>(M[k])[t],Marked_Avail[t],p,1);
    }
    u += 2;
  }
  return Like;
}

// [[Rcpp::export]]
arma::mat Prop_Mat_cup(arma::mat n, double m, int K, int k){
  
  Rcpp::Environment pkg = Rcpp::Environment::namespace_env("extraDistr");
  Rcpp::Function f = pkg["rdunif"];
  SEXP val = f(1,-m,m);
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
int PlusOrMinusOne() {
  return (rand() % 2) * 2 - 1;
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
arma::mat matrix2vector(arma::mat m, const bool byrow=false){
  if (byrow) {
    return m.as_row();
  } else {
    return m.as_col();
  }
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
arma::mat Prop_Mat_all(arma::mat n, int m, int K){
  
  Rcpp::Environment pkg = Rcpp::Environment::namespace_env("extraDistr");
  Rcpp::Function f = pkg["rdunif"];
  int l;
  for(int c=0; c<K; ++c){
    for(int r=0; r<(K-c); ++r){
      SEXP val = f(1,-m,m);
      l = Rcpp::as<int>(val);
      r = (rand()%(K));    
      c = (rand()%(K-r));
      n(r,c) = n(r,c) + (l*PlusOrMinusOne());
    }
  }
  return n;
}





// [[Rcpp::export]]
arma::vec Pois_Prior_K_sum_mat(arma::cube m_cup, arma::mat n_cand, arma::mat m_old, arma::mat O, int w, int K, arma::vec T_K){
  double log_cand=0;
  double log_old=0;
  arma::vec res(2);
  arma::mat n(K+1,K+1);
  n.zeros();
  n = n_cand;
  arma::mat m(K+1,K+1);
  m.zeros();
  m = m_old;
  for(int i = 0; i<K; ++i){
    n += m_cup.slice(i);
    m += m_cup.slice(i);   
  }
  for(int i=0; i<(K); ++i){
    for(int j=0; j<(K-i); ++j){
      log_cand += R::dpois(n(i,j),w*O(i,j),1);
      log_old +=  R::dpois(m(i,j),w*O(i,j),1);
    }
  }
  res[0] = log_cand;
  res[1] = log_old;
  return res;
}


// [[Rcpp::export]]
arma::vec BatchMark_UnmAvail( arma::vec Unmarked,
                              arma::cube m, arma::mat n_0, int K, arma::vec T_K,
                              arma::vec T_K_0){
  arma::vec Unmarked_Avail(sum(T_K));
  Unmarked_Avail = Rep_Avail(IndAvail(n_0,K),T_K);
  int u = -1;
  for(int k=0; k<K; ++k){
    u += T_K[k];
    Unmarked_Avail.subvec(0,u) +=
      Rep_Avail(IndAvail(m.slice(k),K),T_K).subvec(0,u);
  }
  for(int k=0; k<K; ++k){
    u = 0;
    for(int t=accu(T_K_0.subvec(0,k)); t<(accu(T_K_0.subvec(0,k))+T_K[k]); ++t){
      if(t!=accu(T_K_0.subvec(0,k))){
        Unmarked_Avail[t] = Unmarked_Avail[t]-
          accu(Unmarked.subvec(accu(T_K_0.subvec(0,k)),accu(T_K_0.subvec(0,k))+u)) ;
        u += 1;
      }
    }
  }
  return Unmarked_Avail;
}



// [[Rcpp::export]]
arma::mat BatchMark_MrkAvail(arma::vec Unmarked, arma::cube m, int K, arma::vec T_K){
  arma::mat M_avail(K,sum(T_K));
  arma::vec samp(sum(T_K));
  int u = 0;
  for(int k=0; k<K; ++k){
    samp = Rep_Avail(IndAvail(m.slice(k),K),T_K);
    for(int t=(accu(T_K.subvec(0,k))-1); t<sum(T_K); ++t){
      if(t!= (accu(T_K.subvec(0,k))-1)){
        M_avail(k,t) = samp[t];
      }else{
        M_avail(k,t) = accu(Unmarked.subvec(u,t-1));
      }  
    }
    u += T_K[k];
  }
  
  return M_avail;
}


// [[Rcpp::export]]
arma::vec BatchRecup_Avail(arma::vec Unmarked, arma::cube m,
                           int K, arma::vec T_K, arma::vec T_K_0, int k){

  
  arma::vec M_Avail(sum(T_K));
  M_Avail = Rep_Avail(IndAvail(m.slice(k),K),T_K);
  int u;
  u = 1 + accu(T_K_0.subvec(0,k));
  for(int t=0; t<u; ++t){
    M_Avail[t] = 0;
  }
  int idx = accu(T_K_0.subvec(0,k+1));
  for(int t=u; t<accu(T_K_0.subvec(0,k+1)); ++t){
    M_Avail[t] = accu(Unmarked.subvec(accu(T_K_0.subvec(0,k)), t-1));
  }
  
  return M_Avail;
}

// [[Rcpp::export]]
List Generate_Marked_correct(arma::vec Unmarked, arma::cube m ,
                             int K, arma::vec T_K, arma::vec T_K_0, arma::vec p){
  arma::vec M_avail(sum(T_K));
  List M(K);
  arma::vec m_cup(sum(T_K));
  int u = 1;
  for(int k=0; k<K; ++k){
    M_avail = BatchRecup_Avail(Unmarked, m,
                               K, T_K,  T_K_0, k);
    m_cup.zeros();
    for(int t=u; t<sum(T_K); ++t){
      m_cup[t] = R::rbinom(M_avail[t],p[t]);
    }
    u += T_K[k];
    M[k] = m_cup;
  }
  return List::create(_["M"] = M);
  
}

// [[Rcpp::export]]
double BatchMark_Likelihood(List M_recup, arma::vec Unmarked, arma::vec p,
                            arma::cube m, arma::mat n_0, int K, arma::vec T_K,
                            arma::vec T_K_0, int k){
  double Like = 0;
  arma::vec Unmarked_Avail(sum(T_K));
  Unmarked_Avail = BatchMark_UnmAvail( Unmarked,
                                       m, n_0, K,  T_K,
                                       T_K_0);
  
  for(int t=0; t<accu(T_K_0.subvec(0,k+1)); ++t){
    Like += R::dbinom(Unmarked[t],Unmarked_Avail[t],p[t],1);
  }
  
  
  arma::vec M_avail(sum(T_K));
  M_avail.zeros();
  //  //for(int k=0; k<K; ++k){
    M_avail = BatchRecup_Avail(Unmarked,m,K,T_K,T_K_0,k);
    for(int t=accu(T_K_0.subvec(0,k))+1; t<sum(T_K); ++t){
      Like += R::dbinom(as<arma::rowvec>(M_recup[k])[t],M_avail[t],p[t],1);
    }
    //}
  
  return Like;
}

// [[Rcpp::export]]
double Likelihood_Total(List M_recup, arma::vec Unmarked, arma::vec p, arma::cube m,
                        arma::mat n_0, int K, arma::vec T_K, arma::vec T_K_0){
  double Like=0;
  for(int k=0; k<K; ++k){
    Like += BatchMark_Likelihood( M_recup,  Unmarked,  p,
                                  m,  n_0,  K,  T_K,
                                  T_K_0,  k);
  }
  return Like;
}


