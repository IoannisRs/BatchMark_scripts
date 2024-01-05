library(RcppArmadillo)
library(extraDistr)
library(Rcpp)
sourceCpp("BatchMarkRD.cpp")
N = 5000
K = 6
p_true_true = c(0.21, 0.54, 0.48, 0.51, 0.37, 0.26)
p_true = p_true_true
V_true = c(0.41, 0.53, 0.49, 0.57, 0.31, 0.48)
L_true = c(0.34, 0.03, 0.28, 0.23, 0.09, 0.38)
wij = grid_prob_unif_unif_univariate_t_fixed(V_true,L_true,  K)


sim_probs = c()
for(k in 1:(K+1)){
  sim_probs = c(sim_probs,wij$Omega[1:(K-k+2),k])
}
n = N*sim_probs
nij = matrix(0,K+1,K+1)
u = 1
for(i in 1:(K+1)){
  for(j in 1:(K+2-i)){
    nij[j,i] = n[u]
    u = u + 1
  }
}
nij = round(nij)

T_K = c(16,16,16,16,16,16)
T_K_0 = c(0,16,16,16,16,16,16)

p_true = rep(p_true_true, each = 16)

p_runs = list()
w_runs = list()
V_runs = list()
W_runs = list()
for(runs in 1:100){
  Dt_U = Generate_Unmarked( nij,  T_K,  T_K_0,p_true,  K)
  M_recup = array(0,c(K+1,K+1,K))
  for(k in 1:K){
    for(t in 1:T_K[k]){
      M_recup[,,k] =M_recup[,,k] + Dt_U$Cup_List[[k]][,,t]
    }
  }
  
  Dt_M = Generate_Marked_correct(Dt_U$Unmarked, M_recup ,
                                 K,  T_K,  T_K_0, p_true)
  wij_init = grid_prob_unif_unif_univariate_t(diag(0,K+1,K+1),K)
  while(sum(wij_init$V<0.3 | wij_init$V >0.6)>0 | sum(wij_init$L>0.5)>0){
    wij_init = grid_prob_unif_unif_univariate_t(diag(0,K+1,K+1),K)
  }
  n_0_init = matrix(0,K+1,K+1)
  for(i in 1:(K+1)){
    for(j in 1:(K+2-i)){
      n_0_init[i,j] = rpois(1,(sum(Dt_U$Unmarked)*1.5)*wij_init$Omega[i,j]) 
    }
  }
  M_recup_init_trans = array(0,c(7,7,6))
  for(k in 1:K){
    M_recup_init_trans[1:(K+1-k),1:k,k] =  
      rmultinom(1,sum(Dt_U$Cup_List[[k]]),
                as.vector(wij_init$Omega[1:(K+1-k),1:k]))
  }
  
  group_size <- 16
  sums <- tapply(Dt_U$Unmarked, rep(seq_len(ceiling(length(Dt_U$Unmarked) / group_size)),
                                    each = group_size, length.out = length(Dt_U$Unmarked)), sum)
  
  while(mean(Rep_Avail(IndAvail(M_recup_init_trans[,,1],K),T_K) >= Dt_M$M[[1]])!=1 |
        mean(Rep_Avail(IndAvail(M_recup_init_trans[,,2],K),T_K) >= Dt_M$M[[2]])!=1 |
        mean(Rep_Avail(IndAvail(M_recup_init_trans[,,3],K),T_K) >= Dt_M$M[[3]])!=1 | 
        mean(Rep_Avail(IndAvail(M_recup_init_trans[,,4],K),T_K) >= Dt_M$M[[4]])!=1 |
        mean(Rep_Avail(IndAvail(M_recup_init_trans[,,5],K),T_K) >= Dt_M$M[[5]])!=1 |
        mean(Rep_Avail(IndAvail(M_recup_init_trans[,,6],K),T_K) >= Dt_M$M[[6]])!=1 |
        Check_Unmarked(M_recup_init_trans, n_0_init, M_recup_init_trans[,,1],  K,  1-1, sums)!=1 |
        Check_Unmarked(M_recup_init_trans, n_0_init, M_recup_init_trans[,,2],  K,  2-1, sums)!=1 |
        Check_Unmarked(M_recup_init_trans, n_0_init, M_recup_init_trans[,,3],  K,  3-1, sums)!=1 |
        Check_Unmarked(M_recup_init_trans, n_0_init, M_recup_init_trans[,,4],  K,  4-1, sums)!=1 |
        Check_Unmarked(M_recup_init_trans, n_0_init, M_recup_init_trans[,,5],  K,  5-1, sums)!=1 |
        Check_Unmarked(M_recup_init_trans, n_0_init, M_recup_init_trans[,,6],  K,  6-1, sums)!=1){
    wij_init = grid_prob_unif_unif_univariate_t(diag(0,K+1,K+1),K)
    while(sum(wij_init$V<0.3 | wij_init$V >0.6)>0 | sum(wij_init$L>0.5)>0){
      wij_init = grid_prob_unif_unif_univariate_t(diag(0,K+1,K+1),K)
    }
    n_0_init = matrix(0,K+1,K+1)
    for(i in 1:(K+1)){
      for(j in 1:(K+2-i)){
        n_0_init[i,j] = rpois(1,(sum(Dt_U$Unmarked)*1.5)*wij_init$Omega[i,j]) 
      }
    }
    M_recup_init_trans = array(0,c(7,7,6))
    for(k in 1:K){
      M_recup_init_trans[1:(K+1-k),1:k,k] =
        rmultinom(1,sum(Dt_U$Cup_List[[k]]),
                  as.vector(wij_init$Omega[1:(K+1-k),1:k]))
    }
    
  }
  M_recup_init = M_recup_init_trans
  CC = c()
  CC = Dt_M$M[[1]] + Dt_M$M[[2]] + Dt_M$M[[3]] + Dt_M$M[[4]] + Dt_M$M[[5]] + Dt_M$M[[6]]
  CC = CC + Dt_U$Unmarked
  niter = 500000
  indic = seq(2,niter,by=100)
  Sampled_M_iter = list()
  Sampled_M_iter[[1]] =  M_recup_init
  sampled_n_0 = array(0,c(K+1,K+1,niter))
  sampled_n_0[,,1] = n_0_init
  sampled_n_0_iter = array(0,c(K+1,K+1,niter))
  p_sampled = matrix(0,niter,sum(T_K))
  p_sampled[1,] = CC/Rep_Avail(IndAvail(apply(M_recup_init,c(1,2),sum) + n_0_init,K),T_K)
  w_sampled = rep(0,niter)
  w_sampled[1] = sum(n_0_init + apply(M_recup_init,c(1,2),sum))
  Ome_trans = wij_init$Omega
  sampled_n_0_trans = sampled_n_0[,,1]
  p_trans = p_sampled[1,]
  Omega_sampled = array(0,c(K+1,K+1,niter))
  Ome_trans[which(Ome_trans==1,arr.ind = TRUE)] = 0
  Omega_sampled[,,1] = Ome_trans
  a_all = rep(0,niter)
  prop_step = c(25,5,5,5,5,5)
  neg_ind = c()
  V_sampled = matrix(0,niter,K)
  V_sampled[1,] = wij_init$V
  W_sampled = matrix(0,niter,K)
  W_sampled[1,] = wij_init$L
  p_store = matrix(0,length(indic),sum(T_K))
  V_store = matrix(0,length(indic),K)
  W_store = matrix(0,length(indic),K)
  w_store = c()
  store_idx = 1
  accept = matrix(0,niter,K)
  accept_0 = rep(0,niter)
  accept_cup = rep(0,niter)
  for(iter in 2:niter){
    neg_ind = c()
    new_n = list()
    idx_idx = 1
    if(iter != 2){
      new_n_trans = array(0,c(K+1,K+1,K))
      for(idc in 1:K){
        new_n_trans[,,idc] = Prop_Mat_cup(Sampled_M_iter[[iter-1]][,,idc],prop_step[idx_idx],K,idc-1)
        neg_ind[idx_idx] = ifelse(sum(new_n_trans[,,idc]<0)==0,1,0)
        idx_idx = idx_idx + 1
      }
    }else{
      new_n_trans = array(0,c(K+1,K+1,K))
      for(idc in 1:K){
        new_n_trans[,,idc] = Prop_Mat_cup(Sampled_M_iter[[iter-1]][,,idc],prop_step[idx_idx],K,idc-1)
        neg_ind[idx_idx] = ifelse(sum(new_n_trans[,,idc]<0)==0,1,0)
        idx_idx = idx_idx + 1
      }
    }
    
    Sampled_M_iter[[iter]] = Sampled_M_iter[[iter-1]]
    M_trans = Sampled_M_iter[[iter-1]]
    for(k in 1:K){
      if(neg_ind[k]==1 ){
        M_trans[,,k] = new_n_trans[,,k]
        num = 0; den=0
        num = num +
          BatchMark_Likelihood(Dt_M$M,Dt_U$Unmarked,p_sampled[iter-1,],M_trans,
                               sampled_n_0_trans,K,T_K,T_K_0,k-1)+
          dmultinom(as.vector((apply(Sampled_M_iter[[iter]][,,-k],c(1,2),sum)+M_trans[,,k]+sampled_n_0_trans)[1:(K-k+1),1:k]),
                    sum(as.vector((apply(Sampled_M_iter[[iter]][,,-k],c(1,2),sum)+M_trans[,,k]+sampled_n_0_trans)[1:(K-k+1),1:k])),
                    as.vector(Omega_sampled[1:(K-k+1),1:k,iter-1]),log=TRUE)
        den = den +
          BatchMark_Likelihood(Dt_M$M,Dt_U$Unmarked,p_sampled[iter-1,],Sampled_M_iter[[iter]],
                               sampled_n_0_trans,K,T_K,T_K_0,k-1)+dmultinom(as.vector((apply(Sampled_M_iter[[iter]],c(1,2),sum)+sampled_n_0_trans)[1:(K-k+1),1:k]),
                                                                            sum(as.vector((apply(Sampled_M_iter[[iter]],c(1,2),sum)+sampled_n_0_trans)[1:(K-k+1),1:k])),
                                                                            as.vector(Omega_sampled[1:(K-k+1),1:k,iter-1]),log=TRUE)
        if(runif(1,0,1) <= exp(num-den) ){
          Sampled_M_iter[[iter]][,,k]= new_n_trans[,,k]
          accept[iter,k] = 1
        }else{
          Sampled_M_iter[[iter]][,,k] = Sampled_M_iter[[iter-1]][,,k]
          M_trans[,,k] = Sampled_M_iter[[iter-1]][,,k]
        }
      }else{
        Sampled_M_iter[[iter]][,,k] = Sampled_M_iter[[iter-1]][,,k]
      }
    }
    sampled_n_0_trans_new = Prop_Mat(sampled_n_0_trans,1,K)
    if(mean(sampled_n_0_trans_new>=0)==1 ){
      pois_prior = Pois_Prior_K_sum_mat(Sampled_M_iter[[iter]],sampled_n_0_trans_new,sampled_n_0_trans,
                                        Omega_sampled[,,iter-1],w_sampled[iter-1],K,T_K)
      num = 
        Likelihood_Total(Dt_M$M, Dt_U$Unmarked,p_sampled[iter-1,],
                         Sampled_M_iter[[iter]],
                         sampled_n_0_trans_new,K,T_K,T_K_0) +
        pois_prior[1]
      den = 
        Likelihood_Total(Dt_M$M, Dt_U$Unmarked,p_sampled[iter-1,],
                         Sampled_M_iter[[iter]],
                         sampled_n_0_trans,K,T_K,T_K_0) +
        pois_prior[2]
      if(runif(1,0,1) <= exp(num-den) ){
        sampled_n_0_trans = sampled_n_0_trans_new
        accept_0[iter] = 1
      }
    }
    for(c in 1:(K+1)){
      for(r in 1:(K+1)){
        if(c==(K-r+2)){
          sampled_n_0_trans[r,c] = rpois(1,w_sampled[iter-1]*Omega_sampled[r,c,iter-1])
        }
      }
    }
    sampled_n_0[,,iter] = sampled_n_0_trans
    n_tot = sampled_n_0_trans + apply(Sampled_M_iter[[iter]],c(1,2),sum)
    P = grid_prob_unif_unif_univariate_t(n_tot,K)
    P$Omega[which(P$Omega==1,arr.ind = TRUE)] = 0
    Omega_sampled[,,iter] = P$Omega
    W_sampled[iter,] = P$L
    V_sampled[iter,] = P$V
    
    w_sampled[iter] = rgamma(1,sum(n_tot)+0.001,1+0.001)
    NN = Rep_Avail(IndAvail(apply(Sampled_M_iter[[iter]],c(1,2),sum) + sampled_n_0_trans,K),T_K)
    for(t in 1:sum(T_K)){
      p_sampled[iter,t] = rbeta(1,1+CC[t],1+NN[t]-CC[t])
    }
    if(iter %in% indic){
      p_store[store_idx,] = p_sampled[iter,]
      V_store[store_idx,] = V_sampled[iter,]
      W_store[store_idx,] = W_sampled[iter,]
      w_store[store_idx] = w_sampled[iter]
      store_idx = store_idx + 1
    }
  }
  message(runs)
  p_runs[[runs]] = p_store
  w_runs[[runs]] = w_store
  V_runs[[runs]] = V_store
  W_runs[[runs]] = W_store
  saveRDS(p_runs, file = paste("p_runs_BM_sim_alt_unif_bias_16_so",".RDS"))
  saveRDS(w_runs, file = paste("abund_runs_BM_sim_alt_unif_bias_16_so",".RDS")) 
  saveRDS(V_runs, file = paste("V_runs_alt_unif_bias_16_so",".RDS")) 
  saveRDS(W_runs, file = paste("W_runs_alt_unif_bias_16_so",".RDS")) 
  
}

