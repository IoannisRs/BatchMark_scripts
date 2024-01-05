library(extraDistr)
library(Rcpp)
sourceCpp("BatchMarkRemovals.cpp")
K = 11
N = 6000
#entry conditional probabilities
Vj = c(0.13, 0.14, 0.16, 0.14, 0.06, 0.07, 0.04, 0.03, 
        0.03,  0.08, 0.10, 0.02)
#exit conditional probabilities
Vij = c(0.02, 0.11, 0.01, 0.22, 0.12, 0.02, 0.01, 0.01, 
        0.17, 0.16, 0.12, 0.03)
#join entry/exit probabilities
wij = grid_prob_fixed(Vj, Vij, K)
wij$Omega[which(wij$Omega==1)] = 0 
#latent number of entry/exit configurations
nij = matrix(0,K+1,K+1)
for(j in 1:(K+1)){
  for(i in 1:(K+2-j)){
    nij[i,j] = rpois(1,N*wij$Omega[i,j])
  }
}
#number of available individuals per sampling occasion
N_Avail_true_0.1 = c()
N_Avail_true_0.1 = IndAvail(nij,K)
#capture probabilities
p_true = c(0.51, 0.24, 0.21, 0.25, 0.15, 0.15, 0.05, 
           0.16, 0.23, 0.31, 0.13)
#removal probability
l_true = rep(0.1,K)
#Storage of parameters per replication 
{
Vij_0.1_with_list = list()
Vj_0.1_with_list = list()
p_0.1_with_list = list()
U_0.1_with_list = list()
N_0.1_with_list = list()
w_0.1_with_list = list()
Ome_0.1_with_list = list()

Vij_0.1_without_list = list()
Vj_0.1_without_list = list()
p_0.1_without_list = list()
U_0.1_without_list = list()
N_0.1_without_list = list()
w_0.1_without_list = list()
Ome_0.1_without_list = list()
}
#true number of available unmarked individuals 
U_Avail_true_0.1_list = list()

for(runs in 1:100){
  Data = Generate_mult_p(nij,  K,  l_true, p_true)
  Data_M = list()
  for(idx in 1:K){
    m_mat_0 = rmultinom(1,sum(Data$M[[idx]]),
                        rep(1/(length(1:(K+1-idx))*length(1:idx)),(length(1:(K+1-idx))*length(1:idx))))
    m_mat = matrix(0,K+1,K+1,byrow=T)
    t = 1
    for(c in 1:idx){
      for(r in 1:(K+1-idx)){
        m_mat[r,c] = m_mat_0[t] 
        t = t + 1
      }
    }
    Data_M[[idx]] = m_mat
  }
  Data_R = list()
  for(idx in 1:K){
    r_mat_0 = rmultinom(1,sum(Data$R[[idx]]),
                        rep(1/(length(1:(K+1-idx))*length(1:idx)),(length(1:(K+1-idx))*length(1:idx))))
    m_mat = matrix(0,K+1,K+1,byrow=T)
    r_mat = matrix(0,K+1,K+1,byrow=T)
    t = 1
    for(c in 1:idx){
      for(r in 1:(K+1-idx)){
        r_mat[r,c] = r_mat_0[t] 
        t = t + 1
      }
    }
    
    Data_R[[idx]] = r_mat
  }
  nij_init = matrix(0,K+1,K+1)
  for(j in 1:(K+1)){
    for(i in 1:(K+2-j)){
      nij_init[i,j] = rpois(1,4000*(1/((K+1)*(K+2)/2)))
    }
  }
  n_0 =  pmax( nij_init - Reduce('+',Data_R) - Reduce('+',Data_M),0)
  niter = 200000
  {
    #we saved 4000 samples
    indic = seq(2,niter,by=50)
    Sampled_M = list()
    for(i in 1:K){
      Sampled_M[[i]] =  array(0,c(K+1,K+1,length(indic))) 
    }
    for(i in 1:K){
      Sampled_M[[i]][,,1] = Data_M[[i]]
    }
    Sampled_R = list()
    for(i in 1:K){
      Sampled_R[[i]] =  array(0,c(K+1,K+1,length(indic))) 
    }
    for(i in 1:K){
      Sampled_R[[i]][,,1] = Data_R[[i]]
    }
    sampled_n_0 = array(0,c(K+1,K+1,length(indic)))
    sampled_n_0[,,1] = n_0
    sampled_n_0_trans = sampled_n_0[,,1]
    p = matrix(0,length(indic),K)
    p[1,] = (Data$r + apply(Data$Data,2,sum))/
       Available_mult_p(Data_M, Data_R, n_0,  K)
    p_trans = p[1,]
    w = c()
    w[1] = sum(n_0) + sum(Reduce('+',Data_R)) + sum(Reduce('+',Data_M))
    Ome_trans = grid_cell_dir(n_0 + Reduce('+',Data_R) + Reduce('+',Data_M),K)
    u= 1
    prop_m = rep(1,11)
    prop_r = rep(1,11)
    Omega_post = array(0,c(K+1,K+1,length(indic)))
  }
  for(iter in 2:niter){
    new_m = list() 
    old_m = list()
    new_r = list()
    old_r = list()
    if(iter != 2){
      for(idx in 1:K){
        new_r[[idx]] = Prop_Mat_cup( Res_1$n_output[[idx]],prop_r,K,idx-1)  
        while(  (sum(new_r[[idx]]<0) > 0) ){
          new_r[[idx]] = Prop_Mat_cup( Res_1$n_output[[idx]],prop_r,K,idx-1)    
        }
        new_m[[idx]] = Prop_Mat_cup( Res$n_output[[idx]],prop_m,K,idx-1)  
        while( (sum(t(IndAvel(new_m[[idx]],K)) >= Data$Data[idx,]) != K)  |   (sum(new_m[[idx]]<0) > 0) ){
          new_m[[idx]] = Prop_Mat_cup( Res$n_output[[idx]],prop_m,K,idx-1)  
        }
      }
      old_r = Res_1$n_output
      old_m = Res$n_output
    }else{
      for(idx in 1:K){
        new_r[[idx]] = Prop_Mat_cup( Sampled_R[[idx]][,,iter-1],prop_r,K,idx-1) 
        while(  (sum(new_r[[idx]]<0) > 0) ){
          new_r[[idx]] = Prop_Mat_cup( Sampled_R[[idx]][,,iter-1],prop_r,K,idx-1) 
        }
        old_r[[idx]] = Sampled_R[[idx]][,,iter-1]
        new_m[[idx]] = Prop_Mat_cup(Sampled_M[[idx]][,,iter-1],prop_m,K,idx-1) 
        while( (sum(t(IndAvel(new_m[[idx]],K)) >= Data$Data[idx,]) != K)  |   (sum(new_m[[idx]]<0) > 0)    ){
          new_m[[idx]] = Prop_Mat_cup(Sampled_M[[idx]][,,iter-1],prop_m,K,idx-1) 
        }
        old_m[[idx]] = Sampled_M[[idx]][,,iter-1]
      }
    }
    Res =  MH_caps_p(new_m, old_m,  old_r, sampled_n_0_trans, Ome_trans, Data$Unmarked,
                     Data$Data, K,  p_trans, 0.1)
    Res_1 = MH_loss_p(Res$n_output, new_r, old_r, sampled_n_0_trans, Ome_trans, 
                      Data$Unmarked,  K,  p_trans,0.1)
    new_n_miss = Prop_Mat(sampled_n_0_trans,m=1,K)
    sampled_n_0_trans_new  = MH_mis_p(Res$n_output, Res_1$n_output, Data$Unmarked,
                                      new_n_miss, sampled_n_0_trans, Ome_trans, K, w[iter-1], p_trans)
    sampled_n_0_trans = sampled_n_0_trans_new
    n = Reduce('+', Res$n_output) + Reduce('+', Res_1$n_output) +
      sampled_n_0_trans
    w[iter] = rgamma(1,sum(n)+0.001,1+0.001)
    Probs = grid_cell_dir(n,K)
    Ome_trans = Probs
    Nj = Avail_p(Res$n_output, Res_1$n_output, sampled_n_0_trans,  K)
    for(i in 1:K){
      p_trans[i] = rbeta(1,1 + Data$r[i] + sum(Data$Data[,i]),
                         1 + Nj[i] - ( Data$r[i] + sum(Data$Data[,i])) )
    }
    if(iter %in% indic){
      p[u,] = p_trans
      sampled_n_0[,,u] = sampled_n_0_trans
      for(idx in 1:K){
        Sampled_M[[idx]][,,u] = Res$n_output[[idx]]
        Sampled_R[[idx]][,,u] = Res_1$n_output[[idx]]
      }
      Omega_post[,,u] = Ome_trans
      u = u + 1
    }
  }
  M_with_0.1 = Sampled_M
  R_with_0.1 = Sampled_R
  n_0_with_0.1 = sampled_n_0
  p_with_0.1 = p
  abund_with_0.1 = w
  Ome_with_0.1 = Omega_post
  
  p_0.1_with_list[[runs]] = p_with_0.1
  w_0.1_with_list[[runs]] = abund_with_0.1
  Ome_0.1_with_list[[runs]] = Ome_with_0.1
  
  U_Avail_0.1 = matrix(0,length(indic),K)
  for(i in 1:length(indic)){
    U_Avail_0.1[i,] = Unmarked_Available_unif(list(Sampled_M[[1]][,,i],Sampled_M[[2]][,,i],Sampled_M[[3]][,,i],Sampled_M[[4]][,,i],Sampled_M[[5]][,,i],Sampled_M[[6]][,,i],
                                                   Sampled_M[[7]][,,i],Sampled_M[[8]][,,i],Sampled_M[[9]][,,i],Sampled_M[[10]][,,i],Sampled_M[[11]][,,i]
    ),list(Sampled_R[[1]][,,i],Sampled_R[[2]][,,i],Sampled_R[[3]][,,i],Sampled_R[[4]][,,i],Sampled_R[[5]][,,i],
           Sampled_R[[6]][,,i],Sampled_R[[7]][,,i],Sampled_R[[8]][,,i],Sampled_R[[9]][,,i],Sampled_R[[10]][,,i],
           Sampled_R[[11]][,,i]),sampled_n_0[,,i],K)
  }  
  
  U_0.1_with_list[[runs]] = U_Avail_0.1
  
  N_Avail_0.1 = matrix(0,length(indic),K)
  for(i in 1:length(indic)){
    N_Avail_0.1[i,] = IndAvail(Sampled_M[[1]][,,i]+Sampled_M[[2]][,,i]+Sampled_M[[3]][,,i]+
                                Sampled_M[[4]][,,i]+Sampled_M[[5]][,,i]+Sampled_M[[6]][,,i]+
                                Sampled_M[[7]][,,i]+Sampled_M[[8]][,,i]+Sampled_M[[9]][,,i]+
                                Sampled_M[[10]][,,i]+Sampled_M[[11]][,,i] +
                                Sampled_R[[1]][,,i]+Sampled_R[[2]][,,i]+Sampled_R[[3]][,,i]+Sampled_R[[4]][,,i]+Sampled_R[[5]][,,i]+
                                Sampled_R[[6]][,,i]+Sampled_R[[7]][,,i]+Sampled_R[[8]][,,i]+Sampled_R[[9]][,,i]+
                                Sampled_R[[10]][,,i]+
                                Sampled_R[[11]][,,i]+sampled_n_0[,,i],K)
  }
  N_0.1_with_list[[runs]] = N_Avail_0.1
  
  U_Avail_true_0.1 = rep(0,K)
  
  U_Avail_true_0.1 = Unmarked_Available_unif(Data$M,Data$R,Data$n_0,K)

  
  {
    Sampled_M = list()
    for(i in 1:K){
      Sampled_M[[i]] =  array(0,c(K+1,K+1,length(indic))) 
    }
    for(i in 1:K){
      Sampled_M[[i]][,,1] = Data_M[[i]]
    }
    Data_R = list()
    for(i in 1:K){
      Data_R[[i]] = matrix(0,K+1,K+1)
    }
    sampled_n_0 = array(0,c(K+1,K+1,length(indic)))
    sampled_n_0[,,1] = n_0
    p = matrix(0,length(indic),K)
    p[1,] = (apply(Data$Data,2,sum))/
       t(Available_mult_p(Data_M, Data_R, n_0,  K))
    w = c()
    w[1] = sum(n_0) + sum(Reduce('+',Data_M))
    Ome_trans = grid_cell_dir(n_0 + Reduce('+',Data_M),K)
    sampled_n_0_trans = sampled_n_0[,,1]
    p_trans = p[1,]
    u= 1
    Omega_post = array(0,c(K+1,K+1,length(indic)))
  }
  for(iter in 2:niter){
    new_m = list() 
    old_m = list()
    if(iter != 2){
      for(idx in 1:K){
        new_m[[idx]] = Prop_Mat_cup( Res$n_output[[idx]],prop_m,K,idx-1)  
        while( (sum(t(IndAvel(new_m[[idx]],K)) >= Data$Data[idx,]) != K)  |   (sum(new_m[[idx]]<0) > 0)    ){
          new_m[[idx]] = Prop_Mat_cup( Res$n_output[[idx]],prop_m,K,idx-1)  
        }
      }
      old_m = Res$n_output
    }else{
      for(idx in 1:K){
        new_m[[idx]] = Prop_Mat_cup(Sampled_M[[idx]][,,iter-1],prop_m,K,idx-1) 
        while( (sum(t(IndAvel(new_m[[idx]],K)) >= Data$Data[idx,]) != K)  |   (sum(new_m[[idx]]<0) > 0)    ){
          new_m[[idx]] = Prop_Mat_cup(Sampled_M[[idx]][,,iter-1],prop_m,K,idx-1) 
        }
        old_m[[idx]] = Sampled_M[[idx]][,,iter-1]
      }
    }
    Res =  MH_caps_p(new_m, old_m,  Data_R, sampled_n_0_trans, Ome_trans, diag(Data$Data),
                     Data$Data, K,  p_trans,1)
    new_n_miss = Prop_Mat(sampled_n_0_trans,m=1,K)
    sampled_n_0_trans_new  = MH_mis_p(Res$n_output, Data_R, diag(Data$Data),
                                      new_n_miss, sampled_n_0_trans, Ome_trans, K, w[iter-1], p_trans)
    sampled_n_0_trans = sampled_n_0_trans_new
    n = Reduce('+', Res$n_output)  +
      sampled_n_0_trans
    w[iter] = rgamma(1,sum(n)+0.001,1+0.001)
    Probs = grid_cell_dir(n,K)
    Ome_trans = Probs     
    Nj = Avail_p(Res$n_output, Data_R, sampled_n_0_trans,  K)
    for(i in 1:K){
      p_trans[i] = rbeta(1,1  + sum(Data$Data[,i]), 
                         1 + Nj[i] - ( sum(Data$Data[,i])) )
    }
    if(iter %in% indic){
      p[u,] = p_trans
      sampled_n_0[,,u] = sampled_n_0_trans
      for(idx in 1:K){
        Sampled_M[[idx]][,,u] = Res$n_output[[idx]]
      }
      Omega_post[,,u] = Ome_trans
      u = u + 1
    }
  }
  M_without_0.1 = Sampled_M
  n_0_without_0.1 = sampled_n_0
  p_without_0.1 = p
  abund_without_0.1 = w
  Ome_without_0.1 = Omega_post
  
  p_0.1_without_list[[runs]] = p_without_0.1
  w_0.1_without_list[[runs]] = abund_without_0.1
  Ome_0.1_without_list[[runs]] = Ome_without_0.1
  
  U_Avail_without_0.1 = matrix(0,length(indic),K)
  for(i in 1:length(indic)){
    U_Avail_without_0.1[i,] = Unmarked_Available_unif(list(Sampled_M[[1]][,,i],Sampled_M[[2]][,,i],Sampled_M[[3]][,,i],Sampled_M[[4]][,,i],Sampled_M[[5]][,,i],Sampled_M[[6]][,,i],
                                                           Sampled_M[[7]][,,i],Sampled_M[[8]][,,i],Sampled_M[[9]][,,i],Sampled_M[[10]][,,i],Sampled_M[[11]][,,i]
    ),Data_R,sampled_n_0[,,i],K)
  }  
  
  U_0.1_without_list[[runs]] = U_Avail_without_0.1
  
  N_Avail_0.1 = matrix(0,length(indic),K)
  for(i in 1:length(indic)){
    N_Avail_0.1[i,] = IndAvail(Sampled_M[[1]][,,i]+Sampled_M[[2]][,,i]+Sampled_M[[3]][,,i]+
                                Sampled_M[[4]][,,i]+Sampled_M[[5]][,,i]+Sampled_M[[6]][,,i]+
                                Sampled_M[[7]][,,i]+Sampled_M[[8]][,,i]+Sampled_M[[9]][,,i]+
                                Sampled_M[[10]][,,i]+Sampled_M[[11]][,,i] +sampled_n_0[,,i],K)
  }
  N_0.1_without_list[[runs]] = N_Avail_0.1
  U_Avail_true_0.1_list[[runs]] = U_Avail_true_0.1
  message(runs)
  
  saveRDS(w_0.1_with_list, file = paste("abund_0.1_with_list_unif_1",".RDS"))
  saveRDS(N_0.1_with_list, file = paste("N_0.1_with_list_unif_1",".RDS"))
  saveRDS(p_0.1_with_list, file = paste("p_0.1_with_list_unif_1",".RDS"))
  saveRDS(U_0.1_with_list, file = paste("U_0.1_with_list_unif_1",".RDS"))
  saveRDS(Ome_0.1_with_list, file = paste("Ome_0.1_with_list_unif_1",".RDS"))
  
  saveRDS(w_0.1_without_list, file = paste("abund_0.1_without_list_unif_1",".RDS"))
  saveRDS(N_0.1_without_list, file = paste("N_0.1_without_list_unif_1",".RDS"))
  saveRDS(p_0.1_without_list, file = paste("p_0.1_without_list_unif_1",".RDS"))
  saveRDS(U_0.1_without_list, file = paste("U_0.1_without_list_1",".RDS"))
  saveRDS(Ome_0.1_without_list, file = paste("Ome_0.1_without_list_1",".RDS"))
  
  saveRDS(U_Avail_true_0.1, file = paste("U_0.1_true_unif_1",".RDS"))
  saveRDS(IndAvel(nij,K), file = paste("N_0.1_true_unif_1",".RDS"))
  saveRDS(wij, file = paste("w_ij_0.1_true_unif_1",".RDS"))

}


