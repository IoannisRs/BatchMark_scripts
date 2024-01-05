library(gridExtra)
library(ggplot2)
p_with_0.05 = `p_0.05_with_list_unif_1 `
p_without_0.05 = `p_0.05_without_list_unif_1 `
p_with_0.1 = `p_0.1_with_list_unif_1 `
p_without_0.1 = `p_0.1_without_list_unif_1 `
p_with_0.2 = `p_0.2_with_list_unif_1 `
p_without_0.2 = `p_0.2_without_list_unif_1 `


p_true = c(0.51269897, 0.24317193, 0.21823916, 0.25740672, 0.15595798, 0.15334564,
           0.05551861, 0.16449608, 0.23604110, 0.31122140, 0.13177267)

p_PTL_0.05 = matrix(0,100,11)
for(i in 1:100){
  p_PTL_0.05[i,] =  apply(p_with_0.05[[i]][1000:4000,],2,mean) - p_true
}
indc_PTL = list()
for(i in 1:11){
  indc_PTL[[i]] = which(abs(p_PTL_0.05[,i])<1)
}
p_PTWL_0.05 = matrix(0,100,11)
for(i in 1:100){
  p_PTWL_0.05[i,] = apply(p_without_0.05[[i]][1:4000,],2,mean) - p_true
}


####
####
####
p_PTL_0.1 = matrix(0,100,11)
for(i in 1:100){
  p_PTL_0.1[i,] = apply(p_with_0.1[[i]][1000:4000,],2,mean) - p_true
}
indc_PTL = list()
for(i in 1:11){
  indc_PTL[[i]] = which(abs(p_PTL_0.1[,i])<1)
}


p_PTWL_0.1 = matrix(0,100,11)
for(i in 1:100){
  p_PTWL_0.1[i,] = apply(p_without_0.1[[i]][1:4000,],2,mean) - p_true
}


####
####
####
p_PTL_0.2 =  matrix(0,100,11)
for(i in 1:100){
  p_PTL_0.2[i,] = apply(p_with_0.2[[i]][1000:4000,],2,mean) - p_true
}
indc_PTL = list()
for(i in 1:11){
  indc_PTL[[i]] = which(abs(p_PTL_0.2[,i])<1)
}
p_PTWL_0.2 = matrix(0,100,11)
for(i in 1:100){
  p_PTWL_0.2[i,] = apply(p_without_0.2[[i]][1:4000,],2,mean) - p_true
}





df_p = data.frame(x = c(as.vector(p_PTL_0.05[indc_PTL[[1]],1]), as.vector(p_PTL_0.05[indc_PTL[[2]],2]),
                        as.vector(p_PTL_0.05[indc_PTL[[3]],3]), as.vector(p_PTL_0.05[indc_PTL[[4]],4]),
                        as.vector(p_PTL_0.05[indc_PTL[[5]],5]), as.vector(p_PTL_0.05[indc_PTL[[6]],6]),
                        as.vector(p_PTL_0.05[indc_PTL[[7]],7]), as.vector(p_PTL_0.05[indc_PTL[[8]],8]),
                        as.vector(p_PTL_0.05[indc_PTL[[9]],9]), as.vector(p_PTL_0.05[indc_PTL[[10]],10]),
                        as.vector(p_PTL_0.05[indc_PTL[[11]],11]),
                        as.vector(p_PTWL_0.05[indc_PTL[[1]],1]), as.vector(p_PTWL_0.05[indc_PTL[[2]],2]),
                        as.vector(p_PTWL_0.05[indc_PTL[[3]],3]), as.vector(p_PTWL_0.05[indc_PTL[[4]],4]),
                        as.vector(p_PTWL_0.05[indc_PTL[[5]],5]), as.vector(p_PTWL_0.05[indc_PTL[[6]],6]),
                        as.vector(p_PTWL_0.05[indc_PTL[[7]],7]), as.vector(p_PTWL_0.05[indc_PTL[[8]],8]),
                        as.vector(p_PTWL_0.05[indc_PTL[[9]],9]), as.vector(p_PTWL_0.05[indc_PTL[[10]],10]),
                        as.vector(p_PTWL_0.05[indc_PTL[[11]],11])), 
                  y = factor(c(rep(1,length(indc_PTL[[1]])), rep(2,length(indc_PTL[[2]])), 
                               rep(3,length(indc_PTL[[3]])), rep(4,length(indc_PTL[[4]])), 
                               rep(5,length(indc_PTL[[5]])), rep(6,length(indc_PTL[[6]])),
                               rep(7,length(indc_PTL[[7]])), rep(8,length(indc_PTL[[8]])), 
                               rep(9,length(indc_PTL[[9]])), rep(10,length(indc_PTL[[10]])), 
                               rep(11,length(indc_PTL[[11]])),
                               rep(1,length(indc_PTL[[1]])), rep(2,length(indc_PTL[[2]])), 
                               rep(3,length(indc_PTL[[3]])), rep(4,length(indc_PTL[[4]])), 
                               rep(5,length(indc_PTL[[5]])), rep(6,length(indc_PTL[[6]])),
                               rep(7,length(indc_PTL[[7]])), rep(8,length(indc_PTL[[8]])), 
                               rep(9,length(indc_PTL[[9]])), rep(10,length(indc_PTL[[10]])), 
                               rep(11,length(indc_PTL[[11]])))),
                  model = c(rep("PTBM-R",length(indc_PTL[[11]])*11),rep("PTBM",length(indc_PTL[[11]])*11)))  
colnames(df_p) = c("p", "so", "model") 
p_1 = ggplot(df_p, aes(x=so, y=p, fill= model)) + 
  geom_boxplot() +
  scale_fill_manual(values = c("#555555", "grey")) +
  theme(legend.position = "none") + 
  labs(title = "Removal probability 0.05", x="sampling occasions")+
  coord_cartesian(ylim = c(-0.05, 0.05)) +theme(axis.text = element_text(size = 13),
                                                legend.text = element_text(size = 13),
                                                legend.title = element_text(size = 13),
                                                axis.title.x = element_text(size = 13)) 


df_p = data.frame(x = c(as.vector(p_PTL_0.1[indc_PTL[[1]],1]), as.vector(p_PTL_0.1[indc_PTL[[2]],2]),
                        as.vector(p_PTL_0.1[indc_PTL[[3]],3]), as.vector(p_PTL_0.1[indc_PTL[[4]],4]),
                        as.vector(p_PTL_0.1[indc_PTL[[5]],5]), as.vector(p_PTL_0.1[indc_PTL[[6]],6]),
                        as.vector(p_PTL_0.1[indc_PTL[[7]],7]), as.vector(p_PTL_0.1[indc_PTL[[8]],8]),
                        as.vector(p_PTL_0.1[indc_PTL[[9]],9]), as.vector(p_PTL_0.1[indc_PTL[[10]],10]),
                        as.vector(p_PTL_0.1[indc_PTL[[11]],11]),
                        as.vector(p_PTWL_0.1[indc_PTL[[1]],1]), as.vector(p_PTWL_0.1[indc_PTL[[2]],2]),
                        as.vector(p_PTWL_0.1[indc_PTL[[3]],3]), as.vector(p_PTWL_0.1[indc_PTL[[4]],4]),
                        as.vector(p_PTWL_0.1[indc_PTL[[5]],5]), as.vector(p_PTWL_0.1[indc_PTL[[6]],6]),
                        as.vector(p_PTWL_0.1[indc_PTL[[7]],7]), as.vector(p_PTWL_0.1[indc_PTL[[8]],8]),
                        as.vector(p_PTWL_0.1[indc_PTL[[9]],9]), as.vector(p_PTWL_0.1[indc_PTL[[10]],10]),
                        as.vector(p_PTWL_0.1[indc_PTL[[11]],11])), 
                  y = factor(c(rep(1,length(indc_PTL[[1]])), rep(2,length(indc_PTL[[2]])), 
                               rep(3,length(indc_PTL[[3]])), rep(4,length(indc_PTL[[4]])), 
                               rep(5,length(indc_PTL[[5]])), rep(6,length(indc_PTL[[6]])),
                               rep(7,length(indc_PTL[[7]])), rep(8,length(indc_PTL[[8]])), 
                               rep(9,length(indc_PTL[[9]])), rep(10,length(indc_PTL[[10]])), 
                               rep(11,length(indc_PTL[[11]])),
                               rep(1,length(indc_PTL[[1]])), rep(2,length(indc_PTL[[2]])), 
                               rep(3,length(indc_PTL[[3]])), rep(4,length(indc_PTL[[4]])), 
                               rep(5,length(indc_PTL[[5]])), rep(6,length(indc_PTL[[6]])),
                               rep(7,length(indc_PTL[[7]])), rep(8,length(indc_PTL[[8]])), 
                               rep(9,length(indc_PTL[[9]])), rep(10,length(indc_PTL[[10]])), 
                               rep(11,length(indc_PTL[[11]])))),
                  model = c(rep("PTBM-R",length(indc_PTL[[11]])*11),rep("PTBM",length(indc_PTL[[11]])*11)))  
colnames(df_p) = c("p", "so", "model") 
p_2 = ggplot(df_p, aes(x=so, y=p, fill= model)) + 
  geom_boxplot() +
  scale_fill_manual(values = c("#555555", "grey")) +
  theme(legend.position = "none") + 
  labs(title = "Removal probability 0.1", x="sampling occasions",y=NULL) +
  theme(legend.justification = c(1, 0), 
        legend.position = c(0.4, 0.83))+
  coord_cartesian(ylim = c(-0.05, 0.05))+theme(axis.text = element_text(size = 13),
                                               legend.text = element_text(size = 13),
                                               legend.title = element_text(size = 13),
                                               axis.title.x = element_text(size = 13)) 

df_p = data.frame(x = c(as.vector(p_PTL_0.2[indc_PTL[[1]],1]), as.vector(p_PTL_0.2[indc_PTL[[2]],2]),
                        as.vector(p_PTL_0.2[indc_PTL[[3]],3]), as.vector(p_PTL_0.2[indc_PTL[[4]],4]),
                        as.vector(p_PTL_0.2[indc_PTL[[5]],5]), as.vector(p_PTL_0.2[indc_PTL[[6]],6]),
                        as.vector(p_PTL_0.2[indc_PTL[[7]],7]), as.vector(p_PTL_0.2[indc_PTL[[8]],8]),
                        as.vector(p_PTL_0.2[indc_PTL[[9]],9]), as.vector(p_PTL_0.2[indc_PTL[[10]],10]),
                        as.vector(p_PTL_0.2[indc_PTL[[11]],11]),
                        as.vector(p_PTWL_0.2[indc_PTL[[1]],1]), as.vector(p_PTWL_0.2[indc_PTL[[2]],2]),
                        as.vector(p_PTWL_0.2[indc_PTL[[3]],3]), as.vector(p_PTWL_0.2[indc_PTL[[4]],4]),
                        as.vector(p_PTWL_0.2[indc_PTL[[5]],5]), as.vector(p_PTWL_0.2[indc_PTL[[6]],6]),
                        as.vector(p_PTWL_0.2[indc_PTL[[7]],7]), as.vector(p_PTWL_0.2[indc_PTL[[8]],8]),
                        as.vector(p_PTWL_0.2[indc_PTL[[9]],9]), as.vector(p_PTWL_0.2[indc_PTL[[10]],10]),
                        as.vector(p_PTWL_0.2[indc_PTL[[11]],11])), 
                  y = factor(c(rep(1,length(indc_PTL[[1]])), rep(2,length(indc_PTL[[2]])), 
                               rep(3,length(indc_PTL[[3]])), rep(4,length(indc_PTL[[4]])), 
                               rep(5,length(indc_PTL[[5]])), rep(6,length(indc_PTL[[6]])),
                               rep(7,length(indc_PTL[[7]])), rep(8,length(indc_PTL[[8]])), 
                               rep(9,length(indc_PTL[[9]])), rep(10,length(indc_PTL[[10]])), 
                               rep(11,length(indc_PTL[[11]])),
                               rep(1,length(indc_PTL[[1]])), rep(2,length(indc_PTL[[2]])), 
                               rep(3,length(indc_PTL[[3]])), rep(4,length(indc_PTL[[4]])), 
                               rep(5,length(indc_PTL[[5]])), rep(6,length(indc_PTL[[6]])),
                               rep(7,length(indc_PTL[[7]])), rep(8,length(indc_PTL[[8]])), 
                               rep(9,length(indc_PTL[[9]])), rep(10,length(indc_PTL[[10]])), 
                               rep(11,length(indc_PTL[[11]])))),
                  model = c(rep("PTBM-R",length(indc_PTL[[11]])*11),rep("PTBM",length(indc_PTL[[11]])*11)))  
colnames(df_p) = c("p", "so", "model") 
p_3 = ggplot(df_p, aes(x=so, y=p, fill= model)) + 
  geom_boxplot() +
  scale_fill_manual(values = c("#555555", "grey")) +
  theme(legend.position = "none") + 
  labs(title = "Removal probability 0.2", x="sampling occasions",y=NULL)+
  coord_cartesian(ylim = c(-0.05, 0.05))+theme(axis.text = element_text(size = 13),
                                               legend.text = element_text(size = 13),
                                               legend.title = element_text(size = 13),
                                               axis.title.x = element_text(size = 13)) 

grid.arrange(p_1, p_2, p_3, ncol = 3)


####################################
####################################
####################################
####################################

U_with_0.05 = `U_0.05_with_list_unif_1 `
U_without_0.05 = `U_0.05_without_list_1 `
U_with_0.1 = `U_0.1_with_list_unif_1 `
U_without_0.1 = `U_0.1_without_list_1 `
U_with_0.2 = `U_0.2_with_list_unif_1 `
U_without_0.2 = `U_0.2_without_list_1 `

U_0.05_true = `U_0.05_true_unif_1 `
U_PTL_0.05 = matrix(0,100,11)
for(i in 1:100){
  U_PTL_0.05[i,] =  apply(U_with_0.05[[i]][1:4000,],2,mean) - U_0.05_true
}
indc_PTL = list()
for(i in 1:11){
  indc_PTL[[i]] = which(abs(U_PTL_0.05[,i])<10000)
}

U_PTWL_0.05 = matrix(0,100,11)
for(i in 1:100){
  U_PTWL_0.05[i,] = apply(U_without_0.05[[i]][1:4000,],2,mean) - U_0.05_true
}


####
####
####

U_0.1_true = `U_0.1_true_unif_1 `

U_PTL_0.1 = matrix(0,100,11)
for(i in 1:100){
  U_PTL_0.1[i,] = apply(U_with_0.1[[i]][1:4000,],2,mean) - U_0.1_true
}
indc_PTL = list()
for(i in 1:11){
  indc_PTL[[i]] = which(abs(U_PTL_0.1[,i])<10000)
}

U_PTWL_0.1 = matrix(0,100,11)
for(i in 1:100){
  U_PTWL_0.1[i,] = apply(U_without_0.1[[i]][1:4000,],2,mean) - U_0.1_true
}


####
####
####

U_0.2_true = `U_0.2_true_unif_1 `
U_PTL_0.2 =  matrix(0,100,11)
for(i in 1:100){
  U_PTL_0.2[i,] = apply(U_with_0.2[[i]][1:4000,],2,mean) - U_0.2_true
}
indc_PTL = list()
for(i in 1:11){
  indc_PTL[[i]] = which(abs(U_PTL_0.2[,i])<10000)
}

U_PTWL_0.2 = matrix(0,100,11)
for(i in 1:100){
  U_PTWL_0.2[i,] = apply(U_without_0.2[[i]][1:4000,],2,mean) - U_0.2_true
}







df_p = data.frame(x = c(as.vector(U_PTL_0.05[indc_PTL[[1]],1]), as.vector(U_PTL_0.05[indc_PTL[[2]],2]),
                        as.vector(U_PTL_0.05[indc_PTL[[3]],3]), as.vector(U_PTL_0.05[indc_PTL[[4]],4]),
                        as.vector(U_PTL_0.05[indc_PTL[[5]],5]), as.vector(U_PTL_0.05[indc_PTL[[6]],6]),
                        as.vector(U_PTL_0.05[indc_PTL[[7]],7]), as.vector(U_PTL_0.05[indc_PTL[[8]],8]),
                        as.vector(U_PTL_0.05[indc_PTL[[9]],9]), as.vector(U_PTL_0.05[indc_PTL[[10]],10]),
                        as.vector(U_PTL_0.05[indc_PTL[[11]],11]),
                        as.vector(U_PTWL_0.05[indc_PTL[[1]],1]), as.vector(U_PTWL_0.05[indc_PTL[[2]],2]),
                        as.vector(U_PTWL_0.05[indc_PTL[[3]],3]), as.vector(U_PTWL_0.05[indc_PTL[[4]],4]),
                        as.vector(U_PTWL_0.05[indc_PTL[[5]],5]), as.vector(U_PTWL_0.05[indc_PTL[[6]],6]),
                        as.vector(U_PTWL_0.05[indc_PTL[[7]],7]), as.vector(U_PTWL_0.05[indc_PTL[[8]],8]),
                        as.vector(U_PTWL_0.05[indc_PTL[[9]],9]), as.vector(U_PTWL_0.05[indc_PTL[[10]],10]),
                        as.vector(U_PTWL_0.05[indc_PTL[[11]],11])), 
                  y = factor(c(rep(1,length(indc_PTL[[1]])), rep(2,length(indc_PTL[[2]])), 
                               rep(3,length(indc_PTL[[3]])), rep(4,length(indc_PTL[[4]])), 
                               rep(5,length(indc_PTL[[5]])), rep(6,length(indc_PTL[[6]])),
                               rep(7,length(indc_PTL[[7]])), rep(8,length(indc_PTL[[8]])), 
                               rep(9,length(indc_PTL[[9]])), rep(10,length(indc_PTL[[10]])), 
                               rep(11,length(indc_PTL[[11]])),
                               rep(1,length(indc_PTL[[1]])), rep(2,length(indc_PTL[[2]])), 
                               rep(3,length(indc_PTL[[3]])), rep(4,length(indc_PTL[[4]])), 
                               rep(5,length(indc_PTL[[5]])), rep(6,length(indc_PTL[[6]])),
                               rep(7,length(indc_PTL[[7]])), rep(8,length(indc_PTL[[8]])), 
                               rep(9,length(indc_PTL[[9]])), rep(10,length(indc_PTL[[10]])), 
                               rep(11,length(indc_PTL[[11]])))),
                  model = c(rep("PTBM-R",length(indc_PTL[[11]])*11),rep("PTBM",length(indc_PTL[[11]])*11))) 
colnames(df_p) = c("p", "so", "model") 
library(ggplot2)
p_1 = ggplot(df_p, aes(x=so, y=p, fill= model)) + 
  geom_boxplot() +
  scale_fill_manual(values = c("#555555", "grey")) +
  theme(legend.position = "none") + 
  labs(title = "Removal probability 0.05", x="sampling occasions", y = "U")+
  coord_cartesian(ylim = c(-470, 200))+theme(axis.text = element_text(size = 13),
                                             legend.text = element_text(size = 13),
                                             legend.title = element_text(size = 13),
                                             axis.title.x = element_text(size = 13))


df_p = data.frame(x = c(as.vector(U_PTL_0.1[indc_PTL[[1]],1]), as.vector(U_PTL_0.1[indc_PTL[[2]],2]),
                        as.vector(U_PTL_0.1[indc_PTL[[3]],3]), as.vector(U_PTL_0.1[indc_PTL[[4]],4]),
                        as.vector(U_PTL_0.1[indc_PTL[[5]],5]), as.vector(U_PTL_0.1[indc_PTL[[6]],6]),
                        as.vector(U_PTL_0.1[indc_PTL[[7]],7]), as.vector(U_PTL_0.1[indc_PTL[[8]],8]),
                        as.vector(U_PTL_0.1[indc_PTL[[9]],9]), as.vector(U_PTL_0.1[indc_PTL[[10]],10]),
                        as.vector(U_PTL_0.1[indc_PTL[[11]],11]),
                        as.vector(U_PTWL_0.1[indc_PTL[[1]],1]), as.vector(U_PTWL_0.1[indc_PTL[[2]],2]),
                        as.vector(U_PTWL_0.1[indc_PTL[[3]],3]), as.vector(U_PTWL_0.1[indc_PTL[[4]],4]),
                        as.vector(U_PTWL_0.1[indc_PTL[[5]],5]), as.vector(U_PTWL_0.1[indc_PTL[[6]],6]),
                        as.vector(U_PTWL_0.1[indc_PTL[[7]],7]), as.vector(U_PTWL_0.1[indc_PTL[[8]],8]),
                        as.vector(U_PTWL_0.1[indc_PTL[[9]],9]), as.vector(U_PTWL_0.1[indc_PTL[[10]],10]),
                        as.vector(U_PTWL_0.1[indc_PTL[[11]],11])), 
                  y = factor(c(rep(1,length(indc_PTL[[1]])), rep(2,length(indc_PTL[[2]])), 
                               rep(3,length(indc_PTL[[3]])), rep(4,length(indc_PTL[[4]])), 
                               rep(5,length(indc_PTL[[5]])), rep(6,length(indc_PTL[[6]])),
                               rep(7,length(indc_PTL[[7]])), rep(8,length(indc_PTL[[8]])), 
                               rep(9,length(indc_PTL[[9]])), rep(10,length(indc_PTL[[10]])), 
                               rep(11,length(indc_PTL[[11]])),
                               rep(1,length(indc_PTL[[1]])), rep(2,length(indc_PTL[[2]])), 
                               rep(3,length(indc_PTL[[3]])), rep(4,length(indc_PTL[[4]])), 
                               rep(5,length(indc_PTL[[5]])), rep(6,length(indc_PTL[[6]])),
                               rep(7,length(indc_PTL[[7]])), rep(8,length(indc_PTL[[8]])), 
                               rep(9,length(indc_PTL[[9]])), rep(10,length(indc_PTL[[10]])), 
                               rep(11,length(indc_PTL[[11]])))),
                  model = c(rep("PTBM-R",length(indc_PTL[[11]])*11),rep("PTBM",length(indc_PTL[[11]])*11))) 
colnames(df_p) = c("p", "so", "model") 
p_2 = ggplot(df_p, aes(x=so, y=p, fill= model)) + 
  geom_boxplot() +
  scale_fill_manual(values = c("#555555", "grey")) +
  theme(legend.position = "none") + 
  labs(title = "Removal probability 0.1", x="sampling occasions",y=NULL)+
  theme(legend.justification = c(1, 0), 
        legend.position = c(0.95, 0.84))+
  coord_cartesian(ylim = c(-470, 200))+theme(axis.text = element_text(size = 13),
                                             legend.text = element_text(size = 13),
                                             legend.title = element_text(size = 13),
                                             axis.title.x = element_text(size = 13))
p_2


df_p = data.frame(x = c(as.vector(U_PTL_0.2[indc_PTL[[1]],1]), as.vector(U_PTL_0.2[indc_PTL[[2]],2]),
                        as.vector(U_PTL_0.2[indc_PTL[[3]],3]), as.vector(U_PTL_0.2[indc_PTL[[4]],4]),
                        as.vector(U_PTL_0.2[indc_PTL[[5]],5]), as.vector(U_PTL_0.2[indc_PTL[[6]],6]),
                        as.vector(U_PTL_0.2[indc_PTL[[7]],7]), as.vector(U_PTL_0.2[indc_PTL[[8]],8]),
                        as.vector(U_PTL_0.2[indc_PTL[[9]],9]), as.vector(U_PTL_0.2[indc_PTL[[10]],10]),
                        as.vector(U_PTL_0.2[indc_PTL[[11]],11]),
                        as.vector(U_PTWL_0.2[indc_PTL[[1]],1]), as.vector(U_PTWL_0.2[indc_PTL[[2]],2]),
                        as.vector(U_PTWL_0.2[indc_PTL[[3]],3]), as.vector(U_PTWL_0.2[indc_PTL[[4]],4]),
                        as.vector(U_PTWL_0.2[indc_PTL[[5]],5]), as.vector(U_PTWL_0.2[indc_PTL[[6]],6]),
                        as.vector(U_PTWL_0.2[indc_PTL[[7]],7]), as.vector(U_PTWL_0.2[indc_PTL[[8]],8]),
                        as.vector(U_PTWL_0.2[indc_PTL[[9]],9]), as.vector(U_PTWL_0.2[indc_PTL[[10]],10]),
                        as.vector(U_PTWL_0.2[indc_PTL[[11]],11])), 
                  y = factor(c(rep(1,length(indc_PTL[[1]])), rep(2,length(indc_PTL[[2]])), 
                               rep(3,length(indc_PTL[[3]])), rep(4,length(indc_PTL[[4]])), 
                               rep(5,length(indc_PTL[[5]])), rep(6,length(indc_PTL[[6]])),
                               rep(7,length(indc_PTL[[7]])), rep(8,length(indc_PTL[[8]])), 
                               rep(9,length(indc_PTL[[9]])), rep(10,length(indc_PTL[[10]])), 
                               rep(11,length(indc_PTL[[11]])),
                               rep(1,length(indc_PTL[[1]])), rep(2,length(indc_PTL[[2]])), 
                               rep(3,length(indc_PTL[[3]])), rep(4,length(indc_PTL[[4]])), 
                               rep(5,length(indc_PTL[[5]])), rep(6,length(indc_PTL[[6]])),
                               rep(7,length(indc_PTL[[7]])), rep(8,length(indc_PTL[[8]])), 
                               rep(9,length(indc_PTL[[9]])), rep(10,length(indc_PTL[[10]])), 
                               rep(11,length(indc_PTL[[11]])))),
                  model = c(rep("PTBM-R",length(indc_PTL[[11]])*11),rep("PTBM",length(indc_PTL[[11]])*11)))  
colnames(df_p) = c("p", "so", "model") 
p_3 = ggplot(df_p, aes(x=so, y=p, fill= model)) + 
  geom_boxplot() +
  theme(legend.position = "none")  +
  scale_fill_manual(values = c("#555555", "grey")) +
  labs(title = "Removal probability 0.2", x="sampling occasions", y= NULL) +
  coord_cartesian(ylim = c(-470, 200))+theme(axis.text = element_text(size = 13),
                                             legend.text = element_text(size = 13),
                                             legend.title = element_text(size = 13),
                                             axis.title.x = element_text(size = 13))


grid.arrange(p_1, p_2, p_3, ncol = 3)


####################################
####################################
####################################
####################################

N_with_0.05 = `N_0.05_with_list_unif_1 `
N_without_0.05 = `N_0.05_without_list_unif_1 `
N_with_0.1 = `N_0.1_with_list_unif_1 `
N_without_0.1 = `N_0.1_without_list_unif_1 `
N_with_0.2 = `N_0.2_with_list_unif_1 `
N_without_0.2 = `N_0.2_without_list_unif_1 `

N_0.05_true = `N_0.05_true_unif_1 `
N_PTL_0.05 = matrix(0,100,11)
for(i in 1:100){
  N_PTL_0.05[i,] =  apply(N_with_0.05[[i]][1:4000,],2,mean) - N_0.05_true
}
indc_PTL = list()
for(i in 1:11){
  indc_PTL[[i]] = which(abs(N_PTL_0.05[,i])<10000)
}

N_PTWL_0.05 = matrix(0,100,11)
for(i in 1:100){
  N_PTWL_0.05[i,] = apply(N_without_0.05[[i]][1:4000,],2,mean) - N_0.05_true
}


####
####
####



N_0.1_true = `N_0.1_true_unif_1 `

N_PTL_0.1 = matrix(0,100,11)
for(i in 1:100){
  N_PTL_0.1[i,] = apply(N_with_0.1[[i]][1:4000,],2,mean) - N_0.1_true
}
indc_PTL = list()
for(i in 1:11){
  indc_PTL[[i]] = which(abs(N_PTL_0.1[,i])<10000)
}


N_PTWL_0.1 = matrix(0,100,11)
for(i in 1:100){
  N_PTWL_0.1[i,] = apply(N_without_0.1[[i]][1:4000,],2,mean) - N_0.1_true
}


####
####
####

N_0.2_true = `N_0.2_true_unif_1 `
N_PTL_0.2 =  matrix(0,100,11)
for(i in 1:100){
  N_PTL_0.2[i,] = apply(N_with_0.2[[i]][1:4000,],2,mean) - N_0.2_true
}
indc_PTL = list()
for(i in 1:11){
  indc_PTL[[i]] = which(abs(N_PTL_0.2[,i])<10000)
}


N_PTWL_0.2 = matrix(0,100,11)
for(i in 1:100){
  N_PTWL_0.2[i,] = apply(N_without_0.2[[i]][1:4000,],2,mean) - N_0.2_true
}






df_p = data.frame(x = c(as.vector(N_PTL_0.05[indc_PTL[[1]],1]), as.vector(N_PTL_0.05[indc_PTL[[2]],2]),
                        as.vector(N_PTL_0.05[indc_PTL[[3]],3]), as.vector(N_PTL_0.05[indc_PTL[[4]],4]),
                        as.vector(N_PTL_0.05[indc_PTL[[5]],5]), as.vector(N_PTL_0.05[indc_PTL[[6]],6]),
                        as.vector(N_PTL_0.05[indc_PTL[[7]],7]), as.vector(N_PTL_0.05[indc_PTL[[8]],8]),
                        as.vector(N_PTL_0.05[indc_PTL[[9]],8]), as.vector(N_PTL_0.05[indc_PTL[[10]],10]),
                        as.vector(N_PTL_0.05[indc_PTL[[11]],11]),
                        as.vector(N_PTWL_0.05[indc_PTL[[1]],1]), as.vector(N_PTWL_0.05[indc_PTL[[2]],2]),
                        as.vector(N_PTWL_0.05[indc_PTL[[3]],3]), as.vector(N_PTWL_0.05[indc_PTL[[4]],4]),
                        as.vector(N_PTWL_0.05[indc_PTL[[5]],5]), as.vector(N_PTWL_0.05[indc_PTL[[6]],6]),
                        as.vector(N_PTWL_0.05[indc_PTL[[7]],7]), as.vector(N_PTWL_0.05[indc_PTL[[8]],8]),
                        as.vector(N_PTWL_0.05[indc_PTL[[9]],9]), as.vector(N_PTWL_0.05[indc_PTL[[10]],10]),
                        as.vector(N_PTWL_0.05[indc_PTL[[11]],11])), 
                  y = factor(c(rep(1,length(indc_PTL[[1]])), rep(2,length(indc_PTL[[2]])), 
                               rep(3,length(indc_PTL[[3]])), rep(4,length(indc_PTL[[4]])), 
                               rep(5,length(indc_PTL[[5]])), rep(6,length(indc_PTL[[6]])),
                               rep(7,length(indc_PTL[[7]])), rep(8,length(indc_PTL[[8]])), 
                               rep(9,length(indc_PTL[[9]])), rep(10,length(indc_PTL[[10]])), 
                               rep(11,length(indc_PTL[[11]])),
                               rep(1,length(indc_PTL[[1]])), rep(2,length(indc_PTL[[2]])), 
                               rep(3,length(indc_PTL[[3]])), rep(4,length(indc_PTL[[4]])), 
                               rep(5,length(indc_PTL[[5]])), rep(6,length(indc_PTL[[6]])),
                               rep(7,length(indc_PTL[[7]])), rep(8,length(indc_PTL[[8]])), 
                               rep(9,length(indc_PTL[[9]])), rep(10,length(indc_PTL[[10]])), 
                               rep(11,length(indc_PTL[[11]])))),
                  model = cc(rep("PTBM-R",length(indc_PTL[[11]])*11),rep("PTBM",length(indc_PTL[[11]])*11)))
colnames(df_p) = c("p", "so", "model") 
library(ggplot2)
p_1 = ggplot(df_p, aes(x=so, y=p, fill= model)) + 
  geom_boxplot() +
  scale_fill_manual(values = c("#555555", "grey")) +
  theme(legend.position = "none") + 
  labs(title = "Removal probability 0.05", x="sampling occasions", y = "N")+
  coord_cartesian(ylim = c(-650, 240))+theme(axis.text = element_text(size = 13),
                                             legend.text = element_text(size = 13),
                                             legend.title = element_text(size = 13),
                                             axis.title.x = element_text(size = 13))

df_p = data.frame(x = c(as.vector(N_PTL_0.1[indc_PTL[[1]],1]), as.vector(N_PTL_0.1[indc_PTL[[2]],2]),
                        as.vector(N_PTL_0.1[indc_PTL[[3]],3]), as.vector(N_PTL_0.1[indc_PTL[[4]],4]),
                        as.vector(N_PTL_0.1[indc_PTL[[5]],5]), as.vector(N_PTL_0.1[indc_PTL[[6]],6]),
                        as.vector(N_PTL_0.1[indc_PTL[[7]],7]), as.vector(N_PTL_0.1[indc_PTL[[8]],8]),
                        as.vector(N_PTL_0.1[indc_PTL[[9]],8]), as.vector(N_PTL_0.1[indc_PTL[[10]],10]),
                        as.vector(N_PTL_0.1[indc_PTL[[11]],11]),
                        as.vector(N_PTWL_0.1[indc_PTL[[1]],1]), as.vector(N_PTWL_0.1[indc_PTL[[2]],2]),
                        as.vector(N_PTWL_0.1[indc_PTL[[3]],3]), as.vector(N_PTWL_0.1[indc_PTL[[4]],4]),
                        as.vector(N_PTWL_0.1[indc_PTL[[5]],5]), as.vector(N_PTWL_0.1[indc_PTL[[6]],6]),
                        as.vector(N_PTWL_0.1[indc_PTL[[7]],7]), as.vector(N_PTWL_0.1[indc_PTL[[8]],8]),
                        as.vector(N_PTWL_0.1[indc_PTL[[9]],9]), as.vector(N_PTWL_0.1[indc_PTL[[10]],10]),
                        as.vector(N_PTWL_0.1[indc_PTL[[11]],11])), 
                  y = factor(c(rep(1,length(indc_PTL[[1]])), rep(2,length(indc_PTL[[2]])), 
                               rep(3,length(indc_PTL[[3]])), rep(4,length(indc_PTL[[4]])), 
                               rep(5,length(indc_PTL[[5]])), rep(6,length(indc_PTL[[6]])),
                               rep(7,length(indc_PTL[[7]])), rep(8,length(indc_PTL[[8]])), 
                               rep(9,length(indc_PTL[[9]])), rep(10,length(indc_PTL[[10]])), 
                               rep(11,length(indc_PTL[[11]])),
                               rep(1,length(indc_PTL[[1]])), rep(2,length(indc_PTL[[2]])), 
                               rep(3,length(indc_PTL[[3]])), rep(4,length(indc_PTL[[4]])), 
                               rep(5,length(indc_PTL[[5]])), rep(6,length(indc_PTL[[6]])),
                               rep(7,length(indc_PTL[[7]])), rep(8,length(indc_PTL[[8]])), 
                               rep(9,length(indc_PTL[[9]])), rep(10,length(indc_PTL[[10]])), 
                               rep(11,length(indc_PTL[[11]])))),
                  model = c(rep("PTBM-R",length(indc_PTL[[11]])*11),rep("PTBM",length(indc_PTL[[11]])*11))) 
colnames(df_p) = c("p", "so", "model") 
p_2 = ggplot(df_p, aes(x=so, y=p, fill= model)) + 
  geom_boxplot() +
  scale_fill_manual(values = c("#555555", "grey")) +
  theme(legend.position = "none") + 
  labs(title = "Removal probability 0.1", x="sampling occasions",y=NULL) +
  theme(legend.justification = c(1, 0), 
        legend.position = c(0.95, 0.84))+
  coord_cartesian(ylim = c(-650, 240))+theme(axis.text = element_text(size = 13),
                                             legend.text = element_text(size = 13),
                                             legend.title = element_text(size = 13),
                                             axis.title.x = element_text(size = 13))




df_p = data.frame(x = c(as.vector(N_PTL_0.2[indc_PTL[[1]],1]), as.vector(N_PTL_0.2[indc_PTL[[2]],2]),
                        as.vector(N_PTL_0.2[indc_PTL[[3]],3]), as.vector(N_PTL_0.2[indc_PTL[[4]],4]),
                        as.vector(N_PTL_0.2[indc_PTL[[5]],5]), as.vector(N_PTL_0.2[indc_PTL[[6]],6]),
                        as.vector(N_PTL_0.2[indc_PTL[[7]],7]), as.vector(N_PTL_0.2[indc_PTL[[8]],8]),
                        as.vector(N_PTL_0.2[indc_PTL[[9]],8]), as.vector(N_PTL_0.2[indc_PTL[[10]],10]),
                        as.vector(N_PTL_0.2[indc_PTL[[11]],11]),
                        as.vector(N_PTWL_0.2[indc_PTL[[1]],1]), as.vector(N_PTWL_0.2[indc_PTL[[2]],2]),
                        as.vector(N_PTWL_0.2[indc_PTL[[3]],3]), as.vector(N_PTWL_0.2[indc_PTL[[4]],4]),
                        as.vector(N_PTWL_0.2[indc_PTL[[5]],5]), as.vector(N_PTWL_0.2[indc_PTL[[6]],6]),
                        as.vector(N_PTWL_0.2[indc_PTL[[7]],7]), as.vector(N_PTWL_0.2[indc_PTL[[8]],8]),
                        as.vector(N_PTWL_0.2[indc_PTL[[9]],9]), as.vector(N_PTWL_0.2[indc_PTL[[10]],10]),
                        as.vector(N_PTWL_0.2[indc_PTL[[11]],11])), 
                  y = factor(c(rep(1,length(indc_PTL[[1]])), rep(2,length(indc_PTL[[2]])), 
                               rep(3,length(indc_PTL[[3]])), rep(4,length(indc_PTL[[4]])), 
                               rep(5,length(indc_PTL[[5]])), rep(6,length(indc_PTL[[6]])),
                               rep(7,length(indc_PTL[[7]])), rep(8,length(indc_PTL[[8]])), 
                               rep(9,length(indc_PTL[[9]])), rep(10,length(indc_PTL[[10]])), 
                               rep(11,length(indc_PTL[[11]])),
                               rep(1,length(indc_PTL[[1]])), rep(2,length(indc_PTL[[2]])), 
                               rep(3,length(indc_PTL[[3]])), rep(4,length(indc_PTL[[4]])), 
                               rep(5,length(indc_PTL[[5]])), rep(6,length(indc_PTL[[6]])),
                               rep(7,length(indc_PTL[[7]])), rep(8,length(indc_PTL[[8]])), 
                               rep(9,length(indc_PTL[[9]])), rep(10,length(indc_PTL[[10]])), 
                               rep(11,length(indc_PTL[[11]])))),
                  model = c(rep("PTBM-R",length(indc_PTL[[11]])*11),rep("PTBM",length(indc_PTL[[11]])*11)))  
colnames(df_p) = c("p", "so", "model") 
p_3 = ggplot(df_p, aes(x=so, y=p, fill= model)) + 
  geom_boxplot() +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("#555555", "grey")) +
  labs(title = "Removal probability 0.2", x="sampling occasions", y= NULL) +
  coord_cartesian(ylim = c(-650, 240))+theme(axis.text = element_text(size = 13),
                                             legend.text = element_text(size = 13),
                                             legend.title = element_text(size = 13),
                                             axis.title.x = element_text(size = 13))
p_3

grid.arrange(p_1, p_2, p_3, ncol = 3)

####################################
####################################
####################################
####################################

abund_with_0.05 = `abund_0.05_with_list_unif_1 `
abund_without_0.05 = `abund_0.05_without_list_unif_1 `
abund_with_0.1 = `abund_0.1_with_list_unif_1 `
abund_without_0.1 = `abund_0.1_without_list_unif_1 `
abund_with_0.2 = `abund_0.2_with_list_unif_1 `
abund_without_0.2 = `abund_0.2_without_list_unif_1 `

abund_true = 6000
abund_PTL_0.05 = rep(0,100)
for(i in 1:100){
  abund_PTL_0.05[i] =  mean(abund_with_0.05[[i]][1:200000]) - abund_true
}


abund_PTWL_0.05 = rep(0,100)
for(i in 1:100){
  abund_PTWL_0.05[i] = mean(abund_without_0.05[[i]][1:200000]) - abund_true
}


####
####
####


abund_PTL_0.1 = rep(0,100)
for(i in 1:100){
  abund_PTL_0.1[i] = median(abund_with_0.1[[i]][1:200000]) - abund_true
}


abund_PTWL_0.1 = rep(0,100)
for(i in 1:100){
  abund_PTWL_0.1[i] = mean(abund_without_0.1[[i]][1:200000]) - abund_true
}


####
####
####

abund_PTL_0.2 =  rep(0,100)
for(i in 1:100){
  abund_PTL_0.2[i] = mean(abund_with_0.2[[i]][1:200000]) - abund_true
}

abund_PTWL_0.2 = rep(0,100)
for(i in 1:100){
  abund_PTWL_0.2[i] = mean(abund_without_0.2[[i]][1:200000]) - abund_true
}



df_p = data.frame(x = c(as.vector(abund_PTL_0.05),
                        as.vector(abund_PTWL_0.05)),
                  model = c(rep("PTBM-R",100),rep("PTBM",100)))  
colnames(df_p) = c("p", "model") 
library(ggplot2)
p_1 = ggplot(df_p, aes( y=p, fill= model)) + 
  geom_boxplot() +
  scale_fill_manual(values = c("#555555", "grey")) +
  theme(legend.position = "none") + 
  labs(title = "Removal probability 0.05", x=NULL, y = "N") +
  coord_cartesian(ylim = c(-1100,250)) + theme(axis.text = element_text(size = 13),
                                               legend.text = element_text(size = 13),
                                               legend.title = element_text(size = 13))

df_p = data.frame(x = c(as.vector(abund_PTL_0.1),
                        as.vector(abund_PTWL_0.1)), 
                  model = c(rep("PTBM-R",100),rep("PTBM",100)))  
colnames(df_p) = c("p", "model") 
p_2 = ggplot(df_p, aes( y=p, fill= model)) + 
  geom_boxplot() +
  scale_fill_manual(values = c("#555555", "grey")) +
  theme(legend.position = "none") + 
  labs(title = "Removal probability 0.1", x=NULL,y=NULL) +
  theme(legend.justification = c(1, 0), 
        legend.position = c(0.3, 0.81))+
  coord_cartesian(ylim = c(-1100,250)) +theme(axis.text = element_text(size = 13),
                                              legend.text = element_text(size = 13),
                                              legend.title = element_text(size = 13)) 


df_p = data.frame(x = c(as.vector(abund_PTL_0.2),
                        as.vector(abund_PTWL_0.2)), 
                  model = c(rep("PTBM-R",100),rep("PTBM",100)))  
colnames(df_p) = c("p", "model") 
p_3 = ggplot(df_p, aes( y=p, fill= model)) + 
  geom_boxplot() +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("#555555", "grey")) +
  labs(title = "Removal probability 0.2", x=NULL, y= NULL)  +
  coord_cartesian(ylim = c(-1100,250)) + theme(axis.text = element_text(size = 13),
                                               legend.text = element_text(size = 13),
                                               legend.title = element_text(size = 13))

grid.arrange(p_1, p_2, p_3, ncol = 3)


####################################
####################################
####################################

Ome_with_0.05 = `Ome_0.05_with_list_unif_1 `
Ome_without_0.05 = `Ome_0.05_without_list_1 `
Ome_with_0.1 = `Ome_0.1_with_list_unif_1 `
Ome_without_0.1 = `Ome_0.1_without_list_1 `
Ome_with_0.2 = `Ome_0.2_with_list_unif_1 `
Ome_without_0.2 = `Ome_0.2_without_list_1 `

Ome_true_0.05 = `w_ij_0.05_true_unif_1 `
V_true = apply(Ome_true_0.05,2,sum)

Ome_PTL_0.05 = matrix(0,100,12)
for(i in 1:100){
  Ome_PTL_0.05[i,] =  apply(apply(Ome_with_0.05[[i]],c(1,2),mean),2,sum) - V_true
}
indc_PTL = list()
for(i in 1:12){
  indc_PTL[[i]] = which(abs(Ome_PTL_0.05[,i])<1)
}

Ome_PTWL_0.05 = matrix(0,100,12)
for(i in 1:100){
  Ome_PTWL_0.05[i,] = apply(apply(Ome_without_0.05[[i]],c(1,2),mean),2,sum) - V_true
}


####
####
####

Ome_true_0.1 = `w_ij_0.1_true_unif_1 `
V_true = apply(Ome_true_0.1,2,sum)

Ome_PTL_0.1 = matrix(0,100,12)
for(i in 1:100){
  Ome_PTL_0.1[i,] =  apply(apply(Ome_with_0.1[[i]],c(1,2),mean),2,sum) - V_true
}
indc_PTL = list()
for(i in 1:12){
  indc_PTL[[i]] = which(abs(Ome_PTL_0.1[,i])<1)
}

Ome_PTWL_0.1 = matrix(0,100,12)
for(i in 1:100){
  Ome_PTWL_0.1[i,] = apply(apply(Ome_without_0.1[[i]],c(1,2),mean),2,sum) - V_true
}

####
####
####

Ome_true_0.2 = `w_ij_0.2_true_unif_1 `
V_true = apply(Ome_true_0.2,2,sum)

Ome_PTL_0.2 = matrix(0,100,12)
for(i in 1:100){
  Ome_PTL_0.2[i,] =  apply(apply(Ome_with_0.2[[i]],c(1,2),mean),2,sum) - V_true
}
indc_PTL = list()
for(i in 1:12){
  indc_PTL[[i]] = which(abs(Ome_PTL_0.2[,i])<1)

}

Ome_PTWL_0.2 = matrix(0,100,12)
for(i in 1:100){
  Ome_PTWL_0.2[i,] = apply(apply(Ome_without_0.2[[i]],c(1,2),mean),2,sum) - V_true
}



df_p = data.frame(x = c(as.vector(Ome_PTL_0.05[indc_PTL[[1]],1]), as.vector(Ome_PTL_0.05[indc_PTL[[2]],2]),
                        as.vector(Ome_PTL_0.05[indc_PTL[[3]],3]), as.vector(Ome_PTL_0.05[indc_PTL[[4]],4]),
                        as.vector(Ome_PTL_0.05[indc_PTL[[5]],5]), as.vector(Ome_PTL_0.05[indc_PTL[[6]],6]),
                        as.vector(Ome_PTL_0.05[indc_PTL[[7]],7]), as.vector(Ome_PTL_0.05[indc_PTL[[8]],8]),
                        as.vector(Ome_PTL_0.05[indc_PTL[[9]],8]), as.vector(Ome_PTL_0.05[indc_PTL[[10]],10]),
                        as.vector(Ome_PTL_0.05[indc_PTL[[11]],11]), as.vector(Ome_PTL_0.05[indc_PTL[[12]],12]),
                        as.vector(Ome_PTWL_0.05[indc_PTL[[1]],1]), as.vector(Ome_PTWL_0.05[indc_PTL[[2]],2]),
                        as.vector(Ome_PTWL_0.05[indc_PTL[[3]],3]), as.vector(Ome_PTWL_0.05[indc_PTL[[4]],4]),
                        as.vector(Ome_PTWL_0.05[indc_PTL[[5]],5]), as.vector(Ome_PTWL_0.05[indc_PTL[[6]],6]),
                        as.vector(Ome_PTWL_0.05[indc_PTL[[7]],7]), as.vector(Ome_PTWL_0.05[indc_PTL[[8]],8]),
                        as.vector(Ome_PTWL_0.05[indc_PTL[[9]],9]), as.vector(Ome_PTWL_0.05[indc_PTL[[10]],10]),
                        as.vector(Ome_PTWL_0.05[indc_PTL[[11]],11]), as.vector(Ome_PTWL_0.05[indc_PTL[[12]],12])), 
                  y = factor(c(rep(1,length(indc_PTL[[1]])), rep(2,length(indc_PTL[[2]])), 
                               rep(3,length(indc_PTL[[3]])), rep(4,length(indc_PTL[[4]])), 
                               rep(5,length(indc_PTL[[5]])), rep(6,length(indc_PTL[[6]])),
                               rep(7,length(indc_PTL[[7]])), rep(8,length(indc_PTL[[8]])), 
                               rep(9,length(indc_PTL[[9]])), rep(10,length(indc_PTL[[10]])), 
                               rep(11,length(indc_PTL[[11]])), rep(12,length(indc_PTL[[12]])),
                               rep(1,length(indc_PTL[[1]])), rep(2,length(indc_PTL[[2]])), 
                               rep(3,length(indc_PTL[[3]])), rep(4,length(indc_PTL[[4]])), 
                               rep(5,length(indc_PTL[[5]])), rep(6,length(indc_PTL[[6]])),
                               rep(7,length(indc_PTL[[7]])), rep(8,length(indc_PTL[[8]])), 
                               rep(9,length(indc_PTL[[9]])), rep(10,length(indc_PTL[[10]])), 
                               rep(11,length(indc_PTL[[11]])), rep(12,length(indc_PTL[[12]])))),
                  model = c(rep("PTBM-R",length(indc_PTL[[11]])*11),rep("PTBM",length(indc_PTL[[11]])*11)))
colnames(df_p) = c("p", "so", "model") 
library(ggplot2)
p_1 = ggplot(df_p, aes(x=so, y=p, fill= model)) + 
  geom_boxplot() +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("#555555", "grey")) +
  labs(title = "Removal probability 0.05", x="sampling occasions", y= "A") +
  coord_cartesian(ylim = c(-0.02, 0.02))+ theme(axis.text = element_text(size = 13),
                                                legend.text = element_text(size = 13),
                                                legend.title = element_text(size = 13))


df_p = data.frame(x = c(as.vector(Ome_PTL_0.1[indc_PTL[[1]],1]), as.vector(Ome_PTL_0.1[indc_PTL[[2]],2]),
                        as.vector(Ome_PTL_0.1[indc_PTL[[3]],3]), as.vector(Ome_PTL_0.1[indc_PTL[[4]],4]),
                        as.vector(Ome_PTL_0.1[indc_PTL[[5]],5]), as.vector(Ome_PTL_0.1[indc_PTL[[6]],6]),
                        as.vector(Ome_PTL_0.1[indc_PTL[[7]],7]), as.vector(Ome_PTL_0.1[indc_PTL[[8]],8]),
                        as.vector(Ome_PTL_0.1[indc_PTL[[9]],8]), as.vector(Ome_PTL_0.1[indc_PTL[[10]],10]),
                        as.vector(Ome_PTL_0.1[indc_PTL[[11]],11]), as.vector(Ome_PTL_0.1[indc_PTL[[12]],12]),
                        as.vector(Ome_PTWL_0.1[indc_PTL[[1]],1]), as.vector(Ome_PTWL_0.1[indc_PTL[[2]],2]),
                        as.vector(Ome_PTWL_0.1[indc_PTL[[3]],3]), as.vector(Ome_PTWL_0.1[indc_PTL[[4]],4]),
                        as.vector(Ome_PTWL_0.1[indc_PTL[[5]],5]), as.vector(Ome_PTWL_0.1[indc_PTL[[6]],6]),
                        as.vector(Ome_PTWL_0.1[indc_PTL[[7]],7]), as.vector(Ome_PTWL_0.1[indc_PTL[[8]],8]),
                        as.vector(Ome_PTWL_0.1[indc_PTL[[9]],9]), as.vector(Ome_PTWL_0.1[indc_PTL[[10]],10]),
                        as.vector(Ome_PTWL_0.1[indc_PTL[[11]],11]), as.vector(Ome_PTWL_0.1[indc_PTL[[12]],12])), 
                  y = factor(c(rep(1,length(indc_PTL[[1]])), rep(2,length(indc_PTL[[2]])), 
                               rep(3,length(indc_PTL[[3]])), rep(4,length(indc_PTL[[4]])), 
                               rep(5,length(indc_PTL[[5]])), rep(6,length(indc_PTL[[6]])),
                               rep(7,length(indc_PTL[[7]])), rep(8,length(indc_PTL[[8]])), 
                               rep(9,length(indc_PTL[[9]])), rep(10,length(indc_PTL[[10]])), 
                               rep(11,length(indc_PTL[[11]])), rep(12,length(indc_PTL[[12]])),
                               rep(1,length(indc_PTL[[1]])), rep(2,length(indc_PTL[[2]])), 
                               rep(3,length(indc_PTL[[3]])), rep(4,length(indc_PTL[[4]])), 
                               rep(5,length(indc_PTL[[5]])), rep(6,length(indc_PTL[[6]])),
                               rep(7,length(indc_PTL[[7]])), rep(8,length(indc_PTL[[8]])), 
                               rep(9,length(indc_PTL[[9]])), rep(10,length(indc_PTL[[10]])), 
                               rep(11,length(indc_PTL[[11]])), rep(12,length(indc_PTL[[12]])))),
                  model = c(rep("PTBM-R",length(indc_PTL[[11]])*11),rep("PTBM",length(indc_PTL[[11]])*11)))
colnames(df_p) = c("p", "so", "model") 
p_2 = ggplot(df_p, aes(x=so, y=p, fill= model)) + 
  geom_boxplot() +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("#555555", "grey")) +
  labs(title = "Removal probability 0.1", x="sampling occasions", y= NULL) +
  theme(legend.justification = c(1, 0), 
        legend.position = c(0.95, 0.84))+
  coord_cartesian(ylim = c(-0.02, 0.02))+ theme(axis.text = element_text(size = 13),
                                                legend.text = element_text(size = 13),
                                                legend.title = element_text(size = 13))


df_p = data.frame(x = c(as.vector(Ome_PTL_0.2[indc_PTL[[1]],1]), as.vector(Ome_PTL_0.2[indc_PTL[[2]],2]),
                        as.vector(Ome_PTL_0.2[indc_PTL[[3]],3]), as.vector(Ome_PTL_0.2[indc_PTL[[4]],4]),
                        as.vector(Ome_PTL_0.2[indc_PTL[[5]],5]), as.vector(Ome_PTL_0.2[indc_PTL[[6]],6]),
                        as.vector(Ome_PTL_0.2[indc_PTL[[7]],7]), as.vector(Ome_PTL_0.2[indc_PTL[[8]],8]),
                        as.vector(Ome_PTL_0.2[indc_PTL[[9]],9]), as.vector(Ome_PTL_0.2[indc_PTL[[10]],10]),
                        as.vector(Ome_PTL_0.2[indc_PTL[[11]],11]), as.vector(Ome_PTL_0.2[indc_PTL[[12]],12]),
                        as.vector(Ome_PTWL_0.2[indc_PTL[[1]],1]), as.vector(Ome_PTWL_0.2[indc_PTL[[2]],2]),
                        as.vector(Ome_PTWL_0.2[indc_PTL[[3]],3]), as.vector(Ome_PTWL_0.2[indc_PTL[[4]],4]),
                        as.vector(Ome_PTWL_0.2[indc_PTL[[5]],5]), as.vector(Ome_PTWL_0.2[indc_PTL[[6]],6]),
                        as.vector(Ome_PTWL_0.2[indc_PTL[[7]],7]), as.vector(Ome_PTWL_0.2[indc_PTL[[8]],8]),
                        as.vector(Ome_PTWL_0.2[indc_PTL[[9]],9]), as.vector(Ome_PTWL_0.2[indc_PTL[[10]],10]),
                        as.vector(Ome_PTWL_0.2[indc_PTL[[11]],11]), as.vector(Ome_PTWL_0.2[indc_PTL[[12]],12])), 
                  y = factor(c(rep(1,length(indc_PTL[[1]])), rep(2,length(indc_PTL[[2]])), 
                               rep(3,length(indc_PTL[[3]])), rep(4,length(indc_PTL[[4]])), 
                               rep(5,length(indc_PTL[[5]])), rep(6,length(indc_PTL[[6]])),
                               rep(7,length(indc_PTL[[7]])), rep(8,length(indc_PTL[[8]])), 
                               rep(9,length(indc_PTL[[9]])), rep(10,length(indc_PTL[[10]])), 
                               rep(11,length(indc_PTL[[11]])), rep(12,length(indc_PTL[[12]])),
                               rep(1,length(indc_PTL[[1]])), rep(2,length(indc_PTL[[2]])), 
                               rep(3,length(indc_PTL[[3]])), rep(4,length(indc_PTL[[4]])), 
                               rep(5,length(indc_PTL[[5]])), rep(6,length(indc_PTL[[6]])),
                               rep(7,length(indc_PTL[[7]])), rep(8,length(indc_PTL[[8]])), 
                               rep(9,length(indc_PTL[[9]])), rep(10,length(indc_PTL[[10]])), 
                               rep(11,length(indc_PTL[[11]])), rep(12,length(indc_PTL[[12]])))),
                  model = c(rep("PTBM-R",length(indc_PTL[[11]])*11),rep("PTBM",length(indc_PTL[[11]])*11)))
colnames(df_p) = c("p", "so", "model") 
p_3 = ggplot(df_p, aes(x=so, y=p, fill= model)) + 
  geom_boxplot() +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("#555555", "grey")) +
  labs(title = "Removal probability 0.2", x="sampling occasions", y= NULL) +
  coord_cartesian(ylim = c(-0.02, 0.02))+ theme(axis.text = element_text(size = 13),
                                                legend.text = element_text(size = 13),
                                                legend.title = element_text(size = 13))


grid.arrange(p_1, p_2, p_3, ncol = 3)


####################################
####################################
####################################

Ome_with_0.05 = `Ome_0.05_with_list_unif_1 `
Ome_without_0.05 = `Ome_0.05_without_list_1 `
Ome_with_0.1 = `Ome_0.1_with_list_unif_1 `
Ome_without_0.1 = `Ome_0.1_without_list_1 `
Ome_with_0.2 = `Ome_0.2_with_list_unif_1 `
Ome_without_0.2 = `Ome_0.2_without_list_1 `

Ome_true_0.05 = `w_ij_0.05_true_unif_1 `
W_true = apply(Ome_true_0.05,1,sum)

Ome_PTL_0.05 = matrix(0,100,12)
for(i in 1:100){
  Ome_PTL_0.05[i,] =  apply(apply(Ome_with_0.05[[i]],c(1,2),mean),1,sum) - W_true
}
indc_PTL = list()
for(i in 1:12){
  indc_PTL[[i]] = which(abs(Ome_PTL_0.05[,i])<1)

}

Ome_PTWL_0.05 = matrix(0,100,12)
for(i in 1:100){
  Ome_PTWL_0.05[i,] = apply(apply(Ome_without_0.05[[i]],c(1,2),mean),1,sum) - W_true
}


####
####
####
Ome_true_0.1 = `w_ij_0.1_true_unif_1 `
W_true = apply(Ome_true_0.1,1,sum)

Ome_PTL_0.1 = matrix(0,100,12)
for(i in 1:100){
  Ome_PTL_0.1[i,] =  apply(apply(Ome_with_0.1[[i]],c(1,2),mean),1,sum) - W_true
}
indc_PTL = list()
for(i in 1:12){
  indc_PTL[[i]] = which(abs(Ome_PTL_0.1[,i])<1)

}

Ome_PTWL_0.1 = matrix(0,100,12)
for(i in 1:100){
  Ome_PTWL_0.1[i,] = apply(apply(Ome_without_0.1[[i]],c(1,2),mean),1,sum) - W_true
}


####
####
####

Ome_true_0.2 = `w_ij_0.2_true_unif_1 `
W_true = apply(Ome_true_0.2,1,sum)

Ome_PTL_0.2 = matrix(0,100,12)
for(i in 1:100){
  Ome_PTL_0.2[i,] =  apply(apply(Ome_with_0.2[[i]],c(1,2),mean),1,sum) - W_true
}
indc_PTL = list()
for(i in 1:12){
  indc_PTL[[i]] = which(abs(Ome_PTL_0.2[,i])<1)

}

Ome_PTWL_0.2 = matrix(0,100,12)
for(i in 1:100){
  Ome_PTWL_0.2[i,] = apply(apply(Ome_without_0.2[[i]],c(1,2),mean),1,sum) - W_true
}


df_p = data.frame(x = c(as.vector(Ome_PTL_0.05[indc_PTL[[1]],1]), as.vector(Ome_PTL_0.05[indc_PTL[[2]],2]),
                        as.vector(Ome_PTL_0.05[indc_PTL[[3]],3]), as.vector(Ome_PTL_0.05[indc_PTL[[4]],4]),
                        as.vector(Ome_PTL_0.05[indc_PTL[[5]],5]), as.vector(Ome_PTL_0.05[indc_PTL[[6]],6]),
                        as.vector(Ome_PTL_0.05[indc_PTL[[7]],7]), as.vector(Ome_PTL_0.05[indc_PTL[[8]],8]),
                        as.vector(Ome_PTL_0.05[indc_PTL[[9]],9]), as.vector(Ome_PTL_0.05[indc_PTL[[10]],10]),
                        as.vector(Ome_PTL_0.05[indc_PTL[[11]],11]), as.vector(Ome_PTL_0.05[indc_PTL[[12]],12]),
                        as.vector(Ome_PTWL_0.05[indc_PTL[[1]],1]), as.vector(Ome_PTWL_0.05[indc_PTL[[2]],2]),
                        as.vector(Ome_PTWL_0.05[indc_PTL[[3]],3]), as.vector(Ome_PTWL_0.05[indc_PTL[[4]],4]),
                        as.vector(Ome_PTWL_0.05[indc_PTL[[5]],5]), as.vector(Ome_PTWL_0.05[indc_PTL[[6]],6]),
                        as.vector(Ome_PTWL_0.05[indc_PTL[[7]],7]), as.vector(Ome_PTWL_0.05[indc_PTL[[8]],8]),
                        as.vector(Ome_PTWL_0.05[indc_PTL[[9]],9]), as.vector(Ome_PTWL_0.05[indc_PTL[[10]],10]),
                        as.vector(Ome_PTWL_0.05[indc_PTL[[11]],11]), as.vector(Ome_PTWL_0.05[indc_PTL[[12]],12])), 
                  y = factor(c(rep(1,length(indc_PTL[[1]])), rep(2,length(indc_PTL[[2]])), 
                               rep(3,length(indc_PTL[[3]])), rep(4,length(indc_PTL[[4]])), 
                               rep(5,length(indc_PTL[[5]])), rep(6,length(indc_PTL[[6]])),
                               rep(7,length(indc_PTL[[7]])), rep(8,length(indc_PTL[[8]])), 
                               rep(9,length(indc_PTL[[9]])), rep(10,length(indc_PTL[[10]])), 
                               rep(11,length(indc_PTL[[11]])), rep(12,length(indc_PTL[[12]])),
                               rep(1,length(indc_PTL[[1]])), rep(2,length(indc_PTL[[2]])), 
                               rep(3,length(indc_PTL[[3]])), rep(4,length(indc_PTL[[4]])), 
                               rep(5,length(indc_PTL[[5]])), rep(6,length(indc_PTL[[6]])),
                               rep(7,length(indc_PTL[[7]])), rep(8,length(indc_PTL[[8]])), 
                               rep(9,length(indc_PTL[[9]])), rep(10,length(indc_PTL[[10]])), 
                               rep(11,length(indc_PTL[[11]])), rep(12,length(indc_PTL[[12]])))),
                  model = c(rep("PTBM-R",length(indc_PTL[[11]])*11),rep("PTBM",length(indc_PTL[[11]])*11)))
colnames(df_p) = c("p", "so", "model") 
p_1 = ggplot(df_p, aes(x=so, y=p, fill= model)) + 
  geom_boxplot() +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("#555555", "grey")) +
  labs(title = "Removal probability 0.05", x="sampling occasions", y= "A") +
  coord_cartesian(ylim = c(-0.02, 0.02))+ theme(axis.text = element_text(size = 13),
                                                legend.text = element_text(size = 13),
                                                legend.title = element_text(size = 13))


df_p = data.frame(x = c(as.vector(Ome_PTL_0.1[indc_PTL[[1]],1]), as.vector(Ome_PTL_0.1[indc_PTL[[2]],2]),
                        as.vector(Ome_PTL_0.1[indc_PTL[[3]],3]), as.vector(Ome_PTL_0.1[indc_PTL[[4]],4]),
                        as.vector(Ome_PTL_0.1[indc_PTL[[5]],5]), as.vector(Ome_PTL_0.1[indc_PTL[[6]],6]),
                        as.vector(Ome_PTL_0.1[indc_PTL[[7]],7]), as.vector(Ome_PTL_0.1[indc_PTL[[8]],8]),
                        as.vector(Ome_PTL_0.1[indc_PTL[[9]],9]), as.vector(Ome_PTL_0.1[indc_PTL[[10]],10]),
                        as.vector(Ome_PTL_0.1[indc_PTL[[11]],11]), as.vector(Ome_PTL_0.1[indc_PTL[[12]],12]),
                        as.vector(Ome_PTWL_0.1[indc_PTL[[1]],1]), as.vector(Ome_PTWL_0.1[indc_PTL[[2]],2]),
                        as.vector(Ome_PTWL_0.1[indc_PTL[[3]],3]), as.vector(Ome_PTWL_0.1[indc_PTL[[4]],4]),
                        as.vector(Ome_PTWL_0.1[indc_PTL[[5]],5]), as.vector(Ome_PTWL_0.1[indc_PTL[[6]],6]),
                        as.vector(Ome_PTWL_0.1[indc_PTL[[7]],7]), as.vector(Ome_PTWL_0.1[indc_PTL[[8]],8]),
                        as.vector(Ome_PTWL_0.1[indc_PTL[[9]],9]), as.vector(Ome_PTWL_0.1[indc_PTL[[10]],10]),
                        as.vector(Ome_PTWL_0.1[indc_PTL[[11]],11]), as.vector(Ome_PTWL_0.1[indc_PTL[[12]],12])), 
                  y = factor(c(rep(1,length(indc_PTL[[1]])), rep(2,length(indc_PTL[[2]])), 
                               rep(3,length(indc_PTL[[3]])), rep(4,length(indc_PTL[[4]])), 
                               rep(5,length(indc_PTL[[5]])), rep(6,length(indc_PTL[[6]])),
                               rep(7,length(indc_PTL[[7]])), rep(8,length(indc_PTL[[8]])), 
                               rep(9,length(indc_PTL[[9]])), rep(10,length(indc_PTL[[10]])), 
                               rep(11,length(indc_PTL[[11]])), rep(12,length(indc_PTL[[12]])),
                               rep(1,length(indc_PTL[[1]])), rep(2,length(indc_PTL[[2]])), 
                               rep(3,length(indc_PTL[[3]])), rep(4,length(indc_PTL[[4]])), 
                               rep(5,length(indc_PTL[[5]])), rep(6,length(indc_PTL[[6]])),
                               rep(7,length(indc_PTL[[7]])), rep(8,length(indc_PTL[[8]])), 
                               rep(9,length(indc_PTL[[9]])), rep(10,length(indc_PTL[[10]])), 
                               rep(11,length(indc_PTL[[11]])), rep(12,length(indc_PTL[[12]])))),
                  model = c(rep("PTBM-R",length(indc_PTL[[11]])*11),rep("PTBM",length(indc_PTL[[11]])*11)))
colnames(df_p) = c("p", "so", "model") 
p_2 = ggplot(df_p, aes(x=so, y=p, fill= model)) + 
  geom_boxplot() +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("#555555", "grey")) +
  labs(title = "Removal probability 0.1", x="sampling occasions", y= NULL) +
  theme(legend.justification = c(1, 0), 
        legend.position = c(0.95, 0.84))+
  coord_cartesian(ylim = c(-0.02, 0.02))+ theme(axis.text = element_text(size = 13),
                                                legend.text = element_text(size = 13),
                                                legend.title = element_text(size = 13))

df_p = data.frame(x = c(as.vector(Ome_PTL_0.2[indc_PTL[[1]],1]), as.vector(Ome_PTL_0.2[indc_PTL[[2]],2]),
                        as.vector(Ome_PTL_0.2[indc_PTL[[3]],3]), as.vector(Ome_PTL_0.2[indc_PTL[[4]],4]),
                        as.vector(Ome_PTL_0.2[indc_PTL[[5]],5]), as.vector(Ome_PTL_0.2[indc_PTL[[6]],6]),
                        as.vector(Ome_PTL_0.2[indc_PTL[[7]],7]), as.vector(Ome_PTL_0.2[indc_PTL[[8]],8]),
                        as.vector(Ome_PTL_0.2[indc_PTL[[9]],9]), as.vector(Ome_PTL_0.2[indc_PTL[[10]],10]),
                        as.vector(Ome_PTL_0.2[indc_PTL[[11]],11]), as.vector(Ome_PTL_0.2[indc_PTL[[12]],12]),
                        as.vector(Ome_PTWL_0.2[indc_PTL[[1]],1]), as.vector(Ome_PTWL_0.2[indc_PTL[[2]],2]),
                        as.vector(Ome_PTWL_0.2[indc_PTL[[3]],3]), as.vector(Ome_PTWL_0.2[indc_PTL[[4]],4]),
                        as.vector(Ome_PTWL_0.2[indc_PTL[[5]],5]), as.vector(Ome_PTWL_0.2[indc_PTL[[6]],6]),
                        as.vector(Ome_PTWL_0.2[indc_PTL[[7]],7]), as.vector(Ome_PTWL_0.2[indc_PTL[[8]],8]),
                        as.vector(Ome_PTWL_0.2[indc_PTL[[9]],9]), as.vector(Ome_PTWL_0.2[indc_PTL[[10]],10]),
                        as.vector(Ome_PTWL_0.2[indc_PTL[[11]],11]), as.vector(Ome_PTWL_0.2[indc_PTL[[12]],12])), 
                  y = factor(c(rep(1,length(indc_PTL[[1]])), rep(2,length(indc_PTL[[2]])), 
                               rep(3,length(indc_PTL[[3]])), rep(4,length(indc_PTL[[4]])), 
                               rep(5,length(indc_PTL[[5]])), rep(6,length(indc_PTL[[6]])),
                               rep(7,length(indc_PTL[[7]])), rep(8,length(indc_PTL[[8]])), 
                               rep(9,length(indc_PTL[[9]])), rep(10,length(indc_PTL[[10]])), 
                               rep(11,length(indc_PTL[[11]])), rep(12,length(indc_PTL[[12]])),
                               rep(1,length(indc_PTL[[1]])), rep(2,length(indc_PTL[[2]])), 
                               rep(3,length(indc_PTL[[3]])), rep(4,length(indc_PTL[[4]])), 
                               rep(5,length(indc_PTL[[5]])), rep(6,length(indc_PTL[[6]])),
                               rep(7,length(indc_PTL[[7]])), rep(8,length(indc_PTL[[8]])), 
                               rep(9,length(indc_PTL[[9]])), rep(10,length(indc_PTL[[10]])), 
                               rep(11,length(indc_PTL[[11]])), rep(12,length(indc_PTL[[12]])))),
                  model = c(rep("PTBM-R",length(indc_PTL[[11]])*11),rep("PTBM",length(indc_PTL[[11]])*11)))
colnames(df_p) = c("p", "so", "model") 
p_3 = ggplot(df_p, aes(x=so, y=p, fill= model)) + 
  geom_boxplot() +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("#555555", "grey")) +
  labs(title = "Removal probability 0.2", x="sampling occasions", y= NULL) +
  coord_cartesian(ylim = c(-0.02, 0.02))+ theme(axis.text = element_text(size = 13),
                                                legend.text = element_text(size = 13),
                                                legend.title = element_text(size = 13))

grid.arrange(p_1, p_2, p_3, ncol = 3)

