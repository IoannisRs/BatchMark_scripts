library(gridExtra)
library(ggplot2)
V_1 = readRDS("V_runs_alt_unif_bias_1_so .RDS")
V_1_true =  c(0.4114271,
              0.5394092,
              0.4981504,
              0.5770110,
              0.3103370,
              0.4803437)
V_2 = readRDS("V_runs_alt_unif_bias_2_so .RDS")
V_2_true = c(0.4114271,
             0.5394092,
             0.4981504,
             0.5770110,
             0.3103370,
             0.4803437)
V_4 = readRDS("V_runs_alt_unif_bias_4_so .RDS")
V_4_true = c(0.4114271,
             0.5394092,
             0.4981504,
             0.5770110,
             0.3103370,
             0.4803437)
V_6 = readRDS("V_runs_alt_unif_bias_8_so .RDS")
V_6_true =  c(0.4114271,
              0.5394092,
              0.4981504,
              0.5770110,
              0.3103370,
              0.4803437)
V_7 = readRDS("V_runs_alt_unif_bias_16_so .RDS")
V_7_true =  c(0.4114271,
              0.5394092,
              0.4981504,
              0.5770110,
              0.3103370,
              0.4803437)

V_1_RMSE = matrix(0,100,6)
for(i in 1:100){
  for(j in 1000:5000){
    V_1_RMSE[i,] = V_1_RMSE[i,] + (V_1[[i]][j,] - V_1_true)^2
  }
  V_1_RMSE[i,] = sqrt(V_1_RMSE[i,]/4000)
}

V_2_RMSE = matrix(0,100,6)
for(i in 1:100){
  for(j in 1000:5000){
    V_2_RMSE[i,] = V_2_RMSE[i,] + (V_2[[i]][j,] - V_2_true)^2
  }
  V_2_RMSE[i,] = sqrt(V_2_RMSE[i,]/4000)
}

V_4_RMSE = matrix(0,100,6)
for(i in 1:100){
  for(j in 1000:5000){
    V_4_RMSE[i,] = V_4_RMSE[i,] + (V_4[[i]][j,] - V_4_true)^2
  }
  V_4_RMSE[i,] = sqrt(V_4_RMSE[i,]/4000)
}

V_6_RMSE = matrix(0,100,6)
for(i in 1:100){
  for(j in 1000:5000){
    V_6_RMSE[i,] = V_6_RMSE[i,] + (V_6[[i]][j,] - V_6_true)^2
  }
  V_6_RMSE[i,] = sqrt(V_6_RMSE[i,]/4000)
}

V_7_RMSE = matrix(0,100,6)
for(i in 1:100){
  for(j in 1000:5000){
    V_7_RMSE[i,] = V_7_RMSE[i,] + (V_7[[i]][j,] - V_7_true)^2
  }
  V_7_RMSE[i,] = sqrt(V_7_RMSE[i,]/4000)
}



mean((mean(apply(V_1_RMSE,2,median)) - mean(apply(V_2_RMSE,2,median)))/mean(apply(V_1_RMSE,2,median)))
mean((mean(apply(V_2_RMSE,2,median)) - mean(apply(V_4_RMSE,2,median)))/mean(apply(V_2_RMSE,2,median)))
mean((mean(apply(V_4_RMSE,2,median)) - mean(apply(V_6_RMSE,2,median)))/mean(apply(V_4_RMSE,2,median)))
mean((mean(apply(V_6_RMSE,2,median)) -mean(apply(V_7_RMSE,2,median)))/mean(apply(V_6_RMSE,2,median)))

apply(V_7_RMSE,2,mean)




df_p = data.frame(x = c(as.vector(V_1_RMSE),as.vector(V_2_RMSE),as.vector(V_4_RMSE),
                        as.vector(V_6_RMSE),as.vector(V_7_RMSE)), 
                  K = factor(c(rep(1,dim(V_1_RMSE)[1]*6),rep(2,dim(V_2_RMSE)[1]*6),rep(4,dim(V_4_RMSE)[1]*6),
                               rep(8,dim(V_6_RMSE)[1]*6),rep(16,dim(V_7_RMSE)[1]*6))),
                  t = factor(c(rep(1,dim(V_1_RMSE)[1]),rep(2,dim(V_1_RMSE)[1]),rep(3,dim(V_1_RMSE)[1]),
                               rep(4,dim(V_1_RMSE)[1]),rep(5,dim(V_1_RMSE)[1]),rep(6,dim(V_1_RMSE)[1]),
                               rep(1,dim(V_2_RMSE)[1]),rep(2,dim(V_2_RMSE)[1]),rep(3,dim(V_2_RMSE)[1]),
                               rep(4,dim(V_2_RMSE)[1]),rep(5,dim(V_2_RMSE)[1]),rep(6,dim(V_2_RMSE)[1]),
                               rep(1,dim(V_4_RMSE)[1]),rep(2,dim(V_4_RMSE)[1]),rep(3,dim(V_4_RMSE)[1]),
                               rep(4,dim(V_4_RMSE)[1]),rep(5,dim(V_4_RMSE)[1]),rep(6,dim(V_4_RMSE)[1]),
                               rep(1,dim(V_6_RMSE)[1]),rep(2,dim(V_6_RMSE)[1]),rep(3,dim(V_6_RMSE)[1]),
                               rep(4,dim(V_6_RMSE)[1]),rep(5,dim(V_6_RMSE)[1]),rep(6,dim(V_6_RMSE)[1]),
                               rep(1,dim(V_7_RMSE)[1]),rep(2,dim(V_7_RMSE)[1]),rep(3,dim(V_7_RMSE)[1]),
                               rep(4,dim(V_7_RMSE)[1]),rep(5,dim(V_7_RMSE)[1]),rep(6,dim(V_7_RMSE)[1]))))   
library(ggplot2)
p_1 = ggplot(df_p, aes(x=t, y=x, fill= K)) + 
  geom_boxplot() +
  scale_fill_manual(values = c("grey30", "grey45", "grey65", "grey80","grey90")) + 
  labs( x="Primary periods", y = "RMSE V_i") +
  ylim(0, 0.4)+theme(axis.text = element_text(size = 13),
                     legend.text = element_text(size = 13),
                     legend.title = element_text(size = 13),
                     axis.title.x = element_text(size = 13)) 
p_1
V_data_frame = data.frame(rbind(cbind(V_1_mean,t(V_1_CI)),
                                cbind(V_2_mean,t(V_2_CI)), 
                                cbind(V_4_mean,t(V_4_CI)),
                                cbind(V_7_mean,t(V_7_CI))),K = factor(c(rep(1,5),rep(2,5),rep(4,5),
                                                                        rep(15,5))),
                          t = factor(rep(c("1","2","3","4","5"),4)))

library(ggplot2)
ggplot(V_data_frame, aes(x = t, y = V_1_mean)) +
  geom_errorbar(
    aes(ymin = X2.5., ymax = X97.5., color = K),
    position = position_dodge(0.3),width = 0.2
  ) + 
  geom_point(aes(color = K), position = position_dodge(0.3)) +
  ylab("V") + xlab("t") +
  scale_color_manual(values = c("#333333","#555555","#999999","#CCCCCC"))




W_1 = readRDS("W_runs_alt_unif_bias_1_so .RDS")
W_2 = readRDS("W_runs_alt_unif_bias_2_so .RDS")
W_4 = readRDS("W_runs_alt_unif_bias_4_so .RDS")
W_6 = readRDS("W_runs_alt_unif_bias_8_so .RDS")
W_7 = readRDS("W_runs_alt_unif_bias_16_so .RDS")
W_1_true = c(0.34334899,
             0.03194600,
             0.28896252,
             0.23559938,
             0.09347495,
             0.38127813)
W_1_RMSE = matrix(0,100,6)
for(i in 1:100){
  for(j in 2000:5000){
    W_1_RMSE[i,] = W_1_RMSE[i,] + (W_1[[i]][j,] - W_1_true)^2
  }
  W_1_RMSE[i,] = sqrt(W_1_RMSE[i,]/3000)
}

W_2_RMSE = matrix(0,100,6)
for(i in 1:100){
  for(j in 2000:5000){
    W_2_RMSE[i,] = W_2_RMSE[i,] + (W_2[[i]][j,] - W_1_true)^2
  }
  W_2_RMSE[i,] = sqrt(W_2_RMSE[i,]/3000)
}

W_4_RMSE = matrix(0,100,6)
for(i in 1:100){
  for(j in 2000:5000){
    W_4_RMSE[i,] = W_4_RMSE[i,] + (W_4[[i]][j,] - W_1_true)^2
  }
  W_4_RMSE[i,] = sqrt(W_4_RMSE[i,]/3000)
}

W_6_RMSE = matrix(0,100,6)
for(i in 1:100){
  for(j in 2000:5000){
    W_6_RMSE[i,] = W_6_RMSE[i,] + (W_6[[i]][j,] - W_1_true)^2
  }
  W_6_RMSE[i,] = sqrt(W_6_RMSE[i,]/3000)
}

W_7_RMSE = matrix(0,100,6)
for(i in 1:100){
  for(j in 2000:5000){
    W_7_RMSE[i,] = W_7_RMSE[i,] + (W_7[[i]][j,] - W_1_true)^2
  }
  W_7_RMSE[i,] = sqrt(W_7_RMSE[i,]/3000)
}

((mean(apply(W_1_RMSE,2,median)) - mean(apply(W_2_RMSE,2,median)))/mean(apply(W_1_RMSE,2,median)))
((mean(apply(W_2_RMSE,2,median)) - mean(apply(W_4_RMSE,2,median)))/mean(apply(W_2_RMSE,2,median)))
((mean(apply(W_4_RMSE,2,median)) - mean(apply(W_6_RMSE,2,median)))/mean(apply(W_4_RMSE,2,median)))
((mean(apply(W_6_RMSE,2,median)) - mean(apply(W_7_RMSE,2,median)))/mean(apply(W_6_RMSE,2,median)))



u_1 = dim(W_1_RMSE)[1]
u_2 = dim(W_2_RMSE)[1]
u_3 = dim(W_4_RMSE)[1]
u_4 = dim(W_6_RMSE)[1]
u_5 = dim(W_7_RMSE)[1]
df_p = data.frame(x = c(as.vector(W_1_RMSE),as.vector(W_2_RMSE),as.vector(W_4_RMSE),
                        as.vector(W_6_RMSE),as.vector(W_7_RMSE)), 
                  K = factor(c(rep(1,u_1*5),rep(2,u_2*5),rep(4,u_3*5),rep(8,u_4*5),rep(16,u_5*5))),
                  t = factor(c(rep("1->2",u_1),rep("2->3",u_1),rep("3->4",u_1),rep("4->5",u_1),rep("5->6",u_1),
                               rep("1->2",u_2),rep("2->3",u_2),rep("3->4",u_2),rep("4->5",u_2),rep("5->6",u_2),
                               rep("1->2",u_3),rep("2->3",u_3),rep("3->4",u_3),rep("4->5",u_3),rep("5->6",u_3),
                               rep("1->2",u_4),rep("2->3",u_4),rep("3->4",u_4),rep("4->5",u_4),rep("5->6",u_4),
                               rep("1->2",u_5),rep("2->3",u_5),rep("3->4",u_5),rep("4->5",u_5),rep("5->6",u_5))))   
p_1 = ggplot(df_p, aes(x=t, y=x, fill= K)) + 
  geom_boxplot() +
  scale_fill_manual(values = c("grey30", "grey45", "grey65", "grey80","grey90")) + 
  labs( x="Primary periods", y = "RMSE V_i") +theme(axis.text = element_text(size = 13),
                                                    legend.text = element_text(size = 13),
                                                    legend.title = element_text(size = 13),
                                                    axis.title.x = element_text(size = 13))


w_1 = readRDS("abund_runs_BM_sim_alt_unif_bias_1_so .RDS")
w_2 = readRDS("abund_runs_BM_sim_alt_unif_bias_2_so .RDS")
w_4 = readRDS("abund_runs_BM_sim_alt_unif_bias_4_so .RDS")
w_6 = readRDS("abund_runs_BM_sim_alt_unif_bias_8_so .RDS")
w_7 = readRDS("abund_runs_BM_sim_alt_unif_bias_16_so .RDS")


abund_1_RMSE = rep(0,100)
for(i in 1:100){
  for(j in 1000:5000){
    abund_1_RMSE[i] = abund_1_RMSE[i] + (w_1[[i]][j] - 5000)^2
  }
  abund_1_RMSE[i] = sqrt(abund_1_RMSE[i]/4000)
}

abund_2_RMSE = rep(0,100)
for(i in 1:100){
  for(j in 1000:5000){
    abund_2_RMSE[i] = abund_2_RMSE[i] + (w_2[[i]][j] - 5000)^2
  }
  abund_2_RMSE[i] = sqrt(abund_2_RMSE[i]/4000)
}

abund_4_RMSE = rep(0,100) 
for(i in 1:100){
  for(j in 1000:5000){
    abund_4_RMSE[i] = abund_4_RMSE[i] + (w_4[[i]][j] - 5000)^2
  }
  abund_4_RMSE[i] = sqrt(abund_4_RMSE[i]/4000)
}


abund_6_RMSE = rep(0,100) 
for(i in 1:100){
  for(j in 1000:5000){
    abund_6_RMSE[i] = abund_6_RMSE[i] + (w_6[[i]][j] - 5000)^2
  }
  abund_6_RMSE[i] = sqrt(abund_6_RMSE[i]/4000)
}

abund_7_RMSE = rep(0,100) 
for(i in 1:100){
  for(j in 1000:5000){
    abund_7_RMSE[i] = abund_7_RMSE[i] + (w_7[[i]][j] - 5000)^2
  }
  abund_7_RMSE[i] = sqrt(abund_7_RMSE[i]/4000)
}

mean((median(abund_1_RMSE)-median(abund_2_RMSE))/median(abund_1_RMSE))
mean((median(abund_2_RMSE)-median(abund_4_RMSE))/median(abund_2_RMSE))
mean((median(abund_4_RMSE)-median(abund_6_RMSE))/median(abund_4_RMSE))
mean((median(abund_6_RMSE)-median(abund_7_RMSE))/median(abund_6_RMSE))


df_p = data.frame(x = c(as.vector(abund_1_RMSE),as.vector(abund_2_RMSE),as.vector(abund_4_RMSE),
                        as.vector(abund_6_RMSE),as.vector(abund_7_RMSE)), 
                  K = factor(c(rep(1,100),rep(2,100),rep(4,100),rep(8,100),rep(16,100))))   
p_1 = ggplot(df_p, aes( y=x, fill= K)) + 
  geom_boxplot() +
  scale_fill_manual(values = c("grey30", "grey45", "grey65", "grey80","grey90")) + 
  labs( y = "RMSE abund") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +theme(axis.text = element_text(size = 13),
                                                                              legend.text = element_text(size = 13),
                                                                              legend.title = element_text(size = 13),
                                                                              axis.title.x = element_text(size = 13)) 




p_1 = readRDS("p_runs_BM_sim_alt_unif_bias_1_so .RDS")
p_2 = readRDS("p_runs_BM_sim_alt_unif_bias_2_so .RDS")
p_4 = readRDS("p_runs_BM_sim_alt_unif_bias_4_so .RDS")
p_6 = readRDS("p_runs_BM_sim_alt_unif_bias_8_so .RDS")
p_7 = readRDS("p_runs_BM_sim_alt_unif_bias_16_so .RDS")
RMSE_p_1 = matrix(0,100,6)
RMSE_p_2 = matrix(0,100,12)
RMSE_p_4 = matrix(0,100,24)
RMSE_p_6 = matrix(0,100,48)
RMSE_p_7 = matrix(0,100,96)
p_true_1 = c(0.2174925, 0.5419413, 0.4859932, 0.5113957, 0.3730587, 0.2683062)
p_true_2 = rep(p_true_1,each=2)
p_true_4 = rep(p_true_1,each=4)
p_true_6 = rep(p_true_1,each=8)
p_true_7 = rep(p_true_1,each=16)

for(j in 1:100){
  for(i in 1:6){
    RMSE_p_1[j,i] =   mean((p_1[[j]][1000:5000,i] - p_true_1[i])^2)
  }
}
for(j in 1:100){
  for(i in 1:12){
    RMSE_p_2[j,i] =   mean((p_2[[j]][1000:5000,i] - p_true_2[i])^2)
  }
}  
for(j in 1:100){    
  for(i in 1:24){
    RMSE_p_4[j,i] =   mean((p_4[[j]][1000:5000,i] - p_true_4[i])^2)
  }
}  
for(j in 1:100){    
  for(i in 1:48){
    RMSE_p_6[j,i] =   median((p_6[[j]][1000:5000,i] - p_true_6[i])^2)
  }
}
for(j in 1:100){    
  for(i in 1:90){
    RMSE_p_7[j,i] =   median((p_7[[j]][1000:5000,i] - p_true_7[i])^2)
  }
}

(median(apply(RMSE_p_1,2,median))-median(apply(RMSE_p_2,2,median)))/median(apply(RMSE_p_1,2,median))
(median(apply(RMSE_p_2,2,median))-median(apply(RMSE_p_4,2,median)))/median(apply(RMSE_p_2,2,median))
(median(apply(RMSE_p_4,2,median))-median(apply(RMSE_p_6,2,median)))/median(apply(RMSE_p_4,2,median))
(median(apply(RMSE_p_6,2,median))-median(apply(RMSE_p_7,2,median)))/median(apply(RMSE_p_6,2,median))


RMSE_p_1 = sqrt(RMSE_p_1)
df_RMSE = data.frame(y = (as.vector(RMSE_p_1)), 
                     t = factor(c(rep(1,100),rep(2,100),rep(3,100),rep(4,100),
                                  rep(5,100),rep(6,100))) )

pp1 = ggplot(df_RMSE, aes( y = y, x = t,fill=t)) +  # ggplot function
  geom_boxplot(show.legend = FALSE)  + ylab("p") + xlab("") +
  scale_y_continuous(limits = c(0, 0.11)) +
  scale_fill_manual(values = c(rep("grey30",6))) +theme(axis.text = element_text(size = 13),
                                                        legend.text = element_text(size = 13),
                                                        legend.title = element_text(size = 13),
                                                        axis.title.x = element_text(size = 13)) 


x_axis_2 = c()
for(i in 1:6){
  for(j in 1:2){
    x_axis_2 = c(x_axis_2,paste0(i,".",j))
  }
}
RMSE_p_2 = sqrt(RMSE_p_2)
df_RMSE = data.frame(y = (as.vector(RMSE_p_2)), 
                     t = factor(x_axis_2) )

pp2 = ggplot(df_RMSE, aes( y = y, x = t,fill=t)) +  # ggplot function
  geom_boxplot(show.legend = FALSE)  + ylab("")+ xlab("") +
  scale_y_continuous(limits = c(0, 0.11))+
  scale_fill_manual(values = c(rep("grey45",12))) +theme(axis.text = element_text(size = 13),
                                                         legend.text = element_text(size = 13),
                                                         legend.title = element_text(size = 13),
                                                         axis.title.x = element_text(size = 13)) 

x_axis_4 = c()
for(i in 1:6){
  for(j in 1:4){
    x_axis_4 = c(x_axis_4,paste0(i,".",j))
  }
}
RMSE_p_4 = sqrt(RMSE_p_4)
df_RMSE = data.frame(y = (as.vector(RMSE_p_4)), 
                     t = factor(x_axis_4) )

pp3 = ggplot(df_RMSE, aes( y = y, x = t,fill=t)) +  # ggplot function
  geom_boxplot(show.legend = FALSE)  + ylab("p")+ 
  scale_y_continuous(limits = c(0, 0.11))+
  scale_fill_manual(values = c(rep("grey65",24))) + labs(x = NULL)+theme(axis.text = element_text(size = 13),
                                                                         legend.text = element_text(size = 13),
                                                                         legend.title = element_text(size = 13),
                                                                         axis.title.x = element_text(size = 13))

x_axis_6 = c()
for(i in 1:6){
  for(j in 1:8){
    x_axis_6 = c(x_axis_6,paste0(i,".",j))
  }
}
RMSE_abund_6 = sqrt(RMSE_p_6)
df_RMSE = data.frame(y = (as.vector(RMSE_p_6)), 
                     t = factor(x_axis_6) )
levels(df_RMSE[,2])  = x_axis_6
pp5 = ggplot(df_RMSE, aes( y = y, x = t,fill=t)) +  
  geom_boxplot(show.legend = FALSE) + scale_x_discrete(breaks = df_RMSE[c(seq(2,48,by=4)) ,2]) + ylab("p")+ 
  scale_y_continuous(limits = c(0, 0.11))+ ylab("")+
  scale_fill_manual(values = c(rep("grey80",48))) + xlab(" ")+theme(axis.text = element_text(size = 13),
                                                                    legend.text = element_text(size = 13),
                                                                    legend.title = element_text(size = 13),
                                                                    axis.title.x = element_text(size = 13))

x_axis_7 = c()
for(i in 1:6){
  for(j in 1:16){
    x_axis_7 = c(x_axis_7,paste0(i,".",j))
  }
}
RMSE_abund_7 = sqrt(RMSE_p_7)
df_RMSE = data.frame(y = (as.vector(RMSE_p_7)), 
                     t = factor(x_axis_7) )
levels(df_RMSE[,2])  = x_axis_7
pp7 = ggplot(df_RMSE, aes( y = y, x = t,fill=t)) +  # ggplot function
  geom_boxplot(show.legend = FALSE) + scale_x_discrete(breaks = df_RMSE[c(seq(2,96,by=16),seq(1,96,by=16)) ,2]) + ylab("p")+ 
  scale_y_continuous(limits = c(0, 0.11))+ ylab("")+
  scale_fill_manual(values = c(rep("grey90",96))) + xlab("secondary sampling occasions per primary period")+
  theme(axis.text = element_text(size = 13),    legend.text = element_text(size = 13),
                    legend.title = element_text(size = 13),
                  axis.title.x = element_text(size = 13))

grid.arrange(pp1,pp2,pp3,pp5,pp7,nrow=5)


apply(RMSE_p_1,2,mean)
RMSE_2 = c()
u = 0 
for(i in 1:6){
  RMSE_2[i] = sum(apply(RMSE_p_2,2,median)[(1+u):(u+2)])/2
  u = u + 2
}
RMSE_3 = c()
u = 0 
for(i in 1:6){
  RMSE_3[i] = sum(apply(RMSE_p_4,2,median)[(1+u):(u+4)])/4
  u = u + 4
}
RMSE_4 = c()
u = 0 
for(i in 1:6){
  RMSE_4[i] = sum(apply(RMSE_p_15,2,median)[(1+u):(u+15)])/15
  u = u + 15
}

RMSE_1 = apply(RMSE_p_1,2,median)
RMSE_2
RMSE_3
RMSE_4
mean((RMSE_1-RMSE_2)/RMSE_1)
mean((RMSE_2-RMSE_3)/RMSE_2)
mean((RMSE_3-RMSE_4)/RMSE_3)
