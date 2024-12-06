gibbs_cap_recap <- function(n, a, b, lambda, r, vetor_ns){
  aa <- matrix(nrow = n, ncol = 15) #amostra da posteriori
  colnames(aa) <- c(paste0("p_", 1:14), "N")
  while(T){ N1 <- rpois(1, lambda); if(N1 >= r) break}
  aa[1,] <- c(rbeta(14, a, b), N1) #chutes iniciais
  for(linha in 2:n){ #vou amostrar p_1,p_2,...,p_14,N
    pis <- purrr::map(1:14, function(i) rbeta(1, a + vetor_ns[i],
                              b + aa[linha-1, 15] - vetor_ns[i]))
    pis <- purrr::list_c(pis)
    N_i <- r + rpois(1, lambda*prod(1-pis))
    aa[linha, ] <- c(pis, N_i)
  }; return(aa)}
library(purrr)
v = c(2, 4, 2, 4, 4, 2, 3, 4, 2, 3, 5, 10, 2, 5);
table(v); par(mfrow = c(2,2), mar = c(4,4.1,1,2))
set.seed(236202);
aas_unif <- map(1:3, function(x) gibbs_cap_recap(1000, 1, 1, 30, 12, v))
plot(1:1000, aas_unif[[1]][,1], xlab ="Índice da amostra", ylab = "p1", cex.lab = 1.5, cex.axis = 1.5)
for(i in 2:3) points(1:1000, aas_unif[[i]][,1], col = i) #convergiu
plot(1:1000, aas_unif[[1]][,12], xlab ="Índice da amostra", ylab = "p12", cex.lab = 1.5, cex.axis = 1.5)
for(i in 2:3) points(1:1000, aas_unif[[i]][,12], col = i) #convergiu
plot(1:1000, aas_unif[[1]][,11], xlab ="Índice da amostra", ylab = "p11", cex.lab = 1.5, cex.axis = 1.5)
for(i in 2:3) points(1:1000, aas_unif[[i]][,11], col = i) #convergiu
plot(1:1000, aas_unif[[1]][,15], xlab ="Índice da amostra", ylab = "N", cex.lab = 1.5, cex.axis = 1.5)
for(i in 2:3) points(1:1000, aas_unif[[i]][,15], col = i) #convergiu
aa_unif <- gibbs_cap_recap(11000, 1, 1, 30, 12, v)[1001:11000,]
acf(aa_unif[,1], ylab = "ACF p1", cex.lab = 1.5, cex.axis = 1.5); acf(aa_unif[,12], ylab = "ACF p12", cex.lab = 1.5, cex.axis = 1.5) #nao correlacionados
acf(aa_unif[,11], ylab = "ACF p11", cex.lab = 1.5, cex.axis = 1.5); acf(aa_unif[,15], ylab = "ACF N", cex.lab = 1.5, cex.axis = 1.5) #vou pegar d = 2
aa_unif_f <- gibbs_cap_recap(11000, 1, 1, 30, 12, v)[seq(1001,11000,by=2),]
# acf(aa_unif_f[,15]) #acabou a correlacao
par(mfrow = c(2,3), mar = c(4,4.1,1,2))
hist(aa_unif_f[,1], main = "", probability = T, ylab = "Densidade", xlab = "Valor de p1", cex.lab = 1.5, cex.axis = 1.5)
hist(aa_unif_f[,12], main = "", probability = T, ylab = "Densidade", xlab = "Valor de p12", cex.lab = 1.5, cex.axis = 1.5)
hist(aa_unif_f[,11], main = "", probability = T, ylab = "Densidade", xlab = "Valor de p11", cex.lab = 1.5, cex.axis = 1.5)
hist(aa_unif_f[,15], main = "", probability = T, ylab = "Densidade", xlab = "Valor de N", cex.lab = 1.5, cex.axis = 1.5)
plot(jitter(aa_unif_f[,15]), aa_unif_f[,1], ylab = "p1", xlab = "N + ruído", cex.lab = 1.5, cex.axis = 1.5)
plot(jitter(aa_unif_f[,15]), aa_unif_f[,12], ylab = "p12", xlab = "N + ruído", cex.lab = 1.5, cex.axis = 1.5)

cumsum(table(aa_unif_f[,15])/5000)
plot(seq(0,1,length.out=1000), dbeta(seq(0,1,length.out=1000),1,13))
1/14; 1 - pbeta(0.2, 1, 13)
set.seed(236203)
par(mfrow = c(2,2), mar = c(4,4.1,1,2))
aas_b_ab <- map(1:3, function(x) gibbs_cap_recap(1000, 1, 13, 30, 12, v))
plot(1:1000, aas_b_ab[[1]][,1], xlab ="Índice da amostra", ylab = "p1", cex.lab = 1.5, cex.axis = 1.5)
for(i in 2:3) points(1:1000, aas_b_ab[[i]][,1], col = i) #convergiu
plot(1:1000, aas_b_ab[[1]][,12], xlab ="Índice da amostra", ylab = "p12", cex.lab = 1.5, cex.axis = 1.5)
for(i in 2:3) points(1:1000, aas_b_ab[[i]][,12], col = i) #convergiu
plot(1:1000, aas_b_ab[[1]][,11], xlab ="Índice da amostra", ylab = "p11", cex.lab = 1.5, cex.axis = 1.5)
for(i in 2:3) points(1:1000, aas_b_ab[[i]][,11], col = i) #convergiu
plot(1:1000, aas_b_ab[[1]][,15], xlab ="Índice da amostra", ylab = "N", cex.lab = 1.5, cex.axis = 1.5)
for(i in 2:3) points(1:1000, aas_b_ab[[i]][,15], col = i) #convergiu
aa_b_ab <- gibbs_cap_recap(11000, 1, 13, 30, 12, v)[1001:11000,]
acf(aa_b_ab[,1], ylab = "ACF p1", cex.lab = 1.5, cex.axis = 1.5); acf(aa_b_ab[,12], ylab = "ACF p12", cex.lab = 1.5, cex.axis = 1.5) #correlacionados
acf(aa_b_ab[,11], ylab = "ACF p11", cex.lab = 1.5, cex.axis = 1.5); acf(aa_b_ab[,15], ylab = "ACF N", cex.lab = 1.5, cex.axis = 1.5) #vou pegar d = 4
aa_b_ab_f <- gibbs_cap_recap(16000, 1, 13, 30, 12, v)[seq(1001,16000,by=3),]
#acf(aa_b_ab_f[,1]); acf(aa_b_ab_f[,12]); #acabou a correlacao
#acf(aa_b_ab_f[,11]); acf(aa_b_ab_f[,15]) #acabou a correlacao
par(mfrow = c(2,3), mar = c(4,4.1,1,2))

hist(aa_b_ab_f[,1], main = "", probability = T, ylab = "Densidade", xlab = "Valor de p1", cex.lab = 1.5, cex.axis = 1.5)
hist(aa_b_ab_f[,12], main = "", probability = T, ylab = "Densidade", xlab = "Valor de p12", cex.lab = 1.5, cex.axis = 1.5)
hist(aa_b_ab_f[,11], main = "", probability = T, ylab = "Densidade", xlab = "Valor de p11", cex.lab = 1.5, cex.axis = 1.5)
hist(aa_b_ab_f[,15], main = "", probability = T, ylab = "Densidade", xlab = "Valor de N", cex.lab = 1.5, cex.axis = 1.5)
plot(jitter(aa_b_ab_f[,15]), aa_b_ab_f[,1], ylab = "p1", xlab = "N + ruído", cex.lab = 1.5, cex.axis = 1.5)
plot(jitter(aa_b_ab_f[,15]), aa_b_ab_f[,12],  ylab = "p12", xlab = "N + ruído", cex.lab = 1.5, cex.axis = 1.5)
cumsum(table(aa_b_ab_f[,15])/5000)
