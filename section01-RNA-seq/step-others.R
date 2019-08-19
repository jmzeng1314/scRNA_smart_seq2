load(file = '../input.Rdata')
dat[1:4,1:4] ## log2(cpm+1)
# hc.logCPM=hclust(dist(t( dat ))) 
 
load(file = '../input_rpkm.Rdata')
# hc.logRPKM=hclust(dist(t(log(dat+0.1)))) 
# hc.RPKM=hclust(dist(t( dat ))) 
#save(hc.logCPM,hc.logRPKM,hc.RPKM,file = 'hc_3.Rdata')

load(file = 'hc_3.Rdata')
plot(hc.logRPKM,labels = F)
plot(hc.RPKM,labels = F)
plot(hc.logCPM,labels = F)

g1 = cutree(hc.logRPKM, 4);table(g1)
g2 = cutree(hc.RPKM, 4)  ;table(g2)
g3 = cutree(hc.logCPM, 4)  ;table(g3)
table(g1,g3)

library(Rtsne) 
if(T){ 
  load(file = '../input_rpkm.Rdata') 
  dat[1:4,1:4]
  dat_matrix <- t(dat) 
  dat_matrix =log2(dat_matrix+0.01) 
}
if(T){
  load(file = '../input.Rdata')
  dat[1:4,1:4]
  dat_matrix <- t(dat) 
}
set.seed(42)
tsne_out <- Rtsne(dat_matrix,pca=FALSE,
                  perplexity=27,theta=0.5)  
par(mfrow=c(1,3))
plot(tsne_out$Y,col= g1,sub = 'hc.logRPKM')  
plot(tsne_out$Y,col= g2,sub = 'hc.RPKM') 
plot(tsne_out$Y,col= g3,sub = 'hc.logCPM')  
par(mfrow=c(1,1))
















