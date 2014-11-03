s4vdpca
=======

This R package implements methods for sparse principal component analysis as descriped in the manuscript
"Applying Stability Selection to Consistently Estimate Sparse Principal Components in High-Dimensional Molecular Data" submitted to Oxford Bioinformatics.

```{r}
# install package using devtools
# install.packages('devtools')
library(devtools)                  
install_github('mwsill/s4vdpca')
library(s4vdpca)

# generate a simulated data set using the single-covariance spike model 
p <- 5000    # number of variables
n <- 50      # number of observations
alpha <- .6  # spike index 
beta <- .5   # sparsity index 

# generate a population variance covariance matrix
Sigma <- generate_covar(alpha,beta,p)

# extract first eigenvector
z1 <- Sigma[[2]]

# extract variance covariance matrix
Sigma <- Sigma[[1]]

# sample from multivariate normal distribution using Cholesky decomposition
# see ?rmvn in package broman for details
D <- chol(Sigma)
set.seed(02112014)
x <- matrix(rnorm(n * p), ncol = p) %*% D + rep(rep(0,p), rep(n, p))

#show documentation
?s4vdpca
?rspca

# apply S4VDPCA and RSPCA with different penalization functions, all with GIC5 
# parallelization is not yet available on Windows machines
res1 <- s4vdpca(x, center=TRUE, cores=1, ic_type='gic5')
res2 <- rspca(x, center=TRUE, cores=1, ic_type='gic5') #lasso
res3 <- rspca(x, center=TRUE, cores=1, ic_type='gic5', type='scad') #scad 
res4 <- rspca(x, center=TRUE, cores=1, ic_type='gic5', gamv=1) # adaptive lasso

# plot the information criterion
par(mfrow=c(2,2))
plot(res1$ic, xlab='number of selected features', ylab='GIC 5'
,main='S4VDPCA')
abline(v=res1$minic, col='red')
text(y=max(res1$ic,na.rm=T)-1000,x=res1$minic+100,res1$minic,col='red')
plot(res2$ic, xlab='number of selected features', ylab='GIC 5'
,main='RSPCA lasso')
abline(v=res2$minic, col='red')
text(y=max(res2$ic,na.rm=T)-1000,x=res2$minic+100,res2$minic,col='red')
plot(res3$ic, xlab='number of selected features', ylab='GIC 5'
,main='RSPCA scad')
abline(v=res3$minic, col='red')
text(y=max(res3$ic,na.rm=T)-1000,x=res3$minic+100,res3$minic,col='red')
plot(res4$ic, xlab='number of selected features', ylab='GIC 5'
,main='RSPCA adaptive lasso')
abline(v=res4$minic, col='red')
text(y=max(res4$ic,na.rm=T)-1000,x=res4$minic+100,res4$minic,col='red')

# calculate angle between estimated sparse loadings vector and simulated eigenvector
angle(res1$v,z1)
angle(res2$v,z1)
angle(res3$v,z1)
angle(res4$v,z1)

# calculate number of falsely selected features
type1(z1,res1$v)
type1(z1,res2$v)
type1(z1,res3$v)
type1(z1,res4$v)

# apply regular PCA and calculate angle between loadings vector
# and simulated eigenvector
pca <- prcomp(x)
angle(pca$rotation[,1],z1)

# ssvdpca is the original rspca function by Lee et al. 2010
X  <- scale(x, center=TRUE)
system.time(
res5 <- ssvdpca(X) #lasso
)
# optimized code; search for minimal bic
system.time(
res6 <- rspca(X, center=FALSE, cores=1,steps=100, ic_type='bic') #lasso
)
# optimized code; parallelized search for minimal bic, only on Unix machines
system.time(
res7 <- rspca(X, center=FALSE, cores=4,steps=100, ic_type='bic') #lasso
)
# estimated loadings are the same 
all(res5$v==res6$v)
all(res5$v==res7$v)


###### real data application example
# to perform the analysis the follwing Bioconductor packages are needed
#source("http://bioconductor.org/biocLite.R")
#biocLite(c("Biobase","HTSanalyzeR","org.Hs.eg.db","KEGG.db"))

# load packages
library(Biobase)
library(KEGG.db)
library(org.Hs.eg.db)
library(HTSanalyzeR)

# load medulloblastoma example data set 
data(medullo)

# data is stored as ExpressionSet
medullo

# extract gene expression matrix
X <- scale(t(exprs(medullo)),scale=F)
dim(X)

# estimate first sparse PC using S4VDPCA
set.seed(10102014,"L'Ecuyer-CMRG")
res.s4vdpca <- s4vdpca(X,center=FALSE,cores=1,steps=500,B=1000) 
# subtract first sparse rank1 layer
X <- X - res.s4vdpca$d * res.s4vdpca$u %*% t(res.s4vdpca$v) 
# obtain second sparse PC
res.s4vdpca2 <- s4vdpca(X,center=FALSE,cores=1,steps=500,B=1000)

# perform gene set over representation analysis for PC1
t_all <- res.s4vdpca$v
# get entrez gene ids
names(t_all) <- featureData(medullo)$GENE
# define selected genes
t_hits <- as.character(featureData(medullo)$GENE[res.s4vdpca$v!=0])
t_all <- sort(t_all)
KEGG_pathways <- KeggGeneSets(species = "Hs") 
gscs.kegg <- list("KEGG pathways" = KEGG_pathways) 
gsca.kegg <- new("GSCA", listOfGeneSetCollections=gscs.kegg,
                 geneList=t_all, hits=t_hits)
gsca.kegg <- preprocess(gsca.kegg, species="Hs",
                        initialIDs="Entrez.gene",
                        keepMultipleMappings=TRUE,
                        duplicateRemoverMethod="max",
                        orderAbsValue=FALSE)
gsca.kegg <- analyze(gsca.kegg, para=list(pValueCutoff=0.05,
                                          pAdjustMethod ="BH",
                                          nPermutations=5,
                                          minGeneSetSize=20,
                                          exponent=1),
                      verbose=T, doGSEA=TRUE)
gsca.kegg <- appendGSTerms(gsca.kegg,
                           keggGSCs="KEGG pathways")
length(getTopGeneSets(gsca.kegg,"HyperGeo.results",
                      "KEGG pathways",allSig=TRUE)[[1]])
hypPC1 <- gsca.kegg@result$HyperGeo.results$`KEGG pathways`
g1 <- viewEnrichMap(gsca.kegg, resultName="HyperGeo.results", gscs="KEGG pathways",
                    ntop=6, allSig=FALSE, gsNameType="term", displayEdgeLabel=FALSE,
                    layout="layout.kamada.kawai",plot=FALSE)
                   
# perform gene set over representation analysis for PC2
t_all <- res.s4vdpca2$v
# get entrez gene ids
names(t_all) <- featureData(medullo)$GENE
# define selected genes
t_hits <- as.character(featureData(medullo)$GENE[res.s4vdpca2$v!=0])
t_all <- sort(t_all)
KEGG_pathways <- KeggGeneSets(species = "Hs") 
gscs.kegg <- list("KEGG pathways" = KEGG_pathways) 
gsca.kegg <- new("GSCA", listOfGeneSetCollections=gscs.kegg,
                 geneList=t_all, hits=t_hits)
gsca.kegg <- preprocess(gsca.kegg, species="Hs",
                        initialIDs="Entrez.gene",
                        keepMultipleMappings=TRUE,
                        duplicateRemoverMethod="max",
                        orderAbsValue=FALSE)
gsca.kegg <- analyze(gsca.kegg, para=list(pValueCutoff=0.05,
                                          pAdjustMethod ="BH",
                                          nPermutations=5,
                                          minGeneSetSize=20,
                                          exponent=1),
                      verbose=T, doGSEA=TRUE)
gsca.kegg <- appendGSTerms(gsca.kegg,
                           keggGSCs="KEGG pathways")
length(getTopGeneSets(gsca.kegg,"HyperGeo.results",
                      "KEGG pathways",allSig=TRUE)[[1]])
hypPC2 <- gsca.kegg@result$HyperGeo.results$`KEGG pathways`
g2 <- viewEnrichMap(gsca.kegg, resultName="HyperGeo.results", gscs="KEGG pathways",
                    ntop=6, allSig=FALSE, gsNameType="term", displayEdgeLabel=FALSE,
                    layout="layout.kamada.kawai",plot=FALSE)


###### Biplot with dominant marker genes and KEGG pathways

# assign colors to molecular subtypes
cols <- rep("gray",length(medullo$subtype))
cols[medullo$subtype=="group D"] <- "darkgreen"
cols[medullo$subtype=="group C"] <- "orange"
cols[medullo$subtype=="WNT"] <- 'blue' 
cols[medullo$subtype=="SHH"] <- "darkred"
# assign pch type to molecular subtypes
ptype <- rep(1,length(medullo$subtype))
ptype[medullo$subtype=="group D"] <- 17
ptype[medullo$subtype=="group C"] <- 15
ptype[medullo$subtype=="WNT"] <- 18
ptype[medullo$subtype=="SHH"] <- 20
ptype<- as.numeric(ptype)

agi2symbol <- featureData(medullo)$GENE_SYMBOL

alpha <- 0.37
ptx <- res.s4vdpca$u*(res.s4vdpca$d^alpha)
pty <- res.s4vdpca2$u*(res.s4vdpca2$d^alpha)
arrowsx <- res.s4vdpca$v*(res.s4vdpca$d^(1-alpha))
arrowsy <- res.s4vdpca2$v*(res.s4vdpca2$d^(1-alpha))

cand <- c(order(abs(arrowsx))[1:10],order(abs(arrowsy))[1:10])

ord <- order(apply(cbind(abs(arrowsx),abs(arrowsy)),1,max),decreasing=T)

cand <- 1:20

p.cutoff.labels<-rep("",4)
p.cutoff.labels[c(1,4,6,9)]<-c(0,0.01,0.05,1)    

nf <- layout(matrix(c(2,4,1,3),2,2,byrow = TRUE), c(2.5,1.75), c(1.5,2.5), TRUE)
par(mar = c(4,4,.5,.5)) #(bottom, left, top, right)
plot(ptx,pty,type="n",xlim=c(-3.25,3.25),ylim=c(-3.25,3.25),axes=T
,ylab=paste("PC2 /",sum(arrowsy!=0),"genes",sep=" ")
     ,xlab=paste("PC1 /",sum(arrowsx!=0),"genes",sep=" "),main='',lty=2)
arrows(0, 0, x1 = arrowsx, y1=arrowsy,length=0.05,col = "gray")
par(xpd=NA)
for(i in cand){
  arrows(0, 0, x1 = arrowsx[ord[i]], y1=arrowsy[ord[i]],length=0.05,col = "red")
   text(x=arrowsx[ord[i]]+ sign(arrowsx[ord[i]])*
   .2,y=arrowsy[ord[i]]+sign(arrowsy[ord[i]])*.2,
         labels=agi2symbol[ord[i]],col="red",cex=.75)
  }
sec <- which(cols=="darkgreen"|cols=="blue")
points(ptx,pty,col=cols,pch=ptype,cex=1)
points(ptx[sec],pty[sec],col=cols[sec],pch=ptype[sec],cex=1.25)
par(mar = c(0,3,3,3)) #(bottom, left, top, right)
V(g1)$label.cex <- .85
V(g2)$label.cex <- .85
set.seed(1235)
plot(g1)
par(mar = c(0,3,0,3.5)) #(bottom, left, top, right)
set.seed(1235)
plot(g2)
plot(c(-1,1),c(-1,1),type="n",axes=F,ylab="",xlab="")
points(
  x = rep(-.6, 9), 
  y = seq(-.2, (-.2-(0.05*9)), 
          length.out = 9), 
  pch = 15, col = c("#FF0000","#FF1F1F","#FF3F3F",
                    "#FF5F5F","#FF7F7F","#FF9F9F",
                    "#FFBFBF", "#FFDFDF", "#FFFFFF")
)

text(
  x = rep(-.677, 9), 
  y = seq(-.2, (-.2-(0.05*9)), 
          length.out = 9),
  labels = p.cutoff.labels, 
  cex = .7,
  adj=1
)
text(
  x = -.65,
  y = 0.02,
  labels = "Adjusted\np-values",
  cex=.7,
  adj= 0.5,
  font=2
)

legend(x=0,y=-.05,legend=c('WNT','SHH','group C','group D')
       ,col=c('blue','darkred','orange','darkgreen')
       ,border="black",pch=c(18,20,15,17),cex=1,bty="n")
text(
  x = 0.35,
  y = -.0225,
  labels = "subgroup",
  cex=.7,
  adj= 0.5,
  font=2
)
```

