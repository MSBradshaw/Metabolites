suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringr))
library(ggplot2)

g_data = t(read.table("impute_genotype.txt", header=TRUE, sep="\t", row.names = 1, check.names=FALSE,stringsAsFactors=FALSE))
#Phenotype data
pheno = read.table("Leaf-Pheno-2018-MAGIC-GH.txt", header=TRUE, sep="\t", check.names=FALSE,stringsAsFactors=FALSE)
maternal = read.table('maternal-E-F-generations.tsv',header=TRUE, sep="\t", check.names=FALSE,stringsAsFactors=FALSE)

head(g_data)
head(pheno)

pheno <- pheno[pheno$ID!='blank',]
maternal$L <- gsub('-','_',maternal$L)
rownames(pheno) <- pheno$ID

com_id = intersect(colnames(g_data),rownames(pheno))
#Arrange the matrices accordingly
g_data <- g_data[,com_id]
pheno <- pheno[com_id,]
rownames(maternal) <- maternal$L
mat <- maternal[com_id,]$Mother

dim(pheno)
dim(g_data)
dim(mat)

pheno$Bergamotene <- as.numeric(gsub(',','.',pheno$Bergamotene))
pheno$Linalool <- as.numeric(gsub(',','.',pheno$Linalool))
pheno$LeafMass <- as.numeric(gsub(',','.',pheno$LeafMass))

metabolite <- "Bergamotene"

phenotype = as.numeric(pheno[,metabolite])

#normalize the data and 
phenotype <- log(phenotype)
phenotype[is.infinite(phenotype)] <- 0

lod_RIL = apply(g_data, 1, function(x){
  dat = cbind(1,suppressWarnings(as.factor(x)),phenotype)
  dat = dat[complete.cases(dat), ]
  fit1 = .lm.fit(dat[,1:2], dat[,3], tol=1e-12)
  fit0 = .lm.fit(as.matrix(dat[,1]),dat[,3], tol=1e-12)
  rss0 = sum(fit0$residuals^2)
  rss1 = sum(fit1$residuals^2)
  (dim(dat)[1]/2)*(log10(rss0)-log10(rss1))
})

n = 100
perm = sapply(1:n, function(i){
  if(i%%10==0){
    message(paste0("Bootstrapping ",i))
  }
  pheno_sample = sample(phenotype)
  lod_perm = apply(g_data, 1, function(x){
    dat = cbind(1,suppressWarnings(as.factor(x)),pheno_sample)
    dat = dat[complete.cases(dat), ]
    fit1 = .lm.fit(dat[,1:2], dat[,3], tol=1e-12)
    fit0 = .lm.fit(as.matrix(dat[,1]),dat[,3], tol=1e-12)
    rss0 = sum(fit0$residuals^2)
    rss1 = sum(fit1$residuals^2)
    (dim(dat)[1]/2)*(log10(rss0)-log10(rss1))
  })
  lod_perm[which(lod_perm == max(lod_perm,na.rm=TRUE))]
})

threshold = quantile(as.numeric(unlist(perm)), 0.95, na.rm = TRUE)

min_marker = names(which(lod_RIL == max(lod_RIL,na.rm=TRUE)))
effect_data = data.frame(accession = com_id, genotype=as.factor(g_data[min_marker,]),pheno=phenotype)
effect_data = effect_data[complete.cases(effect_data), ]
effect_data$genotype = as.factor(effect_data$genotype)

ggplot(effect_data, aes(x=genotype,y=pheno)) +
  geom_boxplot(width=0.2, outlier.shape=NA) +
  geom_jitter(width=0.1, size=0.5) +
  ggtitle(paste0(min_marker))
d = quantile(as.numeric(unlist(perm)), 0.95, na.rm = TRUE)

# Some preprocessing of the data
marker = as.data.frame(str_split_fixed(rownames(g_data), "_", 3))
marker$V1 = as.data.frame(str_split_fixed(marker$V1, "-", 2))$V1
marker$V1 = factor(marker$V1, levels=c(paste0("chr",1:12),"scaffold"))
marker$V2 = 1:dim(marker)[1]
marker = marker[,-3]
plotdata = cbind(marker, lod_RIL)
colnames(plotdata) = c("chrom","pos","lod")
head(plotdata)

#Make a pretty QTL plot
p <- ggplot(plotdata, aes(pos, lod, group=chrom,color=chrom)) +
  geom_point(alpha=0.7) +
  geom_line(alpha=0.7) +
  facet_wrap(~ chrom, ncol = 13, scales = "free_x", strip.position="bottom") +
  geom_hline(yintercept = threshold, linetype = "dotted", color='red') +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.spacing.x = unit(0,"line"),
        axis.text.x=element_blank(),strip.placement = "outside",legend.position = "none")+
  labs(x = "Chromosome", y = "LOD", color = "", linetype = "")
p
ggsave(paste(metabolite,'-QTL-log.png',sep=''))

#---------------------CO VARIATE--------------------------------------------------------------------------
#-------------------------------CO VARIATE----------------------------------------------------------------
#-----------------------------------------CO VARIATE------------------------------------------------------

#--0000---00---0-----0----0----0000--0----0---00000-0000-------0----0--0----0----0----000-000-0---0-000-------
#-0------0--0---0---0----0-0---0--0--0---0-0----0---00--------0-0---00-0---0-0---0-----0--0----0-0--0---------
#-0------0--0----0-0----0-0-0--0-0---0--0-0-0---0---0--------0-0-0--0-00--0-0-0--0-----0----0---0-----0-------
#--0000---00------0----0-----0-0--0--0-0-----0--0---0000----0-----0-0--0-0-----0-0000-000-000---0---000-------

#Now do the Co-variate analysis
lod_RIL = apply(g_data, 1, function(x){
  dat = cbind(1,suppressWarnings(as.factor(x)),as.factor(mat),phenotype)
  dat = dat[complete.cases(dat), ]
  fit1 = .lm.fit(dat[,1:3], dat[,4], tol=1e-12)
  fit0 = .lm.fit(as.matrix(dat[,c(1,3)]),dat[,4], tol=1e-12)
  rss0 = sum(fit0$residuals^2)
  rss1 = sum(fit1$residuals^2)
  (dim(dat)[1]/2)*(log10(rss0)-log10(rss1))
})

n = 100
perm = sapply(1:n, function(i){
  if(i%%10==0){
    message(paste0("Bootstrapping ",i))
  }
  pheno_sample = sample(phenotype)
  lod_perm = apply(g_data, 1, function(x){
    dat = cbind(1,suppressWarnings(as.factor(x)),as.factor(mat),phenotype)
    dat = dat[complete.cases(dat), ]
    fit1 = .lm.fit(dat[,1:3], dat[,4], tol=1e-12)
    fit0 = .lm.fit(as.matrix(dat[,c(1,3)]),dat[,4], tol=1e-12)
    rss0 = sum(fit0$residuals^2)
    rss1 = sum(fit1$residuals^2)
    (dim(dat)[1]/2)*(log10(rss0)-log10(rss1))
  })
  lod_perm[which(lod_perm == max(lod_perm,na.rm=TRUE))]
})

marker = as.data.frame(str_split_fixed(rownames(g_data), "_", 3))
marker$V1 = as.data.frame(str_split_fixed(marker$V1, "-", 2))$V1
marker$V1 = factor(marker$V1, levels=c(paste0("chr",1:12),"scaffold"))
marker$V2 = 1:dim(marker)[1]
marker = marker[,-3]
plotdata = cbind(marker, lod_RIL)
colnames(plotdata) = c("chrom","pos","lod")
head(plotdata)

q <- ggplot(data=plotdata, mapping=aes(x=pos, y=lod, group=chrom,color=chrom)) +
  geom_point(alpha=0.7) +
  geom_line(alpha=0.7) +
  facet_wrap(~ chrom, ncol = 13, scales = "free_x", strip.position="bottom") +
  geom_hline(yintercept = threshold, linetype = "dotted", color='red') +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.spacing.x = unit(0,"line"),
        axis.text.x=element_blank(),strip.placement = "outside",legend.position = "none")+
  labs(x = "Chromosome", y = "LOD", color = "", linetype = "")
q
ggsave(paste(metabolite,'Cytoplasm-Co-varied-QTL-log.png',sep=''))
