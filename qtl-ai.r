suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringr))
library(ggplot2)

g_data = read.table("RIL_snp.txt", header=TRUE, sep="\t", row.names = 1, check.names=FALSE,stringsAsFactors=FALSE)
#Phenotype data
pheno = read.table("flower-pheno-AI-RIL-2016-GH.txt", header=TRUE, row.names = 1, sep="\t", check.names=FALSE,stringsAsFactors=FALSE)

g_data[is.na(g_data)] <- 3

head(g_data)
head(pheno)

#remove the rows that have empty in them
pheno <- pheno[!grepl('Empty',rownames(pheno)),]

names <- gsub('\\d\\d\\d_','',rownames(pheno))
names <- gsub('_.*','',names)
names <- gsub('\\.\\w\\w\\w','',names)
rownames(pheno) <- names

com_id = intersect(colnames(g_data),rownames(pheno))

#Arrange the matrices accordingly
g_data = g_data[,com_id]
pheno = pheno[com_id,]

#pheno$Bergamotene <- as.numeric(gsub(',','.',pheno$Bergamotene))
pheno$Linalooll <- as.numeric(gsub(',','.',pheno$Linalooll))

metabolite <- "Linalooll"
#(E)-alpha-Bergamotene* TAB Linalooll

phenotype = as.numeric(pheno[,metabolite])

#log normalize the data
phenotype <- log(phenotype)
phenotype[is.infinite(phenotype)] <- 0
phenotype

lod_RIL = apply(g_data, 1, function(x){
  dat = cbind(1,suppressWarnings(as.factor(x)),phenotype)
  #dat = cbind(1,suppressWarnings(as.factor(as.numeric(x)),phenotype)
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
    #dat = cbind(1,suppressWarnings(as.factor(as.numeric(x)),pheno_sample)
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

#use the highest peak found on chr10
n <- length(lod_RIL)
min_marker = names(which(lod_RIL == sort(lod_RIL,partial=n-1)[n-1]))

#min_marker = names(which(lod_RIL == max(lod_RIL,na.rm=TRUE)))
effect_data = data.frame(accession = com_id, genotype=as.factor(g_data[min_marker,]),pheno=phenotype)
effect_data = effect_data[complete.cases(effect_data), ]
effect_data$genotype = as.factor(effect_data$genotype)

ggplot(effect_data, aes(x=genotype,y=pheno)) +
  geom_boxplot(width=0.2, outlier.shape=NA) +
  geom_jitter(width=0.1, size=0.5) +
  ggtitle(paste0(min_marker))
d = quantile(as.numeric(unlist(perm)), 0.95, na.rm = TRUE)

# Some preprocessing of the data
#add chr to the row name of the AI RILs so that their naming matches that of the MAGIC population
rownames(g_data) <- paste('chr',rownames(g_data),sep='') 
marker = as.data.frame(str_split_fixed(rownames(g_data), "_", 3))
marker$V1 = as.data.frame(str_split_fixed(marker$V1, "-", 2))$V1
marker$V1 = factor(marker$V1, levels=c(paste0("chr",1:12),"scaffold"))
marker$V2 = 1:dim(marker)[1]
marker = marker[,-3]
plotdata = cbind(marker, lod_RIL)
colnames(plotdata) = c("chrom","pos","lod")
head(plotdata)

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
ggsave(paste(metabolite,'-QTL-log-AI-Flower.pdf',sep=''))

#---------------------CO VARIATE--------------------------------------------------------------------------
#-------------------------------CO VARIATE----------------------------------------------------------------
#-----------------------------------------CO VARIATE------------------------------------------------------

#--0000---00---0-----0----0----0000--0----0---00000-0000-------0----0--0----0----0----000-000-0---0-000-------
#-0------0--0---0---0----0-0---0--0--0---0-0----0---00--------0-0---00-0---0-0---0-----0--0----0-0--0---------
#-0------0--0----0-0----0-0-0--0-0---0--0-0-0---0---0--------0-0-0--0-00--0-0-0--0-----0----0---0-----0-------
#--0000---00------0----0-----0-0--0--0-0-----0--0---0000----0-----0-0--0-0-----0-0000-000-000---0---000-------

#Now do the Co-variate analysis
lod_RIL = apply(g_data, 1, function(x){
  dat = cbind(1,suppressWarnings(as.factor(x)),suppressWarnings(as.factor(g_data[min_marker,])),phenotype)
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
    dat = cbind(1,suppressWarnings(as.factor(x)),suppressWarnings(as.factor(g_data[min_marker,])),pheno_sample)
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
ggsave(paste(metabolite,'-QTL-log-Coviaried-AI-Flower.pdf',sep=''))
