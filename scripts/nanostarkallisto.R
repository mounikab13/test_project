####nanocount table

setwd("~/projects/nanocount/data")
test_reads<-read.table("/hpcnfs/home/ieo5306/projects/nanocount/data/counts",sep = "\t",header = TRUE)
View(test_reads)
#####star table

setwd("~/projects/rna_star/analysis")


##EWJ_EWX
EWJ_EWX<-read.table("/hpcnfs/home/ieo5306/projects/rna_star/analysis/EWJ_EWXReadsPerGene.out.tab",sep = "\t")
colnames(EWJ_EWX)<-c("gene_id","unstranded_counts","EWJcounts","EWXcounts")

#EWW_EXE
EWW_EXE<-read.table("/hpcnfs/home/ieo5306/projects/rna_star/analysis/EWW_EXEReadsPerGene.out.tab",sep = "\t")
colnames(EWW_EXE)<-c("gene_id","unstranded_counts","EWWcounts","EXEcounts")

##EZF_EZT
EZF_EZT<-read.table("/hpcnfs/home/ieo5306/projects/rna_star/analysis/EZF_EZTReadsPerGene.out.tab",sep = "\t")
colnames(EZF_EZT)<-c("gene_id","unstranded_counts","EZFcounts","EZTcounts")


##EZR_EZH
EZR_EZH<-read.table("/hpcnfs/home/ieo5306/projects/rna_star/analysis/EZR_EZHReadsPerGene.out.tab",sep = "\t")
colnames(EZR_EZH)<-c("gene_id","unstranded_counts","EZRcounts","EZHcounts")


##FCG_FCI
FCG_FCI<-read.table("/hpcnfs/home/ieo5306/projects/rna_star/analysis/FCG_FCIReadsPerGene.out.tab",sep = "\t")
colnames(FCG_FCI)<-c("gene_id","unstranded_counts","FCGcounts","FCIcounts")

##FCH_FCU
FCH_FCU<-read.table("/hpcnfs/home/ieo5306/projects/rna_star/analysis/FCH_FCUReadsPerGene.out.tab",sep = "\t")
colnames(FCH_FCU)<-c("gene_id","unstranded_counts","FCHcounts","FCUcounts")

#RCX_RDF
RCX_RDF<-read.table("/hpcnfs/home/ieo5306/projects/rna_star/analysis/RCX_RDFReadsPerGene.out.tab",sep = "\t")
colnames(RCX_RDF)<-c("gene_id","unstranded_counts","RCXcounts","RDFcounts")

#RCY_RDG
RCY_RDG<-read.table("/hpcnfs/home/ieo5306/projects/rna_star/analysis/RCY_RDGReadsPerGene.out.tab",sep = "\t")
colnames(RCY_RDG)<-c("gene_id","unstranded_counts","RCYcounts","RDGcounts")


######to get transcript ids for gene ids using biomaRt

library(biomaRt)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("biomaRt")
library(biomaRt)
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
listAttributes(human)
attributes<-listAttributes(human)

###reads contain gene id for every transcript id
reads <- getBM(attributes=c('ensembl_gene_id','ensembl_transcript_id'),mart = human)

class(FCG_FCI$gene_id)
##factor
class(reads$ensembl_gene_id)
##character
FCG_FCI$gene_id<-as.character(FCG_FCI$gene_id)
class(FCG_FCI$gene_id)

###trimming star data of unwanted entries

gene <- FCG_FCI$gene_id
head(gene)
unwanted <- c('N_unmapped',"N_multimapping",'N_noFeature','N_ambiguous')
gene<-gene[!gene %in% unwanted]
head(gene)


###extract transcripts only for star encode genes
FCG_FCI_genes<-reads[reads$ensembl_gene_id%in%gene,]
head(FCG_FCI_genes)


#filtering star data

filt_FCG_FCI<-FCG_FCI[FCG_FCI$gene_id%in%gene,]
colnames(filt_FCG_FCI)<-c("gene_id","unstranded_counts","FCGcounts","FCIcounts")
filt_FCG_FCI<-filt_FCG_FCI[c('gene_id','FCIcounts')]

filt_EWJ_EWX<-EWJ_EWX[EWJ_EWX$gene_id%in%gene,]
filt_EWJ_EWX<-filt_EWJ_EWX[c('gene_id','EWXcounts')]


filt_EWW_EXE<-EWW_EXE[EWW_EXE$gene_id%in%gene,]
filt_EWW_EXE<-filt_EWW_EXE[c('gene_id','EXEcounts')]

filt_EZF_EZT<-EZF_EZT[EZF_EZT$gene_id%in%gene,]
filt_EZF_EZT<-filt_EZF_EZT[c('gene_id','EZTcounts')]


filt_EZR_EZH<-EZR_EZH[EZR_EZH$gene_id%in%gene,]
filt_EZR_EZH<-filt_EZR_EZH[c('gene_id','EZHcounts')]

filt_RCX_RDF<-RCX_RDF[RCX_RDF$gene_id%in%gene,]
filt_RCX_RDF<-filt_RCX_RDF[c('gene_id','RCXcounts')]


filt_RCY_RDG<-RCY_RDG[RCY_RDG$gene_id%in%gene,]
filt_RCY_RDG<-filt_RCY_RDG[c('gene_id','RCYcounts')]


filt_FCH_FCU<-FCH_FCU[FCH_FCU$gene_id%in%gene,]
filt_FCH_FCU<-filt_FCH_FCU[c('gene_id','FCUcounts')]

#filtering nanocount table
nanocount<-test_reads[c('transcript_name','raw')]
nanocount$raw<-nanocount$raw * 1000000
library(dplyr)
nanocount<-cbind(nanocount,mutate(nanocount, ENSEMBL_ID=gsub("(ENST[0-9]+)::.+$", "\\1", transcript_name)))
nanocount<-nanocount[c('ENSEMBL_ID','raw')]


#normalization of star data

filt_EWJ_EWX$EWXcounts<-filt_EWJ_EWX$EWXcounts/sum(filt_EWJ_EWX$EWXcounts)
filt_EWW_EXE$EXEcounts<-filt_EWW_EXE$EXEcounts/sum(filt_EWW_EXE$EXEcounts)
filt_EZF_EZT$EZTcounts<-filt_EZF_EZT$EZTcounts/sum(filt_EZF_EZT$EZTcounts)
filt_EZR_EZH$EZHcounts<-filt_EZR_EZH$EZHcounts/sum(filt_EZR_EZH$EZHcounts)
filt_FCG_FCI$FCIcounts<-filt_FCG_FCI$FCIcounts/sum(filt_FCG_FCI$FCIcounts)
filt_FCH_FCU$FCUcounts<-filt_FCH_FCU$FCUcounts/sum(filt_FCH_FCU$FCUcounts)
filt_RCX_RDF$RCXcounts<-filt_RCX_RDF$RCXcounts/sum(filt_RCX_RDF$RCXcounts)
filt_RCY_RDG$RCYcounts<-filt_RCY_RDG$RCYcounts/sum(filt_RCY_RDG$RCYcounts)

#adding all the normalized values

allstar<- merge(filt_EWJ_EWX,filt_EZF_EZT,by='gene_id')
allstar<-merge(allstar,filt_EZR_EZH,by='gene_id')
allstar<-merge(allstar,filt_EWW_EXE,by='gene_id')
allstar<-merge(allstar,filt_FCG_FCI,by='gene_id')
allstar<-merge(allstar,filt_FCH_FCU,by='gene_id')
allstar<-merge(allstar,filt_RCX_RDF,by='gene_id')
allstar<-merge(allstar,filt_RCY_RDG,by='gene_id')


# calculate mean value of all the normalized values

sum<-rowSums(allstar[,2:9])
average<-sum/8
allstar<-cbind(allstar,average)

finalstar<-data.frame(allstar$gene_id,allstar$average)
finalstar$allstar.average<-finalstar$allstar.average*1000000


###gene aggregation of nanocount table
colnames(nanocount)<-c("ensembl_transcript_id","raw")
test_nano<-merge(nanocount,FCG_FCI_genes,by="ensembl_transcript_id")
res <- with(test_nano, tapply(raw, ensembl_gene_id, sum))
res2 <- data.frame(gene=names(res), raw=res, row.names = NULL)

colnames(finalstar)<-c('gene_id','star_reads')
colnames(res2)<-c('gene_id','nanoreads')


finaldata<-merge(finalstar,res2,by='gene_id')

plot(1+finaldata$reads,1+finaldata$nanoreads, log='xy')

library(ggplot2)
ggplot(finaldata, aes(x=star_reads, y=nanoreads)) + geom_point(alpha=0.1) + scale_x_log10() + scale_y_log10() + theme_bw()

#pearson coefficient

cor(finaldata$star_reads,finaldata$nanoreads)
## 0.66835

#spearman coefficients

cor(finaldata$star_reads,finaldata$nanoreads,method="spearman")
## 0.8448692

#kallisto and nanocount quantification

###transcript quantification
setwd("~/projects/kallisto/data/result")
kallisto<-read.csv("abundance.tsv",sep="\t")

kallisto_final<-cbind(mutate(kallisto, target_id=gsub("(ENST[0-9]+)::.+$", "\\1", target_id)))
kal1<-cbind.data.frame(kallisto_final$target_id,kallisto_final$tpm)
colnames(kal1)<-c("ensembl_transcript_id","kallisto_counts")

nano_kallisto<-merge(kal1,nanocount,by="ensembl_transcript_id")
nano_kallisto$kallisto_counts<-nano_kallisto$kallisto_counts*1000000
colnames(nano_kallisto)<-c("ensembl_transcript_id","kallisto_counts","nano_reads")

plot(1+nano_kallisto$kallisto_counts,1+nano_kallisto$raw, log='xy')

ggplot(nano_kallisto, aes(x=kallisto_counts, y=nano_reads)) + geom_point(alpha=0.2) + scale_x_log10() + scale_y_log10() + theme_bw()

cor(nano_kallisto$kallisto_counts,nano_kallisto$nano_reads)
# 0.6589651

cor(nano_kallisto$kallisto_counts,nano_kallisto$nano_reads,method="spearman")
# 0.5104981


##gene aggregation of transcripts

kal2<-merge(kal1,FCG_FCI_genes,by="ensembl_transcript_id")

kal3 <- with(kal2, tapply(kallisto_counts, ensembl_gene_id, sum))
kal4 <- data.frame(gene_id=names(kal3), kallisto_reads=kal3, row.names = NULL)

star_kallisto<-merge(kal4,finalstar,by="gene_id")

ggplot(star_kallisto, aes(x=kallisto_reads, y=star_reads)) + geom_point(alpha=0.2) + scale_x_log10() + scale_y_log10() + theme_bw()

cor(star_kallisto$kallisto_reads,star_kallisto$star_reads)
# 0.4141877

cor(star_kallisto$kallisto_reads,star_kallisto$star_reads,method="spearman")
#0.7772911

plot(2+star_kallisto$kallisto_reads,2+star_kallisto$reads, log='xy')

#nanocount and kallisto

nanokal<-merge(res2,kal4,by="gene_id")
ggplot(nanokal, aes(x=kallisto_reads, y=nanoreads)) + geom_point(alpha=0.2) + scale_x_log10() + scale_y_log10() + theme_bw()


cor(nanokal$nanoreads,nanokal$kallisto_reads)
#0.7690679

cor(nanokal$nanoreads,nanokal$kallisto_reads,method = "spearman")
#0.8680328
