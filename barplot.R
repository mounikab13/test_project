attributes<-c("number of reads in fastq","reads mapped to genome")
nano_values<-c(13361612,16710724)
star_reads<-c(138137622,229692250)
nano_transcripts<-c(33146409)
kallisto_reads<-c(227818)
kallisto<-read.table("abundance.tsv",header=TRUE,row.names = 1)
graph1<-data.frame(nano_values,star_reads)
row.names(graph1)<-attributes
fastq_reads<-c(13361612,138137622)
genome_reads<-c(5336687,229692250)
names(fastq_reads)<-c("Nanocount","STAR")
names(genome_reads)<-c("Nanocount","STAR")
barplot(fastq_reads,main = "Number of reads in fastq")
barplot(genome_reads,main="number of reads mapped to genome")
transcript_reads<-c(9390854,227818)
names(transcript_reads)<-c("Nanocount","KALLISTO")
barplot(transcript_reads,main = "reads mapped to transcriptome")
data<-data.frame(fastq_reads,genome_reads,transcript_reads)

r<-kallisto[kallisto$tpm>0,]
length(r$tpm)

###

library(ggplot2)

nanocount1<-c(13361612,16710724,33146409,60021,18295)
illumina_reads<-c(138137622,229692250,227818,168079,60676)
desc<-c("reads in fastq","mapped to genome","mapped to transcriptome","different transcripts","different genes")
data<-data.frame(desc,nanocount1,illumina_reads)
names(nanocount1)<-desc
names(illumina_reads)<-desc
barplot(illumina_reads)
barplot(nanocount1)
ggplot(data,aes(y=c(nanocount1,illumina_reads),x=c(desc,desc))) + geom_bar(position = "dodge",stat = "identity") + scale_y_log10()
ggplot(data,aes(y=illumina_reads,x=desc)) + geom_bar(position = "dodge",stat = "identity") + scale_y_log10()




condition<-c(rep(c("nanocount","illumina"),5))
value<-c(13361612,138137622,16710724,229692250,33146409,227818,60021,168079,18295,60676)
type<-c(desc,desc)


data1<-data.frame(type,condition,value)
ggplot(data1, aes(fill=condition,y=value, x=type)) + 
  geom_bar(position="dodge", stat="identity") + scale_y_log10()

