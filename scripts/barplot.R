###

library(ggplot2)

nanocount1<-c(13361612,16710724,33146409,60021,18295)
illumina_reads<-c(138137622,229692250,227818,168079,60676)
desc<-c("reads in fastq","mapped to genome","mapped to transcriptome","different transcripts","different genes")
data<-data.frame(desc,nanocount1,illumina_reads)
names(nanocount1)<-desc
names(illumina_reads)<-desc

ggplot(data,aes(y=c(nanocount1,illumina_reads),x=c(desc,desc))) + geom_bar(position = "dodge",stat = "identity") + scale_y_log10()
ggplot(data,aes(y=illumina_reads,x=desc)) + geom_bar(position = "dodge",stat = "identity") + scale_y_log10()

condition<-c(rep(c("nanocount","illumina"),5))
value<-c(13361612,138137622,16710724,229692250,33146409,227818,60021,168079,18295,60676)
type<-c(desc,desc)


data1<-data.frame(type,condition,value)
ggplot(data1, aes(fill=condition,y=value, x=type)) + 
  geom_bar(position="dodge", stat="identity") + scale_y_log10()

