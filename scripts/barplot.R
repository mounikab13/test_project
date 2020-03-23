###

library(ggplot2)

condition<-c(rep(c("nanocount","illumina"),5))
value<-c(13361612,2210201946,5336687,198964151,9390854,583418846,60021,168079,18295,60676)
type<-c("reads in fastq","reads in fastq","mapped to genome","mapped to genome","mapped to transcriptome","mapped to transcriptome","different transcripts","different transcripts","different genes","different genes")



data1<-data.frame(type,condition,value)
ggplot(data1, aes(fill=condition,y=value, x=type)) + 
  geom_bar(position="dodge", stat="identity") + scale_y_log10()


