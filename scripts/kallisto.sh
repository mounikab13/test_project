
## transcriptome fasta obtained from bedtools getfasta 

kallisto index --make-unique --index=human rna.fa

kallisto quant -i human -o /hpcnfs/home/ieo5306/projects/kallisto/data/result ENCFF000EWJ.fastq ENCFF000EWX.fastq ENCFF000EZF.fastq ENCFF000EZT.fastq ENCFF000EZR.fastq ENCFF000EZH.fastq ENCFF000FCG.fastq ENCFF000FCI.fastq ENCFF000FCH.fastq ENCFF000FCU.fastq ENCFF001RCX.fastq ENCFF001RDF.fastq ENCFF001RCY.fastq ENCFF001RDG.fastq 

