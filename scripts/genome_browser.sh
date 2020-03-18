
#!/bin/bash

##nanopore:


samtools view -S -b genome.sam | samtools sort > genome.bam


bedtools bamtobed -split -bed12 -i genome.bam > DNA.bed

bedparse convertChr --assembly hg38 --target ucsc -s DNA.bed > dna.bed

wget -P /hpcnfs/home/ieo5306/projects/genome_analysis/data http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bedToBigBed

chmod +x bedToBigBed

./bedToBigBed dna.bed hg38.chrom.sizes mydna.bb


##STAR:


samtools view -S -b EWJ_EWXAligned.out.sam | samtools sort > /hpcnfs/home/ieo5306/projects/rna_star/browser/EWJ_EWX.bam

samtools view -S -b EWW_EXEAligned.out.sam | samtools sort > /hpcnfs/home/ieo5306/projects/rna_star/browser/EWW_EXE.bam

samtools view -S -b EZF_EZTAligned.out.sam | samtools sort > /hpcnfs/home/ieo5306/projects/rna_star/browser/EZF_EZT.bam

samtools view -S -b EZR_EZHAligned.out.sam | samtools sort > /hpcnfs/home/ieo5306/projects/rna_star/browser/EZR_EZH.bam

samtools view -S -b FCG_FCIAligned.out.sam | samtools sort > /hpcnfs/home/ieo5306/projects/rna_star/browser/FCG_FCI.bam

samtools view -S -b FCH_FCUAligned.out.sam | samtools sort > /hpcnfs/home/ieo5306/projects/rna_star/browser/FCH_FCU.bam

samtools view -S -b RCX_RDFAligned.out.sam | samtools sort > /hpcnfs/home/ieo5306/projects/rna_star/browser/RCX_RDF.bam

samtools view -S -b RCY_RDGAligned.out.sam | samtools sort > /hpcnfs/home/ieo5306/projects/rna_star/browser/RCY_RDG.bam


bedtools bamtobed -split -bed12 -i EWJ_EWX.bam > ewj_ewx.bed

bedtools bamtobed -split -bed12 -i EWW_EXE.bam > eww_exe.bed

bedtools bamtobed -split -bed12 -i EZF_EZT.bam > ezf_ezt.bed

bedtools bamtobed -split -bed12 -i EZR_EZH.bam > ezr_ezh.bed

bedtools bamtobed -split -bed12 -i FCG_FCI.bam > fcg_fci.bed

bedtools bamtobed -split -bed12 -i FCH_FCU.bam > fch_fcu.bed

bedtools bamtobed -split -bed12 -i RCX_RDF.bam > rcx_rdf.bed

bedtools bamtobed -split -bed12 -i RCY_RDG.bam > rcy_rdg.bed

bedparse convertChr --assembly hg38 --target ucsc -s ewj_ewx.bed > EWJ_EWX.bed

bedparse convertChr --assembly hg38 --target ucsc -s eww_exe.bed > EWW_EXE.bed

bedparse convertChr --assembly hg38 --target ucsc -s ezf_ezt.bed > EZF_EZT.bed

bedparse convertChr --assembly hg38 --target ucsc -s ezr_ezh.bed > EZR_EZH.bed

bedparse convertChr --assembly hg38 --target ucsc -s fcg_fci.bed > FCG_FCI.bed

bedparse convertChr --assembly hg38 --target ucsc -s fch_fcu.bed > FCH_FCU.bed

bedparse convertChr --assembly hg38 --target ucsc -s rcx_rdf.bed > RCX_RDF.bed

bedparse convertChr --assembly hg38 --target ucsc -s rcy_rdg.bed > RCY_RDG.bed


wget -P /hpcnfs/home/ieo5306/projects/rna_star/browser http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bedToBigBed

chmod +x bedToBigBed

./bedToBigBed ewj_ewx.bed hg38.chrom.sizes star1.bb

