



## PART 1: read in relevant files

## variants list
data <- read.table('../files/variants.tsv',h=T,stringsAsFactors=F,sep='\t')

## list of candidate and known genes
list <- read.table('../files/genelist.tsv',h=T,stringsAsFactors=F,sep='\t')
data <- subset(data,data$id %in% list$id)

## low quality vairants
low <- read.table('../files/low_quality.tsv',h=T,stringsAsFactors=F,sep='\t')
data <- subset(data,data$snv %nin% low$snv)

freq_db <- strsplit("1000g_all,1000g_chbs,ExAc_EAS,gnomAD_t_EAS_AF,ExAc_all,gnomAD_exome_AF,gnomAD_genome_AF,gnomAD_total_AF,gnomAD_total_MAXPOP,HB.5000(AF)",",")[[1]]



for ( db in freq_db){
	data <- subset(data,db > 0.001)

}

