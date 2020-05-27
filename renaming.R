#install.packages("tidyr")
library("tidyr")

#get fasta sequence names
#get taxonomy file
#isolate accession number from fasta seqnames
#add column with accession numbers to fasta seqname files
#isolate accession number from taxonomy file
#merge taxonomy and fasta seqnames files by accession number column
#final file with only sequence names and associated taxonomy

gene.names<-read.table(file=snakemake@input[["fasta"]])
gene.tax<-read.table(file=snakemake@input[["tax"]], sep="\t")
gene.names.columns<-separate(data = gene.names, col = V1, into = c("accession", "rest"), sep = "\\_")
gene.names$accession<-gene.names.columns$accession
gene.tax.columns<-separate(data = gene.tax, col = V1, into = c("accession", "rest"), sep = "\\.")
gene.taxonomy<-merge(gene.names,gene.tax.columns, by.x = "accession", by.y = "accession", all.x = FALSE)
gene.taxonomy.final<-gene.taxonomy[,c(2,4)]

write.table(gene.taxonomy.final, file=snakemake@output[["final_tax"]], sep="\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
