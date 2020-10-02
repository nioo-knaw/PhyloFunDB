setwd("P:/Project Database/dbtest/nirk_bac/")
library(tidyr)

#names archive - contains name of rep sequence and the list of sequence names associated to it
otulist<-read.table(choose.files(), sep="\t")

#taxonomy file, containing sequence name and taxonomy string
otutaxseqnames<-read.table(choose.files(), sep="\t", header = FALSE)

#convert list of names in table
list_seqs_table<-separate_rows(otulist, V2, sep=",")

#merge both tables by sequence names
expanded_taxonomy<-merge(list_seqs_table,otutaxseqnames, by.x = "V1", by.y = "V1", all.x = TRUE)
taxonomy<-expanded_taxonomy[,c(2,3)]

#extract table to verify
write.table(taxonomy, file="Expanded_taxonomy.txt", sep="\t", quote = FALSE, row.names = FALSE, col.names = FALSE)


        


