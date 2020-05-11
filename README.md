# Database pmoa

Database pmoA construction, pipeline and stuff related

**#before starting**

download db from fungene for framebot

min size aa = 600

hmm coverage = 99

remove duplicated sequences based on sequence

seqkit rmdup nosz.fungene.fasta -s -o nosz.fungene.clean.fasta

count remaining sequences

grep &quot;\&gt;&quot; nosz.fungene.clean.fasta \&gt; wc -l

201 sequences for framebot database

**+++++++++++++++++++++++++++++++++++++++++++++++++++++**

**#Download sequences using Entrez Direct from NCBI - it may take several hours.**

esearch -db nucleotide -query &quot;nosZ[gene]&quot; | efetch -format gpc | xtract -pattern INSDFeature -if INSDFeature\_key -equals CDS -and INSDQualifier\_value -equals nosZ -or INSDQualifier\_value -contains &#39;nitrous-oxide reductase&#39; -element INSDInterval\_accession -element INSDInterval\_from -element INSDInterval\_to | sort -u -k1,1 | uniq | xargs -n 3 sh -c &#39;efetch -db nuccore -id &quot;$0&quot; -seq\_start &quot;$1&quot; -seq\_stop &quot;$2&quot; -format fasta&#39;\&gt;nosZ.fasta

**Plus taxonomy file for formatting later:**

esearch -db nucleotide -query &quot;nosZ[gene]&quot; | efetch -format gpc | xtract -pattern INSDSeq -if INSDFeature\_key -equals CDS -and INSDQualifier\_value -equals nosZ -or INSDQualifier\_value -contains &#39;nitrous-oxide reductase&#39; -element INSDSeq\_accession-version -element INSDSeq\_taxonomy | sort -u -k1,1 | uniq \&gt;tax\_nosz.txt

**Rename\_seqs**

sed -i -e &#39;s/[.]/|/&#39; -e &#39;s/ /|/g&#39; nosz.fasta

cut -d &quot;|&quot; -f 1,3,4 nosz.fasta\&gt;nos.renamed.fasta

sed -i -e &#39;s/[|]/\_/g&#39; -e &#39;s/[.]//g&#39; -e &#39;s/[,]//g&#39; nos.renamed.fasta

in= nosz.fasta

out= nos.renamed.fasta

**#Remove Unverified sequences (it was not necessary for nosz)**

grep UNVERIFIED nos.renamed.fasta\&gt;unverified.accnos

cut -c 2- unverified.accnos \&gt; unverified\_.accnos

#mothur

Mothur \&gt; Remove.seqs(accnos=unverified\_.accnos, fasta=nos.renamed.fasta)

Out= nos.renamed.pick.fasta

25915 seqs

**Remove**  **short and containing ambiguous bases**  **sequences -Mothur**

trim.seqs(fasta=nos.renamed.fasta, minlength=350, maxambig=0, processors=10)

in= nos.renamed.fasta

out = nos.renamed.trim.fasta

16256 seqs

**#Framebot**

java -jar /mnt/nfs/bioinfdata/data\_other/tools/RDPTools/FrameBot.jar framebot -o nosz.framebot -N ./nosz\_fungene/nosz.fungene.clean.fasta nos.renamed.trim.fasta

in= nos.renamed.trim.fasta

db=nosz.fungene.clean.fasta

out=nosz.framebot\_corr\_nucl.fasta

nosz.framebot\_corr\_prot.fasta

16213 corrected seqs

**Align sequences**

mafft --thread 10 --auto nosz.framebot\_corr\_nucl.fasta \&gt;nosz.aligned.fasta

in= nosz.framebot\_corr\_nucl.fasta

out= nosz.aligned.fasta

**Screen/filter alignment**

screen.seqs(fasta=nosz.aligned.fasta, optimize=start-end, criteria=96, processors=10)

in= nosz.aligned.fasta

out= nosz.aligned.good.fasta

**filter alignment**

filter.seqs(fasta=nosz.aligned.good.fasta, vertical=T, trump=., processors=10)

in= nosz.aligned.good.fasta

out=nosz.aligned.good.filter.fasta

**Dereplicate seqs (necessary for chimera.vsearch without reference)**

unique.seqs(fasta=nosz.aligned.good.filter.fasta)

out= nosz.aligned.good.filter.names

nosz.aligned.good.filter.unique.fasta

**chimera detection without reference**

chimera.vsearch(fasta=nosz.aligned.good.filter.unique.fasta, name=nosz.aligned.good.filter.names)

in = nosz.aligned.good.filter.unique.fasta

nosz.aligned.good.filter.names

out= nosz.aligned.good.filter.unique.denovo.vsearch.accnos

**chimera removal**

remove.seqs(accnos=nosz.aligned.good.filter.unique.denovo.vsearch.accnos, fasta=nosz.aligned.good.filter.unique.fasta, name=nosz.aligned.good.filter.names)

in= nosz.aligned.good.filter.unique.denovo.vsearch.accnos

nosz.aligned.good.filter.unique.fasta

nosz.aligned.good.filter.uniquefasta

out=nosz.aligned.good.filter.pick.names

nosz.aligned.good.filter.unique.pick.fasta

14801 seqs

**Screen/filter alignment**

screen.seqs(fasta=nosz.aligned.good.filter.unique.pick.fasta, name=nosz.aligned.good.filter.pick.names, optimize=start-end, criteria=96, processors=10)

in= nosz.aligned.good.filter.unique.pick.fasta

nosz.aligned.good.filter.pick.names

out= nosz.aligned.good.filter.unique.pick.good.fasta

nosz.aligned.good.filter.pick.good.names

filter.seqs(fasta=nosz.aligned.good.filter.unique.pick.good.fasta, vertical=T, trump=., processors=10)

in= nosz.aligned.good.filter.unique.pick.good.fasta

out= nosz.aligned.good.filter.unique.pick.good.filter.fasta

Liu et al. 2019 Cutoff nosZ 89% similarity = 11% dissimilarity 0.11

**Distance matrix**

dist.seqs(fasta=nosz.aligned.good.filter.unique.pick.good.filter.fasta, cutoff=0.3)

in= nosz.aligned.good.filter.unique.pick.good.filter.fasta

out= nosz.aligned.good.filter.unique.pick.good.filter.dist

**Clustering**

cluster(column=nosz.aligned.good.filter.unique.pick.dist, name=nosz.aligned.good.filter.pick.names, method=average, cutoff=0.3)

in= nosz.aligned.good.filter.unique.pick.dist

nosz.aligned.good.filter.pick.names

out= nosz.aligned.good.filter.unique.pick.an.list

**OTU representatives for phylogenetic tree**

get.oturep(column=nosz.aligned.good.filter.unique.pick.good.filter.dist, fasta=nosz.aligned.good.filter.unique.pick.good.filter.fasta, name=nosz.aligned.good.filter.pick.good.names, list=nosz.aligned.good.filter.unique.pick.good.filter.an.list, cutoff=0.11)

in= nosz.aligned.good.filter.unique.pick.good.filter.dist

nosz.aligned.good.filter.unique.pick.good.filter.fasta

nosz.aligned.good.filter.pick.good.names

nosz.aligned.good.filter.unique.pick.good.filter.an.list

out= nosz.aligned.good.filter.unique.pick.good.filter.an.0.11.rep.fasta

nosz.aligned.good.filter.unique.pick.good.filter.an.0.11.rep.names

**Get redundant dataset (final database):**

deunique.seqs(fasta=nosz.aligned.good.filter.unique.pick.good.filter.fasta, name=nosz.aligned.good.filter.pick.good.names)

in= nosz.aligned.good.filter.unique.pick.good.filter.fasta

nosz.aligned.good.filter.pick.good.names

out= nosz.aligned.good.filter.unique.pick.redundant.fasta

**Get the representatives tree:**

iqtree -s nosz.aligned.good.filter.unique.pick.good.filter.an.0.11.rep.fasta -m MFP -alrt 1000 -bb 1000 -nt 10

in= nosz.aligned.good.filter.unique.pick.good.filter.an.0.11.rep.fasta

out=nosz.aligned.good.filter.unique.pick.good.filter.an.0.11.rep.fasta.treefile

**Placement new sequences in the tree**

mafft --thread 10 --auto otu.rep.011\_named\_to\_align.fasta\&gt; otu.rep.011\_aligned.fasta

in= otu.rep.011\_named\_to\_align.fasta

out= otu.rep.011\_aligned.fasta

raxmlHPC -f v -s otu.rep.011\_aligned.fasta -t tree\_newick.txt -m GTRCAT -H -n TEST

in= otu.rep.011\_aligned.fasta

tree\_newick.txt

out= otu.rep.011\_aligned.fasta.treefile