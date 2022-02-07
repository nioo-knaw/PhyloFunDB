# Database pipeline 

Pipeline for specific gene database construction and update.

**Introduction**

The pipeline is based on the following workflow: 

![Db_pipeline](/uploads/e85a83817130c639c7537502834a46f2/Db_pipeline.png)


**Before starting, you can download a framebot database (protein sequences) for your specific gene. However, for some genes, there is no fungene database available, so it is also possible to skip the Framebot step.**

Download reference database from fungene, for framebot: http://fungene.cme.msu.edu/

In fungene, set the parameters that are suitable for your gene, such as:

min size aa = 600

hmm coverage = 99

**To remove duplicated sequences based on protein sequence you can use the following command, with the software seqkit:**

```
conda create -n seqkit
source activate seqkit
conda install seqkit

seqkit rmdup {gene}.fungene.fasta -s -o {gene}.fungene.clean.fasta

```

**Getting started**

**1.** Logon to the place where you will analysis your data, e.g. server

**2.** Create a local copy of the pipeline in a project folder

`git clone https://gitlab.bioinf.nioo.knaw.nl/OhanaC/database-pmoa.git`

**3.** Enter the pipeline folder with: 

 `cd database-pmoa`

**4.** The configuration of the pipeline needs to be set in the file **config.yaml**. Adjust the settings: 

```
gene: gene name
full_name: "protein full name"
minlength: minimum sequence length
cutoff_otu: cut off for OTU clustering (generally found in the literature)
cutoff_dm: cut off for distance matrix (in general, 0.25 is good enough)
framebot_db: false if there is no framebot reference database, otherwise, true
update: false, as you want to create a new dabatase
```
**5.** Check if the config file is correct and which steps will be run

`snakemake -n`

**6.** Run the pipeline. -j specifies the number of threads. Conda is the package manager. Optionally do this in a tmux session.

`snakemake -j 8 --use-conda`

 
_____________________________________________________________________________________________________________________________

# Updating the old pipeline 

Some time after you built your specific gene pipeline, it is possible to update it with the newest sequences uploaded to the NCBI database, setting a date range for downloading new sequences and adding the new OTUs to the reference tree. 

The update pipeline is based on the following workflow: 

![Update_pipeline](/uploads/15a708d2183a2f30acb37ddbff13eca8/Update_pipeline.png)

The most recent sequences will be downloaded, within a date range, processed and added to the initial reference tree.

**Getting started**

**1.** Logon to the place where you will analysis your data, e.g. server

**2.** Create a local copy of the pipeline in a project folder

`git clone https://gitlab.bioinf.nioo.knaw.nl/OhanaC/database-pmoa.git`

**3.** Enter the pipeline folder with: 

 `cd database-pmoa`

**4.** The configuration of the pipeline needs to be set in the file **config.update.yaml**. Adjust the settings: 

```
gene: gene name
full_name: "protein full name"
minlength: minimum sequence length
cutoff_otu: cut off for OTU clustering (generally found in the literature)
cutoff_dm: cut off for distance matrix (in general, 0.25 is good enough)
framebot_db: false if there is no framebot reference database, otherwise, true
update: true, very important
mindate: date after the sequences in the initial database were donwloaded (yyyy/mm/dd)
maxdate: current day (yyyy/mm/dd)
path_to_tree: "path_to_the_reference_tree_of_the_database"
path_to_seqs: "path_to_the_sequences_used_to_build_the_reference_tree_of_the_database" - the new sequences need to be aligned to the sequences in the tree
path_to_db: "path_to_the_fasta_file_of_the_full_database"
path_to_tax: "path_to_the_taxonomy_file_of_the_full_database"

```
**5.** Check if the config file is correct and which steps will be run

`snakemake -n -s Snakefile.update`

**6.** Run the pipeline. -j specifies the number of threads. Conda is the package manager. Optionally do this in a tmux session.

`snakemake -j 8 -s Snakefile.update --use-conda`

___________________________________________________________________________________________________________________________________________________

**Updating the pmoA database**

For updating the pmoA database, there is a specific config file, **config.update.pmoa.yaml**, and a specific Snakefile, **Snakefile.update.pmoa**
The main difference in this file is that the query words to recover sequences belonging to pmoA and other copper monooxygenase genes are already included in the Snakefile.

In order to update the pmoA database, you can follow the steps to update any other gene, then check the date parameters in the config.update.pmoa.yaml file:

```
mindate: 2020/02/06
maxdate: 2020/09/03
```
Also check if the paths to the tree file, to the tree sequences, full database fasta and full taxonomy file are correct. 

Check if everything is correct:

`snakemake -n -s Snakefile.update.pmoa`

Then, run the pmoA update pipeline

`snakemake -j 8 -s Snakefile.update.pmoa --use-conda`

_____________________________________________________________________________________________________________________________________________________

# Refining sequence taxonomy 

After getting your database files, it is still necessary to check the taxonomy/clustering of the sequences in the phylogenetic tree and improve the unassigned/unclassified ones

It is also possible to download metadata from all sequences using the Entrez Direct from NCBI and add that information to the taxonomy string

```
conda create -n entrez
source activate entrez
conda install -c bioconda entrez-direct
```
Example nirK gene:

`esearch -db nucleotide -query "nirK[gene]" | efetch -format gpc | xtract -insd source organism mol_type strain country isolation_source | sort | uniq >metadata_nirk.txt`

After having your tree ready and metadata downloaded (optional):

**1.** Check if there are cultivated representatives in the OTU groups - not always the OTU representative sequence is a cultivated/known organism, the program is not able to distinguish that - check in the file "interm/{gene}.aligned.good.filter.unique.pick.good.filter.an.{cutoff_otu}.rep.names"

**2.** Look at the tree – check if there are defined clades in the literature

**3.** It is possible to match the sequences with their full string taxonomy donwloaded for NCBI using the using VLOOKUP function (Excel)- it makes the checking and formatting of the taxonomy file easier. 

**4.** After checking and refining/correcting the taxonomies of the OTUs (until genus level, work on species/strains in the full taxonomy list), it is necessary to expand the taxonomy to all the sequences in the OTU group (remember that not all the sequences in the db are in the tree, only the OTU representatives) - use the **expand_taxonomy.R** script

**5.** In the full taxonomy file, check whether the cultivated representatives have their correct taxonomy (Excel). Add the species and environmental origin to the last level of the taxonomy

**6.** Some formatting parameters that have to be observed: 
	
- Names and number of sequences in the .fasta and .taxonomy file must be equal
- Formatting will depend on the software to be used – remove all spaces, avoid different characters
	- for mothur, strings should end with “;”
	- for qiime2, strings should end without “;” - also, you have to remove the gaps from the sequences.






