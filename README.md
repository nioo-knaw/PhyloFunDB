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

1. Logon to the place where you will analysis your data, e.g. server

2. Create a local copy of the pipeline in a project folder

`git clone https://gitlab.bioinf.nioo.knaw.nl/OhanaC/database-pmoa.git`

3. Enter the pipeline folder with: 

 `cd database-pmoa`

4. The configuration of the pipeline needs to be set in the file **config.yaml**. Adjust the settings: 

```
gene: gene name
full_name: "protein full name"
minlength: minimum sequence length
cutoff_otu: cut off for OTU clustering (generally found in the literature)
cutoff_dm: cut off for distance matrix (in general, 0.25 is good enough)
framebot_db: false if there is no framebot reference database, otherwise, true
```
5. Check if the config file is correct and which steps will be run

`snakemake -n`

6. Run the pipeline. -j specifies the number of threads. Conda is the package manager. Optionally do this in a tmux session.

`snakemake -j 8 --use-conda`

 
_____________________________________________________________________________________________________________________________

# Update pipeline 

Some time after you built your specific gene pipeline, it is possible to update it with the newest sequences uploaded to the NCBI database, setting a date range for downloading new sequences and adding the new OTUs to the reference tree. 

The update pipeline is based on the following workflow: 

![Update_pipeline](/uploads/15a708d2183a2f30acb37ddbff13eca8/Update_pipeline.png)

The most recent sequences will be downloaded, within a date range, processed and added to the initial reference tree.

**Getting started**

1. Logon to the place where you will analysis your data, e.g. server

2. Create a local copy of the pipeline in a project folder

`git clone https://gitlab.bioinf.nioo.knaw.nl/OhanaC/database-pmoa.git`

3. Enter the pipeline folder with: 

 `cd database-pmoa`

4. The configuration of the pipeline needs to be set in the file **config.update.yaml**. Adjust the settings: 

```
gene: gene name
full_name: "protein full name"
mindate: date after the sequences in the initial database were donwloaded
maxdate: current day
minlength: minimum sequence length
cutoff_otu: cut off for OTU clustering (generally found in the literature)
cutoff_dm: cut off for distance matrix (in general, 0.25 is good enough)
framebot_db: false if there is no framebot reference database, otherwise, true
path_to_tree: "path_to_the_reference_tree_of_the_database"
path_to_seqs: "path_to_the_sequences_used_to_build_the_reference_tree_of_the_database" - the new sequences need to be aligned to the sequences in the tree
path_to_db: "path_to_the_fasta_file_of_the_full_database"
path_to_tax: "path_to_the_taxonomy_file_of_the_full_database"

```
5. Check if the config file is correct and which steps will be run

`snakemake -n`

6. Run the pipeline. -j specifies the number of threads. Conda is the package manager. Optionally do this in a tmux session.

`snakemake -j 8 -s Snakefile.update --use-conda`


