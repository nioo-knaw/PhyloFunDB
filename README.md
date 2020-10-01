# Database pmoa

Database pmoA construction, pipeline and stuff related

**Introduction**

The pipeline is based on the following workflow: 

![Db_pipeline](/uploads/d22d9c524e55c0ce2deb9e8a81b89a6f/Db_pipeline.png)


**Before starting, you can download a framebot database (protein sequences) for your specific gene. However, for some genes, there is no fungene database available, so it is also possible to skip the Framebot step**

Download db from fungene for framebot: http://fungene.cme.msu.edu/

In fungene, set the parameters that are suitable for your gene, such as:

min size aa = 600

hmm coverage = 99

**To remove duplicated sequences based on protein sequence you can use the following command, with the software seqkit:**

``` 
seqkit rmdup nosz.fungene.fasta -s -o nosz.fungene.clean.fasta

```

**Getting started**

1. Logon to the place where you will analysis your data, e.g. server

2. Create a local copy of the pipeline in a project folder

`git clone https://gitlab.bioinf.nioo.knaw.nl/OhanaC/database-pmoa.git`

3. Enter the pipeline folder with: 

 `cd database-pmoa`

4. The configuration of the pipeline needs to be set in the file config.yaml. Adjust the settings.




 








_____________________________________________________________________________________________________________________________

**Update pipeline**

The update pipeline is based on the following workflow: 

![Update_pipeline](/uploads/d86716da1c1bc60667f570fd81d154f9/Update_pipeline.png)

The most recent sequences can be downloaded, within a date range, processed and added to the initial reference tree.


