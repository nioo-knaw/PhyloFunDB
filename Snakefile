configfile: "config.yaml"

rule final:
    input: expand("{gene}.fasta \
                   tax_{gene}.txt".split(), gene=config["gene"])

rule download_ncbi:
    output: 
        "{gene}.fasta"
    conda: 
        "envs/entrez.yaml"
    params:
        gene="{gene}",
        full_name=config["full_name"]
    message: "Retrieving gene sequences from NCBI. Note, this can take a (long) while"
    shell:"""esearch -db nucleotide -query "{params.gene}[gene]" | \
           efetch -format gpc | \
           xtract -pattern INSDFeature -if INSDFeature_key -equals CDS -and INSDQualifier_value -equals {params.gene} -or INSDQualifier_value -contains '{params.full_name}' -element INSDInterval_accession -element INSDInterval_from -element INSDInterval_to | \
           sort -u -k1,1 | \
           uniq | xargs -n 3 sh -c 'efetch -db nuccore -id "$0" -seq_start "$1" -seq_stop "$2" -format fasta' > {output}
           """

rule download_taxonomy:
    output:
        "tax_{gene}.txt"
    conda:
        "envs/entrez.yaml"
    params:
        gene="{gene}",
        full_name=config["full_name"]
    message: "Retrieving taxaonomy from NCBI. Note, this can take a (long) while"
    shell:"""
esearch -db nucleotide -query "{params.gene}[gene]" | \
efetch -format gpc | \
xtract -pattern INSDSeq -if INSDFeature_key -equals CDS -and INSDQualifier_value -equals {params.gene} -or INSDQualifier_value -contains '{params.full_name}' -element INSDSeq_accession-version -element INSDSeq_taxonomy | sort -u -k1,1 | uniq > {output}
"""

