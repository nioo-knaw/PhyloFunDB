configfile: "config.yaml"

rule final:
    input: expand("{gene}.framebot_corr_nucl.fasta, \
                   tax_{gene}.txt".split(), gene=config["gene"], minlength={config['minlength']})

rule download_ncbi:
    output: 
        "{gene}_fasta"
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
    message: "Retrieving taxonomy from NCBI. Note, this can take a (long) while"
    shell:"""
esearch -db nucleotide -query "{params.gene}[gene]" | \
efetch -format gpc | \
xtract -pattern INSDSeq -if INSDFeature_key -equals CDS -and INSDQualifier_value -equals {params.gene} -or INSDQualifier_value -contains '{params.full_name}' -element INSDSeq_accession-version -element INSDSeq_taxonomy | sort -u -k1,1 | uniq > {output}
"""

rule rename_seqs:
    input:
        "{gene}.fasta"
    output:
        "{gene}.renamed.fasta"
    shell:
       """sed -e 's/[.]/-/' -e 's/ /-/g' {input} | \
       stdbuf -o0 cut -d "-" -f 1,4,5| \
       sed -e 's/-/_/g' -e 's/[.]//g' -e 's/[,]//g' > {output}
       """
       
rule get_unverified:
    input:
        "{gene}.renamed.fasta"
    output:
        "{gene}.unverified.accnos"
    shell:
        """grep UNVERIFIED {input} | \
        stdbuf -o0 cut -c 2- > {output}
        """
        
rule remove_unverified:
        input:
            accnos="{gene}.unverified.accnos",
            fasta="{gene}.renamed.fasta"
        output:
            "{gene}.renamed.pick.fasta"
        conda:
            "envs/mothur.yaml"
        shell:
            '''
            mothur "#remove.seqs(accnos={input.accnos}, fasta={input.fasta})"
            '''
            
rule mothur_trim:
        input:
            "{gene}.renamed.pick.fasta"
        output:
            "{gene}.renamed.pick.trim.fasta"
        conda:
            "envs/mothur.yaml"
        params:
            minlength={config['minlength']}
        threads:10
        shell:
            '''
            mothur "#trim.seqs(fasta={input}, minlength={params.minlength}, maxambig=0, processors={threads})"
            '''
rule framebot:
        input:
            fasta="nosZ.renamed.pick.trim.fasta",
            db_framebot="nosZ.fungene.clean.fasta"
        output:
            "nosZ.framebot_corr_nucl.fasta",
        params:
            "nosZ.framebot"
        conda:
            "envs/rdptools.yaml"
        shell:
            "java -jar /mnt/nfs/bioinfdata/home/NIOO/ohanac/.conda/envs/rdptools/share/rdptools-2.0.2-1/FrameBot.jar framebot -o {params} -N {input.db_framebot} {input.fasta}"
       
       
       
       
       