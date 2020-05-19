configfile: "config.yaml"
wildcard_constraints:
   gene = '\w+'

rule final:
    input: expand("{gene}.aligned.good.filter.unique.pick.good.filter.an.0.11.rep.fasta.treefile \
                  {gene}.aligned.good.filter.unique.pick.good.filter.redundant.fasta" \
                   tax_{gene}.txt".split(), gene=config["gene"], minlength={config['minlength']})

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
            fasta="{gene}.renamed.pick.trim.fasta",
            db_framebot="{gene}.fungene.clean.fasta"
        output:
            "{gene}.framebot_corr_nucl.fasta",
        params:
            "{gene}.framebot"
        conda:
            "envs/rdptools.yaml"
        shell:
            "java -jar /mnt/nfs/bioinfdata/home/NIOO/ohanac/.conda/envs/rdptools/share/rdptools-2.0.2-1/FrameBot.jar framebot -o {params} -N {input.db_framebot} {input.fasta}"
       
rule align:
        input:
            expand("{gene}.framebot_corr_nucl.fasta", gene=config["gene"])
        output:
            expand("{gene}.aligned.fasta", gene=config["gene"])
        conda:
            "envs/mafft.yaml"
        threads:10
        shell:
            "mafft --thread {threads} --auto {input} >{output}"
            
rule screen_alignment:
        input:
            expand("{gene}.aligned.fasta", gene=config["gene"])
        output:
            expand("{gene}.aligned.good.fasta", gene=config["gene"])
        conda:
            "envs/mothur.yaml"
        threads:10
        shell:
            '''
            mothur "#screen.seqs(fasta={input}, optimize=start-end, criteria=96, processors={threads})"
            '''
rule filter_alignment:
        input:
            expand("{gene}.aligned.good.fasta", gene=config["gene"])
        output:
            expand("{gene}.aligned.good.filter.fasta", gene=config["gene"])
        conda:
            "envs/mothur.yaml"
        threads:10
        shell:
            '''
            mothur "#filter.seqs(fasta={input}, vertical=T, trump=., processors={threads})"
            '''
rule unique1:
        input:
            expand("{gene}.aligned.good.filter.fasta", gene=config["gene"])
        output:
            expand("{gene}.aligned.good.filter.unique.fasta", gene=config["gene"]),
            expand("{gene}.aligned.good.filter.names", gene=config["gene"])
        conda:
            "envs/mothur.yaml"
        shell:
            '''
            mothur "#unique.seqs(fasta={input})"
            '''            
            
rule chimera_vsearch:
        input:
            fasta=expand("{gene}.aligned.good.filter.unique.fasta", gene=config["gene"]),
            name=expand("{gene}.aligned.good.filter.names", gene=config["gene"])
        output:
            expand("{gene}.aligned.good.filter.unique.denovo.vsearch.accnos", gene=config["gene"])
        conda:
            "envs/mothur.yaml"
        shell:
            '''
            mothur "#chimera.vsearch(fasta={input.fasta}, name={input.name})"
            '''
            
rule chimera_removal:
        input:
            accnos=expand("{gene}.aligned.good.filter.unique.denovo.vsearch.accnos", gene=config["gene"]),
            fasta=expand("{gene}.aligned.good.filter.unique.fasta", gene=config["gene"]),
            name=expand("{gene}.aligned.good.filter.names", gene=config["gene"])
        output:
            expand("{gene}.aligned.good.filter.pick.names", gene=config["gene"]),
            expand("{gene}.aligned.good.filter.unique.pick.fasta", gene=config["gene"])
        conda:
            "envs/mothur.yaml"
        shell:
            '''
            mothur "#remove.seqs(accnos={input.accnos}, fasta={input.fasta}, name={input.name})"
            '''
            
rule screen_alignment2:
        input:
            fasta=expand("{gene}.aligned.good.filter.unique.pick.fasta", gene=config["gene"]),
            name=expand("{gene}.aligned.good.filter.pick.names", gene=config["gene"])
        output:
            expand("{gene}.aligned.good.filter.unique.pick.good.fasta", gene=config["gene"]),
            expand("{gene}.aligned.good.filter.pick.good.names", gene=config["gene"])
        conda:
            "envs/mothur.yaml"
        threads:10
        shell:
            '''
            mothur "#screen.seqs(fasta={input.fasta}, name={input.name}, optimize=start-end, criteria=96, processors={threads})"
            '''
            
rule filter_alignment2:
        input:
            expand("{gene}.aligned.good.filter.unique.pick.good.fasta", gene=config["gene"]),
        output:
            expand("{gene}.aligned.good.filter.unique.pick.good.filter.fasta", gene=config["gene"])
        conda:
            "envs/mothur.yaml"
        threads:10
        shell:
            '''
            mothur "#filter.seqs(fasta={input}, vertical=T, trump=., processors={threads})"
            '''
            
rule distance_matrix:
        input:
            expand("{gene}.aligned.good.filter.unique.pick.good.filter.fasta", gene=config["gene"]),
        output:
            expand("{gene}.aligned.good.filter.unique.pick.good.filter.dist", gene=config["gene"])
        conda:
            "envs/mothur.yaml"
        threads:10
        shell:
            '''
            mothur "#dist.seqs(fasta={input}, cutoff=0.3)"
            '''
            
rule clustering:
        input:
            column=expand("{gene}.aligned.good.filter.unique.pick.good.filter.dist", gene=config["gene"]),
            name=expand("{gene}.aligned.good.filter.pick.good.names", gene=config["gene"])
        output:
            expand("{gene}.aligned.good.filter.unique.pick.good.filter.an.list", gene=config["gene"])
        conda:
            "envs/mothur.yaml"
        shell:
            '''
            mothur "#cluster(column={input.column}, name={input.name}, method=average, cutoff=0.3)"
            '''
rule otu_reps:
        input:
            column=expand("{gene}.aligned.good.filter.unique.pick.good.filter.dist", gene=config["gene"]),
            fasta=expand("{gene}.aligned.good.filter.unique.pick.good.filter.fasta", gene=config["gene"]),
            name=expand("{gene}.aligned.good.filter.pick.good.names", gene=config["gene"]),
            list=expand("{gene}.aligned.good.filter.unique.pick.good.filter.an.list", gene=config["gene"])
        output:
            expand("{gene}.aligned.good.filter.unique.pick.good.filter.an.0.11.rep.fasta", gene=config["gene"]),
            expand("{gene}.aligned.good.filter.unique.pick.good.filter.an.0.11.rep.names", gene=config["gene"])
        conda:
            "envs/mothur.yaml"
        shell:
            '''
            mothur "#get.oturep(column={input.column}, fasta={input.fasta}, name={input.name}, list={input.list}, cutoff={config[cutoff_otu]})"
            '''

rule final_database:
        input:
            fasta=expand("{gene}.aligned.good.filter.unique.pick.good.filter.fasta", gene=config["gene"]),
            name=expand("{gene}.aligned.good.filter.pick.good.names", gene=config["gene"])
        output:
            expand("{gene}.aligned.good.filter.unique.pick.good.filter.redundant.fasta", gene=config["gene"])
        conda:
            "envs/mothur.yaml"
        shell:
            '''
            mothur "#deunique.seqs(fasta={input.fasta}, name={input.name})"
            '''            
rule iqtree:
        input:
           expand("{gene}.aligned.good.filter.unique.pick.good.filter.an.0.11.rep.fasta", gene=config["gene"])
        output:
            expand("{gene}.aligned.good.filter.unique.pick.good.filter.an.0.11.rep.fasta.treefile", gene=config["gene"])
        conda:
            "envs/iqtree.yaml"
        threads:10
        shell:
            "iqtree -s {input}  -m MFP -alrt 1000 -bb 1000 -nt {threads}"