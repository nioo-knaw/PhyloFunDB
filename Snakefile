configfile: "config.yaml"
wildcard_constraints:
   gene = '\w+'

rule final:
    input: expand("results/{gene}.treefile \
                  results/{gene}.aligned.good.filter.unique.pick.good.filter.redundant.fasta \
                   results/{gene}.taxonomy.final.txt".split(), gene=config["gene"])

rule download_ncbi:
    output: 
        "interm/{gene}.fasta"
    conda: 
        "envs/entrez.yaml"
    params:
        gene="{gene}",
        full_name=config["full_name"]
    message: "Retrieving gene sequences from NCBI. Note, this can take a (long) while"
    shell:"""
           esearch -db nucleotide -query "{params.gene}[gene]" | \
           efetch -format gpc | \
           xtract -pattern INSDFeature -if INSDFeature_key -equals CDS -and INSDQualifier_value -equals {params.gene} -or INSDQualifier_value -contains '{params.full_name}' -element INSDInterval_accession -element INSDInterval_from -element INSDInterval_to | \
           sort -u -k1,1 | \
           uniq | xargs -n 3 sh -c 'efetch -db nuccore -id "$0" -seq_start "$1" -seq_stop "$2" -format fasta' > {output}
           """

rule download_taxonomy:
    output:
        "interm/tax_{gene}.txt"
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
        "interm/{gene}.fasta"
    output:
        "interm/{gene}.renamed.fasta"
    shell:
       """sed -e 's/[.]/-/' -e 's/ /-/g' -e 's/_//g' {input} | \
       stdbuf -o0 cut -d "-" -f 1,4,5| \
       sed -e 's/-/_/g' -e 's/[.]//g' -e 's/[,]//g' > {output}
       """
       
rule get_unverified:
    input:
        "interm/{gene}.renamed.fasta"
    output:
        "interm/{gene}.unverified.accnos"
    shell:
        """grep UNVERIFIED {input} | \
        stdbuf -o0 cut -c 2- > {output}
        """
        
rule remove_unverified:
        input:
            accnos="interm/{gene}.unverified.accnos",
            fasta="interm/{gene}.renamed.fasta"
        output:
            "interm/{gene}.renamed.pick.fasta"
        conda:
            "envs/mothur.yaml"
        shell:
            '''
            mothur "#remove.seqs(accnos={input.accnos}, fasta={input.fasta})"
            '''
            
rule mothur_trim:
        input:
            "interm/{gene}.renamed.pick.fasta"
        output:
            "interm/{gene}.renamed.pick.trim.fasta"
        conda:
            "envs/mothur.yaml"
        params:
            minlength=config['minlength']
        threads:10
        shell:
            '''
            mothur "#trim.seqs(fasta={input}, minlength={params.minlength}, maxambig=0, processors={threads})"
            '''
rule framebot:
        input:
            fasta="interm/{gene}.renamed.pick.trim.fasta",
            db_framebot="dbs/{gene}.fungene.clean.fasta"
        output:
            "interm/{gene}.framebot_corr_nucl.fasta"
        params:
            "interm/{gene}.framebot"
        conda:
            "envs/rdptools.yaml"
        shell:
            "FrameBot framebot -o {params} -N {input.db_framebot} {input.fasta}"
       
rule align:
        input:
            "interm/{gene}.framebot_corr_nucl.fasta"
        output:
            "interm/{gene}.aligned.fasta"
        conda:
            "envs/mafft.yaml"
        threads:10
        shell:
            "mafft --thread {threads} --auto {input} >{output}"
            
rule screen_alignment:
        input:
            "interm/{gene}.aligned.fasta"
        output:
            "interm/{gene}.aligned.good.fasta"
        conda:
            "envs/mothur.yaml"
        threads:10
        shell:
            '''
            mothur "#screen.seqs(fasta={input}, optimize=start-end, criteria=96, processors={threads})"
            '''
rule filter_alignment:
        input:
            "interm/{gene}.aligned.good.fasta"
        output:
            "interm/{gene}.aligned.good.filter.fasta"
        conda:
            "envs/mothur.yaml"
        threads:10
        shell:
            '''
            mothur "#filter.seqs(fasta={input}, vertical=T, trump=., processors={threads})"
            '''
rule unique1:
        input:
            "interm/{gene}.aligned.good.filter.fasta"
        output:
            "interm/{gene}.aligned.good.filter.unique.fasta",
            "interm/{gene}.aligned.good.filter.names"
        conda:
            "envs/mothur.yaml"
        shell:
            '''
            mothur "#unique.seqs(fasta={input})"
            '''            
            
rule chimera_vsearch:
        input:
            fasta="interm/{gene}.aligned.good.filter.unique.fasta",
            name="interm/{gene}.aligned.good.filter.names"
        output:
            "interm/{gene}.aligned.good.filter.unique.denovo.vsearch.accnos"
        conda:
            "envs/mothur.yaml"
        shell:
            '''
            mothur "#chimera.vsearch(fasta={input.fasta}, name={input.name})"
            '''
            
rule chimera_removal:
        input:
            accnos="interm/{gene}.aligned.good.filter.unique.denovo.vsearch.accnos",
            fasta="interm/{gene}.aligned.good.filter.unique.fasta",
            name="interm/{gene}.aligned.good.filter.names"
        output:
            "interm/{gene}.aligned.good.filter.pick.names",
            "interm/{gene}.aligned.good.filter.unique.pick.fasta"
        conda:
            "envs/mothur.yaml"
        shell:
            '''
            mothur "#remove.seqs(accnos={input.accnos}, fasta={input.fasta}, name={input.name})"
            '''
            
rule screen_alignment2:
        input:
            fasta="interm/{gene}.aligned.good.filter.unique.pick.fasta",
            name="interm/{gene}.aligned.good.filter.pick.names"
        output:
            "interm/{gene}.aligned.good.filter.unique.pick.good.fasta",
            "interm/{gene}.aligned.good.filter.pick.good.names"
        conda:
            "envs/mothur.yaml"
        threads:10
        shell:
            '''
            mothur "#screen.seqs(fasta={input.fasta}, name={input.name}, optimize=start-end, criteria=96, processors={threads})"
            '''
            
rule filter_alignment2:
        input:
            "interm/{gene}.aligned.good.filter.unique.pick.good.fasta"
        output:
            "interm/{gene}.aligned.good.filter.unique.pick.good.filter.fasta"
        conda:
            "envs/mothur.yaml"
        threads:10
        shell:
            '''
            mothur "#filter.seqs(fasta={input}, vertical=T, trump=., processors={threads})"
            '''
            
rule distance_matrix:
        input:
            "interm/{gene}.aligned.good.filter.unique.pick.good.filter.fasta"
        output:
            "interm/{gene}.aligned.good.filter.unique.pick.good.filter.dist"
        conda:
            "envs/mothur.yaml"
        threads:10
        shell:
            '''
            mothur "#dist.seqs(fasta={input}, cutoff=0.35, processors={threads})"
            '''
            
rule clustering:
        input:
            column="interm/{gene}.aligned.good.filter.unique.pick.good.filter.dist",
            name="interm/{gene}.aligned.good.filter.pick.good.names"
        output:
            "interm/{gene}.aligned.good.filter.unique.pick.good.filter.an.list"
        conda:
            "envs/mothur.yaml"
        shell:
            '''
            mothur "#cluster(column={input.column}, name={input.name}, method=average, cutoff=0.35)"
            '''
rule otu_reps:
        input:
            column="interm/{gene}.aligned.good.filter.unique.pick.good.filter.dist",
            fasta="interm/{gene}.aligned.good.filter.unique.pick.good.filter.fasta",
            name="interm/{gene}.aligned.good.filter.pick.good.names",
            list="interm/{gene}.aligned.good.filter.unique.pick.good.filter.an.list"
        output:
            "interm/{gene}.aligned.good.filter.unique.pick.good.filter.an.0.11.rep.fasta",
            "interm/{gene}.aligned.good.filter.unique.pick.good.filter.an.0.11.rep.names"
        conda:
            "envs/mothur.yaml"
        shell:
            '''
            mothur "#get.oturep(column={input.column}, fasta={input.fasta}, name={input.name}, list={input.list}, cutoff={config[cutoff_otu]})"
            '''

rule final_database:
        input:
            fasta="interm/{gene}.aligned.good.filter.unique.pick.good.filter.fasta",
            name="interm/{gene}.aligned.good.filter.pick.good.names"
        output:
            "results/{gene}.aligned.good.filter.unique.pick.good.filter.redundant.fasta"
        conda:
            "envs/mothur.yaml"
        shell:
            '''
            mothur "#deunique.seqs(fasta={input.fasta}, name={input.name}, outputdir=./results)"
            '''
rule tax_format:
        input:
            "interm/tax_{gene}.txt"
        output:
            "interm/newtax_{gene}.txt"
        shell:
            """
            sed -e 's/_//g' -e 's/ //g' -e 's/$/;/' {input} > {output}
            """

rule get_fasta_names:
        input:
            "results/{gene}.aligned.good.filter.unique.pick.good.filter.redundant.fasta"
        output:
            "interm/{gene}_fasta_final.names.txt"
        shell:
            """
            grep ">" {input} | stdbuf -o0 cut -c 2- > {output}
            """

rule tax_format_final:
        input:
            fasta="interm/{gene}_fasta_final.names.txt",
            tax="interm/newtax_{gene}.txt"
        output:
            final_tax="results/{gene}.taxonomy.final.txt"
        conda:
            "envs/tidyr.yaml"
        script:
            "renaming.R"
          
rule iqtree:
        input:
           "interm/{gene}.aligned.good.filter.unique.pick.good.filter.an.0.11.rep.fasta"
        output:
            "results/{gene}.treefile"
        params:
            "results/{gene}"
        conda:
            "envs/iqtree.yaml"
        threads:10
        shell:
            "iqtree -s {input}  -m MFP -alrt 1000 -bb 1000 -nt {threads} -pre {params}"
