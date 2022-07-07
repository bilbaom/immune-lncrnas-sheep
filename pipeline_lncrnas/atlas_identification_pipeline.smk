import pandas as pd
import glob
import re
import os

### rule complete pipeline ###

rule all:
    input:
        "%s/classification_table.txt" % (config["project-folder"])
        
### setup report ###

report: "%s/report/workflow.rst" % (config["project-folder"])

### rules ###
# custom lncRNA rules #

rule gffcompare_compare_with_annotation:
    """
    Classify all transcripts based on their location relative to the reference annotation.
    """
    input:
        gtf="%s/GTF_merged/merged_STRG.gtf" % (config["project-folder"]),
        annotation=config["annotation"]
    output:
        gtf="%s/gffcompare/merged_STRG_master.annotated.gtf" % (config["project-folder"]),
        list="%s/gffcompare/list.txt" % (config["project-folder"])
    params:
        name="merged_STRG_master",
        dir="%s/gffcompare" % (config["project-folder"])
    log:
        "%s/logs/gffcompare/gffcompare.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/gffcompare.benchmark.tsv" % (config["project-folder"])
    shell:"""
        echo "{input.gtf}" > {output.list}
        gffcompare -i {output.list} -r {input.annotation} -o {params.name} &> {log}
        mkdir -p {params.dir}
        mv merged_STRG_master.* {params.dir}
    """

rule filter_lncrnas:
    """
    Extract lncRNA candidates from gtf file.
    """
    input:
        gtf="%s/gffcompare/merged_STRG_master.annotated.gtf" % (config["project-folder"])
    output:
        gtf="%s/potential_lncrna.gtf" % (config["project-folder"]),
        ids="%s/potentialids.txt" % (config["project-folder"])
    params:
        dir=config["project-folder"],
        scripts=config["script-folder"]
    log:
        "%s/logs/custom_lnc/filter.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/custom_filter.benchmark.tsv" % (config["project-folder"])
    shell:"""
        mkdir -p {params.dir}
        python {params.scripts}/filter2.py {input.gtf} {params.dir} &> {log}
    """

rule get_transcript_sequences:
    """
    GTF to FASTA transcripts
    """
    input:
        gtf="%s/potential_lncrna.gtf" % (config["project-folder"]),
        genome=config["genome"]
    output:
        fa="%s/potential_lncrna.fa" % (config["project-folder"])
    params:
        dir=config["project-folder"]
    log:
        "%s/logs/custom_lnc/get_sequences.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/get_sequences.benchmark.tsv" % (config["project-folder"])
    shell:"""
        gffread {input.gtf} -w {output.fa} -g {input.genome} &> {log}
    """

rule translate_sequences:
    """
    translate sequences in the 3 possible frames (for HMMER)
    """
    input:
        "%s/potential_lncrna.fa" % (config["project-folder"])
    output:
        "%s/potential_lncrna_tr.fa" % (config["project-folder"])
    params:
        scripts=config["script-folder"]
    log:
        "%s/logs/custom_lnc/translate.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/translate.benchmark.tsv" % (config["project-folder"])
    shell:"""
        python {params.scripts}/translate3frames.py {input} {output} &> {log}
    """

rule run_hmmer_press:
    """
    Prepare an HMM database for hmmscan (download of Pfam library needed)
    """
    input:
        config["params"]["hmmer"]["pfam"]
    output:
        "%s/Pfam-A.hmm.h3i" % (config["ref-folder"])
    log:
        "%s/logs/custom_lnc/hmmer_press.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/hmmer_press.benchmark.tsv" % (config["project-folder"])
    shell:"""
        hmmpress {input} &> {log}
    """

rule run_hmmer:
    """
    Run HMMER to find protein domains (download of Pfam library needed)
    """
    input:
        prot="%s/potential_lncrna_tr.fa" % (config["project-folder"]),
        press="%s/Pfam-A.hmm.h3i" % (config["ref-folder"])
    output:
        "%s/hmmer_out.txt" % (config["project-folder"])
    params:
        pfam=config["params"]["hmmer"]["pfam"],
        eval=config["params"]["hmmer"]["eval"]
    log:
        "%s/logs/custom_lnc/hmmer.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/hmmer.benchmark.tsv" % (config["project-folder"])
    shell:"""
        hmmscan --tblout {output} -E {params.eval} {params.pfam} {input.prot} &> {log}
    """

rule run_cpat:
    """
    Run CPAT (calculate coding potential)
    Model needs to be trained before to obtain logit.RData and hexamer.tsv
    """
    input:
        fa="%s/potential_lncrna.fa" % (config["project-folder"]),
        logit=config["params"]["cpat"]["logit-data"],
        hexamer=config["params"]["cpat"]["hexamer-data"]
    output:
        "%s/cpat_results.txt" % (config["project-folder"])
    params:
        name="cpat_results"
    log:
        "%s/logs/custom_lnc/cpat.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/cpat.benchmark.tsv" % (config["project-folder"])
    shell:"""
        cpat.py -g {input.fa} -o {params.name} -x {input.hexamer} -d {input.logit} &> {log}
        cat {params.name}.ORF_prob.best.tsv {params.name}.no_ORF.txt > {output}
    """

rule classify_ncrna_feelnc_02_codingpot:
    """
    lncRNA coding potential (FEELnc).
    """
    input:
        candidates="%s/potential_lncrna.gtf" % (config["project-folder"]),
        known_mrna=config["params"]["feelnc"]["coding-seq"],
        known_lnc=config["params"]["feelnc"]["noncoding-seq"]
    output:
        out="%s/custom_feelnc.lncRNA.gtf" % (config["project-folder"]),
    params:
        genome=config["genome"],
        outf=config["project-folder"],
        outname="custom_feelnc"
    log:
        "%s/logs/custom_lnc/feelnc_02.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/feelnc_02.benchmark.tsv" % (config["project-folder"])
    conda:
        "%s/env_feelnc.yaml" % (config["project-folder"])
    shell:"""
        FEELnc_codpot.pl -i {input.candidates} -a {input.known_mrna} -l {input.known_lnc} -g {params.genome} -o {params.outname} --outdir={params.outf} --keeptmp &> {log}
    """

rule merge_codpot:
    """
    Get a list of trustworthy lncRNA transcripts merging results of HMMER, CPAT and FEELnc
    """
    input:
        hmmer="%s/hmmer_out.txt" % (config["project-folder"]),
        cpat="%s/cpat_results.txt" % (config["project-folder"]),
        feelnc="%s/custom_feelnc.lncRNA.gtf" % (config["project-folder"])
    output:
        ids="%s/lnc_ids.txt" % (config["project-folder"]),
        venn="%s/codpot_venn.png" % (config["project-folder"])
    params:
        scripts=config["script-folder"]
    log:
        "%s/logs/custom_lnc/merge_codplot.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/merge_codplot.benchmark.tsv" % (config["project-folder"])
    shell:"""
        python {params.scripts}/cpvenn.py {input.cpat} {input.feelnc} {input.hmmer} {output.ids} {output.venn} &> {log}
    """

rule classify_lncrnas:
    """
    Classify lncRNAs (python script).
    """
    input:
        lnclist="%s/lnc_ids.txt" % (config["project-folder"]),
        gtf="%s/gffcompare/merged_STRG_master.annotated.gtf" % (config["project-folder"]),
        annotation=config["annotation"]
    output:
        txt="%s/classification_table.txt" % (config["project-folder"])
    params:
        dir=config["project-folder"],
        scripts=config["script-folder"]
    log:
        "%s/logs/custom_lnc/classify.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/custom_classify.benchmark.tsv" % (config["project-folder"])
    shell:"""
        python {params.scripts}/classification_lncrnas.py {input.lnclist} {input.gtf} {input.annotation} {params.dir} &> {log}
    """
