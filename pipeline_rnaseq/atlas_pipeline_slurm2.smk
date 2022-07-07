import pandas as pd
# Open file with accesion numbers for download (modified from ncbi run selector)
# Separate stranded and unstranded samples
samplesst = []
samplesun = []
samples = []
stranddict = {}
with open(config["sample-table"], encoding="utf-8") as f:
    next(f)
    for line in f:
        tabline = line.rsplit(",")
        samples.append(tabline[2])
        if tabline[18] == "Reverse":
            samplesst.append(tabline[2])
            stranddict[tabline[2]] = "--rf-stranded"
        if tabline[18] == "Unstranded":
            samplesun.append(tabline[2])
            stranddict[tabline[2]] = ""

print("This is the project folder: " + config["project-folder"])
print("This is the data folder: " + config["data-folder"])
print("number of reverse stranded samples: " + str(len(samplesst)))
print("number of unstranded samples: " + str(len(samplesun)))
print("total number of samples: " + str(len(samples)))
print(len(stranddict))

### rule complete pipeline ###

rule all:
    input:
        "%s/kallisto/kallisto_counts.tsv" % (config["data-folder"])

### setup report ###

report: "%s/report/workflow.rst" % (config["project-folder"])

### rules ###
# Download samples from SRA for transcriptome assembly #
rule get_ids:
    """ Get every SRA id for each EXP id"""
    input:
        index="%s/genomeindex" % (config["data-folder"]),
        csv=ancient(config["run-table"])
    output:
        "%s/ids/{samples}.txt" % (config["data-folder"])
    params:
        samples="{samples}",
    log:
        "%s/logs/getids_q/getids_q.{samples}.log" % (config["project-folder"])
    threads: 1
    run:
        import time
        import pandas as pd
        runcsv = pd.read_csv(input[1])
        runcsv = runcsv[runcsv["Experiment"]==params[0]]
        runcsv["Run"].to_csv(output[0], header=None, index=None, sep='\t')
        time.sleep(int(params[0][-1]))
        
rule build_transcriptome:
    """
    Build fasta file with transcripts for kallisto pseudoalignment quantification
    """
    input:
        gtf=ancient("%s/GTF_merged/merged_STRG.gtf" % (config["data-folder"])),
        fa=config["genome"]
    output:
        "%s/GTF_merged/merged_STRG.fasta" % (config["data-folder"])
    log:
        "%s/logs/build_transcriptome.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/build_transcriptome.benchmark.tsv" % (config["project-folder"])
    threads: 8
    conda:
        "%s/envs/env_atlas_map.yaml" % (config["data-folder"])
    shell:"""
        gffread {input.gtf} -g {input.fa} -w {output} 2> {log};
    """

# Download files for quantification #

rule download_sra_quant:
    """ Download fastq files from SRA and merge runs"""
    input:
        ids=ancient("%s/ids/{samples}.txt" % (config["data-folder"]))
    output:
        R1="%s/fastq/{samples}_1.fastq.gz" % (config["data-folder"]),
        R2="%s/fastq/{samples}_2.fastq.gz" % (config["data-folder"])
    params:
        samples="{samples}",
        fastqdir="%s/fastq" % (config["data-folder"]),
        datadir=config["data-folder"]
    log:
        "%s/logs/download_sra_quant/download_sra_quant.{samples}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/download_sra_quant/download_sra_quant.{samples}.benchmark.tsv" % (config["project-folder"])
    threads: 8
    conda:
        "%s/envs/env_atlas_sra_down.yaml" % (config["data-folder"])
    priority: 1
    shell:"""
        mkdir -p {params.fastqdir};
        # esearch -db sra -query {params.samples} | efetch --format runinfo | cut -d ',' -f 1 | grep "RR" > {params.samples}.txt;
        # sleep 5
        sort -u {input.ids} | parallel -j 2 "prefetch {{}} -O {params.datadir}/{params.samples}";
        echo "Downloaded {params.samples}";
        fastq-dump --split-3 {params.samples}/*/*.sra -O {params.datadir}/{params.samples} > {log} 2>&1;
        echo "format changed to fastq: {params.samples}";
        cat {params.datadir}/{params.samples}/*_1.fastq > {params.datadir}/{params.samples}/{params.samples}_1.fastq;
        cat {params.datadir}/{params.samples}/*_2.fastq > {params.datadir}/{params.samples}/{params.samples}_2.fastq;
        pigz {params.datadir}/{params.samples}/{params.samples}_1.fastq;
        pigz {params.datadir}/{params.samples}/{params.samples}_2.fastq;
        mv {params.datadir}/{params.samples}/{params.samples}_1.fastq.gz {params.fastqdir}/{params.samples}_1.fastq.gz;
        mv {params.datadir}/{params.samples}/{params.samples}_2.fastq.gz {params.fastqdir}/{params.samples}_2.fastq.gz;

        rm -r {params.datadir}/{params.samples};
    """

rule cutadapt_trim_reads_quant:
    """
    Trim adapters and low quality reads (CUTADAPT).
    """
    input:
        R1=ancient("%s/fastq/{samples}_1.fastq.gz" % (config["data-folder"])),
        R2=ancient("%s/fastq/{samples}_2.fastq.gz" % (config["data-folder"]))
    output:
        R1=temp("%s/fastq_trim/{samples}_1_trim.fastq.gz" % (config["data-folder"])),
        R2=temp("%s/fastq_trim/{samples}_2_trim.fastq.gz" % (config["data-folder"]))
    params:
        phread_score=config["params"]["cutadapt"]["phread_score"],
        adapter_file_R1=config["params"]["cutadapt"]["adapter_R1"],
        adapter_file_R2=config["params"]["cutadapt"]["adapter_R2"],
        pair_filter=config["params"]["cutadapt"]["pair_filter"],
        min_length=config["params"]["cutadapt"]["min_length"],
        outdir="%s/fastq_trim" % (config["data-folder"])
    log:
        "%s/logs/cutadapt_quant/cutadapt_quant.{samples}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/cutadapt_quant/cutadapt_quant.{samples}.benchmark.tsv" % (config["project-folder"])
    threads: 8
    conda:
        "%s/envs/env_atlas_qc.yaml" % (config["data-folder"])
    priority: 2
    shell: """
        mkdir -p {params.outdir}
        cutadapt --cores={threads} -q {params.phread_score} -a file:{params.adapter_file_R1} -A file:{params.adapter_file_R2} -o {output.R1} -p {output.R2} --minimum-length={params.min_length} --pair-filter={params.pair_filter} {input} > {log} 2>&1
        rm {input.R1}
        rm {input.R2}
    """

rule kallisto_index:
    """
    Build transcriptome index for kallisto pseudoalignment quantification
    """
    input:
        "%s/GTF_merged/merged_STRG.fasta" % (config["data-folder"])
    output:
        idx="%s/GTF_merged/merged_STRG.fasta.idx" % (config["data-folder"])
    params:
        idx="merged_STRG.fasta.idx",
        dir="%s/GTF_merged" % (config["data-folder"])
    log:
        "%s/logs/kallistoindex.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/kallistoindex.benchmark.tsv" % (config["project-folder"])
    threads: 8
    conda:
        "%s/envs/env_atlas_quant.yaml" % (config["data-folder"])
    shell:"""
        kallisto index --index={params.idx} {input} 2> {log};
        mv {params.idx} {params.dir}
    """

rule quantify_kallisto:
    """
    Quantify transcript expression with kallisto (--rf-stranded for stranded)
    """
    input:
        R1="%s/fastq_trim/{samples}_1_trim.fastq.gz" % (config["data-folder"]),
        R2="%s/fastq_trim/{samples}_2_trim.fastq.gz" % (config["data-folder"]),
        idx="%s/GTF_merged/merged_STRG.fasta.idx" % (config["data-folder"]),
    output:
        "%s/kallisto/{samples}/abundance.tsv" % (config["data-folder"])
    params:
        dir="%s/kallisto/{samples}" % (config["data-folder"]),
        strand=lambda wcs: stranddict[wcs.samples]
    log:
        "%s/logs/kallisto.{samples}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/kallisto.{samples}.benchmark.tsv" % (config["project-folder"])
    threads: 4
    conda:
        "%s/envs/env_atlas_quant.yaml" % (config["data-folder"])
    priority: 3
    shell:"""
        mkdir -p {params.dir};
        kallisto quant {params.strand} --bias -i {input.idx} -o {params.dir} {input.R1} {input.R2} 2> {log};
    """


rule merge_kallisto_output:
    """
    Merge kallisto output abundance.tsv files into single matrix.
    """
    input:
        expand("%s/kallisto/{samples}/abundance.tsv" % (config["data-folder"]), samples=samples)
    output:
        "%s/kallisto/kallisto_counts.tsv" % (config["data-folder"])
    params:
        smDir="%s/kallisto/" % (config["data-folder"]),
        pyDir=config["reference-folder"]
    log:
        "%s/logs/merge_kallisto.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/merge_kallisto.benchmark.tsv" % (config["project-folder"])
    conda:
        "%s/envs/env_atlas_quant.yaml" % (config["data-folder"])
    shell:"""
        python {params.pyDir}/merge_kallisto.py -f {params.smDir} -o {output} > {log} 2>&1
    """
