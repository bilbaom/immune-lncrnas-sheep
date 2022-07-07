import pandas as pd
import glob
import re
import os

print("This is the project folder: " + config["project-folder"])
print("This is the data folder: " + config["data-folder"])

# Open file with accesion numbers for download (modified from ncbi run selector)
samples = []
with open(config["sample-table"], encoding="utf-8") as f:
    next(f)
    for line in f:
        tabline = line.rsplit(",")
        if tabline[18] == "Reverse":
            samples.append(tabline[2])

print(samples[0:5])
print("number of samples: " + str(len(samples)))

### rule complete pipeline ###

rule all:
    input:
        "%s/GTF_merged/merged_STRG.gtf" % (config["data-folder"]),
        expand("%s/BAM/{samples}_Log.final.out" % (config["data-folder"]), samples=samples)

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
    threads: 1
    run:
        import time
        import pandas as pd
        runcsv = pd.read_csv(input[1])
        runcsv = runcsv[runcsv["Experiment"]==params[0]]
        runcsv["Run"].to_csv(output[0], header=None, index=None, sep='\t')
        time.sleep(int(params[0][-1]))

        
rule download_sra:
    """ Download fastq files from SRA and merge runs"""
    input:
        ids=ancient("%s/ids/{samples}.txt" % (config["data-folder"]))
    output:
        R1="%s/fastq/{samples}_1temp.fastq.gz" % (config["data-folder"]),
        R2="%s/fastq/{samples}_2temp.fastq.gz" % (config["data-folder"])
    params:
        samples="{samples}",
        fastqdir="%s/fastq" % (config["data-folder"]),
        datadir=config["data-folder"]
    log:
        "%s/logs/download_sra/download_sra.{samples}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/download_sra/download_sra.{samples}.benchmark.tsv" % (config["project-folder"])
    threads: 8
    conda:
        "%s/envs/env_atlas_sra_down.yaml" % (config["data-folder"])
    priority: 1
    shell:"""
        mkdir -p {params.fastqdir};
        # esearch -db sra -query {params.samples} | efetch --format runinfo | cut -d ',' -f 1 | grep "RR" > {params.samples}.txt;
        # sleep 5;
        sort -u {input.ids} | parallel -j 2 "prefetch {{}} -O {params.datadir}/{params.samples}";
        echo "Downloaded {params.samples}";
        fastq-dump --split-3 {params.samples}/*/*.sra -O {params.datadir}/{params.samples} > {log} 2>&1;
        echo "format changed to fastq: {params.samples}";
        cat {params.datadir}/{params.samples}/*_1.fastq > {params.datadir}/{params.samples}/{params.samples}_1.fastq;
        cat {params.datadir}/{params.samples}/*_2.fastq > {params.datadir}/{params.samples}/{params.samples}_2.fastq;
        gzip {params.datadir}/{params.samples}/{params.samples}_1.fastq;
        gzip {params.datadir}/{params.samples}/{params.samples}_2.fastq;
        mv {params.datadir}/{params.samples}/{params.samples}_1.fastq.gz {params.fastqdir}/{params.samples}_1temp.fastq.gz;
        mv {params.datadir}/{params.samples}/{params.samples}_2.fastq.gz {params.fastqdir}/{params.samples}_2temp.fastq.gz;

        rm -r {params.datadir}/{params.samples};
    """

# Raw data analysis #
rule cutadapt_trim_reads:
    """
    Trim adapters and low quality reads (CUTADAPT).
    """
    input:
        R1=ancient("%s/fastq/{samples}_1temp.fastq.gz" % (config["data-folder"])),
        R2=ancient("%s/fastq/{samples}_2temp.fastq.gz" % (config["data-folder"]))
    output:
        R1=temp("%s/fastq_trim/{samples}_1_trimtemp.fastq.gz" % (config["data-folder"])),
        R2=temp("%s/fastq_trim/{samples}_2_trimtemp.fastq.gz" % (config["data-folder"]))
    params:
        phread_score=config["params"]["cutadapt"]["phread_score"],
        adapter_file_R1=config["params"]["cutadapt"]["adapter_R1"],
        adapter_file_R2=config["params"]["cutadapt"]["adapter_R2"],
        pair_filter=config["params"]["cutadapt"]["pair_filter"],
        min_length=config["params"]["cutadapt"]["min_length"],
        outdir="%s/fastq_trim" % (config["data-folder"])
    log:
        "%s/logs/cutadapt/cutadapt.{samples}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/cutadapt/cutadapt.{samples}.benchmark.tsv" % (config["project-folder"])
    threads: 12
    conda:
        "%s/envs/env_atlas_qc.yaml" % (config["data-folder"])
    priority: 2
    shell: """
        mkdir -p {params.outdir}
        cutadapt --cores={threads} -q {params.phread_score} -a file:{params.adapter_file_R1} -A file:{params.adapter_file_R2} -o {output.R1} -p {output.R2} --minimum-length={params.min_length} --pair-filter={params.pair_filter} {input} > {log} 2>&1
        rm {input.R1}
        rm {input.R2}
    """


rule build_star_index:
    """ Creates the idex from the genome and annotation necessary for mapping"""
    input:
        fa=config["genome"],
        gtf=config["annotation"]
    output:
        dir="%s/genomeindex" % (config["data-folder"])
    log:
        "%s/logs/starindex/starindex.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/starindex/starindex.benchmark.tsv" % (config["project-folder"])
    threads: 12
    conda:
        "%s/envs/env_atlas_map.yaml" % (config["data-folder"])
    shell:"""
        mkdir -p {output.dir};
        STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir {output.dir} --genomeFastaFiles {input.fa} --sjdbGTFfile {input.gtf} --sjdbOverhang 74 2> {log};
    """

rule star_map:
    """
    Map the samples to the genome (STAR).
    """
    input:
        index="%s/genomeindex" % (config["data-folder"]),
        annotation=config["annotation"],
        R1="%s/fastq_trim/{samples}_1_trimtemp.fastq.gz" % (config["data-folder"]),
        R2="%s/fastq_trim/{samples}_2_trimtemp.fastq.gz" % (config["data-folder"])
    output:
        file=temp("%s/BAM/{samples}.bam" % (config["data-folder"])),
        log="%s/BAM/{samples}_Log.final.out" % (config["data-folder"]),
    params:
        dir=directory("%s/BAM" % (config["data-folder"]))
    log:
        "%s/logs/star_map/star_map.{samples}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/star_map/star_map.{samples}.benchmark.tsv" % (config["project-folder"])
    threads: 12
    conda:
        "%s/envs/env_atlas_map.yaml" % (config["data-folder"])
    priority: 3
    shell:"""
        mkdir -p {params.dir};
        printf \"%s\t%s\t%s\t%s\t%s\t%s\n\" {input.index} {input.annotation} {input.R1} {input.R2} {output} {log} {threads}
        [ ! -d \"{params.dir}\" ] && mkdir {params.dir}
        STAR --genomeDir {input.index} \
            --readFilesIn {input.R1} {input.R2} \
            --readFilesCommand zcat \
            --outFilterType BySJout \
            --outSAMunmapped Within \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMstrandField intronMotif \
            --outSAMattrIHstart 0 \
            --outFilterIntronMotifs RemoveNoncanonical \
            --runThreadN {threads} \
            --quantMode TranscriptomeSAM \
            --outWigType bedGraph \
            --outWigStrand Stranded \
            --alignSoftClipAtReferenceEnds No \
            --outFileNamePrefix {wildcards.samples}_ 2> {log};
        mv {wildcards.samples}_Aligned.sortedByCoord.out.bam {output.file};
        mv {wildcards.samples}_Log.final.out {params.dir};
        rm {wildcards.samples}_Log.progress.out {wildcards.samples}_Signal.UniqueMultiple.str2.out.bg {wildcards.samples}_Signal.Unique.str2.out.bg {wildcards.samples}_Aligned.toTranscriptome.out.bam {wildcards.samples}_Log.out {wildcards.samples}_Signal.UniqueMultiple.str1.out.bg {wildcards.samples}_Signal.Unique.str1.out.bg {wildcards.samples}_SJ.out.tab
    """

rule transcriptome_assembly_stringtie:
    """
    Transcript assembly (StringTie).
    """
    input:
        bam="%s/BAM/{samples}.bam" % (config["data-folder"]),
        annotation=config["annotation"]
    output:
        "%s/Stringtie/{samples}.stringtie.gtf" % (config["data-folder"])
    log:
        "%s/logs/stringtie.{samples}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/stringtie.{samples}.benchmark.tsv" % (config["project-folder"])
    threads: 8
    conda:
        "%s/envs/env_atlas_map.yaml" % (config["data-folder"])
    priority: 4
    shell:"""
        stringtie {input.bam} -p {threads} --rf -G {input.annotation} -v -o {output} 2> {log};
    """

rule merge_samples:
    """
    Merge the gtf files (stringtie merge).
    """
    input:
        gtfs=expand("%s/Stringtie/{samples}.stringtie.gtf" % (config["data-folder"]), samples=samples),
        annotation=config["annotation"]
    output:
        "%s/GTF_merged/merged_STRG.gtf" % (config["data-folder"])
    params:
        tpm=config["params"]["stringtie"]["tpm"]
    log:
        "%s/logs/stringtiemerge.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/stringtiemerge.benchmark.tsv" % (config["project-folder"])
    threads: 8
    conda:
        "%s/envs/env_atlas_map.yaml" % (config["data-folder"])
    shell:"""
        stringtie --merge -G {input.annotation} -F 0 -T {params.tpm} -o {output} {input.gtfs} 2> {log};
    """
