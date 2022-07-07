#Bioprojects:
#Spleen: PRJEB41457 (2 Male samples) Paired-end -> THIS DATA WAS NOT USED, THERE WERE METADATA ERRORS
#Alveolar macrophagues: PRJEB40528 Single-end

## Download data ################
# Download Spleen data (male samples)
for sra in ERR4896399 ERR4896400 ERR4896401 ERR4896402 ERR4896403 ERR4896404 ERR4896405 ERR4896406 ERR4896407 ERR4896408 ERR4896409 ERR4896410
do
prefetch $sra
fastq-dump --split-3 "$sra".sra
rm "$sra".sra
gzip "$sra"_1.fastq
gzip "$sra"_2.fastq
echo "Finished $sra"
done

# Download macrophage data
for sra in ERR4626957 ERR4626958 ERR4626959 ERR4626960 ERR4626961 ERR4626962 ERR4626963 ERR4626964 ERR4626965 ERR4626966 ERR4626967 ERR4626968
do
prefetch $sra
fastq-dump "$sra".sra
rm "$sra".sra
gzip "$sra".fastq
echo "Finished $sra"
done

## BWA ################
# Create index
bwa index -p index_bwa -a bwtsw /home/labo/datuak/atlas/references/Ovis_aries_rambouillet.Oar_rambouillet_v1.0.dna.toplevel.fa


# Map Spleen data
# ERR4896399
for sra in ERR4896400 ERR4896401 ERR4896402 ERR4896403 ERR4896404 ERR4896405 ERR4896406 ERR4896407 ERR4896408 ERR4896409 ERR4896410
do
bwa mem -t 6 index_bwa fastq_spleen/"$sra"_1.fastq.gz fastq_spleen/"$sra"_2.fastq.gz | samtools view -S -h -b -F 1804 -f 2 -q 30 /dev/stdin | samtools sort -n --threads 6 -o fastq_spleen/"$sra".bam
done

# frogak
bwa mem -t 6 index_bwa fastq_spleen/ERR4896399_1.fastq.gz fastq_spleen/ERR4896399_2.fastq.gz | samtools view -S -h -b -F 1804 -f 2 -q 30 /dev/stdin | samtools sort -n --threads 6 -o fastq_spleen/ERR4896399.bam

samtools fixmate --threads 7 -r ERR4896399.bam ERR4896399_fixmate.bam
samtools view -h -b -F 1804 -f 2 ERR4896399_fixmate.bam | samtools sort --threads 7 -o ERR4896399_1.bam

picard MarkDuplicates INPUT=ERR4896399_1.bam OUTPUT=/dev/stdout METRICS_FILE=ERR4896399_marked_dup_metrics.txt VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=false QUIET=true COMPRESSION_LEVEL=0 | samtools view -h -b -F 1804 -f 2 /dev/stdin > ERR4896399_processed.bam

samtools fixmate --threads 7 -r ERR4896404.bam ERR4896404_fixmate.bam
samtools view -h -b -F 1804 -f 2 ERR4896404_fixmate.bam | samtools sort --threads 7 -o ERR4896404_1.bam

picard MarkDuplicates INPUT=ERR4896404_1.bam OUTPUT=/dev/stdout METRICS_FILE=ERR4896404_marked_dup_metrics.txt VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=false QUIET=true COMPRESSION_LEVEL=0 | samtools view -h -b -F 1804 -f 2 /dev/stdin > ERR4896404_processed.bam

macs2 callpeak -t ERR4896399_processed.bam -c ERR4896404_processed.bam -f BAMPE -g 2.62e+9 -n "ERR4896399_processed" -q 0.05

# Map macrophage data
for sra in ERR4626957 ERR4626958 ERR4626959 ERR4626960 ERR4626961 ERR4626962 ERR4626963 ERR4626964 ERR4626965 ERR4626966 ERR4626967 ERR4626968
do
bwa mem -t 6 index_bwa fastq_macr/"$sra".fastq.gz | samtools view -S -h -b -F 1804 -f 2 -q 30 /dev/stdin | samtools sort -n --threads 6 -o fastq_macr/"$sra".bam
done

## Bowtie2 #############
bowtie2-build /home/labo/datuak/atlas/references/Ovis_aries_rambouillet.Oar_rambouillet_v1.0.dna.toplevel.fa bowtie2_index

# Map Spleen data
for sra in ERR4896399 ERR4896400 ERR4896401 ERR4896402 ERR4896403 ERR4896404 ERR4896405 ERR4896406 ERR4896407 ERR4896408 ERR4896409 ERR4896410
do
bowtie2 -p 6 -q --local -x bowtie2_index -1 fastq_spleen/"$sra"_1.fastq.gz -2 fastq_spleen/"$sra"_2.fastq.gz -S bam_spleen/"$sra"_unsorted.sam &> bam_log/"$sra".log
samtools view -S -h -b bam_spleen/"$sra"_unsorted.sam > bam_spleen/"$sra"_unsorted.bam
sambamba sort -t 6 -o bam_spleen/"$sra"_sorted.bam bam_spleen/"$sra"_unsorted.bam
sambamba view -h -t 6 -f bam -F "[XS] == null and not unmapped and not duplicate" -o bam_spleen/"$sra"_filtered.bam bam_spleen/"$sra"_sorted.bam
samtools index bam_spleen/"$sra"_filtered.bam
rm bam_spleen/"$sra"_unsorted.sam
rm bam_spleen/"$sra"_sorted.bam
rm bam_spleen/"$sra"_sorted.bam.bai
done

# Map macrophage data
for sra in ERR4626957 ERR4626958 ERR4626959 ERR4626960 ERR4626961 ERR4626962 ERR4626963 ERR4626964 ERR4626965 ERR4626966 ERR4626967 ERR4626968
do
bowtie2 -p 4 -q --local -x bowtie2_index -U fastq_macr/"$sra".fastq.gz -S bam_macr/"$sra"_unsorted.sam &> bam_log/"$sra".log
samtools view -S -h -b bam_macr/"$sra"_unsorted.sam > bam_macr/"$sra"_unsorted.bam
sambamba sort -t 4 -o bam_macr/"$sra"_sorted.bam bam_macr/"$sra"_unsorted.bam
sambamba view -h -t 4 -f bam -F "[XS] == null and not unmapped and not duplicate" -o bam_macr/"$sra"_filtered.bam bam_macr/"$sra"_sorted.bam
samtools index bam_macr/"$sra"_filtered.bam
rm bam_macr/"$sra"_unsorted.sam
rm bam_macr/"$sra"_sorted.bam
rm bam_macr/"$sra"_sorted.bam.bai
done

## Chromap ##############

# Create indexes
chromap -i -r /home/labo/datuak/atlas/references/Ovis_aries_rambouillet.Oar_rambouillet_v1.0.dna.toplevel.fa -o index_spleen --min-frag-length 75
chromap -i -r /home/labo/datuak/atlas/references/Ovis_aries_rambouillet.Oar_rambouillet_v1.0.dna.toplevel.fa -o index_macr --min-frag-length 50

# Map Spleen data
for sra in ERR4896399 ERR4896400 ERR4896401 ERR4896402 ERR4896403 ERR4896404 ERR4896405 ERR4896406 ERR4896407 ERR4896408 ERR4896409 ERR4896410
do
chromap --preset chip --trim-adapters -t 4 -x index_spleen -r /home/labo/datuak/atlas/references/Ovis_aries_rambouillet.Oar_rambouillet_v1.0.dna.toplevel.fa -1 "$sra"_1.fastq.gz -2 "$sra"_2.fastq.gz -o "$sra".bed
done

# Map macrophage data
for sra in ERR4626957 ERR4626958 ERR4626959 ERR4626960 ERR4626961 ERR4626962 ERR4626963 ERR4626964 ERR4626965 ERR4626966 ERR4626967 ERR4626968
do
chromap --preset chip --trim-adapters -t 4 -x index_macr -r /home/labo/datuak/atlas/references/Ovis_aries_rambouillet.Oar_rambouillet_v1.0.dna.toplevel.fa -1 "$sra".fastq.gz -o "$sra".bed
done

## Quality check of mappings ##############

multiBamSummary bins --bamfiles ERR4896399_filtered.bam ERR4896400_filtered.bam ERR4896401_filtered.bam ERR4896402_filtered.bam ERR4896403_filtered.bam ERR4896404_filtered.bam ERR4896405_filtered.bam ERR4896406_filtered.bam ERR4896407_filtered.bam ERR4896408_filtered.bam ERR4896409_filtered.bam ERR4896410_filtered.bam -o results.npz

plotCorrelation --corData results.npz --corMethod pearson --skipZeros --whatToPlot heatmap -o heatmap_spleen.png --outFileCorMatrix PearsonCorr.tab

multiBamSummary bins --bamfiles ERR4626957_filtered.bam ERR4626958_filtered.bam ERR4626959_filtered.bam ERR4626960_filtered.bam ERR4626961_filtered.bam ERR4626962_filtered.bam ERR4626963_filtered.bam ERR4626964_filtered.bam ERR4626965_filtered.bam ERR4626966_filtered.bam ERR4626967_filtered.bam ERR4626968_filtered.bam -o results.npz

plotCorrelation --corData results.npz --corMethod pearson --skipZeros --whatToPlot heatmap -o heatmap_macr.png --outFileCorMatrix PearsonCorr.tab

## Peak calling MACS2 ###############

# Call peaks spleen animal 1
for bam in ERR4896399 ERR4896400 ERR4896403 ERR4896401 ERR4896402
do
macs2 callpeak -t bam_spleen/"$bam"_filtered.bam -f BAMPE -g 2.62e+9 -n "$bam"_noimput -q 0.05 --outdir peaks_spleen
macs2 callpeak -t bam_spleen/"$bam"_filtered.bam -c bam_spleen/ERR4896404_filtered.bam -f BAMPE -g 2.62e+9 -n "$bam" -q 0.05 --outdir peaks_spleen
macs2 callpeak -t bam_spleen/"$bam"_filtered.bam -f BAMPE -g 2.62e+9 -n "$bam"_noimput -q 0.05 --broad --outdir peaks_spleen
macs2 callpeak -t bam_spleen/"$bam"_filtered.bam -c bam_spleen/ERR4896404_filtered.bam -f BAMPE -g 2.62e+9 -n "$bam" -q 0.05 --broad --outdir peaks_spleen
done

# Call peaks spleen animal 2
for bam in ERR4896405 ERR4896406 ERR4896409 ERR4896407 ERR4896408 
do
macs2 callpeak -t bam_spleen/"$bam"_filtered.bam -f BAMPE -g 2.62e+9 -n "$bam"_noimput -q 0.05 --outdir peaks_spleen
macs2 callpeak -t bam_spleen/"$bam"_filtered.bam -c bam_spleen/ERR4896410_filtered.bam -f BAMPE -g 2.62e+9 -n "$bam" -q 0.05 --outdir peaks_spleen
macs2 callpeak -t bam_spleen/"$bam"_filtered.bam -f BAMPE -g 2.62e+9 -n "$bam"_noimput -q 0.05 --broad --outdir peaks_spleen
macs2 callpeak -t bam_spleen/"$bam"_filtered.bam -c bam_spleen/ERR4896410_filtered.bam -f BAMPE -g 2.62e+9 -n "$bam" -q 0.05 --broad --outdir peaks_spleen
done

# Call peaks macrophage animal 1 
for bam in ERR4626957 ERR4626959 ERR4626962 ERR4626960 ERR4626961
do
macs2 callpeak -t bam_macr/"$bam"_filtered.bam -c bam_macr/ERR4626958_filtered.bam -f BAM -g 2.62e+9 -n "$bam" -q 0.05 --outdir peaks_macr
macs2 callpeak -t bam_macr/"$bam"_filtered.bam -c bam_macr/ERR4626958_filtered.bam -f BAM -g 2.62e+9 -n "$bam" -q 0.05 --broad --outdir peaks_macr
done

# Call peaks macrophage animal 2 
for bam in ERR4626963 ERR4626964 ERR4626965 ERR4626968 ERR4626966 ERR4626967
do
macs2 callpeak -t bam_macr/"$bam"_filtered.bam -c bam_macr/ERR4626964_filtered.bam -f BAM -g 2.62e+9 -n "$bam" -q 0.05 --outdir peaks_macr
macs2 callpeak -t bam_macr/"$bam"_filtered.bam -c bam_macr/ERR4626964_filtered.bam -f BAM -g 2.62e+9 -n "$bam" -q 0.05 --broad --outdir peaks_macr
done

# Get consensus peaks spleen


# Get consensus peaks macrophages
bedtools intersect -a peaks_macr/ERR4626957_peaks.narrowPeak -b peaks_macr/ERR4626963_peaks.narrowPeak > peaks_macr/CTCF_macr.bed
bedtools intersect -a peaks_macr/ERR4626959_peaks.narrowPeak -b peaks_macr/ERR4626965_peaks.narrowPeak > peaks_macr/H3K27ac_macr.bed
bedtools intersect -a peaks_macr/ERR4626960_peaks.broadPeak -b peaks_macr/ERR4626966_peaks.broadPeak > peaks_macr/H3K27me3_macr.bed
bedtools intersect -a peaks_macr/ERR4626961_peaks.broadPeak -b peaks_macr/ERR4626967_peaks.broadPeak > peaks_macr/H3K4me1_macr.bed
bedtools intersect -a peaks_macr/ERR4626962_peaks.narrowPeak -b peaks_macr/ERR4626968_peaks.narrowPeak > peaks_macr/H3K4me3_macr.bed

