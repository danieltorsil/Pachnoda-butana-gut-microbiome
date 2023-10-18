#!/bin/bash

for part in P1 P2 P3 P4
	do
		# Run FastQC on input files ${part}_R1.fq.gz and ${part}_R2.fq.gz with 16 threads,
		# and save the results in the 'fastQCin' directory.
		fastqc -t 16 -o fastQCin ${part}_R1.fq.gz ${part}_R2.fq.gz

		# Use cutadapt to tr${part} adapter sequences from the input files, and save the
		# tr${part}med reads in ${part}_R1_cutadapt.fq.gz and ${part}_R2_cutadapt.fq.gz.
		cutadapt -j 16 -g AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
		-a GATCGGAAGAGCACACGTCTGAACTCCAGTCACGGATGACTATCTCGTATGCCGTCTTCTGCTTG \
		-G AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
		-A GATCGGAAGAGCACACGTCTGAACTCCAGTCACGGATGACTATCTCGTATGCCGTCTTCTGCTTG \
		-o ${part}_R1_cutadapt.fq.gz -p ${part}_R2_cutadapt.fq.gz ${part}_R1_cutadapt.fq.gz ${part}_R2_cutadapt.fq.gz

		# Use bbduk.sh to further tr${part} and filter the reads, and save them in
		# ${part}_R1_bbduk.fq.gz and ${part}_R2_bbduk.fq.gz.
		bbduk.sh threads=16 in=${part}_R1_cutadapt.fq.gz in2=${part}_R2_cutadapt.fq.gz \
		out=${part}_R1_bbduk.fq.gz out2=${part}_R2_bbduk.fq.gz qtr${part}=lr tr${part}q=20 maq=20 minlen=75

		# Align the tr${part}med reads to a reference genome using Bowtie2 with 16 threads,
		# and save the aligned and unaligned reads in ${part}_mapped_and_unmapped.sam.
		bowtie2 -p 16 -x /media/darwin/Externo/DanielTS/BBDD/GRCh38_noalt_as/GRCh38_noalt_as \
		-1 ${part}_R1_bbduk.fq.gz -2 ${part}_R2_bbduk.fq.gz --un-conc-gz filtered_ > ${part}_mapped_and_unmapped.sam

		# Rename the filtered output files from Bowtie2.
		mv filtered_${part}_R1_bbduk.fq.gz.1 ${part}_filtered_R1.fq.gz
		mv filtered_${part}_R2_bbduk.fq.gz.2 ${part}_filtered_R2.fq.gz

		# Remove the intermediate SAM file.
		rm -r ${part}_mapped_and_unmapped.sam

		# Run FastQC on the filtered reads ${part}_filtered_R1.fq.gz and ${part}_filtered_R2.fq.gz
		# with 16 threads, and save the results in the 'fastQCour' directory.
		fastqc -t 16 -o fastQCout ${part}_filtered_R1.fq.gz ${part}_filtered_R2.fq.gz

		# Assemble the filtered reads using Megahit with 16 threads and save the output
		# in the 'megahit/' directory with the prefix 'a_'.
		megahit -t 16 -1 ${part}_filtered_R1.fq.gz -2 ${part}_filtered_R2.fq.gz -o megahit/ --out-prefix a_

		# Run Quast to assess the quality of the assembly and save the results in the 'QUAST/' directory
		# using 16 threads and the contigs from the Megahit assembly.
		quast.py -o QUAST/ -t 16 megahit/${part}.contigs.fa

		# Create directories for Bowtie2 output.
		mkdir bw2
		mkdir bw2/bw2_index

		# Build Bowtie2 index for the Megahit assembly.
		bowtie2-build megahit/${part}.contigs.fa --threads 16 bw2/bw2_index/${part}_assembly_index

		# Align the filtered reads to the Bowtie2 index and save the output as a SAM file.
		bowtie2 --threads 16 -x  bw2/bw2_index/${part}_assembly_index-1 ${part}_filtered_R1.fq.gz -2 ${part}_filtered_R2.fq.gz -S bw2/bw2_index/${part}_assembly_coverage.sam

		# Convert SAM to BAM format.
		samtools view -b bw2/bw2_index/${part}_assembly_coverage.sam --threads 16 > bw2/bw2_index/${part}_assembly_coverage.bam

		# Sort the BAM file.
		samtools sort --threads 16 bw2/bw2_index/${part}_assembly_coverage.bam > bw2/bw2_index/${part}_assembly_coverage_sorted.bam

		# Generate depth information for the assembly.
		jgi_summarize_bam_contig_depths --outputDepth bw2/depth.txt bw2/bw2_index/${part}_assembly_coverage_sorted.bam

		# Remove intermediate SAM and BAM files.
		rm bw2/bw2_index/${part}_assembly_coverage.sam
		rm bw2/bw2_index/${part}_assembly_coverage.bam
		rm bw2/bw2_index/${part}_assembly_coverage_sorted.bam
		
		# Run the Prodigal tool to predict protein-coding genes and save the output as ${part}_protein_sequences.faa
		prodigal -i megahit/${part}.contigs.fa -a ${part}_protein_sequences.faa
		
		# KofamScan of the assembly
		kofam_scan-master/exec_annotation -c ~/kofam_scan-master/a_config.yml -o ${part}_kofamscan_results.txt ${part}_protein_sequences.faa
		
		# Create directories for Metabat2 and MaxBin2 output.
		mkdir metabat2_depth
		mkdir maxbin2_depth

		# Run Metabat2 to bin contigs.
		metabat2 -i megahit/${part}.contigs.fa -o metabat2_depth/bin -t 16 -a bw2/depth.txt

		# Prepare depth information for MaxBin2.
		cut -f1,3 bw2/depth.txt | tail -n+2 > bw2/depth_maxbin.txt

		# Run MaxBin2 to bin contigs.
		run_MaxBin.pl -contig megahit/${part}.contigs.fa -out maxbin2_depth/bin -abund bw2/depth_maxbin.txt

		# Create a directory for DAS_Tool output.
		mkdir dastool

		# Use DAS_Tool to integrate Metabat2 and MaxBin2 binning results.
		Fasta_to_Scaffolds2Bin.sh -i metabat2_depth -e fa > dastool/metabat2.scaffolds2bin.tsv
		Fasta_to_Scaffolds2Bin.sh -i maxbin2_depth -e fasta > dastool/maxbin2.scaffolds2bin.tsv
		DAS_Tool -i dastool/metabat2.scaffolds2bin.tsv,dastool/maxbin2.scaffolds2bin.tsv -l metabat,maxbin -c megahit/${part}.contigs.fa -o dastool/ --search_engine diamond --threads 16 --write_bins

		# Create directories for CheckM output.
		mkdir checkm_metabat2_depth
		mkdir checkm_maxbin2_depth
		mkdir dastool/checkm

		# Run CheckM to assess bin quality.
		checkm lineage_wf -t 16 -f checkm_metabat2_depth/quality_output.txt -x fa metabat2_depth
		checkm lineage_wf -t 16 -f checkm_maxbin2_depth/quality_output.txt -x fasta maxbin2_depth
		checkm lineage_wf -t 16 -f dastool/quality_dastool_output.txt -x fa dastool/*"DASTool_bins" dastool/checkm


		# Calculate genome coverage using 'coverm' tool.
		coverm genome --coupled ${part}_filtered_R1.fq.gz ${part}_filtered_R2.fq.gz --genome-fasta-files dastool/*fa -o output.tsv

		# Classify genomes with 'gtdbtk'.
		gtdbtk classify_wf --genome_dir dastool/ --out_dir MAGS_classification/ -x fa --cpus 16


		# Annotate genomes with 'prokka'.
		for MAG in dastool/*fa
			do
    				# Run the Prodigal tool to predict protein-coding genes and save the output as ${MAG}_protein_sequences.faa
				prodigal -i ${MAG} -a ${MAG}_protein_sequences.faa
				kofam_scan-master/exec_annotation -c ~/kofam_scan-master/a_config.yml -o ${MAG}_kofamscan_results.txt ${MAG}_protein_sequences.faa
		    	done
	echo

