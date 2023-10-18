#!/bin/bash

# Activate the 'qiime2-2020.8' Conda environment
conda activate qiime2-2020.8

# Generate a Manifest file from fastq files
qiime2input.py -i /media/darwin/Proyectos/2022/PACHNODA_01_2022/16S_suppliers/QIIME2/fastq/ > Manifest

# Import the fastq data using the Manifest file
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path Manifest --output-path ./demux-paired-end.qza --input-format PairedEndFastqManifestPhred33

# Generate a visualization summarizing the demultiplexed data
qiime demux summarize --i-data ./demux-paired-end.qza --o-visualization ./demux-paired-end.qzv

# Perform denoising on paired-end sequences using DADA2
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs ./demux-paired-end.qza \
  --verbose \
  --p-trim-left-f 17 \
  --p-trim-left-r 22 \
  --p-trunc-len-f 294 \
  --p-trunc-len-r 278 \
  --p-n-threads 16 \
  --o-table ./table.qza \
  --o-representative-sequences ./rep-seqs.qza \
  --o-denoising-stats ./denoising-stats.qza

# Generate a tabulated summary of denoising statistics
qiime metadata tabulate --m-input-file ./denoising-stats.qza --o-visualization ./stats-dada2.qzv

# Classify sequences against a reference database (SILVA classifier)
qiime feature-classifier classify-sklearn \
  --i-classifier /media/darwin/Externo/DanielTS/BBDD/silva-138-99-nb-classifier.qza \
  --i-reads ./rep-seqs.qza \
  --verbose \
  --p-n-jobs 16 \
  --o-classification ./taxonomy.qza

# Export the feature table
qiime tools export --input-path ./table.qza --output-path ./table
ls table

# Export the taxonomy data
qiime tools export --input-path ./taxonomy.qza --output-path ./taxonomy
ls taxonomy

# Modify the taxonomy file to replace column headers
# Replace 'Feature ID' with '#OTUID'
sed -i -e 's/Feature ID/#OTUID/g' ./taxonomy/taxonomy.tsv

# Replace 'Taxon' with 'taxonomy'
sed -i -e 's/Taxon/taxonomy/g' ./taxonomy/taxonomy.tsv

# Replace 'Confidence' with 'confidence'
sed -i -e 's/Confidence/confidence/g' ./taxonomy/taxonomy.tsv

# Display the modified taxonomy file
head taxonomy/taxonomy.tsv

# Add taxonomy metadata to the feature table
biom add-metadata -i ./table/feature-table.biom -o ./table-with-taxonomy.biom --observation-metadata-fp ./taxonomy/taxonomy.tsv --sc-separated taxonomy

# Convert the feature table to JSON format
biom convert -i ./table-with-taxonomy.biom -o ./table-with-taxonomy-json2.biom --table-type="OTU table" --to-json
