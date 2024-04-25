#!/bin/bash

#load software
export PATH=/programs/STAR-2.7.9a/bin/Linux_x86_64_static:$PATH
export PATH=/programs/salmon-1.4.0/bin:$PATH

# make all of the output directories
# The -p option means mkdir will create the whole path if it does not exist and refrain from complaining if it does exist
mkdir -p ./results/fastqc/
mkdir -p ./results/STAR/
mkdir -p ./results/qualimap/
mkdir -p ./results/salmon/

#specify paths to raw fastq files, genome + index, transcriptome + index, and annotation file
export raw_fq="$1"
export genome="$2"
export transcriptome="$3"
export gtf="$4"

export cores=16

for files in ${raw_fq} do
	samplename=`basename ${files}`
	#Run FastQC and move output to the appropriate folder
	echo "Starting fastQC for $samplename"
	fastqc -o ./results/fastqc/ ${samplename}

	#Run STAR
	echo "Starting STAR for $samplename"
	STAR --runThreadN $cores --genomeDir $genome --readFilesIn ./raw_data/fastq/${samplename}/${samplename}_1.*.gz ./raw_data/fastq/${samplename}/${samplename}_2.*.gz --readFilesCommand zcat --outFileNamePrefix ${samplename} --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes Standard

	#Run Qualimap
	echo "Starting Qualimap for $samplename"
	/programs/qualimap_v2.2.1/qualimap rnaseq \
	-outdir ./results/qualimap/ \
	-a proportional \
	-bam ./results/STAR/${samplename}.bam \
	-gtf $gtf \
	--java-mem-size=8G

	#Run salmon
	echo "Starting Salmon run for $samplename"
	salmon quant -i $transcriptome \
	-p $cores \
	-l A \
	-1 ./raw_data/fastq/${samplename}/${samplename}_1.*.gz \
	-2 ./raw_data/fastq/${samplename}/${samplename}_2.*.gz \
	-o ./results/salmon/ \
	--seqBias \
	--useVBOpt
done 
