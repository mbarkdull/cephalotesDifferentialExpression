#!/bin/bash

#load software
export PATH=/programs/STAR-2.7.9a/bin/Linux_x86_64_static:$PATH
export PATH=/programs/salmon-1.4.0/bin:$PATH

# make all of the output directories
# The -p option means mkdir will create the whole path if it does not exist and refrain from complaining if it does exist
mkdir -p ./02_alignments/fastqc/
mkdir -p ./02_alignments/STAR/
mkdir -p ./02_alignments/qualimap/
mkdir -p ./02_alignments/salmon/

#specify paths to raw fastq files, genome + index, transcriptome + index, and annotation file, along with the numbers of cores to use:
export raw_reads_directory="$1"
export genome="$2"
export transcriptome="$3"
export gtf="$4"
export cores=$5

# Index the transcriptome for salmon:
salmon index -t $transcriptome -i ./indexed_transcriptome

for files in ${raw_reads_directory} do
	samplename=`basename ${files}`
	#Run FastQC and move output to the appropriate folder
	echo "Starting fastQC for $samplename"
#	fastqc -o ./02_alignments/fastqc/ ${samplename}

	#Run STAR
#	echo "Starting STAR for $samplename"
#	STAR --runThreadN $cores --genomeDir $genome --readFilesIn ./raw_data/fastq/${samplename}/${samplename}_1.*.gz ./raw_data/fastq/${samplename}/${samplename}_2.*.gz --readFilesCommand zcat --outFileNamePrefix ${samplename} --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes Standard

	#Run Qualimap
#	echo "Starting Qualimap for $samplename"
#	/programs/qualimap_v2.2.1/qualimap rnaseq \
#	-outdir ./02_alignments/qualimap/ \
#	-a proportional \
#	-bam ./02_alignments/STAR/${samplename}.bam \
#	-gtf $gtf \
#	--java-mem-size=8G

	#Run salmon
#	echo "Starting Salmon run for $samplename"
#	salmon quant -i $transcriptome \
#	-p $cores \
#	-l A \
#	-1 ./raw_data/fastq/${samplename}/${samplename}_1.*.gz \
#	-2 ./raw_data/fastq/${samplename}/${samplename}_2.*.gz \
#	-o ./02_alignments/salmon/ \
#	--seqBias \
#	--useVBOpt
done
