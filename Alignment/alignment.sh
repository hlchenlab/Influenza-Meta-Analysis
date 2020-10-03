#############################
#Title:STAR Alignment
#Written by: Conor Cremin


mkdir out
fastqDIR=$PWD  #Specify directory were fastq files are located.

for file in $(ls $fastqDIR/*.fastq.gz | cut -c35- | rev | cut -c12- | rev | uniq); ## Change to get sample ids i.e. /path/to/fastq/dir/{sample_id}.fastq.gz
do
	/software/STAR/2.6.0c/bin/STAR --readFilesCommand gunzip -c \
				--outFileNamePrefix $out/${file} \
				--runThreadN 2 --genomeDir $PWD/star_idx \
				--outSAMtype NONE --quantMode GeneCounts \
				--readFilesIn ${file}.1.fastq.gz ${file}.2.fastq.gz
done
