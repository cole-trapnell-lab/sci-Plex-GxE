#$ -S /bin/bash
#$ -l mfree=2G
#$ -l h_rt=6:0:0:0
#$ -pe serial 20

module load STAR/2.5.2b

INPUT=$1
INDEX=$2
OUTPUT=$3

STAR --genomeDir $INDEX --genomeLoad Remove

ls $INPUT | grep "[.]fq[.]gz$" | while read FILE; do
    SAMPLE=`basename "$FILE" _trimmed.fq.gz`

    STAR \
        --runThreadN 20 \
        --genomeDir $INDEX \
        --genomeLoad LoadAndKeep \
        --readFilesIn $INPUT/$FILE \
        --readFilesCommand zcat \
        --outFileNamePrefix $OUTPUT/$SAMPLE. \
        --outSAMtype BAM Unsorted \
        --outSAMstrandField intronMotif
done

STAR --genomeDir $INDEX --genomeLoad Remove

