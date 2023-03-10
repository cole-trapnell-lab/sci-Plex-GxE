#
# To follow this pipeline, copy and paste the commands onto the command line
# in a qlogin session with at least 10G of memory (qlogin -l mfree=10G).
#
# Please read the comments carefully.
#

#
# Change this to a directory that you own
# This directory will contain all the outputs of the pipeline
#

WORKING_DIR=


#
# Change this directory to the local copy of the cloned github repo
# https://github.com/cole-trapnell-lab/sci-Plex-GxE
#

PATH_TO_GIT_REPO = 

 
# Download GNU Datamash
# https://www.gnu.org/software/datamash/
DATAMASH_PATH=

# Point pipeline at scripts directory
SCRIPTS_DIR=$PATH_TO_GIT_REPO/process_from_raw/sci-RNA-seq-pipeline-scripts

# Hash and gRNA sample Sheets can be found associated with each experiment
HASH_BARCODES_FILE=$WORKING_DIR/hashSampleSheet.txt
GRNA_BARCODES_FILE=$WORKING_DIR/gRNA_sample_sheet.txt

#
# STAR_INDEX is the path to the STAR genome index you want to align reads with.
# 

# Used for the mapping of all  data
STAR_INDEX=STAR/human/GRCh38-primary-assembly

# GENE_MODEL_DIR is a directory with BED files for various genomic features
# that are used to assign aligned reads to genes.
GENE_MODEL_DIR=path_to/gene-models/human/


# RT_BARCODES_FILE is a two column tab-delimited file where
# column 1 = well id
# column 2 = RT barcode
#
# e.g.
# A01    ACGTACGTAC
# A02    CGTACGTACG
#
# The id for a cell will be P7-well_P5-well_RT-well, e.g. A01_B01_C01
#
# You can optionally just have the barcode in both columns 1 and 2
# instead of using well ids, e.g.
# ACGTACGTAC   ACGTACGTAC
# CGTACGTACG   CGTACGTACG
#

RT_BARCODES_FILE=$SCRIPTS_DIR/RT_384_3lvl
LIG_BARCODES_FILE=$SCRIPTS_DIR/ligation_384_3lvl

# If you have only one biological sample,
# just write the sample name to a file called combinatorial.indexing.key
#
# If you have multiple samples, you can use a combinatorial.indexing.key file
# to tell the pipeline to automatically annotate which sample a cell came from
# based on its barcodes.

SAMPLE_NAME="sciPlexGxE”

echo "$SAMPLE_NAME" >$WORKING_DIR/combinatorial.indexing.key

#
# In various steps of the pipeline,
# one qsub job will be submitted for each BATCH_SIZE PCR wells of the experiment.
#
# If you have >1 PCR plate, you may want to increase BATCH_SIZE
# so you don't take up too many resources.
# The most memory intensive script that runs in batches in this manner
# uses 10 GB per qsub, or 240 GB total for 96 PCR wells and BATCH_SIZE=4
#

BATCH_SIZE=50

#-------------------------------------------------------------------------------
# Extract FASTQ files from BCL files
#-------------------------------------------------------------------------------

cd $WORKING_DIR

mkdir fastq

nohup bcl2fastq -R /path_to_seq_output -o fastq --no-lane-splitting --minimum-trimmed-read-length 0 --mask-short-adapter-reads 0 --barcode-mismatches 1 --with-failed-reads --sample-sheet RNA_SampleSheet.csv &

#-------------------------------------------------------------------------------
# Put read 1 info (RT well, UMI) into read 2 read name
#-------------------------------------------------------------------------------

cd $WORKING_DIR

mkdir combined-fastq
mkdir file-lists-for-r1-info-munging
mkdir put-r1-info-in-r2-logs

ls fastq/ | grep _R1_ | grep -v Undetermined | split -l $BATCH_SIZE -d - file-lists-for-r1-info-munging/

ls file-lists-for-r1-info-munging | while read BATCH; do

    qsub $SCRIPTS_DIR/put-read1-info-in-read2_3Level.sh \
        $WORKING_DIR/fastq                                      \
        $WORKING_DIR/file-lists-for-r1-info-munging/$BATCH      \
        $SCRIPTS_DIR/                                    \
        $RT_BARCODES_FILE                                       \
        $LIG_BARCODES_FILE                                      \
        $WORKING_DIR/combinatorial.indexing.key                 \
        $WORKING_DIR/combined-fastq                             \
        $WORKING_DIR/put-r1-info-in-r2-logs
done

#-------------------------------------------------------------------------------
# Parse Hash Barcodes
#-------------------------------------------------------------------------------

cd $WORKING_DIR

mkdir file-lists-for-trimming
mkdir hashed-fastq
mkdir hashed-logs

ls combined-fastq/ | split -l $BATCH_SIZE -d - file-lists-for-trimming/

ls file-lists-for-r1-info-munging | while read BATCH; do
    qsub $SCRIPTS_DIR/parse_hash.sh \
        $WORKING_DIR/combined-fastq                             \
        $WORKING_DIR/file-lists-for-trimming/$BATCH             \
        $SCRIPTS_DIR/                                    \
        $HASH_BARCODES_FILE                                     \
        $WORKING_DIR/combinatorial.indexing.key                 \
        $WORKING_DIR/hashed-fastq                               \
        $WORKING_DIR/hashed-logs
done

#-------------------------------------------------------------------------------
# Parse Hash Barcodes
#-------------------------------------------------------------------------------

cd $WORKING_DIR

mkdir $WORKING_DIR/hashRDS

zcat hashed-fastq/*.gz \
    | $DATAMASH_PATH -g 2,3,4 count 3 \
    | $DATAMASH_PATH -g 1 sum 4 \
    | awk -v S=$SAMPLE_NAME '{OFS="\t";} {print S, $0}' \
    > $WORKING_DIR/hashReads.per.cell

zcat hashed-fastq/*.gz \
    | uniq \
    | $DATAMASH_PATH -g 2,3,4 count 3 \
    | $DATAMASH_PATH -g 1 sum 4 \
    | awk -v S=$SAMPLE_NAME '{OFS="\t";} {print S, $0}' \
    > $WORKING_DIR/hashUMIs.per.cell

Rscript $SCRIPTS_DIR/knee-plot.R     \
    $WORKING_DIR/hashUMIs.per.cell          \
    $WORKING_DIR/hashed-logs

zcat hashed-fastq/*.gz \
    | uniq \
    | $DATAMASH_PATH -g 1,2,4,5 count 3  \
    > $WORKING_DIR/hashRDS/hashTable.out 

paste $WORKING_DIR/hashUMIs.per.cell  $WORKING_DIR/hashReads.per.cell \
    | cut -f 1,2,6,3 \
    | awk 'BEGIN {OFS="\t";} {dup = 1-($3/$4); print $1,$2,$3,$4,dup;}' \
    > $WORKING_DIR/hashDupRate.txt

Rscript $SCRIPTS_DIR/hashFilter.R            \
    $WORKING_DIR/hashRDS/hashTable.out              \
    $WORKING_DIR/hashRDS/                           

#-------------------------------------------------------------------------------
# Parse sgRNA Barcodes
#-------------------------------------------------------------------------------

cd $WORKING_DIR

mkdir file-lists-for-trimming
mkdir gRNA-fastq
mkdir gRNA-logs

ls combined-fastq/ | split -l $BATCH_SIZE -d - file-lists-for-trimming/

ls file-lists-for-r1-info-munging | while read BATCH; do
    qsub $SCRIPTS_DIR/parse_gRNA.sh \
        $WORKING_DIR/combined-fastq                             \
        $WORKING_DIR/file-lists-for-trimming/$BATCH             \
        $SCRIPTS_DIR                               \
        $GRNA_BARCODES_FILE                                     \
        $WORKING_DIR/combinatorial.indexing.key                 \
        $WORKING_DIR/gRNA-fastq                               	\
        $WORKING_DIR/gRNA-logs
done

#-------------------------------------------------------------------------------
# Summarize gRNA Barcodes
#-------------------------------------------------------------------------------

cd $WORKING_DIR

mkdir $WORKING_DIR/gRNA_summary

zcat gRNA-fastq/*.gz \
    | $DATAMASH_PATH -g 2,3,4 count 3 \
    | $DATAMASH_PATH -g 1 sum 4 \
    | awk -v S=$SAMPLE_NAME '{OFS="\t";} {print S, $0}' \
    > $WORKING_DIR/gRNA_summary/gRNAs.per.cell

zcat gRNA-fastq/*.gz \
    | uniq \
    | $DATAMASH_PATH -g 2,3,4 count 3 \
    | $DATAMASH_PATH -g 1 sum 4 \
    | awk -v S=$SAMPLE_NAME '{OFS="\t";} {print S, $0}' \
    > $WORKING_DIR/gRNA_summary/gRNA_UMIs.per.cell

Rscript $SCRIPTS_DIR/gRNA_knee-plot.R     \
    $WORKING_DIR/gRNA_summary/gRNA_UMIs.per.cell          \
    $WORKING_DIR/gRNA_summary

zcat gRNA-fastq/*.gz \
    | uniq \
    | $DATAMASH_PATH -g 1,2,4,5 count 3  \
    > $WORKING_DIR/gRNA_summary/gRNA_Table.out 

paste $WORKING_DIR/gRNA_summary/gRNA_UMIs.per.cell  $WORKING_DIR/gRNA_summary/gRNAs.per.cell \
    | cut -f 1,2,6,3 \
    | awk 'BEGIN {OFS="\t";} {dup = 1-($3/$4); print $1,$2,$3,$4,dup;}' \
    > $WORKING_DIR/gRNA_summary/gRNA_dup_rate.per.cell

zcat gRNA-fastq/*.gz \
    | $DATAMASH_PATH -g 2,3,4 count 3 \
    > $WORKING_DIR/gRNA_summary/reads_per_UMI.per.cell

#-------------------------------------------------------------------------------
# Trim poly-A tails
#-------------------------------------------------------------------------------

cd $WORKING_DIR

mkdir file-lists-for-trimming
mkdir trimmed-fastq

ls combined-fastq/ | split -l $BATCH_SIZE -d - file-lists-for-trimming/

ls file-lists-for-trimming | while read BATCH; do
    qsub $SCRIPTS_DIR/run-trim-galore.sh     \
        $WORKING_DIR/combined-fastq                         \
        $WORKING_DIR/file-lists-for-trimming/$BATCH         \
        $WORKING_DIR/trimmed-fastq
done

#-------------------------------------------------------------------------------
# Align reads using STAR
#-------------------------------------------------------------------------------

cd $WORKING_DIR

mkdir aligned-reads

qsub $SCRIPTS_DIR/STAR-alignReads-moreMem.sh \
    $WORKING_DIR/trimmed-fastq                              \
    $STAR_INDEX                                             \
    $WORKING_DIR/aligned-reads

#-------------------------------------------------------------------------------
# Filter ambiguously-mapped reads and sort BAM files
# Also count rRNA reads
#-------------------------------------------------------------------------------

cd $WORKING_DIR

mkdir file-lists-for-samtools-sort
mkdir aligned-reads-filtered-sorted
mkdir rRNA-read-counts

ls aligned-reads/ | grep "[.]Aligned[.]out[.]bam$" | split -l $BATCH_SIZE -d - file-lists-for-samtools-sort/

ls file-lists-for-samtools-sort | while read BATCH; do
    qsub $SCRIPTS_DIR/samtools-filter-sort.sh    \
        $WORKING_DIR/aligned-reads                              \
        $WORKING_DIR/file-lists-for-samtools-sort/$BATCH        \
        $WORKING_DIR/aligned-reads-filtered-sorted
done

#
# You can run these jobs in parallel to the samtools-filter-sort.sh jobs
#

ls file-lists-for-samtools-sort | while read BATCH; do
    qsub $SCRIPTS_DIR/count-rRNA-reads.sh        \
        $WORKING_DIR/aligned-reads                              \
        $WORKING_DIR/file-lists-for-samtools-sort/$BATCH        \
        $GENE_MODEL_DIR/latest.rRNA.gene.regions.union.bed      \
        $WORKING_DIR/rRNA-read-counts
done 

#-------------------------------------------------------------------------------
# Split reads in BAM files into BED intervals --- Subbed in Andrew's code here
#-------------------------------------------------------------------------------

cd $WORKING_DIR

mkdir file-lists-for-rmdup
mkdir aligned-reads-rmdup-split-bed

ls aligned-reads-filtered-sorted/ | grep "[.]bam$" | split -l $BATCH_SIZE -d - file-lists-for-rmdup/

ls file-lists-for-rmdup | while read BATCH; do
    qsub $SCRIPTS_DIR/rmdup-and-make-split-bed_AH.sh    \
        $WORKING_DIR/aligned-reads-filtered-sorted                  \
        $WORKING_DIR/file-lists-for-rmdup/$BATCH                    \
        $SCRIPTS_DIR                                       \
        $WORKING_DIR/aligned-reads-rmdup-split-bed
done

#-------------------------------------------------------------------------------
# Assign reads to genes, using the BED files as input
#-------------------------------------------------------------------------------

cd $WORKING_DIR

mkdir file-lists-for-assign-reads-to-genes
mkdir unique-read-to-gene-assignments

ls aligned-reads-rmdup-split-bed/ | grep "[.]bed$" | split -l $BATCH_SIZE -d - file-lists-for-assign-reads-to-genes/

ls file-lists-for-assign-reads-to-genes | while read BATCH; do
    qsub $SCRIPTS_DIR/assign-reads-to-genes.sh       \
        $WORKING_DIR/aligned-reads-rmdup-split-bed                  \
        $WORKING_DIR/file-lists-for-assign-reads-to-genes/$BATCH    \
        $GENE_MODEL_DIR/latest.exons.bed                            \
        $GENE_MODEL_DIR/latest.genes.bed                            \
        $SCRIPTS_DIR/                                               \
        $WORKING_DIR/unique-read-to-gene-assignments
done

#-------------------------------------------------------------------------------
# Compute the duplication rate and proportion of reads that are from rRNA
#-------------------------------------------------------------------------------

cd $WORKING_DIR

mkdir UMI-counts-by-sample
mkdir file-lists-for-UMI-counting

ls aligned-reads-rmdup-split-bed/ | while read FILE; do
    PCR_WELL=`basename $FILE .bed`
    echo "$PCR_WELL"
done \

| split -l $BATCH_SIZE -d - file-lists-for-UMI-counting/

ls file-lists-for-UMI-counting | while read BATCH; do
    qsub $SCRIPTS_DIR/count-UMI-per-sample.sh    \
        $WORKING_DIR/aligned-reads-filtered-sorted              \
        $WORKING_DIR/aligned-reads-rmdup-split-bed              \
        $WORKING_DIR/file-lists-for-UMI-counting/$BATCH         \
        $WORKING_DIR/UMI-counts-by-sample
done

#-------------------------------------------------------------------------------

cat UMI-counts-by-sample/*.UMI.count | sort -k1,1 \
| $DATAMASH_PATH -g 1 sum 2 \
>total.UMI.count.by.sample

cat UMI-counts-by-sample/*.read.count | sort -k1,1 \
| $DATAMASH_PATH -g 1 sum 2 \
>total.read.count.by.sample

mkdir final-output

cat rRNA-read-counts/* | sort -k1,1 | $DATAMASH_PATH -g 1 sum 2 sum 3 \
| join - total.UMI.count.by.sample \
| join - total.read.count.by.sample \
| awk 'BEGIN {
    printf "%-18s    %13s    %8s    %10s    %8s\n",
        "sample", "n.reads", "pct.rRNA", "n.UMI", "dup.rate";
} {
    printf "%-18s    %13d    %7.1f%%    %10d    %7.1f%%\n",
        $1, $3, 100 * $2/$3, $4, 100 * (1 - $4/$5);
}' \
>final-output/rRNA.and.dup.rate.stats

cat final-output/rRNA.and.dup.rate.stats

#-------------------------------------------------------------------------------
# Make knee plots
#-------------------------------------------------------------------------------

cd $WORKING_DIR

module load modules modules-init modules-gs
module load pypy/3.5.6.0

source /path_to/bin/activate

START_TIME=$SECONDS

cat unique-read-to-gene-assignments/* \
	| pypy  $SCRIPTS_DIR/count_UMIs.py --cell - \
	> UMIs.per.cell.barcode_2
ELAPSED_TIME=$(($SECONDS - $START_TIME))

echo "Processed $FILE in $ELAPSED_TIME seconds"

mkdir final-output/knee-plots

Rscript $SCRIPTS_DIR/knee-plot.R            \
    UMIs.per.cell.barcode_2                   \
    $WORKING_DIR/final-output/knee-plots

cd final-output && tar czf ../knee-plots.tar.gz knee-plots/ && cd ..

# Inspect knee-plots.tar.gz
# Decide on a cutoff for the number of UMIs per cell barcode needed
# for a cell barcode to be considered a cell and included in the CellDataSet.

# Once you've decided on a UMI-per-cell cutoff,
# re-run knee-plot.R as follows to make versions of the plots
# that show the cutoff as a horizontal red line on the plot.

cd $WORKING_DIR

UMI_PER_CELL_CUTOFF=500

Rscript $SCRIPTS_DIR/knee-plot.R            \
    UMIs.per.cell.barcode_2                   \
    $WORKING_DIR/final-output/knee-plots    \
    $UMI_PER_CELL_CUTOFF

cd final-output && tar czf ../knee-plots.tar.gz knee-plots/ && cd ..


#-------------------------------------------------------------------------------
# Make the final UMI count matrix
#-------------------------------------------------------------------------------

cd $WORKING_DIR

cp $GENE_MODEL_DIR/latest.gene.annotations final-output/gene.annotations

mkdir file-lists-for-UMI-count-rollup
mkdir UMI-count-rollup

ls unique-read-to-gene-assignments/ | split -l $BATCH_SIZE -d - file-lists-for-UMI-count-rollup/

ls file-lists-for-UMI-count-rollup | while read BATCH; do
    qsub $SCRIPTS_DIR/UMI-count-rollup.sh        \
        $WORKING_DIR/unique-read-to-gene-assignments            \
        $WORKING_DIR/file-lists-for-UMI-count-rollup/$BATCH     \
        $WORKING_DIR/UMI-count-rollup
done

#-------------------------------------------------------------------------------

cd $WORKING_DIR

cat UMI-count-rollup/* | gzip >prelim.UMI.count.rollup.gz

# Make samples.to.exclude file.
# Each line lists a sample name to exclude.
# You can leave it empty.

touch samples.to.exclude

# This uses the UMI_PER_CELL_CUTOFF variable
# defined in the knee plot section.

echo "UMI_PER_CELL_CUTOFF = $UMI_PER_CELL_CUTOFF"

gunzip <prelim.UMI.count.rollup.gz \
| $DATAMASH_PATH -g 1 sum 3 \
| tr '|' '\t' \
| awk -v CUTOFF=$UMI_PER_CELL_CUTOFF '
    ARGIND == 1 {
        exclude[$1] = 1;
    } $3 >= CUTOFF && !($1 in exclude) {
        print $2 "\t" $1; 
    }' samples.to.exclude - \
| sort -k1,1 -S 4G \
>final-output/cell.annotations

gunzip <prelim.UMI.count.rollup.gz \
| tr '|' '\t' \
| awk '{
    if (ARGIND == 1) {
        gene_idx[$1] = FNR;
    } else if (ARGIND == 2) {
        cell_idx[$1] = FNR;
    } else if ($2 in cell_idx) {
        printf "%d\t%d\t%d\n",
            gene_idx[$3], cell_idx[$2], $4;
    }
}' final-output/gene.annotations final-output/cell.annotations - \
>final-output/UMI.count.matrix

Rscript ~/SCRIPTS_DIR/makeCDS.R \
        $WORKING_DIR/final-output/UMI.count.matrix \
        $WORKING_DIR/final-output/gene.annotations \
        $WORKING_DIR/final-output/cell.annotations \
        final-output

#-------------------------------------------------------------------------------
# Using the knee plots, define cells that are background cells and 
# those cells that are intermediate cells
#-------------------------------------------------------------------------------

cd $WORKING_DIR

UMI_PER_CELL_CUTOFF=500

 awk '$3 >= CUTOFF { print $2 }' CUTOFF=$UMI_PER_CELL_CUTOFF UMIs.per.cell.barcode_2 \
    > $WORKING_DIR/hashRDS/real.cells.csv

awk '$3 <= CUTOFF {print $2}' CUTOFF=$UMI_PER_CELL_CUTOFF UMIs.per.cell.barcode_2 \
    > $WORKING_DIR/hashRDS/background.cells.csv

Rscript $SCRIPTS_DIR/chiSq_3lvl.R \
        $WORKING_DIR/hashRDS/background.cells.csv \
        $WORKING_DIR/hashRDS/real.cells.csv \
        $WORKING_DIR/hashed-fastq/ \
        $HASH_BARCODES_FILE \
        $WORKING_DIR/final-output/cds.RDS \
        final-output
#-------------------------------------------------------------------------------

# This is the end of the pipeline.
# The data is now ready to load into R as a Monocle CellDataSet.

# Here is an R function to do that.

load.cds = function(mat.path, gene.annotation.path, cell.annotation.path) {
    df = read.table(
        mat.path,
        col.names = c("gene.idx", "cell.idx", "count"),
        colClasses = c("integer", "integer", "integer"))    

    gene.annotations = read.table(
        gene.annotation.path,
        col.names = c("id", "gene_short_name"),
        colClasses = c("character", "character"))

    cell.annotations = read.table(
        cell.annotation.path,
        col.names = c("cell", "sample"),
        colClasses = c("character", "factor"))

    rownames(gene.annotations) = gene.annotations$id
    rownames(cell.annotations) = cell.annotations$cell

    # add a dummy cell to ensure that all genes are included in the matrix
    # even if a gene isn't expressed in any cell

    df = rbind(df, data.frame(
        gene.idx = c(1, nrow(gene.annotations)),
        cell.idx = rep(nrow(cell.annotations)+1, 2),
        count = c(1, 1)))

    mat = sparseMatrix(i = df$gene.idx, j = df$cell.idx, x = df$count)
    mat = mat[, 1:(ncol(mat)-1)]

    rownames(mat) = gene.annotations$id
    colnames(mat) = cell.annotations$cell

    pd = new("AnnotatedDataFrame", data = cell.annotations)
    fd = new("AnnotatedDataFrame", data = gene.annotations)

    cds = new_cell_data_set(mat, cell_metadata = pd, gene_metadata = fd)

    colData(cds)$n.umi = apply(exprs(cds), 2, sum)

    return(cds)
}

# Example usage:

cds = load.cds(
    "final-output/UMI.count.matrix",
    "final-output/gene.annotations",
    "final-output/cell.annotations")
```