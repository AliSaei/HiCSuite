
## 1. Define global variables for the project


```bash
# Define the location for the project 
PROJECT="./example_project"

# define the location of the raw HiC reads
HIC_READS="path/to/hic/reads"

# set the path to the assembled contigs (refrence genome)
REF="path/to/assembled/contigs"
REF_NAME=$(basename $REF)
```

## 2. QC of raw Hi-C reads


```bash
# define output and log directories
OUT="${PROJECT}/001.fastqc"
LOG="${OUT}/logs"

# create directiores
mkdir -p $OUT
mkdir -p $LOG

# retrieve the path to the raw file names
FILES=`ls /workspace/hrpazs/petunia_genome/HiC/007.trim_fastp/*.fastq.gz`

# run FASTQC on the raw data
module load FastQC
module load MultiQC

command="fastqc $FILES -q -t 2 -o ${OUT} && multiqc $OUT -o $OUT"

# Submit the command to lsf for execution 
bsub -o ${LOG}/FQC.out -e ${LOG}/FQC.err -n 4 -J FASTQC-fp $command

```

## 3. Trimmimg the Hi-C reads


```bash
# define output and log directories
OUT="${PROJECT}/002.trim_fastp"
LOG="${OUT}/logs"

# create directiores
mkdir -p $OUT
mkdir -p $LOG

module load fastp/0.19.4

COMMAND="fastp \
--in1 $RAW/PG_PETUNIA_HiC_HJTG3DSXX_ATGAGC_L004_R1.fastq.gz \
--out1 $OUT/PG_PETUNIA_HiC_HJTG3DSXX_ATGAGC_L004_trimmed_R1.fastq.gz \
--in2 $RAW/PG_PETUNIA_HiC_HJTG3DSXX_ATGAGC_L004_R2.fastq.gz \
--out2 $OUT/PG_PETUNIA_HiC_HJTG3DSXX_ATGAGC_L004_trimmed_R2.fastq.gz \
--trim_tail1 25 \
--trim_tail2 25 \
--disable_trim_poly_g \
--disable_quality_filtering \
--disable_length_filtering \
--thread 16
"

# Submit the command to lsf for execution
bsub -o $LOG/fastp.out -e $LOG/fastp.err -n 16 -J fastp $COMMAND
```

## 4. Aligning Hi-C reads (Phase Genomic) to the assembled contig 
For further information please visit Phase Genomics documentations on aligning Hi-C reads [here](https://phasegenomics.github.io/2019/09/19/hic-alignment-and-qc.html).

### 4.1. Index the assembled contigs (refrence sequence)


```bash
# define output and log directories
OUT="$PROJECT/003.bwa/index"
LOG="$OUT/logs"


# create directiores
mkdir -p $OUT
mkdir -p $LOG

# load bwa
module load bwa/0.7.17

command="bwa index -p ${OUT}/$REF_NAME -a bwtsw $REF"

# Submit the command to LSF job scheduler for execution on the cluster:
bsub -o $LOG/index-bwa.out -e $LOG/index-bwa.err -n 20 -J bwa_index $COMMAND

# or if running localy on one machine just execute the command by running:
$command > $LOG/index_bwa.out 2> $LOG/index_bwa.err
```

### 4.2. Align the Hi-C reads to the indexed refrence sequence


```bash
# define input, output and log directories
IN="${PROJECT}/002.trim_fastp"
OUT="$PROJECT/003.bwa/mem"
LOG="$OUT/logs"

# set the indexed refrence directory
INDEX="$PROJECT/002.bwa/index"

# create directiores
mkdir -p $OUT
mkdir -p $LOG

# get basenames of Hi-C reads and remove trailing read number and file extension
NAME=`basename $(ls ${IN}/*.fastq.gz) | sort -u | uniq | sed 's/_R[1,2].fastq.gz//g'`

module load bwa/0.7.17
module load samtools/1.9

for name in $NAME
    do
        # define read 2 and read 2 of the Hi-C reads
        READ1="${name}_R1.fastq.gz"
        READ1="${name}_R2.fastq.gz"
        
        # run bwa to map the reads then samtools to filter out unmapped reads, unmmaped mate, not primary alignment and 
        # supplementary alignment
        
        command="bwa mem -5SP -t 50 $INDEX $IN/$read1 $IN/$read2 | \
        samtools view -S -h -b -F 2316 > $OUT/$name.bam"
          
        bsub -o $LOG/align_bwa.out -e $LOG/align_bwa.err -n 50 -J bwa_align $command
     done
```

## 5. Filtering the Hi-C reads (Optional)
For further information please visit LACHESIS User's Manual [here](https://github.com/shendurelab/LACHESIS#4-filtering-the-hi-c-reads).


```bash
IN="$PROJECT/003.bwa/mem"


# get basename of bam file
BAM=$(basename $IN | sed 's/.bam//')

for bam in $BAM
    do
        command="./scripts/PreprocessSAMs.pl $IN/ $REF MBOI"
        bsub -o $LOG/PreprocessSAMs.out -e $LOG/PreprocessSAMs.err -J filter $command
    done    
```

## 6. Generating Hi-C contact file


```bash
# generating contact file

command="samtools view $ASSEMBLY_RAGOO_MAP | \
/workspace/hrpazs/bilberry_genome/filter_reads.awk -v isize=0 q=0 | \
/workspace/hrpazs/software/dryhic/inst/src/reads_to_bins_whole.awk -v w=100 > \
$WORKDIR/mapped_ragoo_contact_map_100_qualFilt.txt"
bsub -J convert -o $WORKDIR/mapped_ragoo_convert_to_contactMap_100.out -e $WORKDIR/mapped_ragoo_convert_to_contactMap_100.err $command
```
