<h1>Sequencing Data Alignment</h1>


<h2>Description</h2>
This project consists of a Bash script that walks the user through retrieving a reference genome, downloading sequencing data, aligning reads, and processing output. This script is useful for variant calling and downstream analysis of sequencing data.
<br />


<h2>Languages Used</h2>

- <b>Bash</b> 

<h2>Tools Used </h2>

- <b>Burrows-Wheeler Aligner</b> 
- <b>NCBI Entrez Direct</b>
- <b>SRA Toolkit</b>
- <b>SAMtools</b>

<h2>Script:</h2>


# Stop on error
set -uex

# Define directory name
DIR="refs"

# Create a new directory
mkdir -p "$DIR"

# Define the accession number for the reference genome
ACC="AF086833"

# Define path to reference genome file
REF="$DIR/$ACC.fa"

# Get the genome in GenBank format and extract the source only
bio fetch "$ACC" --format fasta > "$REF"

# Build the index using bwa 
bwa index "$REF"

# Obtain SRA data and save into new file
QUERY="PRJNA257197"
RUNINFO="$DIR/runinfo.csv"
esearch -db sra -query "$QUERY" | efetch -format runinfo > "$RUNINFO"

# Pick a run from the file
RUN="SRR1972739"

# Define the number of reads to subset
NUM_READS="10000"

# Subset the data to 10K reads
fastq-dump -X "$NUM_READS" --split-files "$RUN" -O "$DIR"

# Perform alignment in paired-end mode
OUTPUT_SAM="$DIR/output.sam"
R1="$DIR/${RUN}_1.fastq"
R2="$DIR/${RUN}_2.fastq"
bwa mem "$REF" "$R1" "$R2" > "$OUTPUT_SAM"

# Convert SAM to sorted BAM
OUTPUT_BAM="$DIR/output.bam"
samtools sort -o "$OUTPUT_BAM" "$OUTPUT_SAM"

# Index the BAM file
samtools index "$OUTPUT_BAM"

<!--
 ```diff
- text in red
+ text in green
! text in orange
# text in gray
@@ text in purple (and bold)@@
```
--!>
