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

set -uex

DIR="refs"

mkdir -p "$DIR"

ACC="AF086833"

REF="$DIR/$ACC.fa"

bio fetch "$ACC" --format fasta > "$REF"

bwa index "$REF"

QUERY="PRJNA257197"
RUNINFO="$DIR/runinfo.csv"
esearch -db sra -query "$QUERY" | efetch -format runinfo > "$RUNINFO"

RUN="SRR1972739"

NUM_READS="10000"

fastq-dump -X "$NUM_READS" --split-files "$RUN" -O "$DIR"

OUTPUT_SAM="$DIR/output.sam"
R1="$DIR/${RUN}_1.fastq"
R2="$DIR/${RUN}_2.fastq"
bwa mem "$REF" "$R1" "$R2" > "$OUTPUT_SAM"

OUTPUT_BAM="$DIR/output.bam"
samtools sort -o "$OUTPUT_BAM" "$OUTPUT_SAM"

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
