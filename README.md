# Ivy's SNP Pipeline
<p> This pipeline takes all files in a fastq directory, maps the files, call variants with two methods and lastly will annotate the overlapping calls. Bam files can be provided instead of fastq to skip mapping step</p> <br>

## Included are the following tools and uses <br>
1. Alignment with BWA
2. Sorting by Coordinate
3. Variant Calling with Mutect2 and Varscan2
4. Joint Method Variant Calling (Overlap between methods)
5. Functional Annotation with Funcotater

## Dependencies
This workflow requires the folling before usage
1. samtools
2. bwa
3. gatk
4. java
5. bcftools
6. references + databases (reference assembly, funcotater database, gnomAD germline AF only file)
   
## Configure Variables
This pipeline utilizes snakemake for workflow management and provided is usage analysis <br>

### Specify File Locations
Run This command to configure snakemake run correctly. It provides absolute path of all files and dependencies needed. Adjust to your own needs <br> <br>
The following are optional if installed to system
1. bwa
2. samtools
3. bcftools
<br>
You Must Specify the following
1. fastq_dir
2. gatk executable
3. varscan jar file
4. reference fasta file
5. funcotator database location
<br>
You can configure the number of threads available with the --threads parameter
```
python config.py \
   --fastq_dir /path/to/fastqs \
    --bwa /path/to/bwa \
    --samtools /path/to/samtools \
    --bcftools /path/to/bcftools \
    --gatk /path/to/gatk \
    --varscan /path/to/varscan.jar \
    --reference_fasta /path/to/Homo_sapiens_assembly38.fasta \
    --funcotator_db /path/to/funcotator_dataSources.v1.7.20200521s \
    --output config.yaml \
    --threads num_threads

```
## Run Snakemake Command
Contains automatic log file generation and error handling
```
snakemake -s variant_pipeline.smk
```
