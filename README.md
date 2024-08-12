# Ivy's SNP Pipeline
Included are the following tools and uses <br>
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
   
## Snakemake
This pipeline utilizes snakemake for workflow management and provided is usage analysis <br>
<p> This pipeline takes all files in a fastq directory, maps the files, call variants with two methods and lastly will annotate the overlapping calls. <br>
Bam files can be provided instead of fastq to skip mapping step</p>
### Specify File Locations
optional if installed to PATH: gatk,bwa,samtools,bcftools

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
