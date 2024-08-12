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
This pipeline utilizes snakemake for workflow management and provided is usage analysis

### Specify File Locations
optional if installed to PATH: gatk,bwa,samtools,bcftools
```
python config.py --ref hg38.fasta \
--ref-version hg38 \
--funcotator-db funcotator_dataSources.v1.7.20200521s \
--varscan VarScan.v2.3.9.jar \
--gatk gatk-4.3.0.0/ \
--bcfools bcftools/ \
--samtools samtools/ 
```
### Run Snakemake Command
Contains automatic log file generation and error handling
```
snakemake -s variant_pipeline.smk --config variant_config.
```
