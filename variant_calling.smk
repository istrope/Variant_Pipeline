"""
Author: Ivy Strope
Last Edited: 8/12/24
Contact: ivystrope@gmail.com
"""
#######################################
#     SPECIFY WILDCARD VARIABLES
#######################################
configfile:'config.yaml'
threads=config['threads']
samtools=config['samtools']
bwa=config['bwa']
bcftools=config['bcftools']
gatk=config['gatk']
varscan=config['varscan']
reference=config['reference_fasta']
funcotator_db=config['funcotator_db']

#####################################
#     SPECIFY FINAL OUTPUT FILE
#####################################
fastq_dir=config['fastq_dir']
FILES=glob_wildcards(f'{fastq_dir}/{{sample}}.fastq.gz')
SAMPLES=FILES.sample
rule all:
    input:
        expand("annotated/{sample}_funcotated.vcf",sample=SAMPLES)

#####################################
#      MAP READS WITH BWA
#####################################

rule map_reads:
    input:
        fastq=fastq_dir + "{sample}.fastq.gz",
        ref=reference
    output:
        sorted="{sample}_sorted.bam",
        mapped=temp('{sample}.sam'),
        bam=temp("{sample}.bam")
    params:
        bwa=bwa,
        samtools=samtools,
        threads=threads
    shell:
        """
        header=$(zcat {input.fastq} | head -n 1)
        id=$(echo $header | head -n 1 | cut -f 1-4 -d":" | sed 's/@//' | sed 's/:/_/g')
        sm=$(echo $header | head -n 1 | grep -Eo "[ATGCN]+$")
        
        {params.bwa} mem -t {params.threads} -R $(echo "@RG\tID:$id\tSM:$id"_"$sm\tLB:$id"_"$sm\tPL:ILLUMINA") {input.ref} {input.fastq} > {output.mapped}
        {params.samtools} view -Sb -@ {params.threads} -o {output.bam} {output.mapped}
        {params.samtools} sort -@ {params.threads} -o {output.sorted} {output.bam}
        {params.samtools} index -@ {params.threads} $sorted
        """
##################################
#   RUN VARIANT CALLING
##################################
rule gatk_calling:
    input:
        bam="{sample}_sorted.bam",
        ref=reference
    output:
        vcf="mutect2/{sample}_filt.bcf"
    params:
        gatk=gatk,
        bcftools=bcftools,
        threads=threads
    shell:
        """
        name="{wildcards.sample}"
        unfiltered="mutect2/${name}_unfiltered.vcf"
        filtered="{output.vcf}"
        
        {params.gatk} Mutect2 -R {input.ref} --f1r2-tar-gz mutect2/${name}_f1r2.tar.gz --native-pair-hmm-threads {params.threads} --germline-resource ref/somatic-hg38_af-only-gnomad.hg38.vcf.gz -I {input.bam} -O $unfiltered
        {params.gatk} LearnReadOrientationModel -I mutect2/${name}_f1r2.tar.gz -O mutect2/${name}_read_orientation.tar.gz
        {params.gatk} GetPileupSummaries -I {input.bam} -V ref/somatic-hg38_small_exac_common_3.hg38.vcf.gz -L ref/somatic-hg38_small_exac_common_3.hg38.vcf.gz -O mutect2/${name}_pileup_summaries.table
        {params.gatk} CalculateContamination -I mutect2/${name}_pileup_summaries.table -tumor-segmentation mutect2/${name}_segments.table -O mutect2/${name}_contamination.table
        {params.gatk} FilterMutectCalls -R {input.ref} -V $unfiltered --tumor-segmentation mutect2/${name}_segments.table --contamination-table mutect2/${name}_contamination.table --ob-priors mutect2/${name}_read_orientation.tar.gz -O $filtered
        {params.bcftools} view --threads {params.threads} $filtered -Ob -o {output.vcf}
        """
rule varscan_calling:
    input:
        bam='{sample}_sorted.bam',
        ref=reference
    output:
        vcf=temp('varscan/{sample}.vcf'),
        filtered='varscan/{sample}_filtered.bcf'
    params:
        varscan=varscan,
        samtools=samtools,
        threads=threads
    shell:
        """
        {params.samtools} mpileup -B -f {input.ref} | java -jar {params.varscan} pileup2snp --output-vcf 1 > {output.vcf}
        bcftools view --threads {params.threads} -i 'QUAL>50 && DP>10' -Ob -o {output.filtered} {output.vcf}
        """
#######################################
#   MERGE VARIANT CALLS TOGETHER
#######################################
rule merge_calls:
    input:
        mutect2="mutect2/{sample}_filt.bcf",
        varscan="varscan/{sample}_filtered.bcf"
    params:
        bcftools=bcftools
    output:
        overlap=temp('variants/{sample}_overlapping_variants.vcf'),
        compressed='variants/{sample}_overlapping_variants.bcf'
    shell:
        """
        {params.bcftools} isec -Ov -n=2 -w1 {input.mutect2} {input.varscan} -o {output.overlap}
        {params.bcftools} view -Ob -o {output.compressed} {output.overlap}
        """
######################################
#   PERFORM FUNCTIONAL ANNOTATION
######################################
rule functional_annotation:
    input:
        vcf="{sample}_overlapping_variants.vcf",
        ref=reference,
        funcotator=funcotator_db
    output:
        vcf="annotated/{sample}_funcotated.vcf"
    shell:
        """
        gatk Funcotator --variant {input.vcf} --reference {input.ref} --ref-version hg38 --data-sources-path {input.funcotator} --output {output.vcf} --output-file-format VCF
        """
