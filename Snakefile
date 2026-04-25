########################################
# RULE ALL (final target)
########################################

rule all:
    input:
        "results/climatomics/Monte-Carlo_permutation_robustTF.csv",
        "data/cistromics/ciscross.complete"


########################################
# 1. TRANSCRIPTOMICS
########################################

rule transcriptomics:
    input:
        array_dir="data/transcriptomics/microarray",
        rna_dir="data/transcriptomics/RNA-seq"
    output:
        all="results/transcriptomics/all_DEG_logFC.csv",
        final="results/transcriptomics/robustDEG.csv",
        threshold="results/transcriptomics/binomial_threshold.txt",
        deg=expand(
            "results/transcriptomics/{col}.csv",
            col=["upER","upLR","upVR","downER","downLR","downVR"]
        )
    conda:
        "envs/r_env.yaml"
    script:
        "scripts/1.Transcriptomics.R"

rule prepare_ciscross:
    input:
        deg=expand(
            "results/transcriptomics/{col}.csv",
            col=["upER","upLR","upVR","downER","downLR","downVR"]
        )
    output:
        flag="data/cistromics/ciscross.complete"
    message:
        "Waiting for manual cistromics step (max 30 min)..."
    shell:
        """
        for i in {{1..60}}; do
            if [ -f {output.flag} ]; then
                echo "Found ciscross.complete"
                exit 0
            fi
            echo "Waiting for data/cistromics/ciscross.complete... ($i/60)"
            sleep 30
        done

        echo ""
        echo "Timeout after 5 minutes."
        echo "Please complete manual step:"
        echo "  touch data/cistromics/ciscross.complete"
        exit 1
        """
########################################
# 2. CISTROMICS (TF enrichment)
########################################

rule cistromics:
    input:
        threshold="results/transcriptomics/binomial_threshold.txt",
        count_table="results/transcriptomics/robustDEG.csv",
        ciscross_dir="data/cistromics",
        flag="data/cistromics/ciscross.complete"
    output:
        tf_enrichment="results/cistromics/tf_cis_enrichment.csv"
    conda:
        "envs/r_env.yaml"
    script:
        "scripts/2.Cistromics.R"


########################################
# 3. TF → GENOMIC REGIONS (BED)
########################################

rule genomics_bed:
    input:
        tf_enrichment="results/cistromics/tf_cis_enrichment.csv"
    output:
        tf_bed="data/genomics/tf_regions.bed"
    conda:
        "envs/r_env.yaml"
    script:
        "scripts/3.Extract_TF_regions.R"


########################################
# 4. INDEX VCF (required for bcftools)
########################################

rule index_vcf:
    input:
        vcf="data/genomics/1001genomes_snp-short-indel_only_ACGTN.vcf.gz"
    output:
        tbi="data/genomics/1001genomes_snp-short-indel_only_ACGTN.vcf.gz.tbi"
    conda:
        "envs/bcftools_env.yaml"
    shell:
        "bcftools index -t {input.vcf}"


########################################
# 5. SUBSET VCF BY TF REGIONS
########################################

rule subset_vcf_by_tf:
    input:
        vcf="data/genomics/1001genomes_snp-short-indel_only_ACGTN.vcf.gz",
        bed="data/genomics/tf_regions.bed",
        tbi="data/genomics/1001genomes_snp-short-indel_only_ACGTN.vcf.gz.tbi"
    output:
        vcf="results/genomics/tf_variants.vcf.gz"
    conda:
        "envs/bcftools_env.yaml"
    shell:
        """
        bcftools view -R {input.bed} {input.vcf} | \
        bcftools norm -m-any -Oz -o {output.vcf}
        """


########################################
# 6. ANNOTATE VARIANTS (VEP)
########################################

rule annotate_variants:
    input:
        vcf="results/genomics/tf_variants.vcf.gz"
    output:
        vcf="results/genomics/tf_variants_annotated.vcf"
    conda:
        "envs/vep_env.yaml"
    shell:
        """
        vep -i {input.vcf} \
            --offline \
            --cache --dir_cache resource/vep \
            --cache_version 62 \
            --species arabidopsis_thaliana \
            --vcf \
            --force_overwrite \
            --variant_class \
            --o stdout \
        | filter_vep -filter "IMPACT is HIGH" \
            --o {output.vcf}
        """


########################################
# 7. EXTRACT GENOTYPES (ALT / REF / NA)
########################################

rule extract_variant_genotypes:
    input:
        vcf="results/genomics/tf_variants_annotated.vcf",
        bed="data/genomics/tf_regions.bed"
    output:
        alt="results/genomics/tf_variants_genotype_alt.csv",
        ref="results/genomics/tf_variants_genotype_ref.csv",
        na="results/genomics/tf_variants_genotype_na.csv"
    conda:
        "envs/r_env.yaml"
    script:
        "scripts/4.Extract_variant_genotypes.R"


########################################
# 8. CLIMATOMICS (Monte Carlo test)
########################################

rule climatomics:
    input:
        alt="results/genomics/tf_variants_genotype_alt.csv",
        CHELSA_temperature="data/climatomics/CHELSA_min_temp_ALL_1058_accessions.xlsx"
    output:
        mc="results/climatomics/Monte-Carlo_permutation_robustTF.csv"
    conda:
        "envs/r_env.yaml"
    script:
        "scripts/5.Climatomics.R"
