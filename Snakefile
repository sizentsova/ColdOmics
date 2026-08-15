configfile: "config.yaml"

########################################
# CONFIGURATION
########################################

CIS      = config["cistromics"]
PHASES   = ["ER", "LR", "VR"]
CIS_TAG  = "{}p_{}".format(CIS["upstream"], CIS["fdr"])
CIS_INDEX = "resource/ciscross/index/ciscross_{}_up{}.rds".format(
    CIS["collection"], CIS["upstream"])
CIS_REPORTS = expand(
    "data/cistromics/ciscross_{ph}_" + CIS_TAG + ".txt", ph=PHASES)

# scripts/2.Cistromics.R looks the reports up by their exact file names
# (see its tf_list block), so the config must produce those names.
if CIS_TAG != "1500p_0.05":
    raise ValueError(
        "scripts/2.Cistromics.R expects reports named "
        "ciscross_<PHASE>_1500p_0.05.txt, but this config produces "
        "ciscross_<PHASE>_{}.txt. Either set cistromics.upstream=1500 and "
        "cistromics.fdr=0.05, or update the tf_list block in "
        "scripts/2.Cistromics.R to match.".format(CIS_TAG)
    )

if CIS["source"] not in ("local", "web"):
    raise ValueError("cistromics.source must be 'local' or 'web'")

########################################
# RULE ALL (final target)
########################################

rule all:
    input:
        "results/climatomics/Monte-Carlo_permutation_robustTF.csv",
        CIS_REPORTS,
        "results/regulomics/Connectivity_network_summary.csv",
        "results/integrated/robust_TF_multilayer.csv"

########################################
# 1. TRANSCRIPTOMICS
########################################

rule transcriptomics:
    input:
        array="data/transcriptomics/microarray_DEG.rds",
        rna="data/transcriptomics/rna_DEG.rds"
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

########################################
# 2a. CISTROMICS INPUT (CisCross reports)
########################################

if CIS["source"] == "local":

    # Offline CisCross: runs the vendored engine on the step-1 DEG lists and
    # writes reports in the web service's format. No network access needed.
    rule ciscross_local:
        input:
            upER="results/transcriptomics/upER.csv",
            downER="results/transcriptomics/downER.csv",
            upLR="results/transcriptomics/upLR.csv",
            downLR="results/transcriptomics/downLR.csv",
            upVR="results/transcriptomics/upVR.csv",
            downVR="results/transcriptomics/downVR.csv",
            engine="resource/ciscross/enrichment.R",
            index=CIS_INDEX
        output:
            CIS_REPORTS
        params:
            fdr=CIS["fdr"],
            background=CIS.get("background", "")
        message:
            "Running CisCross locally ({}, {} bp promoters)".format(
                CIS["collection"], CIS["upstream"])
        conda:
            "envs/r_env.yaml"
        script:
            "scripts/2a.CisCross_local.R"

else:

    # Manual route: reports downloaded from the CisCross web service and placed
    # in data/cistromics/ by hand (see README).
    rule prepare_ciscross:
        input:
            deg=expand(
                "results/transcriptomics/{col}.csv",
                col=["upER","upLR","upVR","downER","downLR","downVR"]
            )
        output:
            CIS_REPORTS
        message:
            "Waiting for manual CisCross step (max 30 min)..."
        shell:
            """
            for i in $(seq 1 60); do
                if [ -f "{output[0]}" ] && [ -f "{output[1]}" ] \
                   && [ -f "{output[2]}" ]; then
                    echo "Found all CisCross reports"
                    exit 0
                fi

                echo "Waiting for CisCross results... ($i/60)"
                sleep 30
            done

            echo ""
            echo "Timeout after 30 minutes."
            echo "Place the CisCross reports in data/cistromics/:"
            echo "  {output}"
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
        reports=CIS_REPORTS
    output:
        tf_enrichment="results/cistromics/tf_cis_enrichment.csv"
    conda:
        "envs/r_env.yaml"
    script:
        "scripts/2.Cistromics.R"

########################################
# 3. TF → Network Connectivity
########################################

rule calculate_connectivity:
    input:
        ciscross_dir="data/cistromics",
        reports=CIS_REPORTS
    output:
        connectivity_per_TF="results/regulomics/Connectivity_per_TF.csv",
        connectivity_network="results/regulomics/Connectivity_network_summary.csv"
    conda:
        "envs/r_env.yaml"
    script:
        "scripts/3.Regulomics.R"

########################################
# 4.1. TF → GENOMIC REGIONS (BED)
########################################

rule genomics_bed:
    input:
        tf_enrichment="results/cistromics/tf_cis_enrichment.csv"
    output:
        tf_bed="data/genomics/tf_regions.bed"
    conda:
        "envs/r_env.yaml"
    script:
        "scripts/4.Extract_TF_regions.R"


########################################
# 4.2. INDEX VCF (required for bcftools)
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
# 4.3. SUBSET VCF BY TF REGIONS
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
# 4.4. ANNOTATE VARIANTS (VEP)
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
            --no_stats \
            --o stdout \
        | filter_vep -filter "IMPACT is HIGH" \
            --o {output.vcf}
        """


########################################
# 4.5. EXTRACT GENOTYPES (ALT / REF / NA)
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
# 5. CLIMATOMICS (Monte Carlo test)
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


########################################
# 6. INTEGRATE MULTI-OMICS
########################################

rule integrate:
    input:
        tf_cis="results/cistromics/tf_cis_enrichment.csv",
        tf_reg="results/regulomics/Connectivity_per_TF.csv",
        tf_gen="results/genomics/tf_variants_genotype_alt.csv",
        tf_clima="results/climatomics/Monte-Carlo_permutation_robustTF.csv"
    output:
        tf_int="results/integrated/robust_TF_multilayer.csv"
    conda:
        "envs/r_env.yaml"
    script:
        "scripts/6.Integrate_multi-omics_layers.R"


