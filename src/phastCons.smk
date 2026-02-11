CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
wildcard_constraints:
    chrom = "chr([0-9]+|X|Y)"
rule all:
    input:
        expand("../data/multiple_sequence_alignments/split_t2t_primates_maf_files/{chrom}.maf", chrom=CHROMS),
        expand("../data/genomes/chm13_split_chroms/{chrom}.fa", chrom=CHROMS),
        expand("../data/multiple_sequence_alignments/split_t2t_primates_ss_files/{chrom}.ss", chrom=CHROMS),
        "../data/multiple_sequence_alignments/long_cds_gffs/complete_annotation.longestIsoform.gff3",
        expand("../data/multiple_sequence_alignments/long_cds_gffs/{chrom}.gff3", chrom=CHROMS),
        expand("../data/multiple_sequence_alignments/primate_neutral_model/{chrom}.4d.codons.ss", chrom=CHROMS),
        "../data/multiple_sequence_alignments/primate_neutral_model/all_chroms.4d.codons.ss",
        "../data/multiple_sequence_alignments/primate_neutral_model/all_chroms.4d.sites.ss",
        "../data/multiple_sequence_alignments/primate_neutral_model/t2t_primates.neutral.mod",
        expand("../data/multiple_sequence_alignments/primate_phastcons_wig/{chrom}.phastCons.wig", chrom=CHROMS),
        "../data/multiple_sequence_alignments/primate_phastcons.bw"
rule hal2maf:
    input:
        hal = "../data/multiple_sequence_alignments/8-t2t-apes-2023v2.hal"
    output:
        maf = "../data/multiple_sequence_alignments/split_t2t_primates_maf_files/{chrom}.maf"
    params:
        ref = "hs1",
        targets = "hs1,GCA_028885655.2,GCA_028885625.2,GCA_029281585.2,GCA_029289425.2,GCA_028858775.2,GCA_028878055.2"
    shell:
        """
        hal2maf \
            --noDupes \
            --onlyOrthologs \
            --refSequence {wildcards.chrom} \
            --refGenome {params.ref} \
            --targetGenomes {params.targets} \
            {input.hal} \
            /dev/stdout | sed 's/hs1\.//g' | sed '2d' > {output.maf}
        """
rule split_chroms_fasta:
    input:
        genome = "../data/genomes/chm13.v2.0.fasta",
    output:
        chrom_fasta = "../data/genomes/chm13_split_chroms/{chrom}.fa",
    params:
        split_chroms_dir = "../data/genomes/chm13_split_chroms",
    shell:
        """
        mkdir -p {params.split_chroms_dir}
        samtools faidx {input.genome} {wildcards.chrom} > {output.chrom_fasta}
        """
rule convert_maf_to_ss:
    input:
        chrom_fa = "../data/genomes/chm13_split_chroms/{chrom}.fa",
        chrom_maf = "../data/multiple_sequence_alignments/split_t2t_primates_maf_files/{chrom}.maf",
    output:
        chrom_ss = "../data/multiple_sequence_alignments/split_t2t_primates_ss_files/{chrom}.ss"
    params:
        outdir = "../data/multiple_sequence_alignments/split_t2t_primates_ss_files"
    shell:
        """
        mkdir -p {params.outdir}
        msa_view \
            --in-format MAF \
            --refseq {input.chrom_fa} \
            --out-format SS \
            {input.chrom_maf} > {output.chrom_ss}

        sed -i 's/NAMES = chr[0-9XYM]*/NAMES = hs1/' {output}
        """
rule create_long_isoform_gff:
    input:
        gff3 = "../data/annotation_data/chm13/complete_annotation.gff3"
    output:
        gff3_long = "../data/multiple_sequence_alignments/long_cds_gffs/complete_annotation.longestIsoform.gff3",
    params:
        outdir = "../data/multiple_sequence_alignments/long_cds_gffs"
    shell:
        """
        mkdir -p {params.outdir}
        agat_sp_keep_longest_isoform.pl \
            --gff {input.gff3} \
            -o {output.gff3_long}
        """
rule split_long_gff:
    input:
        gff3_long = "../data/multiple_sequence_alignments/long_cds_gffs/complete_annotation.longestIsoform.gff3"
    output:
        "../data/multiple_sequence_alignments/long_cds_gffs/{chrom}.gff3"
    shell:
        """
        awk '$1 == "{wildcards.chrom}" || /^#/' {input.gff3_long} > {output}
        """
rule extract_4d_codons:
    input:
        maf = "../data/multiple_sequence_alignments/split_t2t_primates_maf_files/{chrom}.maf",
        gff = "../data/multiple_sequence_alignments/long_cds_gffs/{chrom}.gff3"
    output:
        ss_4d = "../data/multiple_sequence_alignments/primate_neutral_model/{chrom}.4d.codons.ss"
    params:
        ref = "hs1"
    shell:
        """
        msa_view \
            --in-format MAF \
            --4d \
            --features {input.gff} \
            {input.maf} > {output.ss_4d}

        sed -i 's/NAMES = chr[0-9XYM]*/NAMES = hs1/' {output}
        """
rule combine_4d_codons:
    input:
        expand("../data/multiple_sequence_alignments/primate_neutral_model/{chrom}.4d.codons.ss", chrom=CHROMS)
    output:
        agg_ss = "../data/multiple_sequence_alignments/primate_neutral_model/all_chroms.4d.codons.ss"
    params:
        targets = "hs1,GCA_028885655,GCA_028885625,GCA_029281585,GCA_029289425,GCA_028858775,GCA_028878055",
        ref = "hs1"
    shell:
        """
        msa_view \
            --aggregate {params.targets} \
            --out-format SS \
            --tuple-size 3 \
            --unordered-ss \
            {input} > {output.agg_ss}
        """
rule extract_4d_sites:
    input:
        ss = "../data/multiple_sequence_alignments/primate_neutral_model/all_chroms.4d.codons.ss"
    output:
        ss = "../data/multiple_sequence_alignments/primate_neutral_model/all_chroms.4d.sites.ss"
    shell:
        """
        msa_view {input.ss} --in-format SS --out-format SS --tuple-size 1 > {output.ss}
        """
rule phylofit_model:
    input:
        agg_ss = "../data/multiple_sequence_alignments/primate_neutral_model/all_chroms.4d.sites.ss",
        tree = "../data/filtered_tree.nh"
    output:
        mod = "../data/multiple_sequence_alignments/primate_neutral_model/t2t_primates.neutral.mod"
    threads: 100
    shell:
        """
        phyloFit \
            --tree {input.tree} \
            --subst-mod REV \
            --EM \
            --precision MED \
            --threads {threads} \
            --seed 42 \
            --out-root ../data/multiple_sequence_alignments/primate_neutral_model/t2t_primates.neutral \
            {input.agg_ss}
        """

rule phastcons:
    input:
        chrom_ss = "../data/multiple_sequence_alignments/split_t2t_primates_ss_files/{chrom}.ss",
        neutral_mod = "../data/multiple_sequence_alignments/primate_neutral_model/t2t_primates.neutral.mod"
    output:
        chrom_wig = "../data/multiple_sequence_alignments/primate_phastcons_wig/{chrom}.phastCons.wig",
        conserved_bed = "../data/multiple_sequence_alignments/primate_phastcons_bed/{chrom}.conserved.phastCons.bed"
    threads: 16
    params:
        outdir_wig = "../data/multiple_sequence_alignments/primate_phastcons_wig",
        outdir_bed = "../data/multiple_sequence_alignments/primate_phastcons_bed"
    shell:
        """
        mkdir -p {params.outdir_wig} {params.outdir_bed}
        phastCons \
            --msa-format SS \
            --target-coverage 0.3 \
            --expected-length 45 \
            --rho 0.3 \
            --viterbi {output.conserved_bed} \
            --threads {threads} \
            {input.chrom_ss} {input.neutral_mod} \
        | sed 's/chrom=(null)/chrom={wildcards.chrom}/g' \
        > {output.chrom_wig}
        """
rule combine_wig_to_bigwig:
    input:
        wigs = expand("../data/multiple_sequence_alignments/primate_phastcons_wig/{chrom}.phastCons.wig", chrom=CHROMS),
        beds = expand("../data/multiple_sequence_alignments/primate_phastcons_bed/{chrom}.conserved.phastCons.bed", chrom=CHROMS),
        chrom_sizes = "../data/annotation_data/chm13/chm13.v2.0.chrom.sizes"
    output:
        bigwig = "../data/multiple_sequence_alignments/primate_phastcons.bw",
        bed = "../data/multiple_sequence_alignments/primate_conserved.bed"
    shell:
        """
        cat {input.wigs} > {output.bigwig}.tmp.wig
        cat {input.beds} > {output.bed}.tmp
        sort -k1,1 -k2,2n {output.bed}.tmp > {output.bed}

        wigToBigWig {output.bigwig}.tmp.wig {input.chrom_sizes} {output.bigwig}

        rm {output.bigwig}.tmp.wig {output.bed}.tmp
        """