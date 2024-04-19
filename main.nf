params.outdir = "out"
params.include_pcr = false

process ancestral_genotypes {
    input:
    tuple val(id), path(ancestral_reconstruction)

    output:
    path "ancestral_genotypes.tsv"

    publishDir "$params.outdir", mode: 'copy'

    script:
    """
    echo "CHROM\tPOS\tREF\tALT\tCTVT\tHT" > ancestral_genotypes.tsv
    bcftools query \
      -s CTVT,HT \
      -f "%CHROM\t%POS\t%REF\t%ALT\t[%GT\t]\n" ${ancestral_reconstruction[0]} \
      | sed 's/\t\$//' \
      >> ancestral_genotypes.tsv
    """
}

process gene_information {
    input:
    path gtf

    output:
    path "gene_information.tsv"

    publishDir "$params.outdir", mode: 'copy'

    script:
    """
    gunzip -c ${gtf} | awk -F'\\t' \
    'BEGIN { OFS="\\t" }
    {
        if (NR == 1) {
            print "CHROM", "FEATURE", "START", "END", "biotype", "geneID", "geneName"
        }
    }
    {
        if (\$3 == "gene") {
            split(\$9, arr, ";")
            for(i in arr) {
                if(arr[i] ~ /gene_id/) {split(arr[i], id, " "); geneId = id[2]}
                if(arr[i] ~ /gene_name/) {split(arr[i], name, " "); geneName = name[2]}
                if(arr[i] ~ /gene_biotype/) {split(arr[i], biotype, " "); geneBiotype = biotype[2]}
            }
            sub(/"/, "", geneId); sub(/"/, "", geneId);  # Remove quotation marks
            sub(/"/, "", geneName); sub(/"/, "", geneName);
            sub(/"/, "", geneBiotype); sub(/"/, "", geneBiotype);
            print \$1, \$3, \$4, \$5, geneBiotype, geneId, geneName
        }
    }' > gene_information.tsv
    """
}

process extract_exome_variants {
    input: 
    path dataset
    path gtf

    output: 
    path "exome_variants.tsv.gz"

    publishDir "$params.outdir", mode: 'copy'

    script:
    """
    exome_variants.R \
      --dataset $dataset \
      --gtf $gtf
    """
}

process allele_counting {
    input:
    tuple val(id), path(alignment), path(variants) 

    output:
    path "${id}.allele_counts.tsv.gz"

    publishDir "$params.outdir", mode: 'copy'

    script:
    if ( params.include_pcr)
        """
        alco \
          --bamfile ${alignment[0]} \
          --locifile $variants \
          -g 10000 \
          -F 2816 \
          > ${id}.allele_counts.tsv
        annotate_germline_status.R \
          --counts ${id}.allele_counts.tsv \
          --variants $variants \
          --output ${id}.allele_counts.tsv.gz
        rm ${id}.allele_counts.tsv
        """
    else
        """
        alco \
          --bamfile ${alignment[0]} \
          --locifile $variants \
          -g 10000 \
          -F 3840 \
          > ${id}.allele_counts.tsv
        annotate_germline_status.R \
          --counts ${id}.allele_counts.tsv \
          --variants $variants \
          --output ${id}.allele_counts.tsv.gz
        rm ${id}.allele_counts.tsv
        """
}

process combine_counts {
    input:
    path(counts)

    output:
    path "combined_counts.tsv.gz"

    publishDir "$params.outdir", mode: 'copy'

    script:
    """
    header=\$(gunzip -c ${counts[0]} | head -n 1)
    printf "%s\\t%s\\n" "samplename" "\$header" > header.tsv
    for f in ${counts}; do
      id=\$(basename \$f .allele_counts.tsv.gz)
      gunzip -c \$f \
        | awk -v id="\$id" '{print id "\\t" \$0}' \
        | tail -n +2 \
        >> combined_counts.tsv
    done
    cat header.tsv combined_counts.tsv > tmp && mv tmp combined_counts.tsv
    rm header.tsv
    gzip combined_counts.tsv
    """
} 

process compute_size_factors {
    input:
    path(counts)

    output:
    path "size_factors.tsv"
    path "size_factors.alt.tsv"

    publishDir "$params.outdir", mode: 'copy'

    script:
    """
    compute_size_factors.R \
      --input ${counts} \
      --output size_factors.tsv
    """
}

process cram_to_bam {
    input:
    tuple val(id), path(cram)

    output:
    tuple val(id), path("${id}.bam*")

    script:
    """
    samtools view -b --threads ${task.cpus} ${cram[0]} > ${id}.bam
    samtools index -@ ${task.cpus} ${id}.bam
    """
}


process subread_counts {
    label 'subread'

    input:
    tuple val(id), path(bams), path(gtf)

    output:
    path "subread_counts.${id}.txt", emit: "counts"
    path "subread_counts.${id}.txt.summary", emit: "summary"

    script:
    if ( params.include_pcr)
        """
        featureCounts \
          -a ${gtf} \
          -F GTF \
          -t exon \
          --countReadPairs \
          -p \
          -o subread_counts.${id}.txt \
          -T ${task.cpus} \
          ${bams[0]}
        """
    else
        """
        featureCounts \
          -a ${gtf} \
          -F GTF \
          -t exon \
          --ignoreDup \
          --countReadPairs \
          -p \
          -o subread_counts.${id}.txt \
          -T ${task.cpus} \
          ${bams[0]}
        """
}

process compute_deseq_size_factors {
    input:
    tuple val(id), path(counts) 

    output:
    tuple val(id), path("size_factors.${id}.tsv")

    publishDir "$params.outdir", mode: 'copy'

    script:
    """
    compute_size_factors_with_deseq.R \
      --counts ${counts} \
      --output size_factors.${id}.tsv
    """
}

process make_base_table {
    input:
    path dataset
    path counts
    path variants
    path ht_regions
    path ancestral_states
    path ctvt_a

    output:
    path "base_table.tsv"

    publishDir "$params.outdir", mode: 'copy'

    script:
    """
    make_base_table.R \
      --dataset $dataset \
      --counts $counts \
      --variants $variants \
      --ht_regions $ht_regions \
      --ancestral_states $ancestral_states \
      --ctvt_a $ctvt_a \
      --output base_table.tsv
    """
}

process plot_expression {
    input:
    tuple val(id), path(size_factors), path(table)

    output:
    path "expression.*.${id}.pdf"
    path "debug_tables_${id}/*.csv"

    publishDir "$params.outdir/expression_plots", mode: 'copy'

    script:
    """
    ht_expression.R \
        --table ${table} \
        --size_factors ${size_factors} \
        --scatterplot_out expression.scatterplot.${id}.pdf \
        --beeswarm_out expression.beeswarm.${id}.pdf \
        --debug_table_out debug_tables_${id}
    """
}

def remove_duplicate_filepair_keys(primary_ch, secondary_ch) {
    // primary_ch and secondary_ch are channels produced by
    // `fromFilePairs`. If any keys are present in both channels,
    // remove the key from the secondary channel.
    // This will avoid duplicating work on a bam file if it has
    // a bai and a csi index.
    return primary_ch.concat(secondary_ch).groupTuple() |
        map { it -> tuple(it[0], it[1][0]) }
}

workflow {
    // INPUT HANDLING
    dataset = Channel.fromPath("$params.dataset", checkIfExists: true)
    gtf = Channel.fromPath("$params.gtf", checkIfExists: true)
    ancestral_reconstruction = Channel.fromFilePairs("$params.ancestral{,.csi}",
        checkIfExists: true)
    ctvt_a = Channel.fromPath("$params.ctvt_a", checkIfExists: true)
    ht_regions = Channel.fromPath("$params.ht_regions", checkIfExists: true)

    // Gather all bam/cram alignment files from the input directory
    bam_bai_inputs   = Channel.fromFilePairs("${params.alignments}/*.bam{,.bai}")
    bam_csi_inputs   = Channel.fromFilePairs("${params.alignments}/*.bam{,.csi}")
    bam_inputs = remove_duplicate_filepair_keys(bam_csi_inputs, bam_bai_inputs)
    cram_inputs = Channel.fromFilePairs("${params.alignments}/*.cram{,.crai}")
    alignments = bam_inputs.concat(cram_inputs)
    // END INPUT HANDLING

    //BEGIN WORKFLOW
    bam_alignments = alignments | cram_to_bam
    gene_info = gene_information(gtf)
    ancestral_genotypes = ancestral_genotypes(ancestral_reconstruction)
    variants = extract_exome_variants(dataset, gtf)
    counts = bam_alignments.combine(variants) | allele_counting
    combined = counts.collect() | combine_counts
    size_factors = combined | compute_size_factors

    // DESeq2 can compute size factors for us
    // but we need to generate fragment counts first using subread
    fragment_counts = subread_counts(bam_alignments.combine(gtf))

    // Build a channel for each size factor computation
    frag_counts = 
        fragment_counts.counts
            .collect()
            .map { tuple('subread', it) }
    size_factors = compute_deseq_size_factors(frag_counts)
    base_table = make_base_table(dataset, combined, variants, ht_regions, ancestral_genotypes, ctvt_a)

    // Plot expression
    size_factors.combine(base_table) | plot_expression
} 
