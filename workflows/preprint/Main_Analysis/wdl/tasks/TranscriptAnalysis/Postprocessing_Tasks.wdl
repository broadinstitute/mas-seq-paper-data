version 1.0

import "../Structs.wdl"

task CreateCountMatrixFromAnnotatedBam {

    meta {
        description : "Creates a count matrix TSV file from the given annotated bam file.  Bam file must contain tags that indicate the gene/transcript (XG), cell barcode (CB), and umi (BX) of the read."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    input {
        File annotated_transcriptome_bam
        String prefix = "umi_tools_group"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size_gb = 20 + 11*ceil(size(annotated_transcriptome_bam, "GB"))

    command <<<
        set -euxo pipefail
        /python_scripts/create_count_matrix_from_annotated_bam.py -b ~{annotated_transcriptome_bam} -o ~{prefix}.tsv
    >>>

    output {
        File count_matrix = "~{prefix}.tsv"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             32,
        disk_gb:            disk_size_gb,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-transcript_utils:0.0.3"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task MergeBarcodeCounts {
    meta {
        description : "Merge all counts for each unique barcode in the given TSV file.  Assumes file is unheadered and have two columns: BARCODE COUNT.  Merging performed by adding all COUNTs for each BARCODE."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    input {
        File barcode_count_tsv
        String prefix = "merged_counts"

        RuntimeAttr? runtime_attr_override
    }

    # 20 gb - baseline storage for safety
    # 1x for the file itself
    # 2x for the results and some wiggle room.
    Int disk_size_gb = 20 + (3 * ceil(size(barcode_count_tsv, "GB")))

    command {
        /python_scripts/merge_barcode_counts.py ~{barcode_count_tsv}
        if [[ "~{prefix}.tsv" != "merged_counts.tsv" ]] ; then
            mv merged_counts.tsv "~{prefix}.tsv"
        fi
    }

    output {
        File merged_counts = "~{prefix}.tsv"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             16,
        disk_gb:            disk_size_gb,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-transcript_utils:0.0.6"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}


task CreateCountMatrixAnndataFromTsv {

    meta {
        description : "Creates a python anndata object from the given countmatrix tsv.  Expects the input to have been generated by CreateCountMatrixFromAnnotatedBam.  The resulting anndata object can be directly read into scanpy for single-cell analysis."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    input {
        File count_matrix_tsv
        File genome_annotation_gtf_file

        Boolean force_anndata_gencode_overwrite = false

        String prefix = "umi_tools_group"

        File? overlap_intervals
        String? overlap_interval_label
        File? gencode_reference_gtf_file

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size_gb = 20 + 4*ceil(size(count_matrix_tsv, "GB")) + 4*ceil(size(genome_annotation_gtf_file, "GB"))

    String overlap_intervals_arg = if defined(overlap_intervals)  then " --overlap-intervals " else ""
    String overlap_interval_label_arg = if defined(overlap_interval_label) then " --overlap-interval-label " else ""
    String gencode_reference_gtf_file_arg = if defined(gencode_reference_gtf_file) then " --gencode-reference-gtf " else ""

    String force_gencode_overwrite_flag = if force_anndata_gencode_overwrite then " --force-overwrite-gencode-overlaps " else ""

    command <<<
        set -euxo pipefail
        /python_scripts/create_count_matrix_anndata_from_tsv.py \
            -t ~{count_matrix_tsv} \
            -g ~{genome_annotation_gtf_file} \
            ~{overlap_intervals_arg}~{default="" sep=" --overlap-intervals " overlap_intervals} \
            ~{overlap_interval_label_arg}~{default="" sep=" --overlap-interval-label " overlap_interval_label} \
            ~{gencode_reference_gtf_file_arg}~{default="" sep=" --gencode-reference-gtf " gencode_reference_gtf_file} \
            ~{force_gencode_overwrite_flag} \
            -o ~{prefix}
    >>>

    output {
        File transcript_gene_count_anndata_h5ad = "~{prefix}_tx_gene_counts_adata.h5ad"
        Array[File] pickles = glob("*.pickle")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             32,
        disk_gb:            disk_size_gb,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-transcript_utils:0.0.7"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task SubsetCountsMatrixByGenes {

    meta {
        description : "Subsets a count matrix TSV file to contain only the transcripts from the given list of genes.  Assumes the count matrix was produced by comparison with Gencode (due to data formatting) and that the table is a TSV with samples as rows and transcripts as columns."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    input {
        File count_matrix_tsv
        Array[String] gene_names
    }

    parameter_meta {
        count_matrix_tsv : "TSV file containing the counts of each transcript expressed in a sample.  (Transcripts in columns.  Samples in rows.  One header row.)"
        gene_names : "Array of gene names for which to keep data from the given count matrix TSV."
    }

    # We're subsetting the file, so we should be able to get away with very little space here.
    # 1x for the file itself
    # 2x for the results and some wiggle room.
    Int disk_size = 3 * ceil(size(count_matrix_tsv, "GB"))

    command {
        /python_scripts/subset_count_matrix_by_gene.py ~{count_matrix_tsv} ~{sep=' ' gene_names}
    }

    output {
        File subset_count_matrix_tsv = "count_matrix_subset_by_gene.tsv"
        File subset_count_matrix_h5ad = "count_matrix_subset_by_gene.h5ad"
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/lr-transcript_utils:0.0.1"
        memory: 16 + " GiB"
        disks: "local-disk " + disk_size + " HDD"
        boot_disk_gb: 10
        preemptible: 0
        cpu: 8
    }
}
