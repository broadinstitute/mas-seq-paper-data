version 1.0

import "tasks/PBUtils.wdl" as PB
import "tasks/Utils.wdl" as Utils
import "tasks/Finalize.wdl" as FF
import "tasks/AlignReads.wdl" as AR
import "tasks/Cartographer.wdl" as CART
import "tasks/ReadsMetrics.wdl" as RM
import "tasks/AlignedMetrics.wdl" as AM
import "tasks/Ten_X_Tool.wdl" as TENX
import "tasks/JupyterNotebooks.wdl" as JUPYTER
import "tasks/Longbow.wdl" as LONGBOW

import "tasks/StringTie2.wdl"

import "tasks/TranscriptAnalysis/UMI_Tools.wdl" as UMI_TOOLS
import "tasks/TranscriptAnalysis/Postprocessing_Tasks.wdl" as TX_POST
import "tasks/TranscriptAnalysis/Preprocessing_Tasks.wdl" as TX_PRE

workflow PB10xMasSeqSingleFlowcellv3 {

    meta {
        description : "This workflow is designed to process data from the MASSeq v2 protocol and produce aligned reads that are ready for downstream analysis (e.g. transcript isoform identification).  It takes in a raw PacBio run folder location on GCS and produces a folder containing the aligned reads and other processed data."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    input {
        String gcs_input_dir
        String gcs_out_root_dir = "gs://broad-dsde-methods-long-reads-outgoing/PB10xMasSeqSingleFlowcellv3"

        File cell_barcode_whitelist = "gs://broad-dsde-methods-long-reads/resources/MASseq_0.0.2/737K-august-2016.txt"

        # NOTE: Reference for un-split CCS reads:
        File ref_fasta =  "gs://broad-dsde-methods-long-reads/resources/references/grch38_noalt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"
        File ref_fasta_index = "gs://broad-dsde-methods-long-reads/resources/references/grch38_noalt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa.fai"
        File ref_fasta_dict = "gs://broad-dsde-methods-long-reads/resources/references/grch38_noalt/GCA_000001405.15_GRCh38_no_alt_analysis_set.dict"

        # NOTE: Reference for array elements:
        File transcriptome_ref_fasta =  "gs://broad-dsde-methods-long-reads/resources/gencode_v37/gencode.v37.pc_transcripts.fa"
        File transcriptome_ref_fasta_index = "gs://broad-dsde-methods-long-reads/resources/gencode_v37/gencode.v37.pc_transcripts.fa.fai"
        File transcriptome_ref_fasta_dict = "gs://broad-dsde-methods-long-reads/resources/gencode_v37/gencode.v37.pc_transcripts.dict"

        File genome_annotation_gtf = "gs://broad-dsde-methods-long-reads/resources/gencode_v37/gencode.v37.primary_assembly.annotation.gtf"

        File jupyter_template_static = "gs://broad-dsde-methods-long-reads/resources/MASseq_0.0.2/MAS-seq_QC_report_template-static.ipynb"
        File workflow_dot_file = "gs://broad-dsde-methods-long-reads/resources/MASseq_0.0.2/PB10xMasSeqArraySingleFlowcellv2.dot"

        File intervals_of_interest = "gs://broad-dsde-methods-long-reads/resources/MASseq_0.0.2/gencode.37.TCR_intervals.tsv"
        String interval_overlap_name = "is_tcr_overlapping"

        File short_read_umis_tsv

        String starcode_extra_params = "--dist 2 --sphere"

        String expanded_cbc_tag = "CR"

        File? illumina_barcoded_bam

        File? cell_barcode_freq_tsv

        # Default here is 0 because ccs uncorrected reads all seem to have RQ = -1.
        # All pathologically long reads also have RQ = -1.
        # This way we preserve the vast majority of the data, even if it has low quality.
        # We can filter it out at later steps.
        Float min_read_quality = 0.0
        Int max_reclamation_length = 60000

        Boolean is_SIRV_data = false
        String mas_seq_model = "mas15v2"

        Int ccs_lev_dist = 2
        Int clr_lev_dist = 3

        Int ccs_umi_padding = 2
        Int clr_umi_padding = 2

        Int ccs_cbc_padding = 3
        Int clr_cbc_padding = 3

        # Add a suffix here for our out directory so we can label runs:
        String out_dir_suffix = ""

        String? sample_name
    }

    parameter_meta {
        gcs_input_dir : "Input folder on GCS in which to search for BAM files to process."
        gcs_out_root_dir : "Root output GCS folder in which to place results of this workflow."

        cell_barcode_whitelist : "Text file containing a whitelist of cell barcodes for the single-cell library prep."

        ref_fasta : "FASTA file containing the reference sequence to which the input data should be aligned before splitting into array elements."
        ref_fasta_index : "FASTA index file for the given ref_fasta file."
        ref_fasta_dict : "Sequence dictionary file for the given ref_fasta file."

        transcriptome_ref_fasta : "FASTA file containing the reference sequence to which the array elements should be aligned."
        transcriptome_ref_fasta_index : "FASTA index file for the given transcriptome_ref_fasta file."
        transcriptome_ref_fasta_dict : "Sequence dictionary file for the given transcriptome_ref_fasta file."

        genome_annotation_gtf : "Gencode GTF file containing genome annotations for the organism under study (usually humans).  This must match the given reference version and transcriptiome reference (usually hg38)."

        jupyter_template_static : "Jupyter notebook / ipynb file containing a template for the QC report which will contain static plots.  This should contain the same information as the jupyter_template_interactive file, but with static images."
        workflow_dot_file : "DOT file containing the representation of this WDL to be included in the QC reports.  This can be generated with womtool."

        intervals_of_interest : "[optional] An interval list file containing intervals to mark in the final anndata object as overlapping the transcripts.  Defaults to a T-cell receptor interval list."
        interval_overlap_name : "[optional] The name of the annotation to add to the final anndata object for the column containing the overlap flag for transcripts that overlap intervals in the given intervals_of_interest file.  Default: is_tcr_overlapping"

        illumina_barcoded_bam : "[optional] Illumina short reads file from a replicate of this same sample.  Used to perform cell barcode corrections."

        min_read_quality : "[optional] Minimum read quality for reads to have to be included in our data (Default: 0.0)."
        max_reclamation_length : "[optional] Maximum length (in bases) that a read can be to attempt to reclaim from CCS rejection (Default: 60000)."

        is_SIRV_data : "[optional] true if and only if the data in this sample are from the SIRV library prep.  false otherwise (Default: false)"
        mas_seq_model : "[optional] built-in mas-seq model to use (Default: mas15)"

        sample_name : "[optional] The name of the sample to associate with the data in this workflow."
    }

    # Create some runtime attributes that will force google to do network transfers really fast:
    RuntimeAttr fast_network_attrs = object {
        cpu_cores:  4,
        mem_gb:     32,
        disk_type:  "LOCAL",
        preemptible_tries:  0
    }
    RuntimeAttr fast_network_attrs_preemptible = object {
        cpu_cores:  4,
        mem_gb:     32,
        disk_type:  "LOCAL",
        preemptible_tries:  1
    }

    # Call our timestamp so we can store outputs without clobbering previous runs:
    call Utils.GetCurrentTimestampString as t_01_WdlExecutionStartTimestamp { input: }

    String outdir = sub(gcs_out_root_dir, "/$", "")

    call PB.FindBams as t_02_FindBams { input: gcs_input_dir = gcs_input_dir }
    call PB.FindZmwStatsJsonGz as t_03_FindZmwStatsJsonGz { input: gcs_input_dir = gcs_input_dir }

    # Check here if we found ccs bams or subread bams:
    Boolean use_subreads = t_02_FindBams.has_subreads
    Array[String] top_level_bam_files = if use_subreads then t_02_FindBams.subread_bams else t_02_FindBams.ccs_bams

    # Make sure we have **EXACTLY** one bam file to run on:
    if (length(top_level_bam_files) != 1) {
        call Utils.FailWithWarning as t_04_WARN1 { input: warning = "Error: Multiple BAM files found.  Cannot continue!" }
    }

    if (use_subreads) {
        call Utils.FailWithWarning as t_05_WARN2 { input: warning = "Error: This workflow now only supports data from the Sequel IIe." }
    }

    # Alias our bam file so we can work with it easier:
    File reads_bam = top_level_bam_files[0]

    call PB.GetPbReadGroupInfo as t_06_GetReadGroupInfo { input: gcs_bam_path = reads_bam }
    call PB.GetRunInfo as t_07_GetRunInfo { input: subread_bam = reads_bam }

    String SM  = select_first([sample_name, t_07_GetRunInfo.run_info["SM"]])
    String PL  = "PACBIO"
    String PU  = t_07_GetRunInfo.run_info["PU"]
    String DT  = t_07_GetRunInfo.run_info["DT"]
    String ID  = PU
    String DS  = t_07_GetRunInfo.run_info["DS"]
    String DIR = SM + "." + ID

    String RG_subreads  = "@RG\\tID:~{ID}.subreads\\tSM:~{SM}\\tPL:~{PL}\\tPU:~{PU}\\tDT:~{DT}"
    String RG_consensus = "@RG\\tID:~{ID}.consensus\\tSM:~{SM}\\tPL:~{PL}\\tPU:~{PU}\\tDT:~{DT}"
    String RG_array_elements = "@RG\\tID:~{ID}.array_elements\\tSM:~{SM}\\tPL:~{PL}\\tPU:~{PU}\\tDT:~{DT}"

    # Check to see if we need to annotate our reads:
    call LONGBOW.CheckForAnnotatedArrayReads as t_08_CheckForAnnotatedReads {
        input:
            bam = reads_bam
    }

    File read_pbi = sub(reads_bam, ".bam$", ".bam.pbi")
    call PB.ShardLongReads as t_09_ShardLongReads {
        input:
            unaligned_bam = reads_bam,
            unaligned_pbi = read_pbi,
            prefix = SM + "_shard",
            num_shards = 300,
    }

    ## No more preemption on this sharding - takes too long otherwise.
    RuntimeAttr disable_preemption_runtime_attrs = object {
        preemptible_tries: 0
    }

    scatter (sharded_reads in t_09_ShardLongReads.unmapped_shards) {

        String fbmrq_prefix = basename(sharded_reads, ".bam")

        # Filter out the kinetics tags from PB files:
        call PB.RemoveKineticsTags as t_10_RemoveKineticsTags {
            input:
                bam = sharded_reads,
                prefix = SM + "_kinetics_removed"
        }

        # Handle setting up the things that we need for further processing of CCS-only reads:
        call PB.FindCCSReport as t_11_FindCCSReport {
            input:
                gcs_input_dir = gcs_input_dir
        }

        # 1 - filter the reads by the minimum read quality:
        call Utils.Bamtools as t_12_FilterS2EByMinReadQuality {
            input:
                bamfile = t_10_RemoveKineticsTags.bam_file,
                prefix = fbmrq_prefix + "_good_reads",
                cmd = "filter",
                args = '-tag "rq":">=' + min_read_quality + '"',
                runtime_attr_override = disable_preemption_runtime_attrs
        }

        # 1.5 - Get the "rejected" reads:
        call Utils.Bamtools as t_13_GetS2ECcsRejectedReads {
            input:
                bamfile = t_10_RemoveKineticsTags.bam_file,
                prefix = fbmrq_prefix + "_rejected_reads",
                cmd = "filter",
                args = '-tag "rq":"<' + min_read_quality + '"',
                runtime_attr_override = disable_preemption_runtime_attrs
        }

        #################################################################################################################
        # 2 - Get reads we can reclaim:
        call Utils.Bamtools as t_14_ExtractS2ECcsReclaimableReads {
            input:
                bamfile = t_10_RemoveKineticsTags.bam_file,
                prefix = fbmrq_prefix + "_reads_for_ccs_reclamation",
                cmd = "filter",
                args = '-tag "rq":"<' + min_read_quality + '" -length "<=' + max_reclamation_length + '"',
                runtime_attr_override = disable_preemption_runtime_attrs
        }

        if ( ! t_08_CheckForAnnotatedReads.bam_has_annotations ) {
            # 3: Longbow annotate ccs reads
            call LONGBOW.Annotate as t_15_AnnotateS2ECCSReads {
                input:
                    reads = t_12_FilterS2EByMinReadQuality.bam_out,
                    model = mas_seq_model
            }
            # 4: Longbow annotate reclaimable reads
            call LONGBOW.Annotate as t_16_AnnotateS2EReclaimableReads {
                input:
                    reads = t_14_ExtractS2ECcsReclaimableReads.bam_out,
                    model = mas_seq_model
            }
        }

        File annotated_S2E_ccs_file = if t_08_CheckForAnnotatedReads.bam_has_annotations then t_12_FilterS2EByMinReadQuality.bam_out else select_first([t_15_AnnotateS2ECCSReads.annotated_bam])
        File annotated_S2E_reclaimable_file = if t_08_CheckForAnnotatedReads.bam_has_annotations then t_14_ExtractS2ECcsReclaimableReads.bam_out else select_first([t_16_AnnotateS2EReclaimableReads.annotated_bam])

        # 5: Longbow filter ccs annotated reads
        call LONGBOW.Filter as t_17_FilterS2ECCSReads {
            input:
                bam = annotated_S2E_ccs_file,
                prefix = SM + "_subshard",
                model = mas_seq_model
        }

        # 6: Longbow filter ccs reclaimable reads
        call LONGBOW.Filter as t_18_FilterS2EReclaimableReads {
            input:
                bam = annotated_S2E_reclaimable_file,
                prefix = SM + "_subshard",
                model = mas_seq_model
        }

        # 7: PBIndex CCS reads
        call PB.PBIndex as t_19_PbIndexS2ELongbowPassedCcsReads {
            input:
                bam = t_17_FilterS2ECCSReads.passed_reads
        }

        call PB.PBIndex as t_20_PbIndexS2ELongbowFailedCcsReads {
            input:
                bam = t_17_FilterS2ECCSReads.failed_reads
        }

        # 8: PBIndex reclaimable reads
        call PB.PBIndex as t_21_PbIndexS2ELongbowPassedReclaimedReads {
            input:
                bam = t_18_FilterS2EReclaimableReads.passed_reads
        }
        call PB.PBIndex as t_22_PbIndexS2ELongbowFailedReclaimableReads {
            input:
                bam = t_18_FilterS2EReclaimableReads.failed_reads
        }

        #####################################################################################################################
        # Now we have CCS and Reclaimed reads.
        #############################################

        # New pipeline steps:
        #     (1) minimap2 CCS reads in hifi mode
        #     (2) minimap2 CLR reads in shitty mode
        #     (3) merge 1 and 2 alignments
        #     <ADD FILTERING HERE>
        #     D2 sphere for cell barcode correction
        #     (4) annotate merged alignments with SQANTI3
        #     (5) cell barcode error correction with starcode (Levenshtein 1)
        #     (6) UMI error correction with umi-tools (Levenshtein 2), grouping by SQANTI3
        #     (7) generate count matrix

        # Shard our CCS reads into smaller problems to do work on array elements:
        call PB.ShardLongReads as t_23_ShardS2ECcsLongbowPassedReads {
            input:
                unaligned_bam = t_17_FilterS2ECCSReads.passed_reads,
                unaligned_pbi = t_19_PbIndexS2ELongbowPassedCcsReads.pbindex,
                prefix = SM + "_ccs_reads_subshard",
                num_shards = 10,
        }
        scatter (s2e_ccs_longbow_passed_shard in t_23_ShardS2ECcsLongbowPassedReads.unmapped_shards) {
            # Segment CCS reads into array elements:
            call LONGBOW.Segment as t_24_SegmentS2ECcsReads {
                input:
                    annotated_reads = s2e_ccs_longbow_passed_shard,
                    prefix = SM + "_ccs_array_elements_subshard",
                    extra_args = "-b",
                    model = mas_seq_model
            }

            # Now remove all -END reads:
            # We do this here to speed up all other calculations.
            call Utils.RemoveMasSeqTruncatedReads as t_25_RemoveMasSeqTruncatedReadsFromCcsReads {
                input:
                    bam_file = t_24_SegmentS2ECcsReads.segmented_bam
            }
        }

        # Merge Filtered CCS reads together:
        call Utils.MergeBams as t_26_MergeCCSArrayElementShards {
            input:
                bams = t_24_SegmentS2ECcsReads.segmented_bam,
                prefix = SM + "_ccs_array_elements_shard"
        }

        # Merge CCS Barcode Conf files:
        call Utils.MergeFiles as t_27_MergeCCSBarcodeConfShards {
            input:
                files_to_merge = t_24_SegmentS2ECcsReads.barcode_conf_file,
                merged_file_name = SM + "_ccs_array_element_barcode_confs_shard.txt"
        }

        # Merge Filtered CCS reads with no ends together:
        call Utils.MergeBams as t_28_MergeCCSArrayElementsNonTruncatedShards {
            input:
                bams = t_25_RemoveMasSeqTruncatedReadsFromCcsReads.bam,
                prefix = SM + "_ccs_array_elements_no_ends_shard"
        }

        #####################

        # Shard our CCS reclaimed reads into smaller problems to do work on array elements:
        call PB.ShardLongReads as t_29_ShardS2ECcsReclaimedReads {
            input:
                unaligned_bam = t_18_FilterS2EReclaimableReads.passed_reads,
                unaligned_pbi = t_21_PbIndexS2ELongbowPassedReclaimedReads.pbindex,
                prefix = SM + "_ccs_reclaimed_reads_subshard",
                num_shards = 10,
        }
        scatter (s2e_ccs_reclaimed_shard in t_29_ShardS2ECcsReclaimedReads.unmapped_shards) {
             # Segment Reclaimed reads into array elements:
            call LONGBOW.Segment as t_30_SegmentS2ECcsReclaimedReads {
                input:
                    annotated_reads = s2e_ccs_reclaimed_shard,
                    prefix = SM + "_ccs_reclaimed_array_elements_subshard",
                    extra_args = "-b",
                    model = mas_seq_model
            }

            # Now remove all -END reads:
            # We do this here to speed up all other calculations.
            call Utils.RemoveMasSeqTruncatedReads as t_31_RemoveMasSeqTruncatedReadsFromCcsReclaimedReads {
                input:
                    bam_file = t_30_SegmentS2ECcsReclaimedReads.segmented_bam
            }
        }

        # Merge Filtered CCS Reclaimed reads together:
        call Utils.MergeBams as t_32_MergeCCSReclaimedArrayElementShards {
            input:
                bams = t_30_SegmentS2ECcsReclaimedReads.segmented_bam,
                prefix = SM + "_ccs_reclaimed_array_elements_shard"
        }

        # Merge CCS Barcode Conf files:
        call Utils.MergeFiles as t_33_MergeCCSReclaimedBarcodeConfShards {
            input:
                files_to_merge = t_30_SegmentS2ECcsReclaimedReads.barcode_conf_file,
                merged_file_name = SM + "_ccs_reclaimed_array_element_barcode_confs_shard.txt"
        }

        # Merge Filtered CCS reads with no ends together:
        call Utils.MergeBams as t_34_MergeCCSReclaimedArrayElementsNonTruncatedShards {
            input:
                bams = t_31_RemoveMasSeqTruncatedReadsFromCcsReclaimedReads.bam,
                prefix = SM + "_ccs_reclaimed_array_elements_no_ends_shard"
        }

        ###############

        # Now align the array elements with their respective alignment presets.
        # NOTE: We use the non-truncated reads because we only want the good stuff.

        # Align CCS reads to the genome:
        call AR.Minimap2 as t_35_AlignCCSArrayElementsToGenome {
            input:
                reads      = [ t_28_MergeCCSArrayElementsNonTruncatedShards.merged_bam ],
                ref_fasta  = ref_fasta,
                map_preset = "splice:hq"
        }

        # Align Reclaimed reads to the genome:
        call AR.Minimap2 as t_36_AlignReclaimedArrayElementsToGenome {
            input:
                reads      = [ t_34_MergeCCSReclaimedArrayElementsNonTruncatedShards.merged_bam ],
                ref_fasta  = ref_fasta,
                map_preset = "splice"
        }

        ##############

        # Now restore the tags to the aligned bam files:
        call TENX.RestoreAnnotationstoAlignedBam as t_37_RestoreAnnotationsToGenomeAlignedCCSBam {
            input:
                annotated_bam_file = t_28_MergeCCSArrayElementsNonTruncatedShards.merged_bam,
                aligned_bam_file = t_35_AlignCCSArrayElementsToGenome.aligned_bam,
                tags_to_ignore = [],
                mem_gb = 8,
        }

        call TENX.RestoreAnnotationstoAlignedBam as t_38_RestoreAnnotationsToGenomeAlignedReclaimedBam {
            input:
                annotated_bam_file = t_34_MergeCCSReclaimedArrayElementsNonTruncatedShards.merged_bam,
                aligned_bam_file = t_36_AlignReclaimedArrayElementsToGenome.aligned_bam,
                tags_to_ignore = [],
                mem_gb = 8,
        }
    }

    # Merge the arrays:
    call Utils.MergeBams as t_39_MergeCCSLongbowAnnotatedArrayReads {
        input:
            bams = annotated_S2E_ccs_file,
            prefix = SM + "_ccs_array_reads_longbow_annotated"
    }
    call PB.PBIndex as t_40_PbIndexMergedCCSLongbowAnnotatedArrayReads { input: bam = t_39_MergeCCSLongbowAnnotatedArrayReads.merged_bam }

    call Utils.MergeBams as t_41_MergeCCSReclaimableLongbowAnnotatedArrayReads {
        input:
            bams = annotated_S2E_reclaimable_file,
            prefix = SM + "_ccs_reclaimable_array_reads_longbow_annotated"
    }
    call PB.PBIndex as t_42_PbIndexMergedCCSReclaimableLongbowAnnotatedArrayReads { input: bam = t_41_MergeCCSReclaimableLongbowAnnotatedArrayReads.merged_bam }

    call Utils.MergeBams as t_43_MergeCCSLongbowPassedArrayReads {
        input:
            bams = t_17_FilterS2ECCSReads.passed_reads,
            prefix = SM + "_ccs_array_reads_longbow_passed"
    }
    call PB.PBIndex as t_44_PbIndexMergedCCSLongbowPassedArrayReads { input: bam = t_43_MergeCCSLongbowPassedArrayReads.merged_bam }

    call Utils.MergeBams as t_45_MergeCCSLongbowFailedArrayReads {
        input:
            bams = t_17_FilterS2ECCSReads.failed_reads,
            prefix = SM + "_ccs_array_reads_longbow_failed"
    }
    call PB.PBIndex as t_46_PbIndexMergedCCSLongbowFailedArrayReads { input: bam = t_45_MergeCCSLongbowFailedArrayReads.merged_bam }

    call Utils.MergeBams as t_47_MergeCCSReclaimedArrayReads {
        input:
            bams = t_18_FilterS2EReclaimableReads.passed_reads,
            prefix = SM + "_ccs_reclaimed_array_reads_longbow_passed"
    }
    call PB.PBIndex as t_48_PbIndexMergedCCSReclaimedArrayReads { input: bam = t_47_MergeCCSReclaimedArrayReads.merged_bam }

    call Utils.MergeBams as t_49_MergeCCSUnreclaimableArrayReads {
        input:
            bams = t_18_FilterS2EReclaimableReads.failed_reads,
            prefix = SM + "_ccs_unreclaimable_array_reads_longbow_failed"
    }
    call PB.PBIndex as t_50_PbIndexMergedCCSUnreclaimableReclaimedArrayReads { input: bam = t_49_MergeCCSUnreclaimableArrayReads.merged_bam }

    call Utils.MergeBams as t_51_MergeLongbowPassedReads {
        input:
#            bams = [t_43_MergeCCSLongbowPassedArrayReads.merged_bam, t_47_MergeCCSReclaimedArrayReads.merged_bam],
            bams = flatten([t_17_FilterS2ECCSReads.passed_reads, t_18_FilterS2EReclaimableReads.passed_reads]),
            prefix = SM + "_longbow_passed_array_reads"
    }
    call PB.PBIndex as t_52_PbIndexMergedLongbowPassingReads { input: bam = t_51_MergeLongbowPassedReads.merged_bam }

    call Utils.MergeBams as t_53_MergeLongbowFailedReads {
        input:
            bams = flatten([t_17_FilterS2ECCSReads.failed_reads, t_18_FilterS2EReclaimableReads.failed_reads]),
            prefix = SM + "_longbow_failed_array_reads"
    }
    call PB.PBIndex as t_54_PbIndexMergedLongbowFailedReads { input: bam = t_53_MergeLongbowFailedReads.merged_bam }

    call Utils.MergeBams as t_55_MergeAllLongbowAnnotatedReads {
        input:
            bams = flatten([annotated_S2E_ccs_file, annotated_S2E_reclaimable_file]),
            prefix = SM + "_all_longbow_annotated_array_reads"
    }
    call PB.PBIndex as t_56_PbIndexMergedAllLongbowAnnotatedReads { input: bam = t_55_MergeAllLongbowAnnotatedReads.merged_bam }


    # Merge Filtered CCS array elements together:
    call Utils.MergeBams as t_57_MergeCCSArrayElements {
        input:
            bams = t_26_MergeCCSArrayElementShards.merged_bam,
            prefix = SM + "_ccs_array_elements"
    }

    # Merge Filtered CCS Reclaimed array elements together:
    call Utils.MergeBams as t_58_MergeCCSReclaimedArrayElements {
        input:
            bams = t_32_MergeCCSReclaimedArrayElementShards.merged_bam,
            prefix = SM + "_ccs_reclaimed_array_elements"
    }

    # Merge Filtered (and end removed) CCS array elements together:
    call Utils.MergeBams as t_59_MergeCCSArrayElementsNonTruncated {
        input:
            bams = t_28_MergeCCSArrayElementsNonTruncatedShards.merged_bam,
            prefix = SM + "_ccs_array_elements_non_truncated"
    }

    # Merge Filtered (and end removed) CCS Reclaimed array elements together:
    call Utils.MergeBams as t_60_MergeCCSReclaimedArrayElementsNonTruncated {
        input:
            bams = t_34_MergeCCSReclaimedArrayElementsNonTruncatedShards.merged_bam,
            prefix = SM + "_ccs_reclaimed_array_elements_non_truncated"
    }

    call Utils.MergeBams as t_61_MergeAllArrayElementsNonTruncated {
        input:
            bams = flatten([t_28_MergeCCSArrayElementsNonTruncatedShards.merged_bam, t_34_MergeCCSReclaimedArrayElementsNonTruncatedShards.merged_bam]),
            prefix = SM + "_array_elements_non_truncated"
    }

    # Merge Aligned CCS array elements together:
    call Utils.MergeBams as t_62_MergeAlignedCCSArrayElements {
        input:
            bams = t_37_RestoreAnnotationsToGenomeAlignedCCSBam.output_bam,
            prefix = SM + "_ccs_array_elements_aligned"
    }

    # Merge Aligned CCS Reclaimed array elements together:
    call Utils.MergeBams as t_63_MergeAlignedCCSReclaimedArrayElements {
        input:
            bams = t_38_RestoreAnnotationsToGenomeAlignedReclaimedBam.output_bam,
            prefix = SM + "_ccs_reclaimed_array_elements_aligned"
    }

    call Utils.MergeBams as t_64_MergeAllAlignedArrayElementsNonTruncated {
        input:
            bams = flatten([t_37_RestoreAnnotationsToGenomeAlignedCCSBam.output_bam, t_38_RestoreAnnotationsToGenomeAlignedReclaimedBam.output_bam]),
            prefix = SM + "_array_elements_non_truncated_aligned"
    }

    # Merge CCS Barcode Conf files:
    call Utils.MergeFiles as t_65_MergeAllCCSBarcodeConfShards {
        input:
            files_to_merge = t_27_MergeCCSBarcodeConfShards.merged_file,
            merged_file_name = SM + "_ccs_array_element_barcode_confs.txt"
    }

    # Merge CCS Barcode Conf files:
    call Utils.MergeFiles as t_66_MergeAllCCSReclaimedBarcodeConfShards {
        input:
            files_to_merge = t_33_MergeCCSReclaimedBarcodeConfShards.merged_file,
            merged_file_name = SM + "_ccs_reclaimed_array_element_barcode_confs.txt"
    }

    ##########################################################################################################################
    #############################################################
    #
    # FILTER THE READS HERE!!!!!!!!!!!!!!!!111!1!11
    #
    ############################################################
    ##########################################################################################################################

#    	• Post-Alignment filters:
#		○ Remove Flags:
#			§ MQ0
#			§ Supplementary
#			§ Secondary
#			§ Unmapped
#		○ Read length
#			§ Keep: <15Kb
#		○ Spread on the reference
#			§ Some fraciton of ref length
#			§ Or max N (splice gap)
#			§ Overlapping multiple genes - remove it.
#				□ Use funcotate segments or similar
#				□ Whitelist genes with overlapping exons
#				□ Bedtools intersection
#					® THOUGH, occasionally it gives different results
#					® Might not be a problem.
#		○ End soft/hard Clipping
#			§ L or R end of read
#               1000?

    RuntimeAttr filterReadsAttrs = object {
        preemptible_tries: 0
    }

    # Remove unmapped, secondary, supplementary, mq0, length > 15kb, end clips > 1kb
    call Utils.FilterMasSeqReadsWithGatk as t_67_AlignmentFilterForCcsArrayElements {
        input:
            bam_file = t_62_MergeAlignedCCSArrayElements.merged_bam,
            bam_index = t_62_MergeAlignedCCSArrayElements.merged_bai,
            prefix = SM + "_CCS_ArrayElements_Annotated_Aligned_PrimaryOnly",
            runtime_attr_override = filterReadsAttrs
    }

    call Utils.FilterMasSeqReadsWithGatk as t_68_AlignmentFilterForReclaimedArrayElements {
        input:
            bam_file = t_63_MergeAlignedCCSReclaimedArrayElements.merged_bam,
            bam_index = t_63_MergeAlignedCCSReclaimedArrayElements.merged_bai,
            prefix = SM + "_Reclaimed_ArrayElements_Annotated_Aligned_PrimaryOnly",
            runtime_attr_override = filterReadsAttrs
    }

    #########################################################################################################################
    ##########################################################################################################################

    # Mehrtash's suggestion for paper
    # 		• run stringtie2
    #		• filter stringtie2 annotations, get rid of super low TPM annotations
    #		polish stringtie2 annotations (e.g. if a transcript is an extension of a GENCODE transcript, propagate the name, if a “novel” stringtie2 gene overlaps with a previously annotated gene > 95%, propagate the name; otherwise, ignore)
    #		run TALON w/ polished stringtie2 annotations
    #       ignore NIC and NNC TALON transcripts (which should be VERY few), only focus on exact matches and ISM

    # Mehrtash's latest suggestion for paper:
    # Running StringTie2 on our sample w/ GENCODE as reference
    # Running TALON on GENCODE "reads" (we need to cook up a bam file from the GENCODE gtf) with a database initialized with the StringTie2 gtf
    # Running TALON on our sample with a darabase initialized with the StringTie2 gtf

    # Merge all alignments together:
    call Utils.MergeBams as t_69_MergeAllAlignedAndFilteredArrayElements {
        input:
            bams = [t_67_AlignmentFilterForCcsArrayElements.bam, t_68_AlignmentFilterForReclaimedArrayElements.bam],
            prefix = SM + "_all_array_elements_aligned_for_txome_discovery"
    }

    call StringTie2.Quantify as t_70_ST2_Quant {
        input:
            aligned_bam = t_69_MergeAllAlignedAndFilteredArrayElements.merged_bam,
            aligned_bai = t_69_MergeAllAlignedAndFilteredArrayElements.merged_bai,
            gtf = genome_annotation_gtf,
            keep_retained_introns = false,
            prefix = SM + "_StringTie2_Quantify",
    }

    call StringTie2.ExtractTranscriptSequences as t_71_ST2_ExtractTranscriptSequences  {
        input:
            ref_fasta = ref_fasta,
            ref_fasta_fai = ref_fasta_index,
            gtf = t_70_ST2_Quant.st_gtf,
            prefix = SM + "_StringTie2_ExtractTranscriptSequences",
    }

    call StringTie2.CompareTranscriptomes as t_72_ST2_CompareTranscriptomes {
        input:
            guide_gtf = genome_annotation_gtf,
            new_gtf = t_70_ST2_Quant.st_gtf,
            prefix = SM + "_StringTie2_CompareTranscriptome",
    }

    ##########################################################################################################################
    ##########################################################################################################################

    # Now we pad our barcodes and correct them:
    Int num_array_element_shards = 50

    # CCS
    call Utils.ShardReads as t_73_ShardS2ECcsArrayElements {
        input:
            bam = t_67_AlignmentFilterForCcsArrayElements.bam,
            bam_index = t_67_AlignmentFilterForCcsArrayElements.bai,
            prefix = SM + "_ccs_array_elements_subshard",
            num_shards = num_array_element_shards,
    }
    scatter (ccsi in range(length(t_73_ShardS2ECcsArrayElements.shards))) {
        File ccs_array_element_shard = t_73_ShardS2ECcsArrayElements.shards[ccsi]

        call Utils.IndexBam as t_74_IndexCcsArrayElementShard {
            input:
                bam = ccs_array_element_shard
        }

        # Now that we've annotated the reads, we can pad the UMIs by a couple of bases to aid in the deduping:
        call LONGBOW.Pad as t_75_LongbowPadCCSArrayElementUMIs {
            input:
                reads = ccs_array_element_shard,
                model = mas_seq_model,
                tag_to_expand = "ZU",
                padding = ccs_umi_padding,
                prefix = SM + "_ccs_array_elements_aligned_annotated_umi_padded_shard_" + ccsi
        }

        call LONGBOW.Pad as t_76_LongbowPadCCSArrayElementCBCs {
            input:
                reads = t_75_LongbowPadCCSArrayElementUMIs.padded_tag_bam,
                model = mas_seq_model,
                tag_to_expand = "CR",
                new_tag_dest = expanded_cbc_tag,
                padding = ccs_cbc_padding,
                prefix = SM + "_ccs_array_elements_aligned_annotated_cbc_padded_shard_" + ccsi
        }

        # Now we should correct our barcodes based on the whitelist:
        call LONGBOW.Correct as t_77_LongbowCorrectCCSCorrectedArrayElementCBCs {
            input:
                reads = t_76_LongbowPadCCSArrayElementCBCs.padded_tag_bam,
                barcode_allow_list = cell_barcode_whitelist,
                barcode_freq_list = cell_barcode_freq_tsv,
                model = mas_seq_model,
                ccs_lev_dist_threshold = ccs_lev_dist,
                clr_lev_dist_threshold = clr_lev_dist,
                prefix = SM + "_ccs_array_elements_aligned_annotated_padded_cbc_corrected_shard_" + ccsi,
                raw_barcode_tag = expanded_cbc_tag,
                corrected_barcode_tag = "CB",
        }

        call TENX.AdjustUmiSequenceWithAdapterAlignment as t_78_AdjustCCSUMIs {
            input:
                bam = t_77_LongbowCorrectCCSCorrectedArrayElementCBCs.corrected_barcodes_bam,
                short_read_umis = short_read_umis_tsv,
                prefix = SM + "_ccs_array_elements_aligned_annotated_padded_cbc_corrected_UMI_adjusted_shard_" + ccsi,
        }
    }
    # Merge Aligned CCS reads together:
    call Utils.MergeBams as t_79_MergeLongbowPaddedCBCCorrectedCCSArrayElements {
        input:
            bams = t_78_AdjustCCSUMIs.output_bam,
            prefix = SM + "_ccs_array_elements_aligned_annotated_padded_CBC_corrected"
    }
    call Utils.MergeBams as t_80_MergeLongbowPaddedCBCUncorrectableCCSArrayElements {
        input:
            bams = t_77_LongbowCorrectCCSCorrectedArrayElementCBCs.uncorrected_barcodes_bam,
            prefix = SM + "_ccs_array_elements_aligned_annotated_padded_CBC_uncorrectable"
    }

    # RECLAIMED
    call Utils.ShardReads as t_81_ShardS2ECcsReclaimedArrayElements {
        input:
            bam = t_68_AlignmentFilterForReclaimedArrayElements.bam,
            bam_index = t_68_AlignmentFilterForReclaimedArrayElements.bai,
            prefix = SM + "_ccs_reclaimed_array_elements_subshard",
            num_shards = num_array_element_shards,
    }
    scatter (cri in range(length(t_81_ShardS2ECcsReclaimedArrayElements.shards))) {
        File ccs_reclaimed_array_element_shard = t_81_ShardS2ECcsReclaimedArrayElements.shards[cri]

        call Utils.IndexBam as t_82_IndexCcsReclaimedArrayElementShard {
            input:
                bam = ccs_reclaimed_array_element_shard
        }

        # Now that we've annotated the reads, we can pad the CBCs and UMIs by a couple of bases to aid in the deduping:
        call LONGBOW.Pad as t_83_LongbowPadCCSReclaimedArrayElementUMIs {
            input:
                reads = ccs_reclaimed_array_element_shard,
                model = mas_seq_model,
                tag_to_expand = "ZU",
                padding = clr_umi_padding,
                prefix = SM + "_ccs_reclaimed_array_elements_aligned_annotated_umi_padded_shard_" + cri
        }

        call LONGBOW.Pad as t_84_LongbowPadCCSReclaimedArrayElementCBCs {
            input:
                reads = t_83_LongbowPadCCSReclaimedArrayElementUMIs.padded_tag_bam,
                model = mas_seq_model,
                tag_to_expand = "CR",
                new_tag_dest = expanded_cbc_tag,
                padding = clr_cbc_padding,
                prefix = SM + "_ccs_reclaimed_array_elements_aligned_annotated_cbc_padded_shard_" + cri
        }

        # Now we should correct our barcodes based on the whitelist:
        call LONGBOW.Correct as t_85_LongbowCorrectCCSReclaimedArrayElementCBCs {
            input:
                reads = t_84_LongbowPadCCSReclaimedArrayElementCBCs.padded_tag_bam,
                barcode_allow_list = cell_barcode_whitelist,
                barcode_freq_list = cell_barcode_freq_tsv,
                model = mas_seq_model,
                ccs_lev_dist_threshold = ccs_lev_dist,
                clr_lev_dist_threshold = clr_lev_dist,
                prefix = SM + "_ccs_reclaimed_array_elements_aligned_annotated_padded_cbc_corrected_shard_" + cri,
                raw_barcode_tag = expanded_cbc_tag,
                corrected_barcode_tag = "CB",
        }

        call TENX.AdjustUmiSequenceWithAdapterAlignment as t_86_AdjustCCSReclaimedUMIs {
            input:
                bam = t_85_LongbowCorrectCCSReclaimedArrayElementCBCs.corrected_barcodes_bam,
                short_read_umis = short_read_umis_tsv,
                prefix = SM + "_ccs_reclaimed_array_elements_aligned_annotated_padded_cbc_corrected_UMI_adjusted_shard_" + cri,
        }
    }
    # Merge Aligned CCS Reclaimed reads together:
    call Utils.MergeBams as t_87_MergeLongbowPaddedCBCCorrectedCCSReclaimedArrayElements {
        input:
            bams = t_86_AdjustCCSReclaimedUMIs.output_bam,
            prefix = SM + "_ccs_reclaimed_array_elements_aligned_annotated_padded_CBC_corrected"
    }
    call Utils.MergeBams as t_88_MergeLongbowPaddedCBCUncorrectableCCSReclaimedArrayElements {
        input:
            bams = t_85_LongbowCorrectCCSReclaimedArrayElementCBCs.uncorrected_barcodes_bam,
            prefix = SM + "_ccs_reclaimed_array_elements_aligned_annotated_padded_CBC_uncorrectable"
    }

    #################
    # Here we restore the original read names to the bam because we're hashing them with Longbow.segment:

    # Restore original read names to CCS reads:
    call TX_PRE.RestoreOriginalReadNames as t_89_RestoreCcsOriginalReadNames {
        input:
            bam = t_79_MergeLongbowPaddedCBCCorrectedCCSArrayElements.merged_bam,
            prefix =  SM + "_CCS_cbc_annotated_array_elements_padded_original_names"
    }

    # Restore original read names to CLR reads:
    call TX_PRE.RestoreOriginalReadNames as t_90_RestoreClrOriginalReadNames {
        input:
            bam = t_87_MergeLongbowPaddedCBCCorrectedCCSReclaimedArrayElements.merged_bam,
            prefix =  SM + "_CLR_cbc_annotated_array_elements_padded_original_names"
    }

    # Merge Aligned CCS and Reclaimed reads together:
    call Utils.MergeBams as t_91_MergeAllAnnotatedArrayElementsWithOriginalNames {
        input:
            bams = [t_89_RestoreCcsOriginalReadNames.bam_out, t_90_RestoreClrOriginalReadNames.bam_out],
            prefix = SM + "_all_cbc_annotated_array_elements_padded_original_names"
    }

    #################
    # Now we have to split the reads again, process them into gff files, run gffcompare and then aggregate the results in a graph

    # We can actually compare the references without needing to scatter:
    call TX_PRE.GffCompare as t_92_GffCompareStringtie2toGencode {
        input:
            gff_ref = t_70_ST2_Quant.st_gtf,
            gff_query = genome_annotation_gtf,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
    }
    call TX_PRE.GffCompare as t_93_GffCompareGencodetoStringtie2 {
        input:
            gff_ref = genome_annotation_gtf,
            gff_query = t_70_ST2_Quant.st_gtf,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
    }

    # Split by contig:
    call TX_PRE.SplitBamByContig as t_94_SplitArrayElementsByContig {
        input:
            bam = t_91_MergeAllAnnotatedArrayElementsWithOriginalNames.merged_bam,
            prefix = SM + "_all_cbc_annotated_array_elements_padded_original_names"
    }

    # For each contig:
    scatter (i in range(length(t_94_SplitArrayElementsByContig.contig_bams))) {

        File contig_bam = t_94_SplitArrayElementsByContig.contig_bams[i]
        String contig_name = t_94_SplitArrayElementsByContig.contig_names[i]

        # Create a GFF file:
        call TX_PRE.ConvertSplicedBamToGff as t_95_ConvertSplicedBamToGff {
            input:
                bam = contig_bam
        }

        # Compare GFF files:
        call TX_PRE.GffCompare as t_96_GffCompareStringtie2toMasSeqReads {
            input:
                gff_ref = t_70_ST2_Quant.st_gtf,
                gff_query = t_95_ConvertSplicedBamToGff.gff,
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index,
        }

        call TX_PRE.GffCompare as t_97_GffCompareGencodetoMasSeqReads {
            input:
                gff_ref = genome_annotation_gtf,
                gff_query = t_95_ConvertSplicedBamToGff.gff,
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index,
        }

        # Create the comparison graph and tsv files:
        call TX_POST.QuantifyGffComparison as t_98_QuantifyGffComparison {
            input:
                genome_gtf = genome_annotation_gtf,
                st2_gencode_refmap = t_92_GffCompareStringtie2toGencode.refmap,
                st2_gencode_tmap = t_92_GffCompareStringtie2toGencode.tmap,
                st2_read_refmap = t_96_GffCompareStringtie2toMasSeqReads.refmap,
                st2_read_tmap = t_96_GffCompareStringtie2toMasSeqReads.tmap,
                gencode_st2_refmap = t_93_GffCompareGencodetoStringtie2.refmap,
                gencode_st2_tmap = t_93_GffCompareGencodetoStringtie2.tmap,
                gencode_read_refmap = t_97_GffCompareGencodetoMasSeqReads.refmap,
                gencode_read_tmap = t_97_GffCompareGencodetoMasSeqReads.tmap,
                prefix = SM + "_all_cbc_annotated_array_elements_padded_" + contig_name
        }
    }

    # Merge our tx equivalance classes assignments and eq classes:
    call TX_POST.CombineEqClassFiles as t_99_CombineEqClassFiles {
        input:
            gene_eq_class_definitions = t_98_QuantifyGffComparison.gene_eq_class_labels_file,
            gene_assignment_files = t_98_QuantifyGffComparison.gene_assignments_file,
            equivalence_class_definitions = t_98_QuantifyGffComparison.tx_equivalence_class_labels_file,
            equivalence_classes = t_98_QuantifyGffComparison.tx_equivalence_class_file,
            prefix = SM + "_all_cbc_annotated_array_elements_padded"
    }

    ############################################################
    # Quantify Transcripts:
    ##########

    # Use old quant method here as a baseline for comparison:
    call TX_POST.CopyEqClassInfoToTag as t_100_CopyEqClassInfoToTag {
        input:
            bam = t_91_MergeAllAnnotatedArrayElementsWithOriginalNames.merged_bam,
            eq_class_file = t_99_CombineEqClassFiles.combined_tx_eq_class_assignments,
            prefix = SM + "_annotated_array_elements_for_quant_with_gene_names"
    }

    call TX_PRE.CorrectUmisWithSetCover as t_101_CorrectUmisWithSetCover {
        input:
            bam = t_100_CopyEqClassInfoToTag.bam_out,
            prefix = SM + "_annotated_array_elements_for_quant_with_gene_names"
    }

    # Because of how we're doing things, we need to pull out the CCS and CCS Reclaimed reads from the output of the
    # set cover correction:
    call Utils.Bamtools as t_102_GetCcsCorrectedReadsWithCorrectedUmis {
        input:
            bamfile = t_101_CorrectUmisWithSetCover.corrected_umi_reads,
            prefix = SM + "_annotated_array_elements_for_quant_with_gene_names.corrected_umis.CCS",
            cmd = "filter",
            args = '-tag "rq":">=' + min_read_quality + '"',
            runtime_attr_override = disable_preemption_runtime_attrs
        }

    call Utils.Bamtools as t_103_GetCcsReclaimedReadsWithCorrectedUmis {
            input:
                bamfile =t_101_CorrectUmisWithSetCover.corrected_umi_reads,
                prefix = SM + "_annotated_array_elements_for_quant_with_gene_names.corrected_umis.CCS_Reclaimed",
                cmd = "filter",
                args = '-tag "rq":"<' + min_read_quality + '"',
                runtime_attr_override = disable_preemption_runtime_attrs
        }

    call UMI_TOOLS.Run_Group as t_104_UMIToolsGroup {
        input:
            aligned_transcriptome_reads = t_101_CorrectUmisWithSetCover.corrected_umi_reads,
            aligned_transcriptome_reads_index = t_101_CorrectUmisWithSetCover.corrected_umi_reads_index,
            do_per_cell = true,
            prefix = SM + "_annotated_array_elements_with_gene_names_with_umi_tools_group_correction"
    }

    # Create CCS count matrix and anndata:
    call TX_POST.CreateCountMatrixFromAnnotatedBam as t_105_CreateCCSCountMatrixFromAnnotatedBam {
        input:
            annotated_transcriptome_bam = t_102_GetCcsCorrectedReadsWithCorrectedUmis.bam_out,
            tx_equivalence_class_assignments = t_99_CombineEqClassFiles.combined_tx_eq_class_assignments,
            umi_tag = "BX",
            prefix = SM + "_ccs_gene_tx_expression_count_matrix"
    }

    call TX_POST.CreateCountMatrixAnndataFromEquivalenceClasses as t_106_CreateCCSCountMatrixAnndataFromEqClasses {
        input:
            count_matrix_tsv = t_105_CreateCCSCountMatrixFromAnnotatedBam.count_matrix,
            genome_annotation_gtf_file = t_70_ST2_Quant.st_gtf,
            gencode_reference_gtf_file = genome_annotation_gtf,
            overlap_intervals = intervals_of_interest,
            overlap_interval_label = interval_overlap_name,
            tx_equivalence_class_assignments = t_99_CombineEqClassFiles.combined_tx_eq_class_assignments,
            tx_equivalence_class_definitions = t_99_CombineEqClassFiles.combined_tx_eq_class_defs,
            gene_equivalence_class_assignments = t_99_CombineEqClassFiles.combined_gene_eq_class_assignments,
            gene_equivalence_class_definitions = t_99_CombineEqClassFiles.combined_gene_eq_class_defs,
            prefix = SM + "_ccs_gene_tx_expression_count_matrix",

            runtime_attr_override = object {mem_gb: 64}
    }

    # Create CLR count matrix and anndata:
    call TX_POST.CreateCountMatrixFromAnnotatedBam as t_107_CreateCLRCountMatrixFromAnnotatedBam {
        input:
            annotated_transcriptome_bam = t_103_GetCcsReclaimedReadsWithCorrectedUmis.bam_out,
            tx_equivalence_class_assignments = t_99_CombineEqClassFiles.combined_tx_eq_class_assignments,
            umi_tag = "BX",
            prefix = SM + "_clr_gene_tx_expression_count_matrix"
    }

    call TX_POST.CreateCountMatrixAnndataFromEquivalenceClasses as t_108_CreateCLRCountMatrixAnndataFromEqClasses {
        input:
            count_matrix_tsv = t_107_CreateCLRCountMatrixFromAnnotatedBam.count_matrix,
            genome_annotation_gtf_file = t_70_ST2_Quant.st_gtf,
            gencode_reference_gtf_file = genome_annotation_gtf,
            overlap_intervals = intervals_of_interest,
            overlap_interval_label = interval_overlap_name,
            tx_equivalence_class_assignments = t_99_CombineEqClassFiles.combined_tx_eq_class_assignments,
            tx_equivalence_class_definitions = t_99_CombineEqClassFiles.combined_tx_eq_class_defs,
            gene_equivalence_class_assignments = t_99_CombineEqClassFiles.combined_gene_eq_class_assignments,
            gene_equivalence_class_definitions = t_99_CombineEqClassFiles.combined_gene_eq_class_defs,
            prefix = SM + "_clr_gene_tx_expression_count_matrix",

            runtime_attr_override = object {mem_gb: 64}
    }

    # Create overall count matrix and anndata:
    call TX_POST.CreateCountMatrixFromAnnotatedBam as t_109_CreateOverallCountMatrixFromAnnotatedBam {
        input:
            annotated_transcriptome_bam = t_101_CorrectUmisWithSetCover.corrected_umi_reads,
            tx_equivalence_class_assignments = t_99_CombineEqClassFiles.combined_tx_eq_class_assignments,
            umi_tag = "BX",
            prefix = SM + "_overall_gene_tx_expression_count_matrix"
    }

    call TX_POST.CreateCountMatrixAnndataFromEquivalenceClasses as t_110_CreateOverallCountMatrixAnndataFromEqClasses {
        input:
            count_matrix_tsv = t_109_CreateOverallCountMatrixFromAnnotatedBam.count_matrix,
            genome_annotation_gtf_file = t_70_ST2_Quant.st_gtf,
            gencode_reference_gtf_file = genome_annotation_gtf,
            overlap_intervals = intervals_of_interest,
            overlap_interval_label = interval_overlap_name,
            tx_equivalence_class_assignments = t_99_CombineEqClassFiles.combined_tx_eq_class_assignments,
            tx_equivalence_class_definitions = t_99_CombineEqClassFiles.combined_tx_eq_class_defs,
            gene_equivalence_class_assignments = t_99_CombineEqClassFiles.combined_gene_eq_class_assignments,
            gene_equivalence_class_definitions = t_99_CombineEqClassFiles.combined_gene_eq_class_defs,
            prefix = SM + "_overall_gene_tx_expression_count_matrix",

            runtime_attr_override = object {mem_gb: 64}
    }


    #################################################
    #     ___   ____      __  __  __ _____ _____ ____  ___ ____ ____
    #    / _ \ / ___|    / / |  \/  | ____|_   _|  _ \|_ _/ ___/ ___|
    #   | | | | |       / /  | |\/| |  _|   | | | |_) || | |   \___ \
    #   | |_| | |___   / /   | |  | | |___  | | |  _ < | | |___ ___) |
    #    \__\_\\____| /_/    |_|  |_|_____| |_| |_| \_\___\____|____/
    #
    #################################################

    call AM.SamtoolsStats as t_111_BaselineArrayElementStats {
        input:
            bam = t_61_MergeAllArrayElementsNonTruncated.merged_bam
    }

    call AM.SamtoolsStats as t_112_AlignedArrayElementStats {
        input:
            bam = t_64_MergeAllAlignedArrayElementsNonTruncated.merged_bam
    }

    call AM.SamtoolsStats as t_113_AlignedFilteredArrayElementStats {
        input:
            bam = t_69_MergeAllAlignedAndFilteredArrayElements.merged_bam
    }

    call AM.SamtoolsStats as t_114_AlignedAnnotatedArrayElementsForQuantStats {
        input:
            bam = t_101_CorrectUmisWithSetCover.corrected_umi_reads
    }

    call LONGBOW.AggregateCorrectLogStats as t_115_AggregateLongbowCorrectStats {
        input:
            longbow_correct_log_files = flatten([t_85_LongbowCorrectCCSReclaimedArrayElementCBCs.log, t_77_LongbowCorrectCCSCorrectedArrayElementCBCs.log]),
            out_name = SM + "_longbow_correct_stats.txt"
    }

    # Get stats on CCS reads:
    call LONGBOW.Stats as t_116_CCS_longbow_stats {
        input:
            reads = t_39_MergeCCSLongbowAnnotatedArrayReads.merged_bam,
            model = mas_seq_model,
            prefix = SM + "_CCS_Corrected",
    }

    # Get stats on Reclaimable reads:
    call LONGBOW.Stats as t_117_Reclaimable_longbow_stats {
        input:
            reads = t_41_MergeCCSReclaimableLongbowAnnotatedArrayReads.merged_bam,
            model = mas_seq_model,
            prefix = SM + "_CCS_Reclaimable",
    }

    # Get stats on Reclaimed reads:
    call LONGBOW.Stats as t_118_Reclaimed_longbow_stats {
        input:
            reads = t_47_MergeCCSReclaimedArrayReads.merged_bam,
            model = mas_seq_model,
            prefix = SM + "_CCS_Reclaimed",
    }

    # Get stats on All Passing reads (overall stats):
    call LONGBOW.Stats as t_119_Passed_longbow_stats {
        input:
            reads = t_51_MergeLongbowPassedReads.merged_bam,
            model = mas_seq_model,
            prefix = SM + "_All_Longbow_Passed",
    }

    # Get stats on All Failed reads (overall stats):
    call LONGBOW.Stats as t_120_Failed_longbow_stats {
        input:
            reads = t_53_MergeLongbowFailedReads.merged_bam,
            model = mas_seq_model,
            prefix = SM + "_All_Longbow_Failed",
    }

    # Get stats on All reads (overall stats):
    call LONGBOW.Stats as t_121_Overall_longbow_stats {
        input:
            reads = t_55_MergeAllLongbowAnnotatedReads.merged_bam,
            model = mas_seq_model,
            prefix = SM + "_Overall",
    }

    ######################################################################
    #             _____ _             _ _
    #            |  ___(_)_ __   __ _| (_)_______
    #            | |_  | | '_ \ / _` | | |_  / _ \
    #            |  _| | | | | | (_| | | |/ /  __/
    #            |_|   |_|_| |_|\__,_|_|_/___\___|
    #
    ######################################################################

    # NOTE: We key all finalization steps on the static report.
    #       This will prevent incomplete runs from being placed in the output folders.

#    File keyfile = t_110_CreateOverallCountMatrixAnndataFromEqClasses.transcript_gene_count_anndata_h5ad

    # This seems to take longer to get to:
    File keyfile = t_104_UMIToolsGroup.output_tsv

    String base_out_dir = outdir + "/" + DIR + out_dir_suffix + "/" + t_01_WdlExecutionStartTimestamp.timestamp_string
    String stats_out_dir = base_out_dir + "/stats"
    String array_element_dir = base_out_dir + "/annotated_array_elements"
    String intermediate_reads_dir = base_out_dir + "/intermediate_reads"

    String meta_files_dir = base_out_dir + "/meta_files"

    String intermediate_array_reads_dir = intermediate_reads_dir + "/array_reads"
    String intermediate_array_elements_dir = intermediate_reads_dir + "/array_elements"

    String quant_dir = base_out_dir + "/quant"

    ##############################################################################################################
    # Finalize gene / tx assignments:
    call FF.FinalizeToDir as t_122_FinalizeEqClasses {
        input:
            files = [
                t_99_CombineEqClassFiles.combined_gene_eq_class_defs,
                t_99_CombineEqClassFiles.combined_gene_eq_class_assignments,
                t_99_CombineEqClassFiles.combined_tx_eq_class_defs,
                t_99_CombineEqClassFiles.combined_tx_eq_class_assignments,
            ],
            outdir = quant_dir + "/eqivalence_classes",
            keyfile = keyfile
    }

    call FF.FinalizeToDir as t_123_FinalizeUmiToolsOutputs {
        input:
            files = [
                t_104_UMIToolsGroup.output_bam,
                t_104_UMIToolsGroup.output_tsv,
            ],
            outdir = quant_dir + "/UMITools",
            keyfile = keyfile
    }

    # CCS:
    call FF.FinalizeToDir as t_124_FinalizeCCSTxAndGeneAssignments {
        input:
            files = [
                t_105_CreateCCSCountMatrixFromAnnotatedBam.count_matrix,
                t_106_CreateCCSCountMatrixAnndataFromEqClasses.transcript_gene_count_anndata_h5ad,
            ],
            outdir = quant_dir + "/CCS",
            keyfile = keyfile
    }

    call FF.FinalizeToDir as t_125_FinalizeCCSRawQuantPickles {
        input:
            files = t_106_CreateCCSCountMatrixAnndataFromEqClasses.pickles,
            outdir = quant_dir + "/CCS",
            keyfile = keyfile
    }

    # CLR:
    call FF.FinalizeToDir as t_126_FinalizeCLRTxAndGeneAssignments {
        input:
            files = [
                t_107_CreateCLRCountMatrixFromAnnotatedBam.count_matrix,
                t_108_CreateCLRCountMatrixAnndataFromEqClasses.transcript_gene_count_anndata_h5ad,
            ],
            outdir = quant_dir + "/CLR",
            keyfile = keyfile
    }

    call FF.FinalizeToDir as t_127_FinalizeCLRRawQuantPickles {
        input:
            files = t_108_CreateCLRCountMatrixAnndataFromEqClasses.pickles,
            outdir = quant_dir + "/CLR",
            keyfile = keyfile
    }

    # Overall:
    call FF.FinalizeToDir as t_128_FinalizeOverallTxAndGeneAssignments {
        input:
            files = [
                t_109_CreateOverallCountMatrixFromAnnotatedBam.count_matrix,
                t_110_CreateOverallCountMatrixAnndataFromEqClasses.transcript_gene_count_anndata_h5ad,
            ],
            outdir = quant_dir + "/Overall",
            keyfile = keyfile
    }

    call FF.FinalizeToDir as t_129_FinalizeOverallRawQuantPickles {
        input:
            files = t_110_CreateOverallCountMatrixAnndataFromEqClasses.pickles,
            outdir = quant_dir + "/Overall",
            keyfile = keyfile
    }

    call FF.FinalizeToDir as t_130_FinalizeRefAndSt2Comparisons {
        input:
            files = [
                t_92_GffCompareStringtie2toGencode.refmap,
                t_92_GffCompareStringtie2toGencode.tmap,
                t_92_GffCompareStringtie2toGencode.tracking,
                t_92_GffCompareStringtie2toGencode.loci,
                t_92_GffCompareStringtie2toGencode.annotated_gtf,
                t_92_GffCompareStringtie2toGencode.stats,
                t_92_GffCompareStringtie2toGencode.log,

                t_93_GffCompareGencodetoStringtie2.refmap,
                t_93_GffCompareGencodetoStringtie2.tmap,
                t_93_GffCompareGencodetoStringtie2.tracking,
                t_93_GffCompareGencodetoStringtie2.loci,
                t_93_GffCompareGencodetoStringtie2.annotated_gtf,
                t_93_GffCompareGencodetoStringtie2.stats,
                t_93_GffCompareGencodetoStringtie2.log,
            ],
            outdir = quant_dir + "/gencode_and_stringtie2",
            keyfile = keyfile
    }

    # Finalize gene / tx assignment by contig:
    # NOTE: According to the scatter/gather documentation in the WDL spec, this will work correctly
    #       (https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md#scatter--gather)
    scatter (i in range(length(t_94_SplitArrayElementsByContig.contig_bams))) {
        String contig = t_94_SplitArrayElementsByContig.contig_names[i]

        call FF.FinalizeToDir as t_131_FinalizeTxAndGeneAssignmentsByContig {
            input:
                files = [
                    t_96_GffCompareStringtie2toMasSeqReads.refmap[i],
                    t_96_GffCompareStringtie2toMasSeqReads.tmap[i],
                    t_96_GffCompareStringtie2toMasSeqReads.tracking[i],
                    t_96_GffCompareStringtie2toMasSeqReads.loci[i],
                    t_96_GffCompareStringtie2toMasSeqReads.annotated_gtf[i],
                    t_96_GffCompareStringtie2toMasSeqReads.stats[i],
                    t_96_GffCompareStringtie2toMasSeqReads.log[i],

                    t_97_GffCompareGencodetoMasSeqReads.refmap[i],
                    t_97_GffCompareGencodetoMasSeqReads.tmap[i],
                    t_97_GffCompareGencodetoMasSeqReads.tracking[i],
                    t_97_GffCompareGencodetoMasSeqReads.loci[i],
                    t_97_GffCompareGencodetoMasSeqReads.annotated_gtf[i],
                    t_97_GffCompareGencodetoMasSeqReads.stats[i],
                    t_97_GffCompareGencodetoMasSeqReads.log[i],

                    t_98_QuantifyGffComparison.gene_assignments_file[i],
                    t_98_QuantifyGffComparison.tx_equivalence_class_labels_file[i],
                    t_98_QuantifyGffComparison.tx_equivalence_class_file[i],
                    t_98_QuantifyGffComparison.graph_gpickle[i],
                ],
                outdir = quant_dir + "/by_contig/" + contig,
                keyfile = keyfile
        }
    }

    ##############################################################################################################
    # Finalize annotated, aligned array elements:
    call FF.FinalizeToDir as t_132_FinalizeIntermediateAnnotatedArrayElements {
        input:
            files = [
                t_57_MergeCCSArrayElements.merged_bam,
                t_57_MergeCCSArrayElements.merged_bai,
                t_58_MergeCCSReclaimedArrayElements.merged_bam,
                t_58_MergeCCSReclaimedArrayElements.merged_bai,
                t_59_MergeCCSArrayElementsNonTruncated.merged_bam,
                t_59_MergeCCSArrayElementsNonTruncated.merged_bai,
                t_60_MergeCCSReclaimedArrayElementsNonTruncated.merged_bam,
                t_60_MergeCCSReclaimedArrayElementsNonTruncated.merged_bai,

                t_61_MergeAllArrayElementsNonTruncated.merged_bam,
                t_61_MergeAllArrayElementsNonTruncated.merged_bai,

                t_64_MergeAllAlignedArrayElementsNonTruncated.merged_bam,
                t_64_MergeAllAlignedArrayElementsNonTruncated.merged_bai,

                t_69_MergeAllAlignedAndFilteredArrayElements.merged_bam,
                t_69_MergeAllAlignedAndFilteredArrayElements.merged_bai,

                t_79_MergeLongbowPaddedCBCCorrectedCCSArrayElements.merged_bam,
                t_79_MergeLongbowPaddedCBCCorrectedCCSArrayElements.merged_bai,
                t_80_MergeLongbowPaddedCBCUncorrectableCCSArrayElements.merged_bam,
                t_80_MergeLongbowPaddedCBCUncorrectableCCSArrayElements.merged_bai,
                t_87_MergeLongbowPaddedCBCCorrectedCCSReclaimedArrayElements.merged_bam,
                t_87_MergeLongbowPaddedCBCCorrectedCCSReclaimedArrayElements.merged_bai,
                t_88_MergeLongbowPaddedCBCUncorrectableCCSReclaimedArrayElements.merged_bam,
                t_88_MergeLongbowPaddedCBCUncorrectableCCSReclaimedArrayElements.merged_bai,

                t_89_RestoreCcsOriginalReadNames.bam_out,
                t_90_RestoreClrOriginalReadNames.bam_out,

                t_91_MergeAllAnnotatedArrayElementsWithOriginalNames.merged_bam,
                t_91_MergeAllAnnotatedArrayElementsWithOriginalNames.merged_bai,

                t_101_CorrectUmisWithSetCover.uncorrected_umi_reads
            ],
            outdir = intermediate_array_elements_dir,
            keyfile = keyfile
    }

    call FF.FinalizeToDir as t_133_FinalizeAnnotatedArrayElements {
        input:
            files = [
                t_101_CorrectUmisWithSetCover.corrected_umi_reads,
                t_101_CorrectUmisWithSetCover.corrected_umi_reads_index,

                t_102_GetCcsCorrectedReadsWithCorrectedUmis.bam_out,
                t_103_GetCcsReclaimedReadsWithCorrectedUmis.bam_out,

                t_69_MergeAllAlignedAndFilteredArrayElements.merged_bam,
                t_69_MergeAllAlignedAndFilteredArrayElements.merged_bai
            ],
            outdir = array_element_dir,
            keyfile = keyfile
    }

    ##############################################################################################################
    # Finalize meta files:
    call FF.FinalizeToDir as t_134_FinalizeMeta {
        input:
            files = [
                cell_barcode_whitelist,
                t_65_MergeAllCCSBarcodeConfShards.merged_file,
                t_66_MergeAllCCSReclaimedBarcodeConfShards.merged_file
            ],
            outdir = meta_files_dir,
            keyfile = keyfile
    }

    call FF.FinalizeToDir as t_135_FinalizeCCSCBCcorrectionLogsToMeta {
        input:
            files = t_77_LongbowCorrectCCSCorrectedArrayElementCBCs.log,
            outdir = meta_files_dir + "/" + "ccs_cbc_correction_logs",
            keyfile = keyfile
    }

    call FF.FinalizeToDir as t_136_FinalizeCCSRejectedCBCcorrectionLogsToMeta {
        input:
            files = t_85_LongbowCorrectCCSReclaimedArrayElementCBCs.log,
            outdir = meta_files_dir + "/" + "ccs_rejected_cbc_correction_logs",
            keyfile = keyfile
    }

    call FF.FinalizeToDir as t_137_FinalizeCCSUmiAdjustmentLogs {
        input:
            files = t_78_AdjustCCSUMIs.log,
            outdir = meta_files_dir + "/" + "umi_adjustment_logs",
            keyfile = keyfile
    }

    call FF.FinalizeToDir as t_138_FinalizeCCSReclaimedUmiAdjustmentLogs {
        input:
            files = t_86_AdjustCCSReclaimedUMIs.log,
            outdir = meta_files_dir + "/" + "umi_adjustment_logs",
            keyfile = keyfile
    }

    ##############################################################################################################
    # Finalize the discovered transcriptome:
    if ( !is_SIRV_data ) {
        call FF.FinalizeToDir as t_139_FinalizeDiscoveredTranscriptome {
            input:
                files = [
                    t_70_ST2_Quant.st_gtf,
                    t_71_ST2_ExtractTranscriptSequences.transcripts_fa,
                    t_71_ST2_ExtractTranscriptSequences.transcripts_fai,
                    t_71_ST2_ExtractTranscriptSequences.transcripts_dict,
                    t_72_ST2_CompareTranscriptomes.annotated_gtf,
                    t_72_ST2_CompareTranscriptomes.loci,
                    t_72_ST2_CompareTranscriptomes.stats,
                    t_72_ST2_CompareTranscriptomes.tracking,
                    t_72_ST2_CompareTranscriptomes.refmap,
                    t_72_ST2_CompareTranscriptomes.tmap,
                ],
                outdir = base_out_dir + "/discovered_transcriptome",
                keyfile = keyfile
        }
    }
    ##############################################################################################################
    # Finalize the intermediate reads files (from raw CCS corrected reads through split array elements)
    call FF.FinalizeToDir as t_140_FinalizeArrayReads {
        input:
            files = [

                t_39_MergeCCSLongbowAnnotatedArrayReads.merged_bam,
                t_39_MergeCCSLongbowAnnotatedArrayReads.merged_bai,
                t_40_PbIndexMergedCCSLongbowAnnotatedArrayReads.pbindex,

                t_41_MergeCCSReclaimableLongbowAnnotatedArrayReads.merged_bam,
                t_41_MergeCCSReclaimableLongbowAnnotatedArrayReads.merged_bai,
                t_42_PbIndexMergedCCSReclaimableLongbowAnnotatedArrayReads.pbindex,

                t_51_MergeLongbowPassedReads.merged_bam,
                t_51_MergeLongbowPassedReads.merged_bai,
                t_52_PbIndexMergedLongbowPassingReads.pbindex,

                t_53_MergeLongbowFailedReads.merged_bam,
                t_53_MergeLongbowFailedReads.merged_bai,
                t_54_PbIndexMergedLongbowFailedReads.pbindex,

                t_55_MergeAllLongbowAnnotatedReads.merged_bam,
                t_55_MergeAllLongbowAnnotatedReads.merged_bai,
                t_56_PbIndexMergedAllLongbowAnnotatedReads.pbindex,

                t_43_MergeCCSLongbowPassedArrayReads.merged_bam,
                t_43_MergeCCSLongbowPassedArrayReads.merged_bai,
                t_44_PbIndexMergedCCSLongbowPassedArrayReads.pbindex,

                t_45_MergeCCSLongbowFailedArrayReads.merged_bam,
                t_45_MergeCCSLongbowFailedArrayReads.merged_bai,
                t_46_PbIndexMergedCCSLongbowFailedArrayReads.pbindex,

                t_47_MergeCCSReclaimedArrayReads.merged_bam,
                t_47_MergeCCSReclaimedArrayReads.merged_bai,
                t_48_PbIndexMergedCCSReclaimedArrayReads.pbindex,

                t_49_MergeCCSUnreclaimableArrayReads.merged_bam,
                t_49_MergeCCSUnreclaimableArrayReads.merged_bai,
                t_50_PbIndexMergedCCSUnreclaimableReclaimedArrayReads.pbindex,

            ],
            outdir = intermediate_array_reads_dir,
            keyfile = keyfile
    }

    ##############################################################################################################
    # Finalize Stats:

    # Write out completion file so in the future we can be 100% sure that this run was good:
    call FF.FinalizeToDir as t_141_FinalizeHighLevelStats {
        input:
            files = [ t_11_FindCCSReport.ccs_report[0], t_115_AggregateLongbowCorrectStats.stats ],
            outdir = stats_out_dir,
            keyfile = keyfile
    }

    call FF.FinalizeToDir as t_142_FinalizeQuantArrayElementStats {
        input:
            files = [
                t_114_AlignedAnnotatedArrayElementsForQuantStats.raw_stats,
                t_114_AlignedAnnotatedArrayElementsForQuantStats.summary_stats,
                t_114_AlignedAnnotatedArrayElementsForQuantStats.first_frag_qual,
                t_114_AlignedAnnotatedArrayElementsForQuantStats.last_frag_qual,
                t_114_AlignedAnnotatedArrayElementsForQuantStats.first_frag_gc_content,
                t_114_AlignedAnnotatedArrayElementsForQuantStats.last_frag_gc_content,
                t_114_AlignedAnnotatedArrayElementsForQuantStats.acgt_content_per_cycle,
                t_114_AlignedAnnotatedArrayElementsForQuantStats.insert_size,
                t_114_AlignedAnnotatedArrayElementsForQuantStats.read_length_dist,
                t_114_AlignedAnnotatedArrayElementsForQuantStats.indel_distribution,
                t_114_AlignedAnnotatedArrayElementsForQuantStats.indels_per_cycle,
                t_114_AlignedAnnotatedArrayElementsForQuantStats.coverage_distribution,
                t_114_AlignedAnnotatedArrayElementsForQuantStats.gc_depth,
            ],
            outdir = stats_out_dir + "/array_elements_for_quant/",
            keyfile = keyfile
    }

    call FF.FinalizeToDir as t_143_FinalizeTxomeDiscoveryArrayElementStats {
        input:
            files = [
                t_113_AlignedFilteredArrayElementStats.raw_stats,
                t_113_AlignedFilteredArrayElementStats.summary_stats,
                t_113_AlignedFilteredArrayElementStats.first_frag_qual,
                t_113_AlignedFilteredArrayElementStats.last_frag_qual,
                t_113_AlignedFilteredArrayElementStats.first_frag_gc_content,
                t_113_AlignedFilteredArrayElementStats.last_frag_gc_content,
                t_113_AlignedFilteredArrayElementStats.acgt_content_per_cycle,
                t_113_AlignedFilteredArrayElementStats.insert_size,
                t_113_AlignedFilteredArrayElementStats.read_length_dist,
                t_113_AlignedFilteredArrayElementStats.indel_distribution,
                t_113_AlignedFilteredArrayElementStats.indels_per_cycle,
                t_113_AlignedFilteredArrayElementStats.coverage_distribution,
                t_113_AlignedFilteredArrayElementStats.gc_depth,
            ],
            outdir = stats_out_dir + "/array_elements_for_transcriptome_discovery/",
            keyfile = keyfile
    }

    call FF.FinalizeToDir as t_144_FinalizeAlignedArrayElementStats {
        input:
            files = [
                t_112_AlignedArrayElementStats.raw_stats,
                t_112_AlignedArrayElementStats.summary_stats,
                t_112_AlignedArrayElementStats.first_frag_qual,
                t_112_AlignedArrayElementStats.last_frag_qual,
                t_112_AlignedArrayElementStats.first_frag_gc_content,
                t_112_AlignedArrayElementStats.last_frag_gc_content,
                t_112_AlignedArrayElementStats.acgt_content_per_cycle,
                t_112_AlignedArrayElementStats.insert_size,
                t_112_AlignedArrayElementStats.read_length_dist,
                t_112_AlignedArrayElementStats.indel_distribution,
                t_112_AlignedArrayElementStats.indels_per_cycle,
                t_112_AlignedArrayElementStats.coverage_distribution,
                t_112_AlignedArrayElementStats.gc_depth,
            ],
            outdir = stats_out_dir + "/aligned_array_elements/",
            keyfile = keyfile
    }

    call FF.FinalizeToDir as t_145_FinalizeBaselineArrayElementStats {
        input:
            files = [
                t_111_BaselineArrayElementStats.raw_stats,
                t_111_BaselineArrayElementStats.summary_stats,
                t_111_BaselineArrayElementStats.first_frag_qual,
                t_111_BaselineArrayElementStats.last_frag_qual,
                t_111_BaselineArrayElementStats.first_frag_gc_content,
                t_111_BaselineArrayElementStats.last_frag_gc_content,
                t_111_BaselineArrayElementStats.acgt_content_per_cycle,
                t_111_BaselineArrayElementStats.insert_size,
                t_111_BaselineArrayElementStats.read_length_dist,
                t_111_BaselineArrayElementStats.indel_distribution,
                t_111_BaselineArrayElementStats.indels_per_cycle,
                t_111_BaselineArrayElementStats.coverage_distribution,
                t_111_BaselineArrayElementStats.gc_depth,
            ],
            outdir = stats_out_dir + "/baseline_array_elements/",
            keyfile = keyfile
    }

    call FF.FinalizeToDir as t_146_FinalizeCCSLongbowStats {
        input:
            files = [
                t_116_CCS_longbow_stats.summary_stats,
                t_116_CCS_longbow_stats.array_length_counts_plot_png,
                t_116_CCS_longbow_stats.array_length_counts_plot_svg,
                t_116_CCS_longbow_stats.ligation_heatmap_nn_png,
                t_116_CCS_longbow_stats.ligation_heatmap_nn_svg,
                t_116_CCS_longbow_stats.ligation_heatmap_png,
                t_116_CCS_longbow_stats.ligation_heatmap_svg,
                t_116_CCS_longbow_stats.ligation_heatmap_nn_reduced_png,
                t_116_CCS_longbow_stats.ligation_heatmap_nn_reduced_svg,
                t_116_CCS_longbow_stats.ligation_heatmap_reduced_png,
                t_116_CCS_longbow_stats.ligation_heatmap_reduced_svg,
            ],
            outdir = stats_out_dir + "/longbow_stats/CCS_Corrected/",
            keyfile = keyfile
    }
    call FF.FinalizeToDir as t_147_FinalizeReclaimableLongbowStats {
        input:
            files = [
                t_117_Reclaimable_longbow_stats.summary_stats,
                t_117_Reclaimable_longbow_stats.array_length_counts_plot_png,
                t_117_Reclaimable_longbow_stats.array_length_counts_plot_svg,
                t_117_Reclaimable_longbow_stats.ligation_heatmap_nn_png,
                t_117_Reclaimable_longbow_stats.ligation_heatmap_nn_svg,
                t_117_Reclaimable_longbow_stats.ligation_heatmap_png,
                t_117_Reclaimable_longbow_stats.ligation_heatmap_svg,
                t_117_Reclaimable_longbow_stats.ligation_heatmap_nn_reduced_png,
                t_117_Reclaimable_longbow_stats.ligation_heatmap_nn_reduced_svg,
                t_117_Reclaimable_longbow_stats.ligation_heatmap_reduced_png,
                t_117_Reclaimable_longbow_stats.ligation_heatmap_reduced_svg,
            ],
            outdir = stats_out_dir + "/longbow_stats/CCS_Reclaimable/",
            keyfile = keyfile
    }
    call FF.FinalizeToDir as t_148_FinalizeReclaimedLongbowStats {
        input:
            files = [
                t_118_Reclaimed_longbow_stats.summary_stats,
                t_118_Reclaimed_longbow_stats.array_length_counts_plot_png,
                t_118_Reclaimed_longbow_stats.array_length_counts_plot_svg,
                t_118_Reclaimed_longbow_stats.ligation_heatmap_nn_png,
                t_118_Reclaimed_longbow_stats.ligation_heatmap_nn_svg,
                t_118_Reclaimed_longbow_stats.ligation_heatmap_png,
                t_118_Reclaimed_longbow_stats.ligation_heatmap_svg,
                t_118_Reclaimed_longbow_stats.ligation_heatmap_nn_reduced_png,
                t_118_Reclaimed_longbow_stats.ligation_heatmap_nn_reduced_svg,
                t_118_Reclaimed_longbow_stats.ligation_heatmap_reduced_png,
                t_118_Reclaimed_longbow_stats.ligation_heatmap_reduced_svg,
            ],
            outdir = stats_out_dir + "/longbow_stats/CCS_Reclaimed/",
            keyfile = keyfile
    }
    call FF.FinalizeToDir as t_149_FinalizeOverallLongbowStats {
        input:
            files = [
                t_121_Overall_longbow_stats.summary_stats,
                t_121_Overall_longbow_stats.array_length_counts_plot_png,
                t_121_Overall_longbow_stats.array_length_counts_plot_svg,
                t_121_Overall_longbow_stats.ligation_heatmap_nn_png,
                t_121_Overall_longbow_stats.ligation_heatmap_nn_svg,
                t_121_Overall_longbow_stats.ligation_heatmap_png,
                t_121_Overall_longbow_stats.ligation_heatmap_svg,
                t_121_Overall_longbow_stats.ligation_heatmap_nn_reduced_png,
                t_121_Overall_longbow_stats.ligation_heatmap_nn_reduced_svg,
                t_121_Overall_longbow_stats.ligation_heatmap_reduced_png,
                t_121_Overall_longbow_stats.ligation_heatmap_reduced_svg,
            ],
            outdir = stats_out_dir + "/longbow_stats/Overall/",
            keyfile = keyfile
    }
    call FF.FinalizeToDir as t_150_FinalizeAllPassedLongbowStats {
        input:
            files = [
                t_119_Passed_longbow_stats.summary_stats,
                t_119_Passed_longbow_stats.array_length_counts_plot_png,
                t_119_Passed_longbow_stats.array_length_counts_plot_svg,
                t_119_Passed_longbow_stats.ligation_heatmap_nn_png,
                t_119_Passed_longbow_stats.ligation_heatmap_nn_svg,
                t_119_Passed_longbow_stats.ligation_heatmap_png,
                t_119_Passed_longbow_stats.ligation_heatmap_svg,
                t_119_Passed_longbow_stats.ligation_heatmap_nn_reduced_png,
                t_119_Passed_longbow_stats.ligation_heatmap_nn_reduced_svg,
                t_119_Passed_longbow_stats.ligation_heatmap_reduced_png,
                t_119_Passed_longbow_stats.ligation_heatmap_reduced_svg,
            ],
            outdir = stats_out_dir + "/longbow_stats/All_Longbow_Passed/",
            keyfile = keyfile
    }
    call FF.FinalizeToDir as t_151_FinalizeAllPassedLongbowStats {
        input:
            files = [
                t_120_Failed_longbow_stats.summary_stats,
                t_120_Failed_longbow_stats.array_length_counts_plot_png,
                t_120_Failed_longbow_stats.array_length_counts_plot_svg,
                t_120_Failed_longbow_stats.ligation_heatmap_nn_png,
                t_120_Failed_longbow_stats.ligation_heatmap_nn_svg,
                t_120_Failed_longbow_stats.ligation_heatmap_png,
                t_120_Failed_longbow_stats.ligation_heatmap_svg,
                t_120_Failed_longbow_stats.ligation_heatmap_nn_reduced_png,
                t_120_Failed_longbow_stats.ligation_heatmap_nn_reduced_svg,
                t_120_Failed_longbow_stats.ligation_heatmap_reduced_png,
                t_120_Failed_longbow_stats.ligation_heatmap_reduced_svg,
            ],
            outdir = stats_out_dir + "/longbow_stats/All_Longbow_Failed/",
            keyfile = keyfile
    }

    ##############################################################################################################
    # Write out completion file so in the future we can be 100% sure that this run was good:
    call FF.WriteCompletionFile as t_152_WriteCompletionFile {
        input:
            outdir = base_out_dir + "/",
            keyfile = keyfile
    }
}
