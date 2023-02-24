import os
import sys
import logging

from typing import NamedTuple, List, Tuple, Dict, Any

import numpy as np
from scipy import sparse as sp
from operator import itemgetter
from collections import defaultdict

import gffutils
import pysam
import umap
from sklearn.cluster import DBSCAN


logger = logging.getLogger()
logger.setLevel(logging.INFO)
log_info = print


class Interval(NamedTuple):
    feature_type: str
    contig: str
    start: int
    end: int
        
        
def get_feature_interval_list(
        gene_db: gffutils.interface.FeatureDB,
        included_feature_types: List[str]) -> List[Interval]:
    
    # fetch .gtf features
    feature_interval_list = []
    for f in gene_db.all_features(featuretype=included_feature_types):
        interval = Interval(f.featuretype, f.chrom, f.start, f.end)
        feature_interval_list.append(interval)
    feature_interval_list = sorted(feature_interval_list, key=lambda interval: interval.start)
    assert len(set(interval.contig for interval in feature_interval_list)) == 1
    
    return feature_interval_list


def fetch_alignments(
        bam_path: str,
        fetch_contig: str,
        fetch_start: int,
        fetch_end: int,
        gene_body_padding: int,
        min_mapping_quality: int,
        max_mapping_quality: int) -> Tuple[List[pysam.AlignedSegment], Dict[str, Any]]:

    # read records
    query_name_to_alignment_list_map = defaultdict(list)
    bam_file = pysam.AlignmentFile(bam_path)
    for record in bam_file.fetch(fetch_contig, fetch_start, fetch_end):
        query_name_to_alignment_list_map[record.qname].append(record)
    bam_file.close()

    num_total_alignments = sum(list(
        len(alignment_list) for alignment_list in query_name_to_alignment_list_map.values()))
    num_total_queries = len(query_name_to_alignment_list_map.keys())

    # parse alignments
    num_queries_with_out_of_gene_body_alignments = 0
    num_queries_with_single_low_quality_alignment = 0
    num_queries_with_single_high_quality_alignment = 0
    num_queries_with_multiple_low_quality_alignments = 0
    num_queries_with_multiple_high_quality_alignments = 0
    num_queries_with_outstanding_alignment = 0
    
    raw_alignment_list = []
    for query_name, alignment_list in query_name_to_alignment_list_map.items():
        
        map_qualities = [
            record.mapping_quality if record.mapping_quality < max_mapping_quality else -np.inf
            for record in alignment_list]
        
        if len(alignment_list) == 1:
            if map_qualities[0] < min_mapping_quality:
                num_queries_with_single_low_quality_alignment += 1
            else:
                num_queries_with_single_high_quality_alignment += 1
                aln = alignment_list[0]
                if aln.reference_start >= (fetch_start - gene_body_padding) and aln.reference_end <= (fetch_end + gene_body_padding):
                    raw_alignment_list.append(aln)
                else:
                    num_queries_with_out_of_gene_body_alignments += 1
        else:
            if all(map_quality < min_mapping_quality for map_quality in map_qualities):
                num_queries_with_multiple_low_quality_alignments += 1
            else:
                if sum(map_quality >= min_mapping_quality for map_quality in map_qualities) > 1:
                    num_queries_with_multiple_high_quality_alignments += 1
                else:
                    num_queries_with_outstanding_alignment += 1            
                aln = alignment_list[np.argmax(map_qualities)]
                if aln.reference_start >= (fetch_start - gene_body_padding) and aln.reference_end <= (fetch_end + gene_body_padding):
                    raw_alignment_list.append(aln)
                else:
                    num_queries_with_out_of_gene_body_alignments += 1

    stats_dict = {
        'num_total_alignments': num_total_alignments,
        'num_total_queries': num_total_queries,
        'num_queries_with_single_low_quality_alignment': num_queries_with_single_low_quality_alignment,
        'num_queries_with_single_high_quality_alignment': num_queries_with_single_high_quality_alignment,
        'num_queries_with_multiple_low_quality_alignments': num_queries_with_multiple_low_quality_alignments,
        'num_queries_with_multiple_high_quality_alignments': num_queries_with_multiple_high_quality_alignments,
        'num_queries_with_outstanding_alignment': num_queries_with_outstanding_alignment,
        'num_queries_with_out_of_gene_body_alignments': num_queries_with_out_of_gene_body_alignments
    }
    
    return raw_alignment_list, stats_dict


def get_overlap_fraction(ref_start: int, ref_end: int, test_start: int, test_end) -> float:
    overlap_length = max(0, min(ref_end, test_end) - max(ref_start, test_start) + 1)
    return overlap_length / (ref_end - ref_start + 1)


def get_vectorized_repr(
        record: pysam.AlignedSegment,
        vectorization_interval_list: List[Interval],
        abs_overlap_quantization_threshold: float,
        rel_overlap_quantization_threshold: float,
        vectorization_strategy: str) -> np.ndarray:
    """Get a vectorized representation of a `pysam.AlignedSegment` with respect to a given interval list.
    
    .. todo:: this a naive O(N*M) algorithm -- beware!
    """
    assert vectorization_strategy in {'relative', 'absolute', 'quantized'}
    
    abs_overlap_vector = np.zeros((len(vectorization_interval_list),))
    rel_overlap_vector = np.zeros((len(vectorization_interval_list),))
    for i, vectorization_interval in enumerate(vectorization_interval_list):
        if vectorization_interval.contig != record.reference_name:
            continue
        for block in record.get_blocks():
            overlap_fraction = get_overlap_fraction(
                ref_start=vectorization_interval.start,
                ref_end=vectorization_interval.end,
                test_start=block[0],
                test_end=block[1])
            vectorization_interval_length = vectorization_interval.end - vectorization_interval.start + 1
            abs_overlap_vector[i] += vectorization_interval_length * overlap_fraction
            rel_overlap_vector[i] += overlap_fraction
    
    if vectorization_strategy == 'relative':
        return np.clip(rel_overlap_vector, 0., 1.)
    elif vectorization_strategy == 'absolute':
        return abs_overlap_vector
    elif vectorization_strategy == 'quantized':
        quantized_vectorization_repr = (abs_overlap_vector > abs_overlap_quantization_threshold) | \
            (rel_overlap_vector > rel_overlap_quantization_threshold)
        return quantized_vectorization_repr.astype(np.int)
    else:
        raise ValueError('Bad vectorization strategy!')

        
def get_binary_reference_coverage(
        record: pysam.AlignedSegment,
        ref_interval_start: int,
        ref_interval_end: int) -> np.ndarray:
    """Get a binary reference coverage of a `pysam.AlignedSegment` in a given genomic interval.
    """

    coverage = np.zeros((ref_interval_end - ref_interval_start + 1,), dtype=np.bool)
    alignment_block_list = record.get_blocks()
    for block in alignment_block_list:
        block_start = min(max(block[0], ref_interval_start), ref_interval_end) - ref_interval_start
        block_end = max(min(block[1], ref_interval_end), ref_interval_start) - ref_interval_start
        coverage[block_start:(block_end + 1)] = 1
    return coverage


def get_binary_reference_coverage_sparse_matrix(
        alignment_list: List[pysam.AlignedSegment],
        ref_interval_start: int,
        ref_interval_end: int) -> sp.csr_matrix:
    """Returns a sparse coverage matrix from a list of `pysam.AlignedSegment`."""
    
    values = []
    row_indices = []
    col_indices = []
    for row_index, alignment in enumerate(alignment_list):
        binary_coverage = get_binary_reference_coverage(
            alignment, ref_interval_start, ref_interval_end)
        nnz_positions = np.nonzero(binary_coverage)[0].tolist()
        if len(nnz_positions) > 0:
            values += [1] * len(nnz_positions)
            row_indices += [row_index] * len(nnz_positions)
            col_indices += nnz_positions

    coverage_matrix_csr_ni = sp.coo_matrix(
        (values, (row_indices, col_indices)),
        shape=(len(alignment_list), ref_interval_end - ref_interval_start + 1)).tocsr()
    
    return coverage_matrix_csr_ni


def get_vectorized_alignments_matrix(
        alignment_list: List[pysam.AlignedSegment],
        vectorization_interval_list: List[Interval],
        abs_overlap_quantization_threshold: float,
        rel_overlap_quantization_threshold: float,
        vectorization_strategy: str) -> np.ndarray:
    """Returns a dense matrix of vectorized alignments."""
    vectorized_alignments_nv = np.zeros((len(alignment_list), len(vectorization_interval_list)))
    for i_record, record in enumerate(alignment_list):
        vectorized_alignments_nv[i_record, :] = \
            get_vectorized_repr(
            record,
            vectorization_interval_list,
            abs_overlap_quantization_threshold,
            rel_overlap_quantization_threshold,
            vectorization_strategy)
    return vectorized_alignments_nv


def get_feature_to_annotation_track_map(
        feature_interval_list: List[Interval],
        included_feature_types: List[str],
        ref_contig: str,
        ref_interval_start: int,
        ref_interval_end: int) -> Dict[str, np.ndarray]:
    feature_to_annotation_track_map = dict()
    for feature_type in included_feature_types:
        feature_to_annotation_track_map[feature_type] = np.zeros((ref_interval_end - ref_interval_start + 1,), dtype=np.int8)
    for feature in feature_interval_list:
        if feature.contig == ref_contig:
            feature_overlap_start = min(max(feature.start, ref_interval_start), ref_interval_end) - ref_interval_start
            feature_overlap_end = max(min(feature.end, ref_interval_end), ref_interval_start) - ref_interval_start
            feature_to_annotation_track_map[feature.feature_type][feature_overlap_start:(feature_overlap_end + 1)] = 1
    return feature_to_annotation_track_map


def get_junctions_and_counts_from_alignments(
        alignment_list: List[pysam.AlignedSegment],
        alignment_block_min_gap: int,
        consider_read_endpoints_as_junctions: bool) -> Tuple[List[int], List[int]]:
    """Counts junctions for a given list of gapped alignments."""
    junction_counter = defaultdict(int)
    for alignment in alignment_list:
        normalized_blocks = merge_alignment_blocks(alignment.get_blocks(), alignment_block_min_gap)
        if len(normalized_blocks) <= 1:
            continue
        for i, block in enumerate(normalized_blocks):
            if i == 0 and not consider_read_endpoints_as_junctions:  # first block
                junction_counter[block[1]] += 1
            elif i == (len(normalized_blocks) - 1) and not consider_read_endpoints_as_junctions:  # last block
                junction_counter[block[0]] += 1
            else:
                junction_counter[block[0]] += 1
                junction_counter[block[1]] += 1
    sorted_junctions_and_counts = sorted(list(junction_counter.items()), key=lambda x: x[0])
    junctions = list(map(itemgetter(0), sorted_junctions_and_counts))
    counts = list(map(itemgetter(1), sorted_junctions_and_counts))
    return junctions, counts


def merge_alignment_blocks(
        raw_alignment_blocks: List[Tuple[int, int]],
        alignment_block_min_gap: int) -> List[Tuple[int, int]]:
    """Merges adjacent gapped alignment blocks if they are close to one another (e.g. due to insertions)."""
    if len(raw_alignment_blocks) <= 1:
        return raw_alignment_blocks
    merged_blocks = []
    prev_start = raw_alignment_blocks[0][0] 
    prev_end = raw_alignment_blocks[0][1]
    for raw_block in raw_alignment_blocks[1:]:
        curr_start, curr_end = raw_block
        assert curr_start >= prev_end, "The blocks must be non-overlapping and coordinate sorted!"
        if curr_start - prev_end < alignment_block_min_gap:
            prev_end = curr_end
            continue
        else:
            merged_blocks.append((prev_start, prev_end))
            prev_start, prev_end = curr_start, curr_end
    merged_blocks.append((prev_start, prev_end))
    return merged_blocks


def merge_junctions(
        junctions: List[int],
        counts: List[int],
        junction_merge_gap: int) -> Tuple[List[int], List[int]]:
    """Merge nearby junctions and add counts.  The final junction location is the mean of all junctions
    in a junction cluster (weighted by their counts and integer rounded).
    
    .. note:: junctions must be position sorted!
    
    """
    merged_junctions = []
    merged_counts = []
    
    curr_junction_cluster_positions = []
    curr_junction_cluster_counts = []
    last_junction = None
    for junction, count in zip(junctions, counts):
        if last_junction is None:
            curr_junction_cluster_positions.append(junction)
            curr_junction_cluster_counts.append(count)
            last_junction = junction
            continue
        else:
            if junction - last_junction <= junction_merge_gap:
                curr_junction_cluster_positions.append(junction)
                curr_junction_cluster_counts.append(count)
                last_junction = junction
                continue
            else:
                total_counts = sum(curr_junction_cluster_counts)
                mean_merged_junction = int(np.round(sum(j * c for j, c in zip(
                    curr_junction_cluster_positions, curr_junction_cluster_counts)) / total_counts))
                merged_junctions.append(mean_merged_junction)
                merged_counts.append(total_counts)
                curr_junction_cluster_positions = [junction]
                curr_junction_cluster_counts = [count]
                last_junction = junction
                
    if len(curr_junction_cluster_positions) > 0:
        total_counts = sum(curr_junction_cluster_counts)
        mean_merged_junction = int(np.round(sum(j * c for j, c in zip(
            curr_junction_cluster_positions, curr_junction_cluster_counts)) / total_counts))
        merged_junctions.append(mean_merged_junction)
        merged_counts.append(total_counts)
        
    return merged_junctions, merged_counts


def filter_junctions(
        junctions: List[int],
        counts: List[int],
        total_reads: int,
        min_rel_junction_read_support: float,
        min_abs_junction_read_support: int) -> Tuple[List[int], List[int]]:
    """Filters out low frequency junctions."""
    min_abs_junction_read_support = min(
        min_abs_junction_read_support,
        total_reads * min_rel_junction_read_support)
    passing_junctions = [
        junctions[idx] for idx in range(len(junctions))
        if counts[idx] >= min_abs_junction_read_support]
    passing_counts = [
        counts[idx] for idx in range(len(junctions))
        if counts[idx] >= min_abs_junction_read_support]
    return passing_junctions, passing_counts

        
def get_alignment_bounds(
        raw_alignment_list: List[pysam.AlignedSegment]) -> Tuple[int, int]:
    """Returns the reference position lower and upper bounds given a list of alignments."""
    reference_start = min(aln.reference_start for aln in raw_alignment_list)
    reference_end = max(aln.reference_end for aln in raw_alignment_list)
    return reference_start, reference_end


def get_vectorization_interval_list_from_features(
        feature_interval_list: List[Interval],
        gene_id: str,
        min_vectorization_interval_length: int) -> List[Interval]:
    # vectorization intervals
    vectorization_interval_list = []
    feature_endpoints = set()
    for interval in feature_interval_list:
        feature_endpoints.add(interval.start)
        feature_endpoints.add(interval.end)
    feature_junctions = list(sorted(feature_endpoints))
    contig = feature_interval_list[0].contig
    for i in range(len(feature_endpoints) - 1):
        start = feature_junctions[i]
        end = feature_junctions[i + 1]
        if end - start + 1 >= min_vectorization_interval_length:        
            vectorization_interval = Interval(f'vectorization_interval_{gene_id}', contig, start, end)
            vectorization_interval_list.append(vectorization_interval)
    return vectorization_interval_list


def get_vectorization_interval_list_from_junctions(
        junctions: List[int],
        reference_start: int,
        reference_end: int,
        reference_contig: str,
        gene_id: str,
        include_end_segments: bool,
        min_vectorization_interval_length: int) -> List[Interval]:
    vectorization_interval_list = []
    if include_end_segments:
        first_segment = Interval(f'vectorization_interval_{gene_id}', reference_contig, reference_start, junctions[0])
        if (first_segment.end - first_segment.start + 1) >= min_vectorization_interval_length:
            vectorization_interval_list.append(first_segment)
    for i in range(len(junctions) - 1):
        segment = Interval(f'vectorization_interval_{gene_id}', reference_contig, junctions[i], junctions[i + 1])
        if (segment.end - segment.start + 1) >= min_vectorization_interval_length:
            vectorization_interval_list.append(segment)
    if include_end_segments:
        last_segment = Interval(f'vectorization_interval_{gene_id}', reference_contig, junctions[-1], reference_end)
        if (last_segment.end - last_segment.start + 1) >= min_vectorization_interval_length:
            vectorization_interval_list.append(last_segment)
    return vectorization_interval_list


## Experimental ##

def get_flank_augmented_binary_array(binary_array_mv: np.ndarray) -> np.ndarray:
    vector_sz = binary_array_mv.shape[-1]
    first_nnz_idx_m = binary_array_mv.argmax(1)
    last_nnz_idx_m = vector_sz - binary_array_mv[:, ::-1].argmax(1) - 1
    flank_augmented_binary_array_mv = binary_array_mv + \
        (np.arange(vector_sz)[None, :] < first_nnz_idx_m[:, None]).astype(binary_array_mv.dtype) + \
        (np.arange(vector_sz)[None, :] > last_nnz_idx_m[:, None]).astype(binary_array_mv.dtype)
    return flank_augmented_binary_array_mv

def get_support(x: np.ndarray) -> np.ndarray:
    x = x > 0.
    first_nnz_idx = x.argmax()
    last_nnz_idx = x.size - x[::-1].argmax()
    x_support = np.zeros((x.size,), dtype=np.bool)
    x_support[first_nnz_idx:last_nnz_idx] = 1
    return x_support

def mutual_support_l1_distance(x, y, metric_alpha):
    x_support = get_support(x)
    y_support = get_support(y)
    mutual_support = x_support & y_support
    distance_in_mutual_support = np.sum(mutual_support * np.abs(x - y))
    distance_raw = np.sum(np.abs(x - y))
    return metric_alpha * distance_in_mutual_support + (1. - metric_alpha) * distance_raw


def get_degradation_aware_distance(
        x: np.ndarray,
        y: np.ndarray,
        segment_lengths: np.ndarray,
        penalty_5p_per_bp: float,
        penalty_3p_per_bp: float,
        penalty_mutual_support_per_bp: float):
    
    support_x = get_support(x)
    support_y = get_support(y)
    support_mutual = support_x & support_y
    support_mutual_0 = x.argmax()
    support_mutual_1 = x.size - x[::-1].argmax()
    support_non_mutual = support_x ^ support_y
    support_5p = support_non_mutual.copy()
    support_3p = support_non_mutual.copy()
    support_5p[support_mutual_0:] = 0
    support_3p[:support_mutual_1] = 0
    
    x_bp = x * segment_lengths
    y_bp = y * segment_lengths
    l1_dist = np.abs(x_bp - y_bp)
    penalty_mutual = penalty_mutual_support_per_bp * np.sum(l1_dist * support_mutual)
    penalty_5p = penalty_5p_per_bp * np.sum(l1_dist * support_5p)
    penalty_3p = penalty_3p_per_bp * np.sum(l1_dist * support_3p)
    
    return penalty_mutual + penalty_5p + penalty_3p


def get_degradation_aware_similarity(
        x: np.ndarray,
        y: np.ndarray,
        segment_lengths: np.ndarray,
        penalty_5p_per_bp: float,
        penalty_3p_per_bp: float,
        penalty_mutual_support_per_bp: float,
        exp_scale: float):
    return np.exp(-exp_scale * get_degradation_aware_distance(
        x, y, segment_lengths, penalty_5p_per_bp, penalty_3p_per_bp, penalty_mutual_support_per_bp))
