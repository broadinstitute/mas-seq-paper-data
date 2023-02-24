import numpy as np
import matplotlib.pylab as plt
from typing import List, Tuple, Dict
from itertools import groupby


def segment_binary_array(arr: np.ndarray) -> List[Tuple]:
    segments = []
    for g in groupby(enumerate(arr), key=lambda t: t[1]):
        value = g[0]
        if int(value) == 1:
            t_list = list(map(lambda t: t[0], g[1]))
            segments.append((min(t_list), max(t_list)))
    return segments


def annotate_axis_with_tracks(
        feature_to_annotation_track_map: Dict[str, np.ndarray],
        included_features: List[str],
        feature_to_color_map: Dict[str, str],
        ax: plt.Axes,
        included_position_mask: np.ndarray,
        alpha_region: float = 0.5,
        alpha_line: float = 1.0,
        lw: float = 1.):
    for feature in included_features:
        segments = segment_binary_array(feature_to_annotation_track_map[feature][included_position_mask])
        color = feature_to_color_map[feature]
        for segment in segments:
            if segment[1] - segment[0] > 0:
                ax.axvspan(segment[0], segment[1], color=color, alpha=alpha_region)
            else:
                ax.axvline(segment[0], color=color, alpha=alpha_line, lw=lw)
                
