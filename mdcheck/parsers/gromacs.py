"""
Parser for GROMACS XVG data files.
"""

from typing import Tuple, List, Dict, Any, Optional
import os
import numpy as np
import pandas as pd


def parse_xvg_file(filepath: str) -> Tuple[np.ndarray, np.ndarray, List[str], Dict[str, Any]]:
    """
    Parses a GROMACS .xvg file, extracting metadata, titles, column legends, and data arrays.

    Parameters
    ----------
    filepath : str
        Path to the .xvg file.

    Returns
    -------
    time_coords : np.ndarray
        1D array of time values (first column).
    data_matrix : np.ndarray
        2D array of observable timeseries (shape: [n_frames, n_observables]).
    column_names : list of str
        Names/labels of the observable columns.
    metadata : dict
        Extracted XVG header metadata (title, xaxis_label, yaxis_label).
    """
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"File not found: {filepath}")
        
    metadata: Dict[str, Any] = {
        "title": "GROMACS Timeseries",
        "xaxis_label": "Time (ps)",
        "yaxis_label": "Value",
        "legends": []
    }
    
    data_lines = []
    legends_dict: Dict[int, str] = {}
    
    with open(filepath, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line_str = line.strip()
            if not line_str:
                continue
                
            if line_str.startswith("#"):
                continue
                
            if line_str.startswith("@"):
                # Header metadata
                tokens = line_str.split()
                if len(tokens) >= 2:
                    if tokens[1] == "title":
                        metadata["title"] = " ".join(tokens[2:]).strip('"')
                    elif tokens[1] == "xaxis" and len(tokens) >= 4 and tokens[2] == "label":
                        metadata["xaxis_label"] = " ".join(tokens[3:]).strip('"')
                    elif tokens[1] == "yaxis" and len(tokens) >= 4 and tokens[2] == "label":
                        metadata["yaxis_label"] = " ".join(tokens[3:]).strip('"')
                    elif tokens[1].startswith("s") and len(tokens) >= 4 and tokens[2] == "legend":
                        try:
                            s_idx = int(tokens[1][1:])
                            leg_text = " ".join(tokens[3:]).strip('"')
                            legends_dict[s_idx] = leg_text
                        except ValueError:
                            pass
                continue
                
            # Numeric data
            try:
                values = [float(x) for x in line_str.split()]
                if values:
                    data_lines.append(values)
            except ValueError:
                continue
                
    if not data_lines:
        raise ValueError(f"No numeric data found in XVG file: {filepath}")
        
    arr = np.array(data_lines, dtype=np.float64)
    time_coords = arr[:, 0]
    
    if arr.shape[1] > 1:
        data_matrix = arr[:, 1:]
    else:
        data_matrix = arr[:, 0:1]
        
    n_cols = data_matrix.shape[1]
    column_names = []
    for c in range(n_cols):
        if c in legends_dict:
            column_names.append(legends_dict[c])
        else:
            base_label = metadata.get("yaxis_label", "Observable")
            column_names.append(f"{base_label}_{c+1}" if n_cols > 1 else base_label)
            
    metadata["legends"] = column_names
    return time_coords, data_matrix, column_names, metadata
