"""
Generic timeseries file loader with automatic format detection.
"""

from typing import Tuple, List, Dict, Any
import os
import numpy as np
import pandas as pd
from mdcheck.parsers.gromacs import parse_xvg_file
from mdcheck.parsers.amber_namd import parse_openmm_log, parse_amber_dat


def load_timeseries_file(filepath: str) -> Tuple[np.ndarray, Dict[str, np.ndarray], Dict[str, Any]]:
    """
    Universally loads a timeseries file with automated format detection (.xvg, .csv, .tsv, .dat, .txt).

    Parameters
    ----------
    filepath : str
        Path to the timeseries file.

    Returns
    -------
    time_coords : np.ndarray
        1D array of time coordinates.
    series_dict : dict
        Mapping of column name -> 1D numpy array.
    metadata : dict
        Metadata regarding the loaded file.
    """
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"Target file not found: {filepath}")
        
    ext = os.path.splitext(filepath)[1].lower()
    
    # 1. GROMACS XVG
    if ext == ".xvg":
        time_coords, data_matrix, col_names, metadata = parse_xvg_file(filepath)
        series_dict = {col_names[i]: data_matrix[:, i] for i in range(len(col_names))}
        return time_coords, series_dict, metadata
        
    # 2. General CSV / TSV / DAT
    # Try reading first few lines to detect delimiter
    with open(filepath, "r", encoding="utf-8", errors="ignore") as f:
        first_line = f.readline()
        
    if first_line.startswith("#") and ("OpenMM" in first_line or "Potential Energy" in first_line):
        time_coords, data_matrix, col_names, metadata = parse_openmm_log(filepath)
        series_dict = {col_names[i]: data_matrix[:, i] for i in range(len(col_names))}
        return time_coords, series_dict, metadata
        
    # Fallback to pandas robust sniffer
    try:
        df = pd.read_csv(filepath, sep=None, engine='python', comment='#')
    except Exception:
        df = pd.read_csv(filepath, delim_whitespace=True, comment='#')
        
    # Strip column names
    df.columns = [str(c).strip() for c in df.columns]
    
    # Look for time column
    time_col = None
    for c in df.columns:
        if any(keyword in c.lower() for keyword in ["time", "step", "frame", "ps", "ns"]):
            time_col = c
            break
            
    if time_col is not None:
        time_coords = pd.to_numeric(df[time_col], errors='coerce').to_numpy(dtype=np.float64)
        obs_df = df.drop(columns=[time_col])
    else:
        time_coords = np.arange(len(df), dtype=np.float64)
        obs_df = df
        
    series_dict = {}
    for col in obs_df.columns:
        numeric_series = pd.to_numeric(obs_df[col], errors='coerce').dropna().to_numpy(dtype=np.float64)
        if len(numeric_series) > 0:
            series_dict[col] = numeric_series
            
    metadata = {
        "title": os.path.basename(filepath),
        "filepath": os.path.abspath(filepath),
        "xaxis_label": time_col if time_col else "Time",
        "n_frames": len(time_coords)
    }
    
    return time_coords, series_dict, metadata
