"""
Parsers for AMBER and OpenMM / NAMD timeseries output files.
"""

from typing import Tuple, List, Dict, Any
import os
import numpy as np
import pandas as pd


def parse_amber_dat(filepath: str) -> Tuple[np.ndarray, np.ndarray, List[str], Dict[str, Any]]:
    """
    Parses AMBER cpptraj output dat files.

    Parameters
    ----------
    filepath : str
        Path to the AMBER .dat file.

    Returns
    -------
    time_coords, data_matrix, column_names, metadata
    """
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"File not found: {filepath}")
        
    df = pd.read_csv(filepath, delim_whitespace=True, comment='#')
    
    time_coords = df.iloc[:, 0].to_numpy(dtype=np.float64)
    data_matrix = df.iloc[:, 1:].to_numpy(dtype=np.float64)
    column_names = list(df.columns[1:])
    
    metadata = {
        "title": "AMBER Data",
        "xaxis_label": str(df.columns[0]),
        "yaxis_label": "Value",
        "legends": column_names
    }
    
    return time_coords, data_matrix, column_names, metadata


def parse_openmm_log(filepath: str) -> Tuple[np.ndarray, np.ndarray, List[str], Dict[str, Any]]:
    """
    Parses OpenMM StateDataReporter CSV/log output files.

    Parameters
    ----------
    filepath : str
        Path to the OpenMM CSV/log file.

    Returns
    -------
    time_coords, data_matrix, column_names, metadata
    """
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"File not found: {filepath}")
        
    # OpenMM CSVs usually have header with quotes like #"Step","Time (ps)","Potential Energy (kJ/mole)"
    df = pd.read_csv(filepath)
    df.columns = [c.strip('#" ') for c in df.columns]
    
    # Identify time or step column
    time_col = None
    for c in df.columns:
        if "time" in c.lower() or "step" in c.lower() or "frame" in c.lower():
            time_col = c
            break
            
    if time_col is not None:
        time_coords = df[time_col].to_numpy(dtype=np.float64)
        obs_df = df.drop(columns=[time_col])
    else:
        time_coords = np.arange(len(df), dtype=np.float64)
        obs_df = df
        
    data_matrix = obs_df.to_numpy(dtype=np.float64)
    column_names = list(obs_df.columns)
    
    metadata = {
        "title": "OpenMM State Data",
        "xaxis_label": time_col if time_col else "Frame",
        "yaxis_label": "Value",
        "legends": column_names
    }
    
    return time_coords, data_matrix, column_names, metadata
