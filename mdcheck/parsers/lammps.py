"""
Parser for LAMMPS log files (thermo output blocks).
"""

from typing import Tuple, List, Dict, Any
import os
import numpy as np


def _is_number(token: str) -> bool:
    try:
        float(token)
        return True
    except ValueError:
        return False


def parse_lammps_log(filepath: str) -> Tuple[np.ndarray, np.ndarray, List[str], Dict[str, Any]]:
    """
    Parses thermodynamic output from a LAMMPS log file (e.g. log.lammps).

    Every `run` prints a header line starting with "Step" followed by numeric rows and ends
    with "Loop time of ...". Consecutive runs with the same thermo columns are concatenated;
    the first row of a continuation run repeats the last step of the previous one and is dropped.

    Parameters
    ----------
    filepath : str
        Path to the LAMMPS log file.

    Returns
    -------
    time_coords, data_matrix, column_names, metadata
        time_coords is the "Time" column when present (LAMMPS time units), otherwise "Step".
    """
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"File not found: {filepath}")

    blocks: List[Tuple[List[str], List[List[float]]]] = []
    header: List[str] = []
    rows: List[List[float]] = []
    in_block = False

    with open(filepath, "r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            tokens = line.split()
            if not tokens:
                continue
            if tokens[0] == "Step" and not in_block:
                header, rows, in_block = tokens, [], True
                continue
            if in_block:
                if len(tokens) == len(header) and all(_is_number(t) for t in tokens):
                    rows.append([float(t) for t in tokens])
                    continue
                # "Loop time of ...", warnings or anything else closes the block
                if rows:
                    blocks.append((header, rows))
                in_block = False
        if in_block and rows:
            blocks.append((header, rows))

    if not blocks:
        raise ValueError(f"No thermo output blocks (header starting with 'Step') found in {filepath}")

    ref_header = blocks[0][0]
    merged: List[List[float]] = []
    n_runs = 0
    for hdr, blk in blocks:
        if hdr != ref_header:
            continue
        if merged and blk and blk[0][0] == merged[-1][0]:
            blk = blk[1:]
        merged.extend(blk)
        n_runs += 1

    data = np.asarray(merged, dtype=np.float64)
    time_name = "Time" if "Time" in ref_header else "Step"
    time_idx = ref_header.index(time_name)
    obs_idx = [i for i, name in enumerate(ref_header) if name not in ("Step", "Time")]

    time_coords = data[:, time_idx]
    data_matrix = data[:, obs_idx]
    column_names = [ref_header[i] for i in obs_idx]

    metadata = {
        "title": "LAMMPS thermo output",
        "xaxis_label": time_name,
        "yaxis_label": "Value",
        "legends": column_names,
        "n_runs_merged": n_runs,
        "n_blocks_found": len(blocks),
    }
    return time_coords, data_matrix, column_names, metadata
