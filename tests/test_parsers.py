"""
Tests for XVG and CSV/TSV timeseries parsers.
"""

import os
import tempfile
import numpy as np
import pytest
from mdcheck.parsers.gromacs import parse_xvg_file
from mdcheck.parsers.generic import load_timeseries_file


def test_parse_xvg_file():
    xvg_content = """# This is a GROMACS comment
@ title "RMSD of Backbone"
@ xaxis label "Time (ps)"
@ yaxis label "RMSD (nm)"
@ s0 legend "Backbone"
@ s1 legend "Binding Site"
0.0   0.12   0.08
10.0  0.15   0.09
20.0  0.18   0.11
30.0  0.20   0.12
"""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".xvg", delete=False) as f:
        f.write(xvg_content)
        f_path = f.name
        
    try:
        t, data, col_names, meta = parse_xvg_file(f_path)
        assert len(t) == 4
        assert t[1] == 10.0
        assert data.shape == (4, 2)
        assert col_names == ["Backbone", "Binding Site"]
        assert meta["title"] == "RMSD of Backbone"
        
        # Test generic loader
        t2, s_dict, meta2 = load_timeseries_file(f_path)
        assert "Backbone" in s_dict
        assert "Binding Site" in s_dict
        assert len(s_dict["Backbone"]) == 4
    finally:
        if os.path.exists(f_path):
            os.remove(f_path)


def test_parse_csv_file():
    csv_content = """Time_ns,RMSD_A,Rg_nm
0.0,1.2,1.5
1.0,1.4,1.52
2.0,1.5,1.51
3.0,1.6,1.49
"""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
        f.write(csv_content)
        f_path = f.name
        
    try:
        t, s_dict, meta = load_timeseries_file(f_path)
        assert len(t) == 4
        assert "RMSD_A" in s_dict
        assert "Rg_nm" in s_dict
        assert s_dict["RMSD_A"][0] == 1.2
    finally:
        if os.path.exists(f_path):
            os.remove(f_path)
