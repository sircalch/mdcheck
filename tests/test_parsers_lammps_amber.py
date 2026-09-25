"""
Tests for the LAMMPS log parser and the AMBER cpptraj header handling.
"""
import numpy as np

from mdcheck.parsers import parse_amber_dat, parse_lammps_log, load_timeseries_file

LAMMPS_LOG = """LAMMPS (2 Aug 2023 - Update 3)
units real
Per MPI rank memory allocation (min/avg/max) = 12.3 | 12.3 | 12.3 Mbytes
   Step          Time          Temp          PotEng         Press
         0   0              300.0         -1000.5         1.0
       100   100            301.2         -1001.0         2.5
       200   200            299.1         -1002.3         -1.5
Loop time of 1.234 on 1 procs for 200 steps with 1000 atoms

WARNING: something harmless
   Step          Time          Temp          PotEng         Press
       200   200            299.1         -1002.3         -1.5
       300   300            300.4         -1001.7         0.5
Loop time of 0.6 on 1 procs for 100 steps with 1000 atoms
"""


def test_lammps_log_merges_runs_and_drops_duplicate_step(tmp_path):
    f = tmp_path / "log.lammps"
    f.write_text(LAMMPS_LOG)
    t, data, names, meta = parse_lammps_log(str(f))
    assert names == ["Temp", "PotEng", "Press"]
    assert t.tolist() == [0.0, 100.0, 200.0, 300.0]
    assert np.allclose(data[:, 0], [300.0, 301.2, 299.1, 300.4])
    assert meta["n_runs_merged"] == 2


def test_lammps_log_autodetected_by_generic_loader(tmp_path):
    f = tmp_path / "log.lammps"
    f.write_text(LAMMPS_LOG)
    t, series, meta = load_timeseries_file(str(f))
    assert set(series) == {"Temp", "PotEng", "Press"}
    assert len(t) == 4


def test_amber_cpptraj_header_kept(tmp_path):
    f = tmp_path / "rmsd.dat"
    f.write_text("#Frame  RMSD  Rg\n1  1.2  15.1\n2  1.3  15.0\n3  1.25  15.2\n")
    t, data, names, _ = parse_amber_dat(str(f))
    assert names == ["RMSD", "Rg"]
    assert t.tolist() == [1.0, 2.0, 3.0]
    assert data.shape == (3, 2)
