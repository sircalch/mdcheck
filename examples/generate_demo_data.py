"""
Generates sample GROMACS .xvg and .csv trajectory files for testing and tutorials.
"""

import os
import numpy as np
import pandas as pd


def generate_sample_files(output_dir: str = "sample_data"):
    os.makedirs(output_dir, exist_ok=True)
    n_frames = 1500
    time_coords = np.linspace(0.0, 75.0, n_frames)  # 75 ns
    
    # 1. GROMACS RMSD XVG
    rng1 = np.random.default_rng(101)
    noise1 = np.zeros(n_frames)
    phi = 0.82
    for i in range(1, n_frames):
        noise1[i] = phi * noise1[i-1] + np.sqrt(1 - phi**2) * rng1.normal(0, 0.015)
    rmsd_nm = 0.12 * (1.0 - np.exp(-time_coords / 6.0)) + 0.10 + noise1
    
    xvg_path = os.path.join(output_dir, "rmsd.xvg")
    with open(xvg_path, "w", encoding="utf-8") as f:
        f.write("# GROMACS XVG generated for MDCheck\n")
        f.write('@ title "Backbone RMSD"\n')
        f.write('@ xaxis label "Time (ns)"\n')
        f.write('@ yaxis label "RMSD (nm)"\n')
        f.write('@ s0 legend "Backbone"\n')
        for t, val in zip(time_coords, rmsd_nm):
            f.write(f"{t:10.3f} {val:10.4f}\n")
            
    # 2. Multi-column CSV with 3 Replicas
    def make_replica_series(seed):
        rng = np.random.default_rng(seed)
        noise = np.zeros(n_frames)
        for i in range(1, n_frames):
            noise[i] = 0.80 * noise[i-1] + np.sqrt(1 - 0.80**2) * rng.normal(0, 0.02)
        return 1.45 + 0.05 * np.exp(-time_coords / 5.0) + noise

    df_replicas = pd.DataFrame({
        "Time_ns": time_coords,
        "Replica_1_Rg": make_replica_series(201),
        "Replica_2_Rg": make_replica_series(202),
        "Replica_3_Rg": make_replica_series(203)
    })
    csv_path = os.path.join(output_dir, "multi_replica_gyrate.csv")
    df_replicas.to_csv(csv_path, index=False)
    
    print(f"Sample data generated in {os.path.abspath(output_dir)}:")
    print(f" - {xvg_path}")
    print(f" - {csv_path}")


if __name__ == "__main__":
    generate_sample_files()
