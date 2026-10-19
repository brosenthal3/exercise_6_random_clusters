# Modified from https://github.com/sta426hs2026/omni-pipelines-py

import os
import argparse
import json
import pandas as pd
import scanpy as sc
import numpy as np

# Parse command line arguments
parser = argparse.ArgumentParser(description="Random clustering as basline to benchmarking single-cell analysis pipelines")
parser.add_argument(
    "--output_dir", "-o", type=str, required=True, help="Output directory"
)
parser.add_argument("--name", "-n", type=str, required=True, help="Dataset name")
parser.add_argument(
    "--data.h5ad", dest="data_h5ad", type=str, required=True, help="Input h5ad file"
)
parser.add_argument(
    "--data.clusters_truth_num",
    dest="clusters_truth_num_path",
    type=str,
    required=True,
    help="input true number of clusters path",
)
parser.add_argument("--seed", type=int, default=42, help="Random seed for reproducibility")
parser.add_argument("--d_cluster", type=int, default=0, help="Additional clusters to add")

args, _ = parser.parse_known_args()

# time object to store time involved (in seconds) in each step
timings: dict[str, float | None] = {
    "gpu_load": None,
    "find_mit_gene": None,
    "filter": None,
    "normalization": None,
    "hvg": None,
    "scaling": None,
    "pca": None,
    "t_sne": None,
    "umap": None,
    "louvain": None,
    "leiden": None,
}

clustering_info: dict[str, dict[str, float | int | None]] = {
    "resolutions": {
        "louvain": None,
        "leiden": None,
    },
    "num_runs": {
        "louvain": None,
        "leiden": None,
    },
}

# Read data
adata = sc.read_h5ad(args.data_h5ad)
adata.var_names_make_unique()
with open(args.clusters_truth_num_path, "r") as f:
    clusters_truth_num = int(f.read())
n_cluster = clusters_truth_num + args.d_cluster

# numpy rng
rng = np.random.default_rng(seed=args.seed)
adata.obs["louvain"] = rng.integers(low=0, high=n_cluster, size=adata.n_obs).astype(str)
adata.obs["leiden"] = rng.integers(low=0, high=n_cluster, size=adata.n_obs).astype(str)

# Save timings as JSON
with open(os.path.join(args.output_dir, f"{args.name}.timings.json"), "w") as f:
    json.dump(timings, f, indent=2)

# Save clustering metadata as JSON
with open(os.path.join(args.output_dir, f"{args.name}.clustering_info.json"), "w") as f:
    json.dump(clustering_info, f, indent=2)

# Save cluster assignments as TSV
adata.obs[["louvain", "leiden"]].reset_index(names="cell_id").to_csv(
    os.path.join(args.output_dir, f"{args.name}.clusters.tsv"),
    sep="\t",
    index=False,
)

# Possibly need to also save HVG and PCA, unsure how omnibenchmark handles this...