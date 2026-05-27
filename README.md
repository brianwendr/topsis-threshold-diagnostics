# Replication Package

Manuscript: **A Three-Level Diagnostic Hierarchy for TOPSIS Threshold Representability on Ordered Alternatives**

This package regenerates the numerical evidence reported in Section 7 of the manuscript and checks the main computational claims used by Corollaries 1-3 and Proposition 1.

## Contents

- `code/run_replication.py`: standard-library Python script that regenerates all outputs.
- `outputs/table1_base_topsis.csv`: base TOPSIS example, corresponding to manuscript Table 1.
- `outputs/table2_modified_topsis.csv`: modified TOPSIS example, corresponding to manuscript Table 2.
- `outputs/table3_LHL_obstruction.csv`: L,H,L obstruction example, corresponding to manuscript Table 3.
- `outputs/k_level_thresholds.csv`: ordered interval thresholds for reference levels 0.35 and 0.70.
- `outputs/stability_radius_check.csv`: matrix-level stability radius check.
- `outputs/alignment_score_check.csv`: Euclidean alignment score identity check.
- `outputs/vikor_type_loss_scores.csv`: VIKOR-type ideal-point loss score check.
- `outputs/theorem_checks.json`: compact Boolean and threshold verification trace.
- `logs/trace.log`: text trace of the replication run.
- `replication_manifest.csv`: file-by-file manifest.

## Requirements

Python 3.10 or later is recommended. The script uses only the Python standard library. No external packages are required.

## Reproduction command

From the package root, run:

```bash
python code/run_replication.py
```

The command overwrites the CSV/JSON outputs in `outputs/` and writes a trace to `logs/trace.log`.

## Expected headline checks

- Base example threshold at `c0 = 0.55`: `T = 5`.
- Modified example threshold at `c0 = 0.55`: `T = 5`.
- Base and modified closeness sequences are nondecreasing.
- Modified matrix violates the matrix-level monotonicity condition but still satisfies the distance-ratio condition.
- K-level thresholds for `c1 = 0.35`, `c2 = 0.70`: `T1 = 3`, `T2 = 6`.
- Base stability margin: `mu(X) = 1`, so the entrywise SC stability radius is `mu(X)/2 = 0.5`.
- Alignment-score expansion and direct squared-distance calculation agree up to machine precision.

## Notes

The package intentionally uses a small deterministic instance so that reviewers can inspect every number manually if desired. The calculations follow the manuscript definitions: vector normalization, weighted normalized scores, TOPSIS ideal points, p-norm distances, and the closeness coefficient.
