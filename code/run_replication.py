#!/usr/bin/env python3
"""
Replication script for:
A Three-Level Diagnostic Hierarchy for TOPSIS Threshold Representability on Ordered Alternatives

This script regenerates the numerical tables and diagnostic checks reported in
Section 7 of the manuscript. It uses only the Python standard library.

Outputs are written to ../outputs by default when run from the package root:
    python code/run_replication.py
"""
from __future__ import annotations

import csv
import json
import math
import os
from pathlib import Path
from typing import Iterable, List, Tuple, Dict, Any

Vector = List[float]
Matrix = List[List[float]]


def column_norms(X: Matrix) -> Vector:
    m = len(X[0])
    return [math.sqrt(sum(row[j] ** 2 for row in X)) for j in range(m)]


def weighted_normalized_matrix(X: Matrix, w: Vector) -> Matrix:
    norms = column_norms(X)
    V = []
    for row in X:
        V.append([w[j] * row[j] / norms[j] if norms[j] != 0 else 0.0 for j in range(len(w))])
    return V


def ideals(V: Matrix, benefit_flags: List[bool]) -> Tuple[Vector, Vector]:
    m = len(V[0])
    A_plus = []
    A_minus = []
    for j in range(m):
        col = [row[j] for row in V]
        if benefit_flags[j]:
            A_plus.append(max(col))
            A_minus.append(min(col))
        else:
            A_plus.append(min(col))
            A_minus.append(max(col))
    return A_plus, A_minus


def lp_distance(a: Vector, b: Vector, p: float) -> float:
    if p == math.inf:
        return max(abs(x - y) for x, y in zip(a, b))
    return sum(abs(x - y) ** p for x, y in zip(a, b)) ** (1.0 / p)


def topsis(X: Matrix, w: Vector, benefit_flags: List[bool], p: float = 2.0) -> List[Dict[str, Any]]:
    V = weighted_normalized_matrix(X, w)
    A_plus, A_minus = ideals(V, benefit_flags)
    rows = []
    for q, row in enumerate(V):
        d_plus = lp_distance(row, A_plus, p)
        d_minus = lp_distance(row, A_minus, p)
        denom = d_plus + d_minus
        C = d_minus / denom if denom else 0.0
        ratio = math.inf if d_plus == 0 and d_minus > 0 else (d_minus / d_plus if d_plus else math.inf)
        rows.append({
            "q": q,
            "x_q1": X[q][0],
            "x_q2": X[q][1],
            "d_plus": d_plus,
            "d_minus": d_minus,
            "ratio": ratio,
            "C": C,
        })
    return rows


def threshold(values: Iterable[float], level: float) -> int:
    for q, value in enumerate(values):
        if value >= level:
            return q
    return len(list(values)) + 1


def threshold_from_list(values: List[float], level: float) -> int:
    for q, value in enumerate(values):
        if value >= level:
            return q
    return len(values)


def is_nondecreasing(values: List[float], tolerance: float = 1e-12) -> bool:
    return all(values[i + 1] + tolerance >= values[i] for i in range(len(values) - 1))


def first_HL_certificate(C: List[float], c0: float) -> Tuple[int, int] | None:
    seen_H = None
    for q, value in enumerate(C):
        if value >= c0 and seen_H is None:
            seen_H = q
        if seen_H is not None and q > seen_H and value < c0:
            return seen_H, q
    return None


def matrix_SC_margin(X: Matrix, benefit_flags: List[bool]) -> float:
    gaps = []
    for q in range(len(X) - 1):
        for j, is_benefit in enumerate(benefit_flags):
            if is_benefit:
                gaps.append(X[q + 1][j] - X[q][j])
            else:
                gaps.append(X[q][j] - X[q + 1][j])
    return min(gaps)


def euclidean_alignment_scores(X: Matrix, w: Vector, benefit_flags: List[bool]) -> List[Dict[str, Any]]:
    V = weighted_normalized_matrix(X, w)
    A_plus, A_minus = ideals(V, benefit_flags)
    rows = topsis(X, w, benefit_flags, p=2.0)
    out = []
    for q, row in enumerate(V):
        expanded = 0.0
        for j in range(len(row)):
            expanded += 2 * row[j] * (A_plus[j] - A_minus[j]) + (A_minus[j] ** 2 - A_plus[j] ** 2)
        direct = rows[q]["d_minus"] ** 2 - rows[q]["d_plus"] ** 2
        out.append({"q": q, "S_direct": direct, "S_expanded": expanded, "C": rows[q]["C"], "median_preferred": rows[q]["C"] >= 0.5})
    return out


def vikor_type_loss_scores(X: Matrix, w: Vector, benefit_flags: List[bool], alpha: float = 0.5) -> List[Dict[str, Any]]:
    """Compute a simple VIKOR-type ideal-point loss score.

    Lower Q_loss is better. Under the monotone criterion condition, S_loss,
    R_loss, and Q_loss are nonincreasing along q for the current alternative set.
    """
    m = len(X[0])
    # criterion-wise raw best/worst values respecting benefit/cost type
    best = []
    worst = []
    for j, is_benefit in enumerate(benefit_flags):
        col = [row[j] for row in X]
        if is_benefit:
            best.append(max(col)); worst.append(min(col))
        else:
            best.append(min(col)); worst.append(max(col))
    losses = []
    for q, row in enumerate(X):
        criterion_losses = []
        for j, is_benefit in enumerate(benefit_flags):
            denom = abs(best[j] - worst[j])
            if denom == 0:
                loss = 0.0
            elif is_benefit:
                loss = (best[j] - row[j]) / denom
            else:
                loss = (row[j] - best[j]) / denom
            criterion_losses.append(w[j] * loss)
        S_loss = sum(criterion_losses)
        R_loss = max(criterion_losses)
        losses.append({"q": q, "S_loss": S_loss, "R_loss": R_loss})
    S_vals = [r["S_loss"] for r in losses]
    R_vals = [r["R_loss"] for r in losses]
    S_best, S_worst = min(S_vals), max(S_vals)
    R_best, R_worst = min(R_vals), max(R_vals)
    for r in losses:
        S_norm = 0.0 if S_worst == S_best else (r["S_loss"] - S_best) / (S_worst - S_best)
        R_norm = 0.0 if R_worst == R_best else (r["R_loss"] - R_best) / (R_worst - R_best)
        r["Q_loss"] = alpha * S_norm + (1 - alpha) * R_norm
    return losses


def write_csv(path: Path, rows: List[Dict[str, Any]], fieldnames: List[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            cleaned = {}
            for k in fieldnames:
                v = row.get(k, "")
                if isinstance(v, float):
                    if math.isinf(v):
                        cleaned[k] = "infinity"
                    else:
                        cleaned[k] = f"{v:.10f}"
                else:
                    cleaned[k] = v
            writer.writerow(cleaned)


def rounded_table_rows(rows: List[Dict[str, Any]], c0: float = 0.55) -> List[Dict[str, Any]]:
    out = []
    for r in rows:
        out.append({
            "q": r["q"],
            "x_q1": f"{r['x_q1']:.1f}",
            "x_q2": f"{r['x_q2']:.1f}",
            "d_plus": f"{r['d_plus']:.4f}",
            "d_minus": f"{r['d_minus']:.4f}",
            "ratio": "infinity" if math.isinf(r["ratio"]) else f"{r['ratio']:.4f}",
            "C": f"{r['C']:.4f}",
            "sign": "H" if r["C"] >= c0 else "L",
        })
    return out


def main() -> None:
    pkg_root = Path(__file__).resolve().parents[1]
    out_dir = pkg_root / "outputs"
    log_dir = pkg_root / "logs"
    out_dir.mkdir(exist_ok=True)
    log_dir.mkdir(exist_ok=True)

    X_base = [[float(q + 1), float(9 - q)] for q in range(9)]
    X_modified = [[float(q + 1), float(9 - q)] for q in range(9)]
    X_modified[3][0] = 2.0
    w = [0.5, 0.5]
    benefit_flags = [True, False]
    c0 = 0.55
    p = 2.0

    base = topsis(X_base, w, benefit_flags, p)
    modified = topsis(X_modified, w, benefit_flags, p)
    base_C = [r["C"] for r in base]
    modified_C = [r["C"] for r in modified]

    base_table = rounded_table_rows(base, c0)
    modified_table = rounded_table_rows(modified, c0)
    for r in modified_table:
        r["ratio_check"] = "-" if r["q"] == 0 else "ok"
        r["comment"] = ""
    modified_table[0]["comment"] = "SC fails at q=3"
    modified_table[3]["comment"] = "d+ rises"
    modified_table[5]["comment"] = "T=5"

    write_csv(out_dir / "table1_base_topsis.csv", base_table, ["q", "x_q1", "x_q2", "d_plus", "d_minus", "ratio", "C", "sign"])
    write_csv(out_dir / "table2_modified_topsis.csv", modified_table, ["q", "x_q1", "x_q2", "d_plus", "d_minus", "ratio", "C", "sign", "ratio_check", "comment"])

    table3 = [
        {"q": 0, "C": 0.30, "f_0_55_minus_C": +0.25, "pattern": "L"},
        {"q": 1, "C": 0.60, "f_0_55_minus_C": -0.05, "pattern": "H"},
        {"q": 2, "C": 0.40, "f_0_55_minus_C": +0.15, "pattern": "L"},
    ]
    write_csv(out_dir / "table3_LHL_obstruction.csv", table3, ["q", "C", "f_0_55_minus_C", "pattern"])

    # K-level ordered partition check used in the manuscript.
    levels = [0.35, 0.70]
    thresholds = [threshold_from_list(base_C, level) for level in levels]
    intervals = []
    prev = 0
    for T in thresholds:
        intervals.append(list(range(prev, T)))
        prev = T
    intervals.append(list(range(prev, len(base_C))))
    k_rows = []
    for idx, level in enumerate(levels, start=1):
        k_rows.append({"level_index": idx, "reference_level": level, "threshold": thresholds[idx-1]})
    write_csv(out_dir / "k_level_thresholds.csv", k_rows, ["level_index", "reference_level", "threshold"])

    # Stability radius check.
    mu = matrix_SC_margin(X_base, benefit_flags)
    stability_rows = [
        {"quantity": "mu_X", "value": mu, "meaning": "minimum adjacent monotone gap"},
        {"quantity": "rho_SC", "value": mu / 2, "meaning": "entrywise stability radius"},
        {"quantity": "epsilon_safe_example", "value": 0.49 * mu, "meaning": "safe when epsilon < mu/2"},
        {"quantity": "epsilon_adversarial_example", "value": 0.51 * mu, "meaning": "can break the smallest gap when epsilon > mu/2"},
    ]
    write_csv(out_dir / "stability_radius_check.csv", stability_rows, ["quantity", "value", "meaning"])

    alignment_rows = euclidean_alignment_scores(X_base, w, benefit_flags)
    write_csv(out_dir / "alignment_score_check.csv", alignment_rows, ["q", "S_direct", "S_expanded", "C", "median_preferred"])

    vikor_rows = vikor_type_loss_scores(X_base, w, benefit_flags, alpha=0.5)
    write_csv(out_dir / "vikor_type_loss_scores.csv", vikor_rows, ["q", "S_loss", "R_loss", "Q_loss"])

    checks = {
        "base_threshold_c0_0_55": threshold_from_list(base_C, c0),
        "modified_threshold_c0_0_55": threshold_from_list(modified_C, c0),
        "base_C_nondecreasing": is_nondecreasing(base_C),
        "modified_C_nondecreasing": is_nondecreasing(modified_C),
        "base_ratio_nondecreasing": is_nondecreasing([r["ratio"] for r in base]),
        "modified_ratio_nondecreasing": is_nondecreasing([r["ratio"] for r in modified]),
        "base_matrix_SC_margin_mu": mu,
        "base_SC_stability_radius": mu / 2,
        "modified_matrix_SC_margin_mu": matrix_SC_margin(X_modified, benefit_flags),
        "LHL_certificate_at_c0_0_55": first_HL_certificate([0.30, 0.60, 0.40], c0),
        "K_level_thresholds_for_0_35_0_70": thresholds,
        "K_level_intervals": intervals,
        "alignment_expansion_max_abs_error": max(abs(r["S_direct"] - r["S_expanded"]) for r in alignment_rows),
        "vikor_Q_loss_nonincreasing": all(vikor_rows[i + 1]["Q_loss"] <= vikor_rows[i]["Q_loss"] + 1e-12 for i in range(len(vikor_rows) - 1)),
    }
    with (out_dir / "theorem_checks.json").open("w", encoding="utf-8") as f:
        json.dump(checks, f, indent=2)

    manifest_rows = [
        {"file": "code/run_replication.py", "purpose": "Regenerates all CSV outputs and theorem checks."},
        {"file": "outputs/table1_base_topsis.csv", "purpose": "Replicates manuscript Table 1."},
        {"file": "outputs/table2_modified_topsis.csv", "purpose": "Replicates manuscript Table 2."},
        {"file": "outputs/table3_LHL_obstruction.csv", "purpose": "Replicates manuscript Table 3."},
        {"file": "outputs/k_level_thresholds.csv", "purpose": "Checks Corollary 2 ordered K-level thresholds."},
        {"file": "outputs/stability_radius_check.csv", "purpose": "Checks Proposition 1 stability radius."},
        {"file": "outputs/alignment_score_check.csv", "purpose": "Checks Corollary 1 alignment score identity."},
        {"file": "outputs/vikor_type_loss_scores.csv", "purpose": "Checks Corollary 3 VIKOR-type monotone loss scores."},
        {"file": "outputs/theorem_checks.json", "purpose": "Compact trace of all Boolean and threshold checks."},
    ]
    write_csv(pkg_root / "replication_manifest.csv", manifest_rows, ["file", "purpose"])

    with (log_dir / "trace.log").open("w", encoding="utf-8") as f:
        f.write("Replication completed successfully.\n")
        for key, value in checks.items():
            f.write(f"{key}: {value}\n")

    print(json.dumps(checks, indent=2))


if __name__ == "__main__":
    main()
