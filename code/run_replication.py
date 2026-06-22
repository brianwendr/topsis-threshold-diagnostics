#!/usr/bin/env python3
"""
Replication script for:
A tolerance-weighted cutoff-validity audit for ordered TOPSIS and VIKOR decisions

The script regenerates the theorem-check examples and the open-data geometry-profile audit used in the manuscript.
The core verification outputs use only the Python standard library. If matplotlib is available, the script also regenerates Fig. 3 from tolerance_sensitivity.csv.
"""
from __future__ import annotations
import csv, json, math, os
from pathlib import Path
from typing import Any, Dict, Iterable, List, Tuple

Vector = List[float]
Matrix = List[List[float]]

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"
LOG = ROOT / "logs" / "trace.log"
DATA = ROOT / "data"
OUT.mkdir(exist_ok=True)
(DATA).mkdir(exist_ok=True)
(ROOT / "logs").mkdir(exist_ok=True)


def column_norms(X: Matrix) -> Vector:
    return [math.sqrt(sum(row[j] ** 2 for row in X)) for j in range(len(X[0]))]


def weighted_normalized_matrix(X: Matrix, w: Vector) -> Matrix:
    norms = column_norms(X)
    return [[w[j] * row[j] / norms[j] if norms[j] else 0.0 for j in range(len(w))] for row in X]


def ideals(V: Matrix, benefit_flags: List[bool]) -> Tuple[Vector, Vector]:
    A_plus, A_minus = [], []
    for j in range(len(V[0])):
        col = [row[j] for row in V]
        if benefit_flags[j]:
            A_plus.append(max(col)); A_minus.append(min(col))
        else:
            A_plus.append(min(col)); A_minus.append(max(col))
    return A_plus, A_minus


def lp_distance(a: Vector, b: Vector, p: float) -> float:
    if p == math.inf:
        return max(abs(x-y) for x, y in zip(a,b))
    return sum(abs(x-y)**p for x, y in zip(a,b))**(1.0/p)


def topsis(X: Matrix, w: Vector, benefit_flags: List[bool], p: float=2.0) -> List[Dict[str, Any]]:
    V = weighted_normalized_matrix(X, w)
    A_plus, A_minus = ideals(V, benefit_flags)
    rows = []
    for q, row in enumerate(V):
        d_plus = lp_distance(row, A_plus, p)
        d_minus = lp_distance(row, A_minus, p)
        denom = d_plus + d_minus
        C = d_minus / denom if denom else 0.0
        ratio = math.inf if d_plus == 0 and d_minus > 0 else (d_minus / d_plus if d_plus else math.inf)
        rows.append({"q": q, "d_plus": d_plus, "d_minus": d_minus, "ratio": ratio, "C": C})
    return rows


def vikor_like_loss(X: Matrix, w: Vector, benefit_flags: List[bool]) -> List[Dict[str, float]]:
    V = weighted_normalized_matrix(X, w)
    A_plus, A_minus = ideals(V, benefit_flags)
    ranges = [abs(A_plus[j]-A_minus[j]) or 1.0 for j in range(len(w))]
    losses = []
    for q,row in enumerate(V):
        gaps = [abs(row[j]-A_plus[j]) / ranges[j] for j in range(len(w))]
        S = sum(w[j]*gaps[j] for j in range(len(w)))
        R = max(w[j]*gaps[j] for j in range(len(w)))
        # normalized VIKOR-style compromise loss with v=0.5
        losses.append({"q": q, "S_loss": S, "R_loss": R, "Q_loss": 0.5*S + 0.5*R})
    return losses


def is_nondecreasing(vals: Iterable[float], tol: float=1e-12) -> bool:
    vals = list(vals)
    return all(vals[i] <= vals[i+1] + tol for i in range(len(vals)-1))


def first_threshold(vals: Iterable[float], level: float) -> int:
    for q, v in enumerate(vals):
        if v >= level:
            return q
    return len(list(vals))


def first_HL_certificate(vals: List[float], level: float):
    for i in range(len(vals)):
        if vals[i] >= level:
            for j in range(i+1, len(vals)):
                if vals[j] < level:
                    return (i,j)
    return None


def cutoff_validity_spectrum(vals: List[float]) -> List[Dict[str, Any]]:
    """Return compressed invalid cutoff intervals (s_j, max_{i<j} s_i].

    A cutoff c in one of these intervals admits a high-low obstruction;
    cutoffs outside their union are thresholdable.
    """
    out = []
    if not vals:
        return out
    prefix_max = vals[0]
    prefix_arg = 0
    for j in range(1, len(vals)):
        if vals[j] < prefix_max - 1e-12:
            out.append({
                "right_index": j,
                "left_witness_index": prefix_arg,
                "invalid_cutoff_open_lower": vals[j],
                "invalid_cutoff_closed_upper": prefix_max,
                "inversion_depth": prefix_max - vals[j],
            })
        if vals[j] > prefix_max + 1e-12:
            prefix_max = vals[j]
            prefix_arg = j
    return out

def max_inversion_depth(vals: List[float]) -> float:
    spec = cutoff_validity_spectrum(vals)
    return max((r["inversion_depth"] for r in spec), default=0.0)


def monotone_repair_radius(vals: List[float]) -> float:
    """Minimum L-infinity distance from vals to a nondecreasing sequence.

    The closed form is half the maximum inversion drop.
    """
    return 0.5 * max_inversion_depth(vals)


def weighted_monotone_repair_radius(vals: List[float], tau: List[float]) -> float:
    """Minimum weighted sup-distance to a nondecreasing sequence.

    The error is max_q |s_q-y_q|/tau_q with tau_q > 0.
    The closed form is max_{i<j} (s_i-s_j)_+/(tau_i+tau_j).
    """
    best = 0.0
    for i in range(len(vals)):
        for j in range(i+1, len(vals)):
            depth = vals[i] - vals[j]
            if depth > best * (tau[i] + tau[j]) + 1e-12:
                best = max(best, depth / (tau[i] + tau[j]))
    return best


def weighted_worst_inversion_pair(vals: List[float], tau: List[float]):
    best = (0.0, None, None)
    for i in range(len(vals)):
        for j in range(i+1, len(vals)):
            score = max(vals[i] - vals[j], 0.0) / (tau[i] + tau[j])
            if score > best[0] + 1e-12:
                best = (score, i, j)
    return best


def worst_inversion_pair(vals: List[float]):
    best = (0.0, None, None)
    for i in range(len(vals)):
        for j in range(i+1, len(vals)):
            depth = vals[i] - vals[j]
            if depth > best[0] + 1e-12:
                best = (depth, i, j)
    return best


def repair_radius_summary(label: str, vals: List[float]) -> Dict[str, Any]:
    depth, i, j = worst_inversion_pair(vals)
    return {
        "case": label,
        "max_inversion_depth": depth,
        "repair_radius": 0.5 * depth,
        "worst_pair": "" if i is None else f"({i},{j})",
        "all_cutoffs_valid": depth <= 1e-12,
    }


def weighted_repair_radius_summary(label: str, vals: List[float], tau: List[float], tau_rule: str="tau_q = 1 + 0.05 q") -> Dict[str, Any]:
    radius, i, j = weighted_worst_inversion_pair(vals, tau)
    return {
        "case": label,
        "weighted_repair_radius": radius,
        "weighted_worst_pair": "" if i is None else f"({i},{j})",
        "tau_rule": tau_rule,
        "all_cutoffs_valid": radius <= 1e-12,
    }


def write_csv(path: Path, rows: List[Dict[str, Any]], fieldnames: List[str] | None = None):
    if not rows:
        return
    if fieldnames is None:
        fieldnames = list(rows[0].keys())
    with path.open('w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for r in rows:
            writer.writerow({k: fmt(r.get(k, "")) for k in fieldnames})


def fmt(x: Any) -> Any:
    if isinstance(x, float):
        if math.isinf(x):
            return "infinity"
        return f"{x:.10f}"
    return x


def synthetic_examples():
    base = [[q+1, 9-q] for q in range(9)]
    modified = [row[:] for row in base]
    modified[3][0] = 2.0
    w = [0.5, 0.5]
    benefits = [True, False]
    base_rows = topsis(base, w, benefits, 2.0)
    mod_rows = topsis(modified, w, benefits, 2.0)
    for r in base_rows:
        r.update({"x_q1": base[r['q']][0], "x_q2": base[r['q']][1], "sign": "H" if r['C'] >= 0.55 else "L"})
    for r in mod_rows:
        q = r['q']
        r.update({"x_q1": modified[q][0], "x_q2": modified[q][1], "sign": "H" if r['C'] >= 0.55 else "L"})
        r["ratio_check"] = "-" if q == 0 else ("ok" if mod_rows[q-1]['ratio'] <= r['ratio'] + 1e-12 else "fail")
        r["comment"] = "SC fails at q=3" if q == 0 else ("d+ rises" if q == 3 else ("T=5" if q == 5 else ""))
    lhl = [{"q":0,"C":0.30,"f_0_55_minus_C":0.25,"pattern":"L"},
           {"q":1,"C":0.60,"f_0_55_minus_C":-0.05,"pattern":"H"},
           {"q":2,"C":0.40,"f_0_55_minus_C":0.15,"pattern":"L"}]
    write_csv(OUT/'table1_base_topsis.csv', base_rows, ["q","x_q1","x_q2","d_plus","d_minus","ratio","C","sign"])
    write_csv(OUT/'table2_modified_topsis.csv', mod_rows, ["q","x_q1","x_q2","d_plus","d_minus","ratio","C","sign","ratio_check","comment"])
    write_csv(OUT/'table3_LHL_obstruction.csv', lhl)
    write_csv(OUT/'cutoff_validity_spectrum_base.csv', cutoff_validity_spectrum([r['C'] for r in base_rows]))
    write_csv(OUT/'cutoff_validity_spectrum_modified.csv', cutoff_validity_spectrum([r['C'] for r in mod_rows]))
    write_csv(OUT/'cutoff_validity_spectrum_LHL.csv', cutoff_validity_spectrum([r['C'] for r in lhl]))
    repair_rows = [
        repair_radius_summary('Synthetic base', [r['C'] for r in base_rows]),
        repair_radius_summary('Modified matrix', [r['C'] for r in mod_rows]),
        repair_radius_summary('L,H,L obstruction', [r['C'] for r in lhl]),
    ]
    write_csv(OUT/'repair_radius_synthetic.csv', repair_rows, ['case','max_inversion_depth','repair_radius','worst_pair','all_cutoffs_valid'])
    tau = [1.0 + 0.05*q for q in range(len(base_rows))]
    # For the low-high-low example, tolerances are derived from triangular fuzzy scores
    # (a_q, s_q, b_q) via tau_q = (b_q-a_q)/2.
    fuzzy_lhl = [
        {"q":0, "a_q":0.25, "s_q":0.30, "b_q":0.35, "tau_q":0.05},
        {"q":1, "a_q":0.54, "s_q":0.60, "b_q":0.66, "tau_q":0.06},
        {"q":2, "a_q":0.33, "s_q":0.40, "b_q":0.47, "tau_q":0.07},
    ]
    write_csv(OUT/'triangular_fuzzy_lhl_tolerances.csv', fuzzy_lhl, ['q','a_q','s_q','b_q','tau_q'])
    tau_lhl = [r['tau_q'] for r in fuzzy_lhl]
    weighted_rows = [
        weighted_repair_radius_summary('Synthetic base', [r['C'] for r in base_rows], tau),
        weighted_repair_radius_summary('Modified matrix', [r['C'] for r in mod_rows], tau),
        weighted_repair_radius_summary('L,H,L obstruction', [r['C'] for r in lhl], tau_lhl, 'tau_q=(b_q-a_q)/2 from triangular fuzzy scores'),
    ]
    write_csv(OUT/'weighted_repair_radius_synthetic.csv', weighted_rows, ['case','weighted_repair_radius','weighted_worst_pair','tau_rule','all_cutoffs_valid'])
    k_rows = [{"level_index":1,"reference_level":0.35,"threshold":first_threshold([r['C'] for r in base_rows],0.35)},
              {"level_index":2,"reference_level":0.70,"threshold":first_threshold([r['C'] for r in base_rows],0.70)}]
    write_csv(OUT/'k_level_thresholds.csv', k_rows)
    # stability margin for base example
    gaps = []
    for q in range(len(base)-1):
        gaps.append(base[q+1][0] - base[q][0])
        gaps.append(base[q][1] - base[q+1][1])
    mu = min(gaps)
    stability = [{"quantity":"mu_X","value":mu,"meaning":"minimum adjacent monotone gap"},
                 {"quantity":"rho_SC","value":mu/2,"meaning":"entrywise stability radius"}]
    write_csv(OUT/'stability_radius_check.csv', stability)
    losses = vikor_like_loss(base, w, benefits)
    write_csv(OUT/'vikor_type_loss_scores.csv', losses)
    checks = {
        "base_threshold_c0_0_55": first_threshold([r['C'] for r in base_rows], 0.55),
        "modified_threshold_c0_0_55": first_threshold([r['C'] for r in mod_rows], 0.55),
        "base_C_nondecreasing": is_nondecreasing([r['C'] for r in base_rows]),
        "base_invalid_cutoff_spectrum_empty": len(cutoff_validity_spectrum([r['C'] for r in base_rows])) == 0,
        "modified_invalid_cutoff_spectrum_empty": len(cutoff_validity_spectrum([r['C'] for r in mod_rows])) == 0,
        "LHL_max_inversion_depth": max_inversion_depth([r['C'] for r in lhl]),
        "LHL_monotone_repair_radius": monotone_repair_radius([r['C'] for r in lhl]),
        "LHL_weighted_repair_radius_from_triangular_fuzzy_tolerances": weighted_monotone_repair_radius([r['C'] for r in lhl], [0.05, 0.06, 0.07]),
        "modified_C_nondecreasing": is_nondecreasing([r['C'] for r in mod_rows]),
        "base_ratio_nondecreasing": is_nondecreasing([r['ratio'] for r in base_rows]),
        "modified_ratio_nondecreasing": is_nondecreasing([r['ratio'] for r in mod_rows]),
        "base_matrix_SC_margin_mu": mu,
        "base_SC_stability_radius": mu/2,
        "LHL_certificate_at_c0_0_55": first_HL_certificate([r['C'] for r in lhl], 0.55),
        "K_level_thresholds_for_0_35_0_70": [r['threshold'] for r in k_rows],
        "vikor_Q_loss_nonincreasing": all(losses[i]['Q_loss'] >= losses[i+1]['Q_loss'] - 1e-12 for i in range(len(losses)-1))
    }
    return checks


def uci_geometry_profiles():
    # The 12 geometry profiles in the public ENB2012 design grid. The full UCI file expands
    # these by orientation, glazing area, and glazing distribution to 768 observations.
    # Ordered here from lower compactness to higher compactness, so q represents a higher
    # geometry-efficiency grade.
    profiles = [
        (0.62,808.5,367.5,220.5,3.5),
        (0.64,784.0,343.0,220.5,3.5),
        (0.66,759.5,318.5,220.5,3.5),
        (0.69,735.0,294.0,220.5,3.5),
        (0.71,710.5,269.5,220.5,3.5),
        (0.74,686.0,245.0,220.5,3.5),
        (0.76,661.5,416.5,122.5,7.0),
        (0.79,637.0,343.0,147.0,7.0),
        (0.82,612.5,318.5,147.0,7.0),
        (0.86,588.0,294.0,147.0,7.0),
        (0.90,563.5,318.5,122.5,7.0),
        (0.98,514.5,294.0,110.25,7.0),
    ]
    rows=[]
    for q,p in enumerate(profiles):
        rows.append({"q":q,"relative_compactness":p[0],"surface_area":p[1],"wall_area":p[2],"roof_area":p[3],"overall_height":p[4]})
    write_csv(DATA/'uci_geometry_profiles.csv', rows)
    return profiles


def open_data_case():
    X = [list(row) for row in uci_geometry_profiles()]
    w = [0.2]*5
    benefits = [True, False, False, False, False]
    summary = []
    all_scores = []
    for p_name, p in [("p1",1.0),("p2",2.0),("pinf",math.inf)]:
        rows = topsis(X,w,benefits,p)
        C = [r['C'] for r in rows]
        cert = first_HL_certificate(C, 0.55)
        threshold = first_threshold(C,0.55)
        accepted_count = sum(1 for x in C if x >= 0.55)
        if accepted_count == 0:
            interp = "valid by empty-set sentinel; no accepted alternative at c0=0.55"
        elif cert is None and not is_nondecreasing(C):
            interp = "valid at c0=0.55 but not globally monotone"
        elif is_nondecreasing(C):
            interp = "thresholdable for all levels"
        else:
            interp = "not thresholdable at c0=0.55"
        summary.append({
            "norm": p_name,
            "C_nondecreasing": is_nondecreasing(C),
            "threshold_at_0_55": threshold,
            "accepted_count_at_0_55": accepted_count,
            "threshold_meaning": "empty acceptance set" if threshold == len(C) else ("only q=%d accepted" % threshold if accepted_count == 1 else "suffix starts at q=%d" % threshold),
            "HL_certificate_at_0_55": "" if cert is None else f"({cert[0]},{cert[1]})",
            "interpretation": interp
        })
        for r in rows:
            all_scores.append({"norm":p_name,"q":r['q'],"C":r['C'],"ratio":r['ratio'],"sign_c0_0_55":"H" if r['C'] >= 0.55 else "L"})
    write_csv(OUT/'open_data_geometry_scores.csv', all_scores, ["norm","q","C","ratio","sign_c0_0_55"])
    wide_scores = []
    for q in range(12):
        row = {"q": q}
        for p_name in ["p1", "p2", "pinf"]:
            row[p_name] = next(r["C"] for r in all_scores if r["norm"] == p_name and r["q"] == q)
        wide_scores.append(row)
    write_csv(OUT/'norm_specific_open_data_scores_wide.csv', wide_scores, ["q","p1","p2","pinf"])
    # compressed invalid cutoff spectrum for each norm
    spec_rows = []
    for p_name in ["p1","p2","pinf"]:
        vals = [r["C"] for r in all_scores if r["norm"] == p_name]
        for row in cutoff_validity_spectrum(vals):
            row = dict(row)
            row["norm"] = p_name
            spec_rows.append(row)
    write_csv(OUT/'open_data_cutoff_validity_spectrum.csv', spec_rows, ["norm","right_index","left_witness_index","invalid_cutoff_open_lower","invalid_cutoff_closed_upper","inversion_depth"])
    write_csv(OUT/'open_data_diagnostic_summary.csv', summary, ["norm","C_nondecreasing","threshold_at_0_55","accepted_count_at_0_55","threshold_meaning","HL_certificate_at_0_55","interpretation"])
    repair_rows = []
    score_by_norm = {}
    for p_name in ["p1", "p2", "pinf"]:
        vals = [r["C"] for r in all_scores if r["norm"] == p_name]
        score_by_norm[p_name] = vals
        repair_rows.append(repair_radius_summary(f"UCI geometry profile {p_name}", vals))
    write_csv(OUT/'repair_radius_open_data.csv', repair_rows, ['case','max_inversion_depth','repair_radius','worst_pair','all_cutoffs_valid'])

    # Full fuzzy-score instantiation for all twelve UCI geometry profiles.
    # We take the p=2 TOPSIS score as the modal score and use the envelope across
    # p=1, p=2, and p=infinity as the triangular support. Thus
    # tau_q=(b_q-a_q)/2 is fully determined by the norm-sensitivity envelope.
    fuzzy_rows = []
    tau = []
    for q in range(12):
        vals = [score_by_norm["p1"][q], score_by_norm["p2"][q], score_by_norm["pinf"][q]]
        a = min(vals)
        s = score_by_norm["p2"][q]
        b = max(vals)
        t = max((b-a)/2.0, 1e-9)
        tau.append(t)
        fuzzy_rows.append({"q": q, "a_q": a, "s_q_modal_p2": s, "b_q": b, "tau_q": t})
    write_csv(OUT/'triangular_fuzzy_open_data_scores.csv', fuzzy_rows, ['q','a_q','s_q_modal_p2','b_q','tau_q'])

    weighted_rows = []
    for p_name in ["p1", "p2", "pinf"]:
        vals = score_by_norm[p_name]
        weighted_rows.append(weighted_repair_radius_summary(
            f"UCI geometry profile {p_name}",
            vals,
            tau,
            "tau_q=(b_q-a_q)/2 from the 12-profile triangular fuzzy score envelope"
        ))
    write_csv(OUT/'weighted_repair_radius_open_data.csv', weighted_rows, ['case','weighted_repair_radius','weighted_worst_pair','tau_rule','all_cutoffs_valid'])

    # Tolerance-scale scenarios based on the fuzzy-score tolerances.
    sens_rows = []
    for p_name in ["p1", "p2", "pinf"]:
        vals = score_by_norm[p_name]
        row = {"norm": p_name}
        for scale in [0.5, 1.0, 2.0]:
            row[f"R_tau_scale_{scale}"] = weighted_monotone_repair_radius(vals, [scale*x for x in tau])
        sens_rows.append(row)
    write_csv(OUT/'tolerance_sensitivity.csv', sens_rows, ['norm','R_tau_scale_0.5','R_tau_scale_1.0','R_tau_scale_2.0'])
    maybe_plot_weighted_repair_sensitivity(sens_rows)

    scenario_rows = [
        {"scenario":"Uniform tolerance", "tolerance_source":"same support width for all alternatives", "interpretation":"equal score reliability", "expected_effect":"baseline repair comparison"},
        {"scenario":"Fuzzy norm-envelope tolerance", "tolerance_source":"triangular score envelope over p=1,2,infinity", "interpretation":"heterogeneous norm-sensitivity spread", "expected_effect":"main full fuzzy-score audit"},
        {"scenario":"Conservative tolerance", "tolerance_source":"half of fuzzy norm-envelope support", "interpretation":"stricter score reliability requirement", "expected_effect":"larger R_tau values"},
        {"scenario":"Liberal tolerance", "tolerance_source":"twice the fuzzy norm-envelope support", "interpretation":"more permissive score reliability band", "expected_effect":"smaller R_tau values"},
    ]
    write_csv(OUT/'tolerance_scenario_summary.csv', scenario_rows, ['scenario','tolerance_source','interpretation','expected_effect'])
    return summary


def maybe_plot_weighted_repair_sensitivity(rows: List[Dict[str, Any]]) -> None:
    """Regenerate Fig. 3 from the same R_tau values reported in tolerance_sensitivity.csv.

    The plot is optional: all numerical checks remain available even when matplotlib
    is not installed.  The y-axis is the tolerance-weighted repair radius R_tau,
    not the uniform repair radius.
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return

    fig_dir = ROOT / "figures"
    fig_dir.mkdir(exist_ok=True)
    scales = [0.5, 1.0, 2.0]
    labels = {"p1": "p=1", "p2": "p=2", "pinf": "p=∞"}
    linestyles = {"p1": "-", "p2": "--", "pinf": ":"}
    markers = {"p1": "o", "p2": "s", "pinf": "^"}

    fig, ax = plt.subplots(figsize=(6.6, 4.1), dpi=300)
    for row in rows:
        norm = row["norm"]
        y = [float(row[f"R_tau_scale_{scale}"]) for scale in scales]
        ax.plot(scales, y, linestyle=linestyles.get(norm, "-"),
                marker=markers.get(norm, "o"), linewidth=1.8,
                markersize=5.2, label=labels.get(norm, norm))
    ymax = max(float(row[f"R_tau_scale_{scale}"]) for row in rows for scale in scales)
    ax.set_xlabel("Tolerance scale γ", fontsize=10)
    ax.set_ylabel(r"Tolerance-weighted repair radius $R_\tau$", fontsize=10)
    ax.set_xticks(scales)
    ax.set_ylim(0, ymax * 1.12)
    ax.grid(True, which="major", linestyle=":", linewidth=0.7, alpha=0.7)
    ax.legend(frameon=False, fontsize=9, loc="upper right")
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    ax.tick_params(labelsize=9)
    fig.tight_layout(pad=1.1)
    fig.savefig(fig_dir / "figure3_weighted_repair_sensitivity.png", dpi=300)
    plt.close(fig)


def main():
    checks = synthetic_examples()
    open_summary = open_data_case()
    checks["open_data_geometry_summary"] = open_summary
    # repair-radius outputs are available in repair_radius_synthetic.csv and repair_radius_open_data.csv
    with (OUT/'theorem_and_open_data_checks.json').open('w', encoding='utf-8') as f:
        json.dump(checks, f, indent=2)
    with LOG.open('w', encoding='utf-8') as f:
        f.write('Replication completed successfully.\n')
        for k,v in checks.items():
            f.write(f'{k}: {v}\n')

if __name__ == '__main__':
    main()
