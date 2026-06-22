# Online Resource 1: Replication package for Soft Computing submission

Manuscript: A tolerance-weighted cutoff-validity audit for ordered TOPSIS and VIKOR decisions

## Contents

- `code/run_replication.py`: Python script used to regenerate all outputs.
- `data/uci_geometry_profiles.csv`: ordered geometry-profile abstraction derived from the UCI Energy Efficiency design grid for the open-data illustration.
- `outputs/`: generated CSV tables, theorem checks, cutoff-validity spectra, stability-radius checks, uniform and tolerance-weighted repair-radius outputs, triangular fuzzy tolerance inputs, and tolerance-sensitivity outputs.
- `figures/`: Figure 1 audit pipeline, Figure 2 open-data TOPSIS score curves, and Figure 3 weighted repair sensitivity.
- `logs/trace.log`: representative trace log.

## How to run

```bash
python code/run_replication.py
```

The script writes generated CSV files into `outputs/` and logs the main validation checks.

## Main outputs used by the manuscript

- `table1_base_topsis.csv`
- `table2_modified_topsis.csv`
- `table3_LHL_obstruction.csv`
- `cutoff_validity_spectrum_LHL.csv`
- `repair_radius_synthetic.csv`
- `weighted_repair_radius_synthetic.csv`
- `triangular_fuzzy_lhl_tolerances.csv`
- `triangular_fuzzy_open_data_scores.csv`
- `open_data_geometry_scores.csv`
- `norm_specific_open_data_scores_wide.csv`
- `open_data_cutoff_validity_spectrum.csv`
- `repair_radius_open_data.csv`
- `weighted_repair_radius_open_data.csv`
- `tolerance_sensitivity.csv`
- `tolerance_scenario_summary.csv`

## Data note

No proprietary data are used. The open-data illustration uses a public, non-human energy-efficiency design dataset through an ordered geometry-profile abstraction.

## Fuzzy-score tolerance note

For the low-high-low theorem check and the twelve-profile UCI illustration, tolerance vectors are instantiated from triangular fuzzy scores (a_q, s_q, b_q) by tau_q=(b_q-a_q)/2. In the UCI illustration, the p=2 TOPSIS score is used as the modal value and the support is the norm envelope across p=1, p=2, and p=infinity. This gives a full fuzzy-score source for the tolerance-weighted repair radius used in the manuscript, not an arbitrary linear tolerance sequence.


## Norm-sensitivity envelope rationale

In the UCI geometry-profile illustration, the triangular fuzzy score is constructed from the envelope of TOPSIS scores obtained under p=1, p=2, and p=infinity. This treats norm-sensitivity as representational uncertainty: when an analyst has not committed to a single distance metric, the spread across plausible norms is used as a data-driven proxy for score ambiguity, analogous to using multiple expert estimates to define the support of a triangular fuzzy number.

## Weighted repair worst-pair trace

The file `outputs/weighted_repair_radius_open_data.csv` reports the pair attaining the global tolerance-weighted repair radius for each norm. In the current replication run, the worst pair is `(4,7)` for p=1, p=2, and p=infinity. These traces support the R_tau values reported in the open-data audit table of the manuscript.

## Threshold sentinel convention

For an ordered set Q={0,...,N}, the threshold T=N+1 denotes the empty acceptance set. In the open-data p=infinity illustration, T=12 therefore means that no alternative reaches c0=0.55.


## Manuscript-level norm trace

The file `outputs/norm_specific_open_data_scores_wide.csv` gives the p=1, p=2, and p=infinity TOPSIS closeness sequences in wide format for all twelve UCI geometry profiles. This is the manuscript-level trace table used to verify that the tolerance-weighted worst pair `(4,7)` in `outputs/weighted_repair_radius_open_data.csv` is not inferred from the fuzzy envelope alone but from each norm-specific score sequence.


## Version note

Version v8.4 aligns with the revised manuscript by emphasizing the post-score audit role, fuzzy-score tolerance instantiation, norm-sensitivity envelope, and the regenerated CSV/trace files used in Tables 3-8 and Figs. 1-3. It also corrects Fig. 3 so that its y-axis uses the same tolerance-weighted repair radius R_tau values reported in Table 7 and tolerance_sensitivity.csv.


Figure regeneration note. The numerical CSV outputs are produced with the Python standard library. If matplotlib is available, run_replication.py also regenerates figures/figure3_weighted_repair_sensitivity.png directly from tolerance_sensitivity.csv.
