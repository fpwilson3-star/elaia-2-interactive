"""Precompute statistics for the ELAIA-2 interactive paper.

Reads the deidentified dataset and writes a JSON file per table/figure into
`data/`. Re-run whenever analysis logic changes.

Usage:  python scripts/precompute.py
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.stats.contingency_tables import StratifiedTable

ROOT = Path(__file__).resolve().parent.parent
CSV_PATH = ROOT / "files" / "deidentified_elaia2_data_missing_NaN3 (1).csv"
DATA_DIR = ROOT / "data"

# Drug-MOI panels shown inside Fig 1 and Fig 2 — fixed regardless of cohort.
MOI_PANELS = {
    "overall": ("Overall", lambda d: pd.Series(True, index=d.index)),
    "nsaid": ("NSAID-exposed", lambda d: d["on_nsaid"] == 1),
    "raasi": ("RAASi-exposed", lambda d: d["on_acearb"] == 1),
    "ppi": ("PPI-exposed", lambda d: d["on_ppi"] == 1),
}


def get_subgroup_filters(df: pd.DataFrame) -> dict:
    """Pre-specified cohort filters for interactive restriction.

    Matches the subgroups used in the published Fig 3 plus the three MOI strata.
    Returns an ordered dict: key → {label, group, filter}.
    """
    creat = df["baseline_creat"].dropna()
    q1, q2 = np.percentile(creat, [33.333, 66.667])
    return {
        "overall":    {"label": "Overall (full cohort)",        "group": "Full cohort",
                       "filter": lambda d: pd.Series(True, index=d.index)},
        "nsaid":      {"label": "NSAID-exposed",                "group": "MOI exposure",
                       "filter": lambda d: d["on_nsaid"] == 1},
        "raasi":      {"label": "RAASi-exposed",                "group": "MOI exposure",
                       "filter": lambda d: d["on_acearb"] == 1},
        "ppi":        {"label": "PPI-exposed",                  "group": "MOI exposure",
                       "filter": lambda d: d["on_ppi"] == 1},
        "non_icu":    {"label": "Non-ICU at randomization",     "group": "ICU status",
                       "filter": lambda d: d["icu_at_rand"] == 0},
        "icu":        {"label": "ICU at randomization",         "group": "ICU status",
                       "filter": lambda d: d["icu_at_rand"] == 1},
        "surgical":   {"label": "Surgical admission",           "group": "Admitting service",
                       "filter": lambda d: d["admit_medical"] == 0},
        "medical":    {"label": "Medical admission",            "group": "Admitting service",
                       "filter": lambda d: d["admit_medical"] == 1},
        "creat_low":  {"label": f"Baseline creatinine ≤ {q1:.2f} mg/dL",
                       "group": "Baseline creatinine tertile",
                       "filter": lambda d: d["baseline_creat"] <= q1},
        "creat_mid":  {"label": f"Baseline creatinine {q1:.2f}–{q2:.2f} mg/dL",
                       "group": "Baseline creatinine tertile",
                       "filter": lambda d: (d["baseline_creat"] > q1) & (d["baseline_creat"] <= q2)},
        "creat_high": {"label": f"Baseline creatinine > {q2:.2f} mg/dL",
                       "group": "Baseline creatinine tertile",
                       "filter": lambda d: d["baseline_creat"] > q2},
        "hosp_1":     {"label": "Hospital 1",                   "group": "Study hospital",
                       "filter": lambda d: d["hospital"] == 1},
        "hosp_2":     {"label": "Hospital 2",                   "group": "Study hospital",
                       "filter": lambda d: d["hospital"] == 2},
        "hosp_3":     {"label": "Hospital 3",                   "group": "Study hospital",
                       "filter": lambda d: d["hospital"] == 3},
        "hosp_4":     {"label": "Hospital 4",                   "group": "Study hospital",
                       "filter": lambda d: d["hospital"] == 4},
    }


def load() -> pd.DataFrame:
    df = pd.read_csv(CSV_PATH)
    assert len(df) == 5060, f"expected 5060 rows, got {len(df)}"
    assert df["alert"].sum() == 2532, "alert arm should have 2532"
    return df


# ---------- formatting helpers ----------

def fmt_median_iqr(series: pd.Series, decimals: int = 0) -> str:
    s = series.dropna()
    if s.empty:
        return "—"
    q1, med, q3 = np.percentile(s, [25, 50, 75])
    if decimals == 0:
        return f"{med:.0f} ({q1:.0f}, {q3:.0f})"
    return f"{med:.{decimals}f} ({q1:.{decimals}f}, {q3:.{decimals}f})"


def fmt_count_pct(series: pd.Series, total: int) -> str:
    n = int(series.sum())
    pct = 100 * n / total if total else 0
    if pct < 1:
        return f"{n} ({pct:.2f}%)"
    if pct < 10:
        return f"{n} ({pct:.1f}%)"
    return f"{n} ({pct:.0f}%)"


def row_continuous(label, col, a, u, decimals=0):
    return {
        "label": label, "type": "continuous",
        "alert": fmt_median_iqr(a[col], decimals),
        "usual": fmt_median_iqr(u[col], decimals),
    }


def row_categorical(label, col, a, u, value=1):
    return {
        "label": label, "type": "categorical",
        "alert": fmt_count_pct(a[col] == value, len(a)),
        "usual": fmt_count_pct(u[col] == value, len(u)),
    }


# ---------- statistics ----------

def wilson_ci(k: int, n: int, z: float = 1.959963984540054) -> tuple[float, float]:
    """95% Wilson score CI for a proportion."""
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    denom = 1 + z * z / n
    center = (p + z * z / (2 * n)) / denom
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denom
    return (max(0.0, center - half), min(1.0, center + half))


def mh_relative_risk(strata: list[tuple[int, int, int, int]]):
    """Mantel-Haenszel RR with 95% CI (Greenland-Robins log variance).

    Each stratum: (a_ev, a_n, u_ev, u_n). Strata with N==0 or m1==0 are skipped.
    Returns (rr, lo, hi) or Nones if undefined (numerator or denominator zero).
    """
    R_sum = 0.0  # Σ a_k n_{2k} / N_k
    S_sum = 0.0  # Σ b_k n_{1k} / N_k
    P_sum = 0.0  # Σ (n_{1k} n_{2k} m_{1k} - a_k b_k N_k) / N_k^2
    for a_ev, a_n, u_ev, u_n in strata:
        N = a_n + u_n
        m1 = a_ev + u_ev
        if N == 0 or m1 == 0:
            continue
        R_sum += a_ev * u_n / N
        S_sum += u_ev * a_n / N
        P_sum += (a_n * u_n * m1 - a_ev * u_ev * N) / (N * N)
    if R_sum == 0 or S_sum == 0:
        return (None, None, None)
    rr = R_sum / S_sum
    var_log = P_sum / (R_sum * S_sum)
    if var_log <= 0:
        return (rr, None, None)
    half = 1.959963984540054 * math.sqrt(var_log)
    return (rr, math.exp(math.log(rr) - half), math.exp(math.log(rr) + half))


def cmh_chi_square(strata: list[tuple[int, int, int, int]]) -> tuple[float, float] | tuple[None, None]:
    """Cochran-Mantel-Haenszel chi-square (1 df), no continuity correction.

    Each stratum: (a_ev, a_n, u_ev, u_n) = (alert events, alert total, usual events, usual total).
    Returns (chi2_stat, p_value). Skips strata with N==0 or no variance.
    """
    numer = 0.0
    var = 0.0
    for a_ev, a_n, u_ev, u_n in strata:
        N = a_n + u_n
        m1 = a_ev + u_ev
        if N < 2 or m1 == 0 or m1 == N:
            continue
        E = a_n * m1 / N
        V = a_n * u_n * m1 * (N - m1) / (N * N * (N - 1))
        numer += a_ev - E
        var += V
    if var <= 0:
        return (None, None)
    stat = numer * numer / var
    p = 1.0 - chi2.cdf(stat, df=1)
    return (stat, p)


def fmt_p(p: float | None) -> str:
    if p is None:
        return "—"
    if p < 0.0001:
        return f"{p:.2e}"
    if p < 0.01:
        return f"{p:.3f}"
    if p < 0.1:
        return f"{p:.2f}"
    return f"{p:.2f}"


# ---------- Table 1 ----------

def build_table1_for(a, u):
    demographics = [
        row_continuous("Age (years)", "age", a, u),
        row_categorical("Female sex", "sex", a, u),
    ]
    location = [
        row_categorical("Medical admission", "admit_medical", a, u),
        row_categorical("Patient in the ICU", "icu_at_rand", a, u),
        row_categorical("Patient in the emergency department", "er_at_rand", a, u),
        row_categorical("Patient in the ward", "ward_at_rand", a, u),
        row_categorical("Hospital 1", "hospital", a, u, value=1),
        row_categorical("Hospital 2", "hospital", a, u, value=2),
        row_categorical("Hospital 3", "hospital", a, u, value=3),
        row_categorical("Hospital 4", "hospital", a, u, value=4),
    ]
    comorbidities = [
        row_categorical("Chronic kidney disease", "ckd_pmhx", a, u),
        row_categorical("Congestive heart failure", "chf_elx_pmhx", a, u),
        row_categorical("COPD", "pulm_disease_elx_pmhx", a, u),
        row_categorical("Diabetes mellitus", "diabetes_elx_pmhx", a, u),
        row_categorical("Hypertension", "htn_elx_pmhx", a, u),
        row_categorical("Malignancy", "malignancy_elx_pmhx", a, u),
        row_categorical("Liver disease", "liverdisease_elx_pmhx", a, u),
    ]
    aki_stage = [
        row_categorical("Stage 0", "aki_stage_at_rand", a, u, value=0),
        row_categorical("Stage 1", "aki_stage_at_rand", a, u, value=1),
        row_categorical("Stage 2", "aki_stage_at_rand", a, u, value=2),
        row_categorical("Stage 3", "aki_stage_at_rand", a, u, value=3),
    ]
    labs = [
        row_continuous("eGFR on admission (ml/min/1.73m²)", "egfr_at_admit", a, u),
        row_continuous("Creatinine at randomization (mg/dL)", "creat_at_rand", a, u, decimals=1),
        row_continuous("Creatinine at admission (mg/dL)", "admit_creatinine", a, u, decimals=1),
        row_continuous("Sodium (meq/L)", "sodium_at_rand", a, u),
        row_continuous("Potassium (meq/L)", "potassium_at_rand", a, u, decimals=1),
        row_continuous("Chloride (meq/L)", "chloride_at_rand", a, u),
        row_continuous("Bicarbonate (meq/L)", "bicarbonate_at_rand", a, u),
        row_continuous("Anion gap (meq/L)", "aniongap_at_rand", a, u),
        row_continuous("Blood urea nitrogen (mg/dL)", "bun_at_rand", a, u),
        row_continuous("White blood cell count (×1000/μL)", "wbcc_at_rand", a, u, decimals=1),
        row_continuous("Hemoglobin (g/dL)", "hemoglobin_at_rand", a, u, decimals=1),
        row_continuous("Platelet count (×1000/μL)", "plateletcount_at_rand", a, u),
        row_continuous("Modified SOFA score", "sofa_at_rand", a, u),
    ]
    exposures = [
        row_categorical("Contrast in prior 72h", "prior_cs72", a, u),
        row_categorical("RAASi", "on_acearb", a, u),
        row_categorical("NSAID", "on_nsaid", a, u),
        row_categorical("PPI", "on_ppi", a, u),
        row_categorical("One medication of interest", "num_med", a, u, value=1),
        row_categorical("Two medications of interest", "num_med", a, u, value=2),
        row_categorical("Three medications of interest", "num_med", a, u, value=3),
    ]
    return [
        {"title": "Demographics", "rows": demographics},
        {"title": "Hospital location", "rows": location},
        {"title": "Comorbidities", "rows": comorbidities},
        {"title": "AKI stage at randomization", "rows": aki_stage},
        {"title": "Laboratory values", "rows": labs},
        {"title": "Exposures prior to AKI", "rows": exposures},
    ]


def build_table1(df, filters):
    out = {}
    for key, spec in filters.items():
        sub = df[spec["filter"](df)]
        a = sub[sub["alert"] == 1]
        u = sub[sub["alert"] == 0]
        out[key] = {
            "label": spec["label"],
            "group": spec["group"],
            "n_alert": int(len(a)),
            "n_usual": int(len(u)),
            "sections": build_table1_for(a, u),
        }
    return {
        "subgroups": out,
        "footnote": (
            "Data are median (interquartile range) or n (%). "
            "Race/ethnicity and depression are not shown because they were removed "
            "during deidentification of the source dataset."
        ),
    }


# ---------- Fig 1 & Fig 2 (binary outcomes: two arms × stratified by hospital) ----------

def strata_by_hospital(sub: pd.DataFrame, outcome: str) -> list[tuple[int, int, int, int]]:
    strata = []
    for _, g in sub.groupby("hospital"):
        a = g[g["alert"] == 1]
        u = g[g["alert"] == 0]
        y = g[outcome].fillna(0)
        a_ev = int(y[g["alert"] == 1].sum())
        u_ev = int(y[g["alert"] == 0].sum())
        strata.append((a_ev, len(a), u_ev, len(u)))
    return strata


def panel_binary(sub: pd.DataFrame, outcome: str) -> dict:
    """Return arm-level summary + CMH stats for one subgroup/outcome."""
    y = sub[outcome].fillna(0).astype(int)
    a = sub[sub["alert"] == 1]
    u = sub[sub["alert"] == 0]
    a_ev = int(y[sub["alert"] == 1].sum())
    u_ev = int(y[sub["alert"] == 0].sum())
    a_p = a_ev / len(a) if len(a) else float("nan")
    u_p = u_ev / len(u) if len(u) else float("nan")
    a_lo, a_hi = wilson_ci(a_ev, len(a))
    u_lo, u_hi = wilson_ci(u_ev, len(u))
    strata = strata_by_hospital(sub, outcome)
    rr, rr_lo, rr_hi = mh_relative_risk(strata)
    stat, p = cmh_chi_square(strata)
    return {
        "alert": {"n": len(a), "events": a_ev, "prop": a_p, "ci_lo": a_lo, "ci_hi": a_hi},
        "usual": {"n": len(u), "events": u_ev, "prop": u_p, "ci_lo": u_lo, "ci_hi": u_hi},
        "rr": rr, "rr_lo": rr_lo, "rr_hi": rr_hi,
        "cmh_p": p,
    }


def build_fig1(df: pd.DataFrame, filters: dict) -> dict:
    """Composite outcome rates by MOI subgroup, precomputed for each cohort restriction."""
    subgroups = {}
    for ckey, spec in filters.items():
        cohort = df[spec["filter"](df)]
        panels = []
        for pkey, (plabel, pfilt) in MOI_PANELS.items():
            sub = cohort[pfilt(cohort)]
            if len(sub) == 0:
                continue
            panels.append({
                "key": pkey, "label": plabel,
                "data": panel_binary(sub, "composite_outcome"),
            })
        subgroups[ckey] = {
            "label": spec["label"],
            "group": spec["group"],
            "n_alert": int((cohort["alert"] == 1).sum()),
            "n_usual": int((cohort["alert"] == 0).sum()),
            "panels": panels,
        }
    return {
        "title": "Composite outcome: progression of AKI, dialysis, or death within 14 days",
        "outcome_label": "Composite outcome rate",
        "subgroups": subgroups,
        "footnote": (
            "Bars are 95% Wilson confidence intervals. P-values are "
            "Cochran-Mantel-Haenszel chi-square tests stratified by study hospital. "
            "The cohort selector restricts the denominator; panels are the published MOI strata."
        ),
    }


def build_fig2(df: pd.DataFrame, filters: dict) -> dict:
    """MOI discontinuation rate within 24h, precomputed for each cohort restriction."""
    df = df.copy()
    df["_any_stopped24"] = (
        ((df["on_nsaid"] == 1) & (df["nsaid_stopped24"] == 1))
        | ((df["on_acearb"] == 1) & (df["acearb_stopped24"] == 1))
        | ((df["on_ppi"] == 1) & (df["ppi_stopped24"] == 1))
    ).astype(int)

    subgroups = {}
    for ckey, spec in filters.items():
        cohort = df[spec["filter"](df)]
        panels = [{
            "key": "any", "label": "Any MOI",
            "data": panel_binary(cohort, "_any_stopped24"),
        }]
        for pkey, col, label in [
            ("nsaid", "nsaid_stopped24", "NSAID (among NSAID-exposed)"),
            ("raasi", "acearb_stopped24", "RAASi (among RAASi-exposed)"),
            ("ppi",   "ppi_stopped24",    "PPI (among PPI-exposed)"),
        ]:
            exposed_col = {"nsaid": "on_nsaid", "raasi": "on_acearb", "ppi": "on_ppi"}[pkey]
            sub = cohort[cohort[exposed_col] == 1]
            if len(sub) == 0:
                continue
            panels.append({"key": pkey, "label": label, "data": panel_binary(sub, col)})
        subgroups[ckey] = {
            "label": spec["label"],
            "group": spec["group"],
            "n_alert": int((cohort["alert"] == 1).sum()),
            "n_usual": int((cohort["alert"] == 0).sum()),
            "panels": panels,
        }
    return {
        "title": "Medication-of-interest discontinuation within 24 hours",
        "outcome_label": "Discontinuation rate",
        "subgroups": subgroups,
        "footnote": (
            "Bars are 95% Wilson confidence intervals. P-values are "
            "Cochran-Mantel-Haenszel chi-square tests stratified by study hospital. "
            "'Any MOI' counts a patient as discontinued if they stopped any MOI "
            "they were receiving at randomization."
        ),
    }


# ---------- Fig 3 (subgroup forest plot) ----------

def logistic_interaction_p(df: pd.DataFrame, outcome: str, subgroup_col: str) -> float | None:
    """Wald p-value for alert×subgroup interaction in logistic regression adjusted for hospital.

    For multi-level subgroups, returns joint Wald p for all interaction coefficients.
    Matches the paper's approach of 'subgroup-by-alert interaction in a binomial regression
    model adjusted for hospital'.
    """
    d = df[[outcome, "alert", subgroup_col, "hospital"]].dropna().copy()
    d[outcome] = d[outcome].astype(int)
    d[subgroup_col] = d[subgroup_col].astype("category")
    formula = f"{outcome} ~ alert * C({subgroup_col}) + C(hospital)"
    try:
        m = smf.glm(formula, data=d, family=sm.families.Binomial()).fit(disp=0)
    except Exception:
        return None
    # Joint Wald test over all rows matching the interaction term
    interaction_idx = [i for i, name in enumerate(m.params.index) if name.startswith(f"alert:C({subgroup_col})")]
    if not interaction_idx:
        return None
    R = np.zeros((len(interaction_idx), len(m.params)))
    for r, i in enumerate(interaction_idx):
        R[r, i] = 1.0
    try:
        wald = m.wald_test(R, scalar=True)
        return float(wald.pvalue)
    except Exception:
        return None


def breslow_day_p(df: pd.DataFrame, outcome: str, strata_col: str) -> float | None:
    """Breslow-Day test of OR homogeneity across levels of strata_col."""
    tables = []
    for _, g in df.groupby(strata_col):
        y = g[outcome].fillna(0).astype(int)
        a = int(y[g["alert"] == 1].sum())
        b = int((g["alert"] == 1).sum() - a)
        c = int(y[g["alert"] == 0].sum())
        d = int((g["alert"] == 0).sum() - c)
        if min(a + b, c + d, a + c, b + d) == 0:
            continue
        tables.append(np.array([[a, b], [c, d]]))
    if len(tables) < 2:
        return None
    try:
        st = StratifiedTable(tables)
        return float(st.test_equal_odds().pvalue)
    except Exception:
        return None


def panel_binary_strata(sub: pd.DataFrame, outcome: str) -> dict:
    """Like panel_binary but returns only the essentials for forest plot rendering."""
    y = sub[outcome].fillna(0).astype(int)
    a_mask = sub["alert"] == 1
    u_mask = sub["alert"] == 0
    a_n, u_n = int(a_mask.sum()), int(u_mask.sum())
    a_ev, u_ev = int(y[a_mask].sum()), int(y[u_mask].sum())
    strata = strata_by_hospital(sub, outcome)
    rr, lo, hi = mh_relative_risk(strata)
    return {
        "n_alert": a_n, "n_usual": u_n,
        "events_alert": a_ev, "events_usual": u_ev,
        "p_alert": a_ev / a_n if a_n else None,
        "p_usual": u_ev / u_n if u_n else None,
        "rr": rr, "rr_lo": lo, "rr_hi": hi,
    }


def build_fig3(df: pd.DataFrame) -> dict:
    """Subgroup forest plot: effect of alert on composite outcome across pre-specified subgroups."""
    d = df.copy()

    # Baseline creatinine tertiles (drop patients with missing baseline_creat)
    creat = d["baseline_creat"].dropna()
    q1, q2 = np.percentile(creat, [33.333, 66.667])
    d["creat_tertile"] = pd.cut(
        d["baseline_creat"],
        bins=[-np.inf, q1, q2, np.inf],
        labels=["low", "mid", "high"],
    )

    subgroups = [
        {
            "name": "ICU at randomization",
            "col": "icu_at_rand",
            "levels": [("Non-ICU", 0), ("ICU", 1)],
        },
        {
            "name": "Admitting service",
            "col": "admit_medical",
            "levels": [("Surgical", 0), ("Medical", 1)],
        },
        {
            "name": f"Baseline creatinine (tertiles: ≤{q1:.2f}, ≤{q2:.2f}, >{q2:.2f})",
            "col": "creat_tertile",
            "levels": [("Low", "low"), ("Mid", "mid"), ("High", "high")],
        },
        {
            "name": "Study hospital",
            "col": "hospital",
            "levels": [(f"Hospital {h}", h) for h in [1, 2, 3, 4]],
            "interaction_method": "breslow-day",
        },
    ]

    out_subgroups = []
    for sg in subgroups:
        col = sg["col"]
        data_available = d[d[col].notna()]
        rows = []
        for label, value in sg["levels"]:
            sub = data_available[data_available[col] == value]
            rows.append({"label": label, **panel_binary_strata(sub, "composite_outcome")})

        if sg.get("interaction_method") == "breslow-day":
            p = breslow_day_p(data_available, "composite_outcome", col)
            method = "Breslow-Day test for OR homogeneity"
        else:
            p = logistic_interaction_p(data_available, "composite_outcome", col)
            method = "alert × subgroup interaction (logistic regression adjusted for hospital)"

        out_subgroups.append({
            "name": sg["name"],
            "interaction_p": p,
            "interaction_method": method,
            "rows": rows,
        })

    return {
        "title": "Effect of alert on composite outcome across pre-specified subgroups",
        "subgroups": out_subgroups,
        "footnote": (
            "Squares are Mantel-Haenszel relative risks (stratified by hospital, except for "
            "the hospital subgroup itself); bars are 95% confidence intervals. "
            "Interaction p-values test whether the alert effect differs across levels of the subgroup."
        ),
    }


# ---------- Table 3 (safety outcomes) ----------

def newcombe_diff_ci(a_ev: int, a_n: int, u_ev: int, u_n: int):
    """Newcombe (score) 95% CI for difference of proportions (alert - usual), as fractions."""
    if a_n == 0 or u_n == 0:
        return (None, None, None)
    p1, p2 = a_ev / a_n, u_ev / u_n
    diff = p1 - p2
    l1, h1 = wilson_ci(a_ev, a_n)
    l2, h2 = wilson_ci(u_ev, u_n)
    lo = diff - math.sqrt((p1 - l1) ** 2 + (h2 - p2) ** 2)
    hi = diff + math.sqrt((h1 - p1) ** 2 + (p2 - l2) ** 2)
    return (diff, lo, hi)


def hodges_lehmann(x: np.ndarray, y: np.ndarray) -> float:
    """HL estimator: median of pairwise differences (x_i - y_j)."""
    x = np.asarray(x)
    y = np.asarray(y)
    return float(np.median(np.subtract.outer(x, y)))


def hl_bootstrap_ci(x: np.ndarray, y: np.ndarray, B: int = 1000, seed: int = 42):
    rng = np.random.default_rng(seed)
    x = np.asarray(x)
    y = np.asarray(y)
    estimates = np.empty(B)
    for b in range(B):
        xb = rng.choice(x, len(x), replace=True)
        yb = rng.choice(y, len(y), replace=True)
        estimates[b] = hodges_lehmann(xb, yb)
    return (float(np.percentile(estimates, 2.5)), float(np.percentile(estimates, 97.5)))


def safety_binary(label: str, col: str, a: pd.DataFrame, u: pd.DataFrame) -> dict:
    a_ev = int(a[col].fillna(0).sum())
    u_ev = int(u[col].fillna(0).sum())
    diff, lo, hi = newcombe_diff_ci(a_ev, len(a), u_ev, len(u))
    diff_text = "—"
    if diff is not None:
        diff_text = f"{100*diff:+.1f} ({100*lo:+.1f}, {100*hi:+.1f})"
    return {
        "label": label, "type": "binary",
        "alert": fmt_count_pct(a[col].fillna(0).astype(bool), len(a)),
        "usual": fmt_count_pct(u[col].fillna(0).astype(bool), len(u)),
        "diff": diff_text,
    }


def safety_continuous(label: str, col: str, a: pd.DataFrame, u: pd.DataFrame, decimals=1) -> dict:
    ax = a[col].dropna().values
    uy = u[col].dropna().values
    if len(ax) == 0 or len(uy) == 0:
        diff_text = "—"
    else:
        hl = hodges_lehmann(ax, uy)
        # Subsample for bootstrap if very large (>2000) to keep precompute snappy.
        if len(ax) * len(uy) > 2_000_000:
            rng = np.random.default_rng(0)
            ax_s = rng.choice(ax, 2000, replace=False) if len(ax) > 2000 else ax
            uy_s = rng.choice(uy, 2000, replace=False) if len(uy) > 2000 else uy
            lo, hi = hl_bootstrap_ci(ax_s, uy_s, B=500)
        else:
            lo, hi = hl_bootstrap_ci(ax, uy, B=1000)
        diff_text = f"{hl:+.{decimals}f} ({lo:+.{decimals}f}, {hi:+.{decimals}f})"
    return {
        "label": label, "type": "continuous",
        "alert": fmt_median_iqr(a[col], decimals),
        "usual": fmt_median_iqr(u[col], decimals),
        "diff": diff_text,
    }


def build_table3(df: pd.DataFrame) -> dict:
    """Safety outcomes by MOI subgroup — matches the paper's Table 3 structure."""
    panels = []

    nsaid = df[df["on_nsaid"] == 1]
    a, u = nsaid[nsaid["alert"] == 1], nsaid[nsaid["alert"] == 0]
    panels.append({
        "label": "NSAID subgroup",
        "n_alert": len(a), "n_usual": len(u),
        "rows": [
            safety_binary("Opioid prescription", "opioid", a, u),
            safety_continuous("Max pain score", "max_pain_score14", a, u, decimals=1),
        ],
    })

    raasi = df[df["on_acearb"] == 1]
    a, u = raasi[raasi["alert"] == 1], raasi[raasi["alert"] == 0]
    panels.append({
        "label": "RAASi subgroup",
        "n_alert": len(a), "n_usual": len(u),
        "rows": [
            safety_continuous("Max systolic BP (mmHg)", "systolic_max14", a, u, decimals=1),
            safety_continuous("Max diastolic BP (mmHg)", "diastolic_max14", a, u, decimals=1),
            safety_binary("Mechanical ventilation", "ventorder", a, u),
        ],
    })

    ppi = df[df["on_ppi"] == 1]
    a, u = ppi[ppi["alert"] == 1], ppi[ppi["alert"] == 0]
    panels.append({
        "label": "PPI subgroup",
        "n_alert": len(a), "n_usual": len(u),
        "rows": [
            safety_binary("PRBC transfusion", "transfuserbc", a, u),
            safety_continuous("Minimum hemoglobin (g/dL)", "min_hemoglobin", a, u, decimals=1),
            safety_continuous("Max pain score", "max_pain_score14", a, u, decimals=1),
        ],
    })

    return {
        "panels": panels,
        "footnote": (
            "Values are n (%) or median (IQR). Difference column shows the alert-minus-usual-care "
            "contrast with 95% CI (Newcombe method for proportions, Hodges-Lehmann estimate with "
            "bootstrap CI for medians). Each panel is restricted to patients receiving that MOI at "
            "randomization, so no further subgroup filter is applied."
        ),
    }


# ---------- Table 2 (secondary outcomes) ----------

def row_binary_outcome(label: str, col: str, sub: pd.DataFrame) -> dict:
    r = panel_binary(sub, col)
    a, u = r["alert"], r["usual"]
    rr_text = "—"
    if r["rr"] is not None:
        rr_text = f"{r['rr']:.2f} ({r['rr_lo']:.2f}–{r['rr_hi']:.2f})"
    return {
        "label": label,
        "type": "binary",
        "alert": fmt_count_pct(pd.Series([True] * a["events"] + [False] * (a["n"] - a["events"])), a["n"]),
        "usual": fmt_count_pct(pd.Series([True] * u["events"] + [False] * (u["n"] - u["events"])), u["n"]),
        "rr": rr_text,
        "p": fmt_p(r["cmh_p"]),
    }


def row_continuous_outcome(label: str, col: str, a: pd.DataFrame, u: pd.DataFrame, decimals=1) -> dict:
    return {
        "label": label,
        "type": "continuous",
        "alert": fmt_median_iqr(a[col], decimals),
        "usual": fmt_median_iqr(u[col], decimals),
        "rr": "—",
        "p": "—",
    }


def build_table2(df: pd.DataFrame, filters: dict) -> dict:
    binary_rows = [
        ("Progression of AKI", "aki_progression14"),
        ("Dialysis", "dialysis14"),
        ("Death", "death14"),
        ("Inpatient mortality", "inpatient_mortality"),
        ("Inpatient dialysis", "inpatient_dialysis"),
        ("Progression to stage 2 AKI", "progression2"),
        ("Progression to stage 3 AKI", "progression3"),
        ("30-day readmission", "readmit30"),
        ("Inpatient kidney consult", "renalconsult14"),
    ]
    continuous_rows = [
        ("Duration of AKI (days)", "duration_of_aki", 1),
        ("Length of stay post-randomization (days)", "los_since_alert", 1),
    ]

    out = {}
    for key, spec in filters.items():
        sub = df[spec["filter"](df)]
        a = sub[sub["alert"] == 1]
        u = sub[sub["alert"] == 0]
        rows = [row_binary_outcome(lab, col, sub) for lab, col in binary_rows]
        rows += [row_continuous_outcome(lab, col, a, u, dec) for lab, col, dec in continuous_rows]
        out[key] = {
            "label": spec["label"],
            "group": spec["group"],
            "n_alert": int(len(a)),
            "n_usual": int(len(u)),
            "rows": rows,
        }
    return {
        "subgroups": out,
        "footnote": (
            "Binary outcomes: counts (percentage), relative risk with 95% CI (log method), "
            "p-value from Cochran-Mantel-Haenszel chi-square stratified by hospital. "
            "Continuous outcomes: median (IQR); p-values and RRs are not applicable. "
            "Hospitalization cost is omitted because it is not included in the deidentified dataset."
        ),
    }


# ---------- driver ----------

def main():
    df = load()
    DATA_DIR.mkdir(exist_ok=True)

    filters = get_subgroup_filters(df)

    (DATA_DIR / "table1.json").write_text(json.dumps(build_table1(df, filters), indent=2))
    print(f"wrote {DATA_DIR / 'table1.json'}")

    fig1 = build_fig1(df, filters)
    (DATA_DIR / "fig1.json").write_text(json.dumps(fig1, indent=2))
    print(f"wrote {DATA_DIR / 'fig1.json'}")

    fig2 = build_fig2(df, filters)
    (DATA_DIR / "fig2.json").write_text(json.dumps(fig2, indent=2))
    print(f"wrote {DATA_DIR / 'fig2.json'}")

    (DATA_DIR / "table2.json").write_text(json.dumps(build_table2(df, filters), indent=2))
    print(f"wrote {DATA_DIR / 'table2.json'}")

    fig3 = build_fig3(df)
    (DATA_DIR / "fig3.json").write_text(json.dumps(fig3, indent=2))
    print(f"wrote {DATA_DIR / 'fig3.json'}")

    (DATA_DIR / "table3.json").write_text(json.dumps(build_table3(df), indent=2))
    print(f"wrote {DATA_DIR / 'table3.json'}")

    # sanity checks against published numbers
    overall = fig1["subgroups"]["overall"]["panels"][0]["data"]
    print(
        f"  overall composite: {overall['alert']['events']}/{overall['alert']['n']} "
        f"({100*overall['alert']['prop']:.1f}%) alert vs "
        f"{overall['usual']['events']}/{overall['usual']['n']} "
        f"({100*overall['usual']['prop']:.1f}%) usual, CMH p={overall['cmh_p']:.3f}"
    )
    ppi = fig1["subgroups"]["overall"]["panels"][3]["data"]
    print(
        f"  PPI composite:     {ppi['alert']['events']}/{ppi['alert']['n']} "
        f"({100*ppi['alert']['prop']:.1f}%) alert vs "
        f"{ppi['usual']['events']}/{ppi['usual']['n']} "
        f"({100*ppi['usual']['prop']:.1f}%) usual, RR={ppi['rr']:.2f}, CMH p={ppi['cmh_p']:.3f}"
    )
    any_moi = fig2["subgroups"]["overall"]["panels"][0]["data"]
    print(
        f"  any MOI stop 24h:  {any_moi['alert']['events']}/{any_moi['alert']['n']} "
        f"({100*any_moi['alert']['prop']:.1f}%) alert vs "
        f"{any_moi['usual']['events']}/{any_moi['usual']['n']} "
        f"({100*any_moi['usual']['prop']:.1f}%) usual, CMH p={any_moi['cmh_p']:.4f}"
    )
    print("  fig3 interaction p-values:")
    for sg in fig3["subgroups"]:
        p = sg["interaction_p"]
        name = sg["name"].encode("ascii", "replace").decode()
        print(f"    {name}: p={p:.3f}" if p is not None else f"    {name}: p=NA")


if __name__ == "__main__":
    main()
