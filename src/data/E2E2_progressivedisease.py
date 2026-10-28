"""Dyslipoproteinemia analysis in ApoE2/E2 homozygotes (UK Biobank).

Compares lipid parameters between ApoE2/E2 homozygotes ("ε2ε2") and those
with no ε2 allele ("No ε2") -- ApoE2 heterozygotes are excluded, since their
lipid profile is intermediate and not the comparison of interest -- across
three triglyceride (TG) categories.

Produces the manuscript's Tables 2-4, Figures 1-2, and Supplemental Figure 1
from the raw UK Biobank extract (the manuscript's Table 1 comes from a
different source and isn't produced here), plus two supporting outputs not
published in the manuscript: Table4_Normality (the normality diagnostics
behind Table 4) and Table5_EquationComparison (a comparison of the two
esVLDLC equations).

Derived lipid metrics
----------------------
- msLDLC / msVLDLC: modified Sampson equation estimates of LDL-C and VLDL-C.
- esVLDLC: enhanced Sampson VLDL-C estimate, using whichever of two
  equations matches the individual's own genotype. ε2ε2 individuals have
  unusually cholesterol-rich VLDL particles, which breaks the general
  population equation, hence the separate genotype-specific equations.
- esVLDLC1 / esVLDLC2: the same two genotype-specific equations, each
  applied to every individual regardless of their genotype, so the two
  equations can be compared directly (see Table5_EquationComparison)
  rather than only ever seeing the genotype-matched value.
- eApoB: novel NIH equation estimating ApoB.

Usage
-----
    python E2E2_progressivedisease.py

Run from anywhere -- paths are resolved relative to this file, not the
current working directory. Expects the raw data at
``data/raw/final_E2E2_dataset.xlsx`` (see README for how to obtain it) and
writes all tables/figures to ``reports/``.
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as stats
import seaborn as sns

# Windows consoles default to a legacy encoding that can't print the "ε"
# characters used throughout this script's output; force UTF-8 so running
# `python E2E2_progressivedisease.py` doesn't crash partway through.
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

# -------------------------------------------------------------------------
# Paths (resolved relative to this file so the script works regardless of
# the caller's working directory)
# -------------------------------------------------------------------------
REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_PATH = REPO_ROOT / "data" / "raw" / "final_E2E2_dataset.xlsx"
REPORTS_DIR = REPO_ROOT / "reports"

# -------------------------------------------------------------------------
# Shared constants
# -------------------------------------------------------------------------
GENOTYPE_COL = "ε2 Status"
EPSILON2_HOMOZYGOTE = "ε2ε2"
NO_EPSILON2 = "No ε2"

TG_CATEGORY_BINS = [0, 180, 350, np.inf]
TG_CATEGORY_LABELS = ["TG <180 mg/dL", "TG 180-349 mg/dL", "TG ≥350 mg/dL"]

# Column/ratio names use short internal abbreviations (e.g. "msLDLC"). This
# maps them to publication-style labels for figure axes and table
# headers/indices, without renaming the underlying DataFrame columns.
DISPLAY_LABEL_REPLACEMENTS = [
    ("esVLDLC", "VLDL-C"),
    ("msLDLC", "LDL-C"),
    ("NonHDLC", "Non-HDL-C"),
    ("HDLC", "HDL-C"),
]


def display_label(name):
    """Map an internal column/ratio name to a publication-style label."""
    label = name
    for old, new in DISPLAY_LABEL_REPLACEMENTS:
        label = label.replace(old, new)
    return label


def format_p_value(p):
    """Format a p-value for display: blank if missing, "<0.0001" if tiny."""
    if pd.isna(p):
        return ""
    if p < 0.0001:
        return "<0.0001"
    return f"{p:.4f}"


# -------------------------------------------------------------------------
# Data loading and feature engineering
# -------------------------------------------------------------------------
def load_data(path):
    """Load the raw UK Biobank extract and keep only the columns used here."""
    df = pd.read_excel(path)
    df = df[
        [
            "apoE genotype Combined",
            "diabetes by lab or doctor",
            "BMI",
            "Sex",
            "Age",
            "Current smoker",
            "Apo B",
            "HDLC",
            "TG",
            "NonHDLC",
            "dLDLC",
        ]
    ]
    df = df.rename(columns={"Apo B": "ApoB", "dLDLC": "direct LDLC"})
    return df


def add_genotype_columns(df):
    """Flag ApoE2 heterozygotes and classify each row's ε2 Status.

    "E2 Heterozygote" is 1 for one copy of the E2 allele (E2E2 individuals
    have two and are not heterozygotes). GENOTYPE_COL distinguishes ε2ε2
    homozygotes, lone-ε2 heterozygotes, and everyone else; heterozygotes are
    dropped later since their profile is intermediate and not of interest.
    """
    has_e2 = df["apoE genotype Combined"].str.contains("E2")
    has_e2e2 = df["apoE genotype Combined"].str.contains("E2E2")

    df["E2 Heterozygote"] = np.where(has_e2 & ~has_e2e2, 1, 0)
    df[GENOTYPE_COL] = np.select(
        condlist=[has_e2e2, has_e2],  # two copies checked before one copy
        choicelist=[EPSILON2_HOMOZYGOTE, "ε2"],
        default=NO_EPSILON2,
    )
    return df


def clean_cohort(df):
    """Drop ApoE2 heterozygotes and keep only the columns used downstream."""
    df = df[df["E2 Heterozygote"] == 0]
    return df[
        [
            GENOTYPE_COL,
            "ApoB",
            "HDLC",
            "TG",
            "NonHDLC",
            "direct LDLC",
            "Age",
            "Sex",
            "BMI",
            "Current smoker",
            "diabetes by lab or doctor",
        ]
    ].copy()


def engineer_features(df):
    """Add derived LDL-C/VLDL-C/ApoB estimates, ratios, and the TG category."""
    df["msLDLC"] = (
        df["NonHDLC"]
        - (df["TG"] / 8.37)
        - (df["NonHDLC"] * df["TG"] / 2640)
        + (df["TG"] ** 2 / 17400)
    )
    df["msVLDLC"] = (
        (df["TG"] / 8.37)
        + (df["NonHDLC"] * df["TG"] / 2640)
        - (df["TG"] ** 2 / 17400)
    )

    # Two alternative equations for VLDL-C, calibrated for different genotypes.
    no_e2_equation = (
        (df["NonHDLC"] / 3.81)
        - (df["HDLC"] / 8.93)
        + (df["TG"] / 7.73)
        + (df["NonHDLC"] * df["TG"] / 2050)
        - (df["TG"] ** 2 / 13300)
        - (df["ApoB"] / 2.49)
        - (df["ApoB"] * df["TG"] / 3550)
        + 7.46
    )
    e2e2_equation = (
        (0.001 * df["NonHDLC"] * df["TG"])
        - (0.526 * df["ApoB"])
        + (0.318 * df["TG"])
        - (0.0004 * df["TG"] ** 2)
        + 38.025
    )

    # esVLDLC: each individual's own genotype-matched equation.
    esvldlc = np.select(
        condlist=[df[GENOTYPE_COL] == NO_EPSILON2, df[GENOTYPE_COL] == EPSILON2_HOMOZYGOTE],
        choicelist=[no_e2_equation, e2e2_equation],
        default=np.nan,
    )
    df["esVLDLC"] = np.where(esvldlc < 0, 0, esvldlc)

    # esVLDLC1/esVLDLC2: both equations applied to every row, for comparison.
    df["esVLDLC1"] = np.where(no_e2_equation < 0, 0, no_e2_equation)
    df["esVLDLC2"] = np.where(e2e2_equation < 0, 0, e2e2_equation)

    df["eApoB"] = (
        (df["NonHDLC"] / 1.90)
        - (df["TG"] / 22.5)
        + (2.75 * (df["NonHDLC"] / df["HDLC"]))
        + (df["TG"] * df["NonHDLC"] / 12200)
        - (df["TG"] ** 2 / 50400)
        + 13.0
    )

    df["esVLDLC/TG"] = df["esVLDLC"] / df["TG"]
    df["esVLDLC/msLDLC"] = df["esVLDLC"] / df["msLDLC"]
    df["msLDLC/dLDLC"] = df["msLDLC"] / df["direct LDLC"]
    df["ApoB/eApoB"] = df["ApoB"] / df["eApoB"]
    df["NonHDLC/ApoB"] = df["NonHDLC"] / df["ApoB"]
    df["TG/ApoB"] = df["TG"] / df["ApoB"]

    df["TG_category"] = pd.cut(
        df["TG"], bins=TG_CATEGORY_BINS, labels=TG_CATEGORY_LABELS
    )
    return df


# -------------------------------------------------------------------------
# Console diagnostics: skewness/kurtosis and Shapiro-Wilk, by genotype
# -------------------------------------------------------------------------
NORMALITY_DIAGNOSTIC_COLUMNS = [
    "Age",
    "ApoB",
    "HDLC",
    "TG",
    "NonHDLC",
    "msLDLC",
    "BMI",
    "esVLDLC",
    "esVLDLC/TG",
    "esVLDLC/msLDLC",
    "ApoB/eApoB",
    "NonHDLC/ApoB",
    "TG/ApoB",
]


def print_normality_diagnostics(df, columns):
    """Print skewness/kurtosis and a Shapiro-Wilk test for each column."""
    for col in columns:
        group1 = df[df[GENOTYPE_COL] == EPSILON2_HOMOZYGOTE][col].dropna()
        group2 = df[df[GENOTYPE_COL] == NO_EPSILON2][col].dropna()

        skew1, kurt1 = stats.skew(group1), stats.kurtosis(group1)
        skew2, kurt2 = stats.skew(group2), stats.kurtosis(group2)

        print(f"--- Results for **{col}** ---")
        print("**Descriptive Statistics (Skew/Kurtosis):**")
        print(
            f" Group 1 ({EPSILON2_HOMOZYGOTE}): Skew={skew1:.4f}, "
            f"Kurtosis (Fisher)={kurt1:.4f}"
        )
        print(
            f" Group 2 ({NO_EPSILON2}): Skew={skew2:.4f}, "
            f"Kurtosis (Fisher)={kurt2:.4f}"
        )

        w_stat1, p_value1 = stats.shapiro(group1)
        w_stat2, p_value2 = stats.shapiro(group2)
        print("**Shapiro-Wilk Test:**")
        print(f" Group 1: W={w_stat1:.4f}, p={p_value1:.4f}")
        print(f" Group 2: W={w_stat2:.4f}, p={p_value2:.4f}\n")


# -------------------------------------------------------------------------
# Table 2: categorical variables by genotype (manuscript Table 2)
# -------------------------------------------------------------------------
def build_table2(df):
    """Counts/percentages for categorical variables by genotype, with chi-squared p-values."""
    cols_to_summarize = ["Sex", "Current smoker", "diabetes by lab or doctor"]

    def summarize_group(group):
        results = {}
        for col in cols_to_summarize:
            if col == "Sex":
                count = (group[col] == "F").sum()
            else:
                count = (group[col] == 1).sum()
            results[(col, "Count")] = count
            results[(col, "Percentage")] = f"{(count / len(group)) * 100:.0f}%"
        return pd.Series(results)

    summary_df = df.groupby(GENOTYPE_COL).apply(summarize_group).T

    p_values = {}
    for col in cols_to_summarize:
        contingency_table = pd.crosstab(df[col], df[GENOTYPE_COL])
        _, p, _, _ = stats.chi2_contingency(contingency_table)
        p_values[col] = p

    p_value_series = pd.Series(
        [p_values[col] for col in cols_to_summarize],
        index=pd.MultiIndex.from_product([cols_to_summarize, ["p-value"]]),
    )

    summary_df_with_p = pd.concat(
        [summary_df, p_value_series.to_frame("p-value")], axis=0
    )
    return summary_df_with_p.sort_index(level=0)


# -------------------------------------------------------------------------
# Table 3: numerical variables by genotype (manuscript Table 3)
# -------------------------------------------------------------------------
TABLE3_NUMERIC_COLUMNS = [
    "Age",
    "ApoB",
    "HDLC",
    "TG",
    "NonHDLC",
    "msLDLC",
    "BMI",
    "esVLDLC",
    "esVLDLC/TG",
    "esVLDLC/msLDLC",
    "ApoB/eApoB",
    "NonHDLC/ApoB",
    "TG/ApoB",
]
TABLE3_NON_PARAMETRIC_COLUMNS = ["TG", "esVLDLC", "esVLDLC/TG", "esVLDLC/msLDLC", "TG/ApoB"]


def build_table3(df):
    """Descriptive stats (mean/SD etc.) for numeric variables by genotype, with p-values."""
    table3 = (
        df[[GENOTYPE_COL] + TABLE3_NUMERIC_COLUMNS].groupby(GENOTYPE_COL).describe().T
    )

    genotypes = df[GENOTYPE_COL].unique()
    p_values = []
    for var in table3.index.get_level_values(0).unique():
        group1 = df[df[GENOTYPE_COL] == genotypes[0]][var].dropna()
        group2 = df[df[GENOTYPE_COL] == genotypes[1]][var].dropna()

        if var in TABLE3_NON_PARAMETRIC_COLUMNS:
            _, p_val = stats.mannwhitneyu(group1, group2, alternative="two-sided")
        else:
            _, p_val = stats.ttest_ind(group1, group2, equal_var=False)

        p_values.extend([p_val] * len(table3.loc[var]))

    table3["p-value"] = p_values
    return table3.rename(index=display_label, level=0)


# -------------------------------------------------------------------------
# Table4_Normality (unpublished): Shapiro-Wilk diagnostics behind Table 4
# -------------------------------------------------------------------------
TABLE4_NORMALITY_METRICS = [
    "ApoB",
    "TG",
    "NonHDLC",
    "HDLC",
    "esVLDLC",
    "esVLDLC/TG",
    "esVLDLC/msLDLC",
    "ApoB/eApoB",
    "NonHDLC/ApoB",
    "TG/ApoB",
]


def build_table4_normality(df):
    """Shapiro-Wilk normality diagnostics for key metrics, by TG category and genotype."""
    results = []
    for metric in TABLE4_NORMALITY_METRICS:
        row = {"Metric": metric}
        for tg_cat in TG_CATEGORY_LABELS:
            df_group_tg = df[df["TG_category"] == tg_cat]
            for genotype in df[GENOTYPE_COL].unique():
                data = df_group_tg[df_group_tg[GENOTYPE_COL] == genotype][metric].dropna()
                col_prefix = f"({tg_cat}, {genotype})"

                row[f"{col_prefix}_Skewness"] = stats.skew(data)
                row[f"{col_prefix}_Kurtosis"] = stats.kurtosis(data, fisher=True)
                try:
                    w_stat, p_value = stats.shapiro(data)
                except Exception:
                    w_stat, p_value = np.nan, np.nan
                row[f"{col_prefix}_W"] = w_stat
                row[f"{col_prefix}_p"] = p_value
        results.append(row)

    df_results = pd.DataFrame(results).set_index("Metric")

    column_tuples = []
    for tg_cat in TG_CATEGORY_LABELS:
        for genotype_name in [EPSILON2_HOMOZYGOTE, NO_EPSILON2]:
            column_tuples.append((tg_cat, f"{genotype_name} Skewness"))
            column_tuples.append((tg_cat, f"{genotype_name} Kurtosis"))
            column_tuples.append((tg_cat, f"{genotype_name} W"))
            column_tuples.append((tg_cat, f"{genotype_name} p"))
    df_results.columns = pd.MultiIndex.from_tuples(
        column_tuples, names=["TG Group", "Statistic"]
    )

    df_normality = df_results.copy()
    for col in df_normality.columns:
        statistic_name = col[1]
        if statistic_name.endswith("p"):
            df_normality[col] = df_normality[col].apply(
                lambda p: (
                    f"{p:.4f}"
                    if pd.notna(p) and p >= 0.0001
                    else "<0.0001" if pd.notna(p) and p < 0.0001 else ""
                )
            )
        elif (
            statistic_name.endswith("W")
            or "Skewness" in statistic_name
            or "Kurtosis" in statistic_name
        ):
            df_normality[col] = df_normality[col].round(3).astype(str).replace("nan", "")

    df_normality.index = df_normality.index.map(display_label)
    return df_normality


# -------------------------------------------------------------------------
# Table 4: median (IQR) by TG category and genotype (manuscript Table 4)
# -------------------------------------------------------------------------
TABLE4_METRICS = [
    "TG",
    "ApoB",
    "NonHDLC",
    "HDLC",
    "esVLDLC",
    "esVLDLC/TG",
    "esVLDLC/msLDLC",
    "msLDLC/dLDLC",
    "ApoB/eApoB",
    "NonHDLC/ApoB",
    "TG/ApoB",
]


def _median_iqr_and_pvalue(df_group, metric):
    """Median/IQR per genotype plus the Mann-Whitney U p-value, for one metric."""
    stats_dict = {}
    e2e2_data = df_group[df_group[GENOTYPE_COL] == EPSILON2_HOMOZYGOTE][metric].dropna()
    not_e2e2_data = df_group[df_group[GENOTYPE_COL] == NO_EPSILON2][metric].dropna()

    for group_name, data in [(EPSILON2_HOMOZYGOTE, e2e2_data), (NO_EPSILON2, not_e2e2_data)]:
        stats_dict[("median", group_name)] = data.median()
        q1, q3 = data.quantile([0.25, 0.75]).tolist()
        stats_dict[("IQR1", group_name)] = q1
        stats_dict[("IQR3", group_name)] = q3

    if len(e2e2_data) < 3 or len(not_e2e2_data) < 3:
        p_value = np.nan
    else:
        _, p_value = stats.mannwhitneyu(e2e2_data, not_e2e2_data, alternative="two-sided")
    stats_dict[("median", "p-value")] = p_value
    return stats_dict


def build_table4(df):
    """Median (IQR) for key metrics by TG category and genotype, with p-values."""
    row_index_tuples = [
        ("", "count"),
        *[(m, "median") for m in TABLE4_METRICS],
        *[(m, "IQR1") for m in TABLE4_METRICS],
        *[(m, "IQR3") for m in TABLE4_METRICS],
    ]
    final_row_index = pd.MultiIndex.from_tuples(row_index_tuples, names=["Metric", "Statistic"])

    genotype_cols = [EPSILON2_HOMOZYGOTE, NO_EPSILON2, "p-value"]
    column_index_tuples = [
        (tg_cat, col) for tg_cat in TG_CATEGORY_LABELS for col in genotype_cols
    ]
    final_column_index = pd.MultiIndex.from_tuples(
        column_index_tuples, names=["TG Group", "Genotype"]
    )

    final_data = {key: {} for key in final_row_index}

    count_row_key = ("", "count")
    for tg_cat in TG_CATEGORY_LABELS:
        df_group = df[df["TG_category"] == tg_cat]
        final_data[count_row_key][(tg_cat, EPSILON2_HOMOZYGOTE)] = (
            df_group[GENOTYPE_COL] == EPSILON2_HOMOZYGOTE
        ).sum()
        final_data[count_row_key][(tg_cat, NO_EPSILON2)] = (
            df_group[GENOTYPE_COL] == NO_EPSILON2
        ).sum()
        final_data[count_row_key][(tg_cat, "p-value")] = np.nan

    for tg_cat in TG_CATEGORY_LABELS:
        df_group = df[df["TG_category"] == tg_cat]
        for metric in TABLE4_METRICS:
            for (stat, col_name), value in _median_iqr_and_pvalue(df_group, metric).items():
                if (metric, stat) in final_data:
                    final_data[(metric, stat)][(tg_cat, col_name)] = value

    df_final = pd.DataFrame(index=final_row_index, columns=final_column_index)
    for (metric, stat), col_data in final_data.items():
        for (tg_group, genotype_col), value in col_data.items():
            df_final.loc[(metric, stat), (tg_group, genotype_col)] = value
    df_final = df_final.fillna(np.nan)

    table4 = df_final.copy()
    table4 = table4.apply(pd.to_numeric, errors="coerce")

    # Keep the count row pinned at the top; sort everything else by Metric.
    count_row = table4.loc[("", "count")].to_frame().T
    data_rows = table4.drop(index=[("", "count")])
    data_rows_sorted = data_rows.sort_index(level="Metric", kind="stable")
    table4 = pd.concat([count_row, data_rows_sorted])

    data_cols = [col for col in table4.columns if col[1] in [EPSILON2_HOMOZYGOTE, NO_EPSILON2]]
    p_value_cols = [col for col in table4.columns if col[1] == "p-value"]
    count_index = ("", "count")

    table4[data_cols] = table4[data_cols].round(2)
    if count_index in table4.index:
        table4.loc[count_index, data_cols] = (
            table4.loc[count_index, data_cols].astype(float).round(0).astype(int)
        )
        table4.loc[count_index, p_value_cols] = np.nan

    table4[p_value_cols] = table4[p_value_cols].map(format_p_value)
    table4 = table4.fillna("")
    return table4.rename(index=display_label, level=0)


# -------------------------------------------------------------------------
# Table5_EquationComparison (unpublished): esVLDLC equation comparison
# -------------------------------------------------------------------------
TABLE5_COLUMNS = ["msLDLC", "msVLDLC", "esVLDLC1", "esVLDLC2"]


def build_table5_equation_comparison(df):
    """Mean/SD comparison of msLDLC, msVLDLC, esVLDLC1, esVLDLC2 by TG category and genotype."""
    return df.groupby(["TG_category", GENOTYPE_COL])[TABLE5_COLUMNS].describe()


# -------------------------------------------------------------------------
# Supplemental Figure 1: distributions of TG/ApoB/NonHDLC-derived metrics
# -------------------------------------------------------------------------
DISTRIBUTION_PARAMETERS = ["TG", "ApoB", "NonHDLC", "NonHDLC/ApoB", "ApoB/eApoB", "esVLDLC/TG"]
DISTRIBUTION_SUBPLOT_LABELS = ["A", "B", "C", "D", "E", "F"]


def make_figure_distributions(df):
    """6-panel step-histogram grid comparing TG/ApoB/NonHDLC-derived metrics by genotype."""
    fig, axes = plt.subplots(3, 2, figsize=(22, 24))

    for ax, label, parameter in zip(
        axes.flatten(), DISTRIBUTION_SUBPLOT_LABELS, DISTRIBUTION_PARAMETERS
    ):
        sns.histplot(
            data=df,
            x=parameter,
            hue=GENOTYPE_COL,
            element="step",
            stat="density",
            common_norm=False,
            ax=ax,
        )
        ax.set_title(label, fontsize=14, pad=8, loc="left", fontweight="bold")
        ax.set_xlabel(display_label(parameter))
        ax.set_ylabel("")

    fig.tight_layout()
    fig.subplots_adjust(hspace=0.16)
    return fig


# -------------------------------------------------------------------------
# Figure 1: boxplots of key lipid ratios by TG category and genotype
# -------------------------------------------------------------------------
BOXPLOT_Y_COLUMNS = ["esVLDLC/TG", "NonHDLC/ApoB", "ApoB/eApoB", "esVLDLC/msLDLC", "TG/ApoB", "ApoB"]
BOXPLOT_SUBPLOT_LABELS = ["A", "B", "C", "D", "E", "F"]
BOXPLOT_CATEGORY_TICK_LABELS = ["<180", "180-349", "≥350"]


def make_figure_boxplots(df):
    """3x2 grid of boxplots for key lipid ratios by TG category and genotype."""
    fig, axes = plt.subplots(nrows=3, ncols=2, figsize=(12, 12), sharex=True)

    for i, (ax, y_col) in enumerate(zip(axes.flatten(), BOXPLOT_Y_COLUMNS)):
        sns.boxplot(
            data=df,
            x="TG_category",
            y=y_col,
            hue=GENOTYPE_COL,
            ax=ax,
            showfliers=False,
        )
        ax.set_title(
            BOXPLOT_SUBPLOT_LABELS[i], fontsize=14, pad=8, loc="left", fontweight="bold"
        )
        ax.set_ylabel(display_label(y_col), fontsize=14)
        ax.set_xlabel("")
        ax.set_xticklabels(BOXPLOT_CATEGORY_TICK_LABELS)
        ax.tick_params(axis="x", labelsize=12)

        if i in (4, 5):
            ax.set_xlabel("TG Group (mg/dL)", fontsize=14, labelpad=12)

        if i == 1:
            # Keep a legend only in subplot B; it applies to all subplots.
            handles, labels = ax.get_legend_handles_labels()
            ax.legend(handles, labels, title="", loc="upper left")
        else:
            ax.get_legend().remove()

    fig.tight_layout()
    fig.subplots_adjust(hspace=0.18, wspace=0.18)
    return fig


# -------------------------------------------------------------------------
# Figure 2: per-TG-category histograms, standardised x-axes
# -------------------------------------------------------------------------
TG_GRID_PARAMETERS = ["ApoB/eApoB", "NonHDLC/ApoB", "esVLDLC/msLDLC", "ApoB"]
TG_GRID_ROW_LABELS = ["A", "B", "C", "D"]


def make_figure_tg_category_grid(df):
    """4x3 grid of histograms (rows=metric, columns=TG category) with standardised x-axes.

    Rows A (ApoB/eApoB) and B (NonHDLC/ApoB) use the TG>=350 column's x-axis
    for all three columns; row C (esVLDLC/msLDLC) takes its lower limit from
    the TG<180 column and its upper limit from the TG 180-349 column; row D
    (ApoB) is left on auto-scaling since its columns already line up.
    """
    fig, axes = plt.subplots(4, 3, figsize=(22, 24), sharey="row")
    original_xlim = {}

    for i, parameter in enumerate(TG_GRID_PARAMETERS):
        for j, category in enumerate(TG_CATEGORY_LABELS):
            ax = axes[i, j]
            subset = df[df["TG_category"] == category]
            show_legend = i == 0 and j == 2

            sns.histplot(
                data=subset,
                x=parameter,
                hue=GENOTYPE_COL,
                element="step",
                stat="density",
                common_norm=False,
                ax=ax,
                legend=show_legend,
            )

            original_xlim[i, j] = ax.get_xlim()

            ax.set_title("")
            ax.set_xlabel("")
            ax.set_ylabel("")

            if i == 0:
                ax.set_title(category, fontsize=18, fontweight="bold", pad=10)

            if j == 0:
                ax.set_title(
                    TG_GRID_ROW_LABELS[i], fontsize=14, pad=8, loc="left", fontweight="bold"
                )
                ax.set_ylabel(
                    display_label(parameter), fontsize=18, fontweight="bold", labelpad=10
                )

            if show_legend:
                ax.get_legend().set_title("Genotype")

    for i in range(3):
        if i in (0, 1):
            xlim = original_xlim[i, 2]
        else:
            xlim = (original_xlim[i, 0][0], original_xlim[i, 1][1])
        for j in range(3):
            axes[i, j].set_xlim(xlim)

    fig.tight_layout()
    return fig


# -------------------------------------------------------------------------
# Main
# -------------------------------------------------------------------------
def main():
    REPORTS_DIR.mkdir(parents=True, exist_ok=True)

    df = load_data(DATA_PATH)
    df = add_genotype_columns(df)
    df = clean_cohort(df)
    df = engineer_features(df)

    plt.rcParams.update({"font.size": 18})
    fig_distributions = make_figure_distributions(df)
    fig_distributions.savefig(
        REPORTS_DIR / "SupplementalFigure1_distributions.png", dpi=300, bbox_inches="tight"
    )
    plt.show()

    print_normality_diagnostics(df, NORMALITY_DIAGNOSTIC_COLUMNS)

    build_table2(df).to_csv(REPORTS_DIR / "Table2_Categorical.csv")
    build_table3(df).to_csv(REPORTS_DIR / "Table3_Numerical.csv")
    build_table4_normality(df).to_csv(REPORTS_DIR / "Table4_Normality.csv")
    build_table4(df).to_csv(REPORTS_DIR / "Table4.csv")

    table5 = build_table5_equation_comparison(df)
    table5.to_csv(REPORTS_DIR / "Table5_EquationComparison.csv")
    table5.to_excel(REPORTS_DIR / "Table5_EquationComparison.xlsx")

    plt.rcParams.update({"font.size": 14})
    fig_boxplots = make_figure_boxplots(df)
    fig_boxplots.savefig(REPORTS_DIR / "Figure1_boxplots.png", dpi=300, bbox_inches="tight")
    plt.show()

    fig_tg_grid = make_figure_tg_category_grid(df)
    fig_tg_grid.savefig(
        REPORTS_DIR / "Figure2_tg_category_histograms.png", dpi=300, bbox_inches="tight"
    )
    plt.show()


if __name__ == "__main__":
    main()
