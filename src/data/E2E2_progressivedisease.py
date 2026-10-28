import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as stats

# -------------------------------------------------------
# Read Excel files
# -------------------------------------------------------
UKB_E2df = pd.read_excel("../../data/raw/final_E2E2_dataset.xlsx")
UKB_E2df.head()
UKB_E2df["BP meds at baseline"].value_counts()
UKB_E2df["Age"].describe()
UKB_E2df = UKB_E2df[
    [
        "apoE genotype Combined",
        "E2E2 FINAL",
        "diabetes by lab or doctor",
        "BMI",
        "Sex",
        "Age",
        "Current smoker",
        "Apo B",
        "HDLC",
        "TC",
        "TG",
        "NonHDLC",
        "dLDLC",
    ]
]
UKB_E2df = UKB_E2df.rename(
    columns={"E2E2 FINAL": "ApoE2/E2", "Apo B": "ApoB", "dLDLC": "direct LDLC"}
)
UKB_E2df["ApoE2/E2"] = UKB_E2df["ApoE2/E2"].replace(
    {"E2E2": "ε2ε2", "Other": "Not ε2ε2"}
)

# -------------------------------------------------------
# Create ApoE2 Heterozygote and ε2 Status columns
# -------------------------------------------------------
## Create ApoE2 Heterozygote column
# 1. Define the conditions
# Contains 'E2' AND does NOT contain 'E2E2'
mask = (UKB_E2df["apoE genotype Combined"].str.contains("E2")) & (
    ~UKB_E2df["apoE genotype Combined"].str.contains("E2E2")
)

# 2. Assign 1 if True, 0 if False
UKB_E2df["E2 Heterozygote"] = np.where(mask, 1, 0)
UKB_E2df["E2 Heterozygote"].value_counts()
## Create ε2 Status column
# 1. Define your conditions
conditions = [
    UKB_E2df["apoE genotype Combined"].str.contains("E2E2"),  # Two copies
    UKB_E2df["apoE genotype Combined"].str.contains(
        "E2"
    ),  # One copy (E2E2 is already caught above)
]

# 2. Define the corresponding labels
choices = ["ε2ε2", "ε2"]

# 3. Apply the logic (default is 'No ε2')
UKB_E2df["ε2 Status"] = np.select(conditions, choices, default="No ε2")
UKB_E2df["ε2 Status"].value_counts()

# -------------------------------------------------------
# Clean Dataframe
# -------------------------------------------------------
## Remove heterozygotes
UKB_E2df = UKB_E2df[UKB_E2df["E2 Heterozygote"] == 0]

UKB_E2_df = UKB_E2df[
    [
        "ε2 Status",
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
]


# -------------------------------------------------------
# Feature Creation
# -------------------------------------------------------
## Create msLDLC
UKB_E2_df["msLDLC"] = (
    UKB_E2_df["NonHDLC"]
    - (UKB_E2_df["TG"] / 8.37)
    - (UKB_E2_df["NonHDLC"] * UKB_E2_df["TG"] / 2640)
    + (UKB_E2_df["TG"] ** 2 / 17400)
)

## Create msVLDLC
UKB_E2_df["msVLDLC"] = (
    (UKB_E2_df["TG"] / 8.37)
    + (UKB_E2_df["NonHDLC"] * UKB_E2_df["TG"] / 2640)
    - (UKB_E2_df["TG"] ** 2 / 17400)
)


# Create esVLDLC using alternative equations for non-ApoE2/E2 individuals and ApoE2/E2 individuals:
UKB_E2_df["esVLDLC"] = np.select(
    condlist=[UKB_E2_df["ε2 Status"] == "No ε2", UKB_E2_df["ε2 Status"] == "ε2ε2"],
    choicelist=[
        (UKB_E2_df["NonHDLC"] / 3.81)
        - (UKB_E2_df["HDLC"] / 8.93)
        + (UKB_E2_df["TG"] / 7.73)
        + (UKB_E2_df["NonHDLC"] * UKB_E2_df["TG"] / 2050)
        - (UKB_E2_df["TG"] ** 2 / 13300)
        - (UKB_E2_df["ApoB"] / 2.49)
        - (UKB_E2_df["ApoB"] * UKB_E2_df["TG"] / 3550)
        + 7.46,
        (0.001 * UKB_E2_df["NonHDLC"] * UKB_E2_df["TG"])
        - (0.526 * UKB_E2_df["ApoB"])
        + (0.318 * UKB_E2_df["TG"])
        - (0.0004 * UKB_E2_df["TG"] ** 2)
        + 38.025,
    ],
    default=np.nan,
)


UKB_E2_df["esVLDLC"] = np.where(UKB_E2_df["esVLDLC"] < 0, 0, UKB_E2_df["esVLDLC"])

# Create esVLDLC1 and esVLDLC2: the same two equations above, but each applied
# to every row (not just the genotype they were derived for), so they can be
# compared directly instead of only taking the genotype-matched value.
UKB_E2_df["esVLDLC1"] = (  # equation used for "No ε2"
    (UKB_E2_df["NonHDLC"] / 3.81)
    - (UKB_E2_df["HDLC"] / 8.93)
    + (UKB_E2_df["TG"] / 7.73)
    + (UKB_E2_df["NonHDLC"] * UKB_E2_df["TG"] / 2050)
    - (UKB_E2_df["TG"] ** 2 / 13300)
    - (UKB_E2_df["ApoB"] / 2.49)
    - (UKB_E2_df["ApoB"] * UKB_E2_df["TG"] / 3550)
    + 7.46
)
UKB_E2_df["esVLDLC1"] = np.where(UKB_E2_df["esVLDLC1"] < 0, 0, UKB_E2_df["esVLDLC1"])

UKB_E2_df["esVLDLC2"] = (  # equation used for "ε2ε2"
    (0.001 * UKB_E2_df["NonHDLC"] * UKB_E2_df["TG"])
    - (0.526 * UKB_E2_df["ApoB"])
    + (0.318 * UKB_E2_df["TG"])
    - (0.0004 * UKB_E2_df["TG"] ** 2)
    + 38.025
)
UKB_E2_df["esVLDLC2"] = np.where(UKB_E2_df["esVLDLC2"] < 0, 0, UKB_E2_df["esVLDLC2"])

# Create eApoB
UKB_E2_df["eApoB"] = (
    (UKB_E2_df["NonHDLC"] / 1.90)
    - (UKB_E2_df["TG"] / 22.5)
    + (2.75 * (UKB_E2_df["NonHDLC"] / UKB_E2_df["HDLC"]))
    + (UKB_E2_df["TG"] * UKB_E2_df["NonHDLC"] / 12200)
    - (UKB_E2_df["TG"] ** 2 / 50400)
    + 13.0
)


# Create ratios
UKB_E2_df["esVLDLC/TG"] = UKB_E2_df["esVLDLC"] / UKB_E2_df["TG"]
UKB_E2_df["esVLDLC/msLDLC"] = UKB_E2_df["esVLDLC"] / UKB_E2_df["msLDLC"]
UKB_E2_df["msLDLC/dLDLC"] = UKB_E2_df["msLDLC"] / UKB_E2_df["direct LDLC"]
UKB_E2_df["ApoB/eApoB"] = UKB_E2_df["ApoB"] / UKB_E2_df["eApoB"]
UKB_E2_df["nonHDLC/ApoB"] = UKB_E2_df["NonHDLC"] / UKB_E2_df["ApoB"]
UKB_E2_df["TG/ApoB"] = UKB_E2_df["TG"] / UKB_E2_df["ApoB"]


# -------------------------------------------------------
# Display label helper
# -------------------------------------------------------
# Column/variable names use short internal abbreviations (e.g. "msLDLC").
# This maps them to publication-style labels for figure axes and table
# headers/indices without renaming the underlying DataFrame columns.
DISPLAY_LABEL_REPLACEMENTS = [
    ("esVLDLC", "VLDL-C"),
    ("msLDLC", "LDL-C"),
    ("NonHDLC", "Non-HDL-C"),
    ("nonHDLC", "Non-HDL-C"),
    ("HDLC", "HDL-C"),
]


def display_label(name):
    label = name
    for old, new in DISPLAY_LABEL_REPLACEMENTS:
        label = label.replace(old, new)
    return label


# -------------------------------------------------------
# Explore Data - Table 1
# -------------------------------------------------------

# Table 1 Normality Tests
## Histograms to visualize distributions

"""for col in [
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
    "nonHDLC/ApoB",
    "TG/ApoB",
]:
    plt.figure(figsize=(10, 6))
    sns.histplot(
        data=UKB_E2_df,
        x=col,
        hue="ε2 Status",
        stat="density",
        common_norm=False,
        multiple="layer",
        kde=True,
        bins=30,
        alpha=0.4,
    )
    plt.title(f"Distribution of {col} by ApoE2/E2 Status")
    plt.xlabel(col)
    plt.ylabel("Density")
    plt.legend(title="Genotype", labels=["ε2ε2", "No ε2"])
    plt.show()"""


# -------------------------------------------------------
# Create histplots of "TG", "ApoB", "NonHDLC", "nonHDLC/ApoB", "ApoB/eApoB", "esVLDLC/TG" by ApoE2/E2 status
# -------------------------------------------------------
plt.rcParams.update({"font.size": 18})
LEGEND_TITLE = "Genotype"
PARAMETERS = ["TG", "ApoB", "NonHDLC", "nonHDLC/ApoB", "ApoB/eApoB", "esVLDLC/TG"]
# Create 4x3 grid (4 rows for parameters, 3 columns for categories)
fig, axes = plt.subplots(3, 2, figsize=(22, 24))
all_axes = axes.flatten()


SUBPLOT_LABELS = ["A", "B", "C", "D", "E", "F"]

# --- Plotting Loop ---
for i, (ax, y_col) in enumerate(zip(all_axes, PARAMETERS)):
    # Create the boxplot without outliers
    sns.histplot(
        data=UKB_E2_df,
        x=PARAMETERS[i],
        hue="ε2 Status",
        element="step",
        stat="density",
        common_norm=False,
        ax=ax,
    )

    ax.set_title(SUBPLOT_LABELS[i], fontsize=14, pad=8, loc="left", fontweight="bold")
    ax.set_xlabel(display_label(y_col))
    ax.set_ylabel("")

plt.tight_layout()
plt.subplots_adjust(hspace=0.16)
plt.show()


## Check for Skewness and Kurtosis and Perform Shapiro-Wilk test for normality
for col in [
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
    "nonHDLC/ApoB",
    "TG/ApoB",
]:
    # --- Prepare Data ---
    group1 = UKB_E2_df[UKB_E2_df["ε2 Status"] == "ε2ε2"][col].dropna()
    group2 = UKB_E2_df[UKB_E2_df["ε2 Status"] == "No ε2"][col].dropna()

    # --- Skewness and Kurtosis ---
    skew1 = stats.skew(group1)
    kurt1 = stats.kurtosis(group1)  # Fisher's kurtosis (Normal = 0)
    skew2 = stats.skew(group2)
    kurt2 = stats.kurtosis(group2)

    # --- Print Results ---
    print(f"--- Results for **{col}** ---")

    print(f"**Descriptive Statistics (Skew/Kurtosis):**")
    print(f" Group 1 (ε2ε2): Skew={skew1:.4f}, Kurtosis (Fisher)={kurt1:.4f}")
    print(f" Group 2 (No ε2): Skew={skew2:.4f}, Kurtosis (Fisher)={kurt2:.4f}")

    # --- Shapiro-Wilk Test (from your original code) ---
    w_stat1, p_value1 = stats.shapiro(group1)
    w_stat2, p_value2 = stats.shapiro(group2)
    print(f"**Shapiro-Wilk Test:**")
    print(f" Group 1: W={w_stat1:.4f}, p={p_value1:.4f}")
    print(f" Group 2: W={w_stat2:.4f}, p={p_value2:.4f}\n")

### -------------------------------------------------------
### Create Table 1b - Numerical Variables
### -------------------------------------------------------

# Define non-parametric columns
non_parametric = ["TG", "esVLDLC", "esVLDLC/TG", "esVLDLC/msLDLC", "TG/ApoB"]

# Create the base table
table1 = (
    UKB_E2_df[
        [
            "Age",
            "ApoB",
            "HDLC",
            "TG",
            "NonHDLC",
            "msLDLC",
            "BMI",
            "ε2 Status",
            "esVLDLC",
            "esVLDLC/TG",
            "esVLDLC/msLDLC",
            "ApoB/eApoB",
            "nonHDLC/ApoB",
            "TG/ApoB",
        ]
    ]
    .groupby("ε2 Status")
    .describe()
    .T
)

# Calculate p-values
p_values = []
variables = [col for col in table1.index.get_level_values(0).unique()]

for var in variables:
    # Get data for each group
    group1 = UKB_E2_df[UKB_E2_df["ε2 Status"] == UKB_E2_df["ε2 Status"].unique()[0]][
        var
    ].dropna()
    group2 = UKB_E2_df[UKB_E2_df["ε2 Status"] == UKB_E2_df["ε2 Status"].unique()[1]][
        var
    ].dropna()

    # Choose appropriate test
    if var in non_parametric:
        # Mann-Whitney U test for non-parametric
        _, p_val = stats.mannwhitneyu(group1, group2, alternative="two-sided")
    else:
        # Independent t-test for parametric
        _, p_val = stats.ttest_ind(group1, group2, equal_var=False)

    # Add p-value for each statistic row (mean, std, min, etc.)
    num_stats = len(table1.loc[var])
    p_values.extend([p_val] * num_stats)

# Add p-values to table
table1["p-value"] = p_values

table1 = table1.rename(index=display_label, level=0)
table1.to_csv("../../reports/Table1b_Numerical.csv")

### -------------------------------------------------------
### Create Table 1a - Categorical Variables
### -------------------------------------------------------
# 1. Define the columns we want to summarize
cols_to_summary = ["Sex", "Current smoker", "diabetes by lab or doctor"]


# 2. Create a helper function to get counts and percentages
def summarize_group(group):
    results = {}
    for col in cols_to_summary:
        # Define what counts as a "positive" case for each column
        if col == "Sex":
            count = (group[col] == "F").sum()
        else:
            # Using str comparison if your data is cast to string, or int if not
            count = (group[col] == 1).sum()

        percentage = (count / len(group)) * 100
        results[(col, "Count")] = count
        results[(col, "Percentage")] = f"{percentage:.0f}%"
    return pd.Series(results)


# 3. Group by ApoE2/E2 and apply the function
summary_df = UKB_E2_df.groupby("ε2 Status").apply(summarize_group).T

### Perform chi-squared tests to compare proportions between ApoE2/E2 and not ApoE2/E2 groups
Categorical_Variable = []
Chi2_statistic = []
P_value = []
for col in cols_to_summary:
    contingency_table = pd.crosstab(UKB_E2_df[col], UKB_E2_df["ε2 Status"])
    chi2, p, df, expected = stats.chi2_contingency(contingency_table)
    Categorical_Variable.append(col)
    Chi2_statistic.append(chi2)
    P_value.append(p)

stats2 = pd.DataFrame(
    {
        "Variable": Categorical_Variable,
        "chi2-statistic": Chi2_statistic,
        "p-value": P_value,
    }
)

# 4. Add p-values to the summary table
# Create a MultiIndex for the p-value column to match summary_df structure
p_value_series = pd.Series(
    [
        stats2.loc[stats2["Variable"] == col, "p-value"].values[0]
        for col in cols_to_summary
    ],
    index=pd.MultiIndex.from_product([cols_to_summary, ["p-value"]]),
)

# Concatenate with the summary dataframe
summary_df_with_p = pd.concat([summary_df, p_value_series.to_frame("p-value")], axis=0)

# Sort the index to group each variable's statistics together
summary_df_with_p = summary_df_with_p.sort_index(level=0)

# 5. Save the summary dataframe with p-values to CSV
summary_df_with_p.to_csv("../../reports/Table1a_Categorical.csv")

# -------------------------------------------------------
# Table 2
# -------------------------------------------------------
# -------------------------------------------------------
# Create TG categories
# -------------------------------------------------------
UKB_E2_df["TG_category"] = pd.cut(
    UKB_E2_df["TG"],
    bins=[0, 180, 350, np.inf],
    labels=["TG <180 mg/dL", "TG 180-349 mg/dL", "TG ≥350 mg/dL"],
)
UKB_E2_df["TG_category"].value_counts()


### Test for normality using Shapiro-Wilk test
# ==============================================================================
# 1. Define Metrics and Groups (Adjust this list if your column names differ)
# ==============================================================================
metrics = [
    "ApoB",
    "TG",
    "NonHDLC",
    "HDLC",
    "esVLDLC",
    "esVLDLC/TG",
    "esVLDLC/msLDLC",
    "ApoB/eApoB",
    "nonHDLC/ApoB",
    "TG/ApoB",
]
tg_categories = ["TG <180 mg/dL", "TG 180-349 mg/dL", "TG ≥350 mg/dL"]


# ==============================================================================
# 3. Perform Shapiro-Wilk Test
# ==============================================================================
results = []

for metric in metrics:
    row = {"Metric": metric}

    for tg_cat in tg_categories:
        # Filter by TG Category
        df_group_tg = UKB_E2_df[UKB_E2_df["TG_category"] == tg_cat]

        for genotype in UKB_E2_df["ε2 Status"].unique():
            # Filter by Genotype (1 or 0)
            data = df_group_tg[df_group_tg["ε2 Status"] == genotype][metric].dropna()

            # Subgroup label for the result table
            col_prefix = f"({tg_cat}, {genotype})"

            # --- Skewness and Kurtosis Calculation ---
            # Skewness: measures asymmetry of the distribution
            skewness = stats.skew(data)
            # Kurtosis: measures the 'tailedness' of the distribution.
            # Use 'fisher' for standard comparison (Normal distribution has kurtosis of 0)
            kurtosis = stats.kurtosis(data, fisher=True)

            # --- Shapiro-Wilk Test (assuming len(data) >= 3 is always true or handled by try/except) ---
            try:
                w_stat, p_value = stats.shapiro(data)
            except Exception:
                w_stat, p_value = np.nan, np.nan

            # --- Store Results ---
            row[f"{col_prefix}_Skewness"] = skewness
            row[f"{col_prefix}_Kurtosis"] = kurtosis
            row[f"{col_prefix}_W"] = w_stat
            row[f"{col_prefix}_p"] = p_value

    results.append(row)

# Convert results list to DataFrame
df_results = pd.DataFrame(results).set_index("Metric")

# Create the MultiIndex for the columns for a clear structure
column_tuples = []
for tg_cat in tg_categories:
    for genotype_name in ["ε2ε2", "No ε2"]:
        # Order the statistics for better readability: Skewness, Kurtosis, W, p
        column_tuples.append((tg_cat, f"{genotype_name} Skewness"))
        column_tuples.append((tg_cat, f"{genotype_name} Kurtosis"))
        column_tuples.append((tg_cat, f"{genotype_name} W"))
        column_tuples.append((tg_cat, f"{genotype_name} p"))

df_results.columns = pd.MultiIndex.from_tuples(
    column_tuples, names=["TG Group", "Statistic"]
)

# Rename the final DataFrame for clarity on its contents
df_normality_stats = df_results


# ==============================================================================
# 4. Apply Formatting
# ==============================================================================

# Format p-values and other statistics
df_normality = df_normality_stats.copy()

for col in df_normality.columns:
    statistic_name = col[1]

    if statistic_name.endswith("p"):
        # --- Format p-values (ending in 'p') ---
        df_normality[col] = df_normality[col].apply(
            lambda p: (
                f"{p:.4f}"
                if pd.notna(p) and p >= 0.0001
                else "<0.0001" if pd.notna(p) and p < 0.0001 else ""
            )
        )
    # Check for W, Skewness, or Kurtosis
    elif (
        statistic_name.endswith("W")
        or "Skewness" in statistic_name
        or "Kurtosis" in statistic_name
    ):
        # --- Format W-statistic, Skewness, and Kurtosis to 3 decimal places ---
        df_normality[col] = df_normality[col].round(3).astype(str).replace("nan", "")

df_normality.index = df_normality.index.map(display_label)
df_normality.to_csv("../../reports/Table2_Normality_Stats.csv")

### Create Table 2
# ==============================================================================
# 1. Define the necessary columns and statistics (Non-Parametric)
# ==============================================================================
all_metrics = [
    "TG",
    "ApoB",
    "NonHDLC",
    "HDLC",
    "esVLDLC",
    "esVLDLC/TG",
    "esVLDLC/msLDLC",
    "msLDLC/dLDLC",
    "ApoB/eApoB",
    "nonHDLC/ApoB",
    "TG/ApoB",
]

# The final row MultiIndex structure
row_index_tuples = [
    ("", "count"),
    *[(m, "median") for m in all_metrics],
    *[(m, "IQR1") for m in all_metrics],
    *[(m, "IQR3") for m in all_metrics],
]

# NOTE: We keep this structure for initial construction, then re-sort later
final_row_index = pd.MultiIndex.from_tuples(
    row_index_tuples, names=["Metric", "Statistic"]
)

# The final column MultiIndex structure
tg_categories = ["TG <180 mg/dL", "TG 180-349 mg/dL", "TG ≥350 mg/dL"]
genotype_cols = ["ε2ε2", "No ε2", "p-value"]
column_index_tuples = [
    (tg_cat, col) for tg_cat in tg_categories for col in genotype_cols
]
final_column_index = pd.MultiIndex.from_tuples(
    column_index_tuples, names=["TG Group", "Genotype"]
)


# ==============================================================================
# 3. Function to calculate statistics and p-value (Mann-Whitney U)
# ==============================================================================
def calculate_stats_and_pvalue(df_group, metric):
    """
    Calculates summary statistics (Median/IQR) and the p-value (Mann-Whitney U test)
    for a given metric within a TG group.
    """
    stats_dict = {}

    # Split the data by genotype
    e2e2_data = df_group[df_group["ε2 Status"] == "ε2ε2"][metric].dropna()
    not_e2e2_data = df_group[df_group["ε2 Status"] == "No ε2"][metric].dropna()

    # Calculate statistics (Median and IQR for all)
    for group_name, data in [("ε2ε2", e2e2_data), ("No ε2", not_e2e2_data)]:
        stats_dict[("median", group_name)] = data.median()
        q1, q3 = data.quantile([0.25, 0.75]).tolist()
        stats_dict[("IQR1", group_name)] = q1
        stats_dict[("IQR3", group_name)] = q3

    # Statistical Testing Logic: Mann-Whitney U test for all
    if len(e2e2_data) < 3 or len(not_e2e2_data) < 3:
        p_value = np.nan
    else:
        _, p_value = stats.mannwhitneyu(
            e2e2_data, not_e2e2_data, alternative="two-sided"
        )

    # Store p-value: it applies to the 'median' row for all variables
    stats_dict[("median", "p-value")] = p_value

    return stats_dict


# ==============================================================================
# 4. Generate the Final DataFrame
# ==============================================================================
final_data = {}
for metric_stat_tuple in final_row_index:
    final_data[metric_stat_tuple] = {}

# 4.1 Handle the special 'count' row first
count_row_key = ("", "count")
for tg_cat in tg_categories:
    df_group = UKB_E2_df[UKB_E2_df["TG_category"] == tg_cat]
    count_e2e2 = len(df_group[df_group["ε2 Status"] == "ε2ε2"])
    count_not_e2e2 = len(df_group[df_group["ε2 Status"] == "No ε2"])

    final_data[count_row_key][(tg_cat, "ε2ε2")] = count_e2e2
    final_data[count_row_key][(tg_cat, "No ε2")] = count_not_e2e2
    final_data[count_row_key][(tg_cat, "p-value")] = np.nan

# 4.2 Loop through all metrics (now treated uniformly)
for tg_cat in tg_categories:
    df_group = UKB_E2_df[UKB_E2_df["TG_category"] == tg_cat]

    for metric in all_metrics:
        stats_results = calculate_stats_and_pvalue(df_group, metric)
        for (stat, col_name), value in stats_results.items():
            if (metric, stat) in final_data:
                final_data[(metric, stat)][(tg_cat, col_name)] = value


# 5. Restructure the final data into the DataFrame
df_final = pd.DataFrame(index=final_row_index, columns=final_column_index)

for (metric, stat), col_data in final_data.items():
    for (tg_group, genotype_col), value in col_data.items():
        df_final.loc[(metric, stat), (tg_group, genotype_col)] = value

df_final = df_final.fillna(np.nan)


# ==============================================================================
# 6. Apply Formatting, Clean Up, RENAME, and SORT for Display
# ==============================================================================
# 1. Rename df_display to Table2_df
Table2_df = df_final.copy()
Table2_df = Table2_df.apply(pd.to_numeric, errors="ignore")

# 2. Sort by 'Metric' to group median/IQR rows together
# Note: The count row (Metric="") must be manually handled to stay at the top.
count_row = Table2_df.loc[("", "count")].to_frame().T
data_rows = Table2_df.drop(index=[("", "count")])

# Sort the remaining data rows by the outer index level ('Metric')
data_rows_sorted = data_rows.sort_index(level="Metric", kind="stable")

# Recombine the data frame with the count row at the top
Table2_df = pd.concat([count_row, data_rows_sorted])


# 3. Standard Formatting
data_cols = [col for col in Table2_df.columns if col[1] in ["ε2ε2", "No ε2"]]
p_value_cols = [col for col in Table2_df.columns if col[1] == "p-value"]
count_index = ("", "count")

# Format data columns to 2 decimal places
Table2_df[data_cols] = Table2_df[data_cols].round(2)

# Format the 'count' row as integers and clear its p-value
if count_index in Table2_df.index:
    Table2_df.loc[count_index, data_cols] = (
        Table2_df.loc[count_index, data_cols].astype(float).round(0).astype(int)
    )
    Table2_df.loc[count_index, p_value_cols] = np.nan


# Function to format p-values
def format_p_value(p):
    if pd.isna(p):
        return ""
    elif p < 0.0001:
        return "<0.0001"
    else:
        return f"{p:.4f}"


Table2_df[p_value_cols] = Table2_df[p_value_cols].applymap(format_p_value)

# Replace any remaining NaNs with an empty string
Table2_df = Table2_df.fillna("")

Table2_df = Table2_df.rename(index=display_label, level=0)
Table2_df.to_csv("../../reports/Table2.csv")


# -------------------------------------------------------
# Table 3: msLDLC, msVLDLC, esVLDLC1, esVLDLC2 by TG category and ε2 Status
# -------------------------------------------------------
table3_cols = ["msLDLC", "msVLDLC", "esVLDLC1", "esVLDLC2"]

table3 = UKB_E2_df.groupby(["TG_category", "ε2 Status"])[table3_cols].describe()

table3.to_csv("../../reports/Table3_LDLC_VLDLC_Comparison.csv")
table3.to_excel("../../reports/Table3_LDLC_VLDLC_Comparison.xlsx")


# -------------------------------------------------------
# Figure 1: Boxplots of esVLDLC/TG, esVLDLC/msLDLC, pApoB/mApoB and nonHDLC/ApoB  by TG category and ApoE2/E2 status
# -------------------------------------------------------
# --- Configuration & Setup ---
subplot_titles = [
    r"A",
    r"B",
    r"C",
    r"D",
    r"E",
    r"F",
]
category_labels = ["<180", "180-349", "≥350"]
y_labels = [
    r"esVLDLC/TG",
    r"nonHDLC/ApoB",
    r"ApoB/eApoB",
    r"esVLDLC/msLDLC",
    r"TG/ApoB",
    r"ApoB",
]

# Create a figure with a 3x2 grid of subplots
fig, axes = plt.subplots(nrows=3, ncols=2, figsize=(12, 12), sharex=True)
all_axes = axes.flatten()


# --- Plotting Loop ---
for i, (ax, y_col) in enumerate(zip(all_axes, y_labels)):
    # Create the boxplot without outliers
    sns.boxplot(
        data=UKB_E2_df,
        x="TG_category",
        y=y_col,
        hue="ε2 Status",
        ax=ax,
        showfliers=False,
    )

    # Set axis properties
    ax.set_title(subplot_titles[i], fontsize=14, pad=8, loc="left", fontweight="bold")
    ax.set_ylabel(display_label(y_labels[i]), fontsize=14)
    ax.set_xlabel("")
    ax.set_xticklabels(category_labels)
    ax.tick_params(axis="x", labelsize=12)

    # Add X label only to bottom subplots
    if i in [4, 5]:
        ax.set_xlabel("TG Group (mg/dL)", fontsize=14, labelpad=12)

    # Handle legend
    if i == 1:
        # Keep legend only in subplot B
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles, labels, title="", loc="upper left")
    else:
        # Remove legend from other subplots
        ax.get_legend().remove()

# Adjust layout and show plot
plt.tight_layout()
plt.subplots_adjust(hspace=0.18, wspace=0.18)
plt.show()


plt.savefig("ukb_boxplots.png", dpi=300, bbox_inches="tight")


# -------------------------------------------------------
# For each TG category, create histplots of Non-HDLC/ApoB, pApoB/mApoB, esVLDLC/msLDLC and ApoB by ApoE2/E2 status
# -------------------------------------------------------
plt.rcParams.update({"font.size": 14})
LEGEND_TITLE = "Genotype"
CATEGORIES = ["TG <180 mg/dL", "TG 180-349 mg/dL", "TG ≥350 mg/dL"]
PARAMETERS = ["ApoB/eApoB", "nonHDLC/ApoB", "esVLDLC/msLDLC", "ApoB"]

# Create 4x3 grid (4 rows for parameters, 3 columns for categories)
fig, axes = plt.subplots(4, 3, figsize=(22, 24), sharey="row")

original_xlim = {}  # (i, j) -> each subplot's auto-scaled x-axis limits

for i in range(4):  # rows (parameters)
    for j in range(3):  # columns (categories)
        ax = axes[i, j]
        parameter = PARAMETERS[i]
        category = CATEGORIES[j]

        subset = UKB_E2_df[(UKB_E2_df["TG_category"] == category)]

        # Determine if the legend should be shown (top-right subplot)
        show_legend_bool = i == 0 and j == 2

        sns.histplot(
            data=subset,
            x=parameter,
            hue="ε2 Status",
            element="step",
            stat="density",
            common_norm=False,
            ax=ax,
            legend=show_legend_bool,
        )

        original_xlim[i, j] = ax.get_xlim()

        # Remove individual subplot titles
        ax.set_title("")
        ax.set_xlabel("")
        ax.set_ylabel("")

        # Add column titles (only for top row)
        if i == 0:
            ax.set_title(category, fontsize=18, fontweight="bold", pad=10)

        # Add row titles (only for first column)
        if j == 0:
            ax.set_title(
                ["A", "B", "C", "D"][i],
                fontsize=14,
                pad=8,
                loc="left",
                fontweight="bold",
            )
            ax.set_ylabel(
                display_label(parameter), fontsize=18, fontweight="bold", labelpad=10
            )

        # Legend customization
        if show_legend_bool:
            legend = ax.get_legend()
            legend.set_title(LEGEND_TITLE)

# Standardise x-axes across categories (columns) within a row, using each
# subplot's own auto-scaled limits captured above: rows A (ApoB/eApoB) and
# B (nonHDLC/ApoB) take the 3rd column's (TG>=350) limits for all columns;
# row C (esVLDLC/msLDLC) takes the lower limit from column 1 (TG<180) and the
# upper limit from column 2 (TG 180-349); row D (ApoB) is left as-is.
for i in range(3):
    if i in (0, 1):
        xlim = original_xlim[i, 2]
    else:
        xlim = (original_xlim[i, 0][0], original_xlim[i, 1][1])
    for j in range(3):
        axes[i, j].set_xlim(xlim)

plt.tight_layout()
plt.show()
