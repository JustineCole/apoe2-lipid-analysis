# Analysis of Lipoprotein Parameters in ApoE2/E2 Individuals

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)

## Overview

This repository contains the analysis code for the manuscript **"Dyslipoproteinemia in Apolipoprotein E2 Homozygotes"** investigating lipid parameter differences between ApoE2/E2 homozygotes and individuals with no ε2 allele, across different triglyceride categories. ApoE2 heterozygotes are excluded from the comparison.

## Citation

If you use this code, please cite our paper:

```bibtex
@article{Cole2026apoe2,
  title={Dyslipoproteinemia in Apolipoprotein E2 Homozygotes},
  author={Justine Cole, Maureen Sampson, Claire Auger, André J. Tremblay, Patrick Couture, Yeganeh Mansourian, Anna Wolska, Rafael Zubiran, James D. Otvos, Ana Cenarro, Fernando Civeira, Alan T. Remaley, Allan D. Sniderman},
  journal={Journal Name},
  year={2026},
  doi={your-doi}
}
```

## Repository Structure

```
.
├── data/
│   ├── raw/              # Raw data files (not included - see Data Availability)
│   └── processed/        # Generated intermediate files
├── reports/              # Generated tables and figures (not included - produced by the script)
├── src/
│   └── data/
│       └── E2E2_progressivedisease.py   # Full analysis pipeline: data prep, tables, and figures
├── .gitignore
├── LICENSE
├── README.md
└── requirements.txt      # Python dependencies
```

## Installation

### Prerequisites
- Python 3.8 or higher

### Using pip

```bash
# Clone the repository
git clone https://github.com/JustineCole/apoe2-lipid-analysis.git
cd apoe2-lipid-analysis

# Create virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt
```

## Data Availability

Due to privacy restrictions, the raw UK Biobank data cannot be shared publicly. To reproduce this analysis:

1. **Apply for UK Biobank data access**: Visit [ukbiobank.ac.uk](https://www.ukbiobank.ac.uk/)
2. **Required data fields**, one row per participant:

   | Column | Contents |
   |---|---|
   | `apoE genotype Combined` | Raw ApoE genotype call (e.g. `E2E2`, `E3E3`, `E3E4`), used to derive ε2 Status |
   | `E2E2 FINAL` | ApoE2/E2 classification used alongside the genotype call (`E2E2` vs `Other`) |
   | `diabetes by lab or doctor` | Diabetes diagnosis flag (1 = yes, 0 = no) |
   | `BMI` | Body mass index (kg/m²) |
   | `Sex` | Sex (`F`/`M`) |
   | `Age` | Age in years |
   | `Current smoker` | Current smoking status flag (1 = yes, 0 = no) |
   | `Apo B` | Apolipoprotein B concentration (mg/dL) |
   | `HDLC` | HDL cholesterol concentration (mg/dL) |
   | `TC` | Total cholesterol concentration (mg/dL) |
   | `TG` | Triglyceride concentration (mg/dL) |
   | `NonHDLC` | Non-HDL cholesterol concentration (mg/dL) |
   | `dLDLC` | Directly measured LDL cholesterol concentration (mg/dL) |

3. **Place the data file** in `data/raw/`:
   - `final_E2E2_dataset.xlsx`: Processed UK Biobank dataset

## Usage

Run the script from anywhere - paths are resolved relative to the script's
own location, not the working directory:

```bash
python src/data/E2E2_progressivedisease.py
```

Each table and figure is also displayed interactively as it's produced, and
written to `reports/`. Table and figure numbers match the manuscript (whose
Table 1 comes from a different source and isn't produced here):
- `Table2_Categorical.csv`, `Table3_Numerical.csv`: Descriptive statistics by genotype (manuscript Tables 2-3)
- `Table4.csv`: Stratified analysis by TG category (manuscript Table 4)
- `Figure1_boxplots.png`, `Figure2_tg_category_histograms.png`: The manuscript's figures
- `SupplementalFigure1_distributions.png`: Exploratory distributions by genotype (published as a supplemental figure)

Plus two supporting outputs considered during analysis but not published in
the manuscript:
- `Table4_Normality.csv`: The normality diagnostics behind Table 4
- `Table5_EquationComparison.csv` / `.xlsx`: A comparison of the two esVLDLC equations (`esVLDLC1`/`esVLDLC2`) by TG category and genotype

## Methods

### Derived Lipid Metrics

1. **msLDLC / msVLDLC** (modified Sampson LDL-C / VLDL-C): Estimated LDL-C and VLDL-C using the modified Sampson equation
2. **esVLDLC** (Estimated VLDL-C): Genotype-specific VLDL-C estimation using the enhanced Sampson equation, with `esVLDLC1`/`esVLDLC2` showing each genotype-specific equation applied to every individual for direct comparison
3. **eApoB** (Estimated Apo B): Calculated apolipoprotein B using a novel NIH equation

### Statistical Analysis

- **Normality testing**: Shapiro-Wilk test with visual inspection
- **Parametric comparisons**: Independent t-tests (Welch's)
- **Non-parametric comparisons**: Mann-Whitney U test
- **Categorical comparisons**: Chi-squared test

### Triglyceride Categories

- Low: TG <180 mg/dL
- Moderate: TG 180-349 mg/dL
- High: TG ≥350 mg/dL

## Dependencies

- pandas >= 1.3.0
- numpy >= 1.21.0
- scipy >= 1.7.0
- matplotlib >= 3.4.0
- seaborn >= 0.11.0

See `requirements.txt` for the complete list.

## Contributing

This code is provided for reproducibility of published results. If you find issues:
1. Open an issue describing the problem
2. Include error messages and your environment details

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact

- **Author**: Justine Cole
- **Email**: justine.cole@gmail.com
- **Institution**: McGill University Health Centre

## Acknowledgments

- UK Biobank participants and team
- This research was supported by an unconditional grant from the Mike and Valeria Rosenbloom Foundation as well as the Intramural Research Program of the NHLBI (HL006275) at the National Institutes of Health

## Manuscript Status

**Submitted**: Currently under review.

---

**Last Updated**: October 2026
