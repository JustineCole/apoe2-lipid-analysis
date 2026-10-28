# Analysis of Lipoprotein Parameters in ApoE2/E2 Individuals

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)

## Overview

This repository contains the analysis code for the manuscript **"Dyslipidemia in Apolipoprotein E2 Homozygotes"** investigating lipid parameter differences between ApoE2/E2 and non-ApoE2/E2 individuals across different triglyceride categories.

## Citation

If you use this code, please cite our paper:

```bibtex
@article{Cole2026apoe2,
  title={Dyslipidemia in Apolipoprotein E2 Homozygotes},
  author={Justine Cole, Maureen Sampson, Ana Cenarro, Yeganeh Mansourian, Anna Wolska, Claire Auger, Jim Otvos, Fernando Civeira, André Tremblay, Patrick Couture, Rafael Zubiran, Alan T. Remaley, Allan D. Sniderman},
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
2. **Required data fields**: `apoE genotype Combined`, `E2E2 FINAL`, `diabetes by lab or doctor`, `BMI`, `Sex`, `Age`, `Current smoker`, `Apo B`, `HDLC`, `TC`, `TG`, `NonHDLC`, `dLDLC`
3. **Place the data file** in `data/raw/`:
   - `final_E2E2_dataset.xlsx`: Processed UK Biobank dataset

## Usage

Run the full analysis from the repository root:

```bash
cd src/data
python E2E2_progressivedisease.py
```

This produces, in `reports/`:
- `Table1a_Categorical.csv`, `Table1b_Numerical.csv`: Descriptive statistics by genotype
- `Table2.csv`, `Table2_Normality_Stats.csv`: Stratified analysis by TG category
- `Table3_LDLC_VLDLC_Comparison.csv`: Comparison of LDL-C/VLDL-C estimation equations by TG category and genotype

and displays the histogram and boxplot figures used in the manuscript.

## Methods

### Derived Lipid Metrics

1. **msLDLC / msVLDLC** (modified Sampson LDL-C / VLDL-C): Estimated LDL-C and VLDL-C using the modified Sampson equation
2. **esVLDLC** (Estimated VLDL-C): Genotype-specific VLDL-C estimation using the enhanced Sampson equation, with `esVLDLC1`/`esVLDLC2` showing each genotype-specific equation applied to every individual for direct comparison
3. **eApoB** (Estimated Apo B): Calculated apolipoprotein B using the Sampson equation

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

**Published**: [Journal link and DOI when published]

---

**Last Updated**: October 2026
