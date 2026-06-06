# Heart Rate Detection from Ballistocardiogram (BCG) Signals
### A Julia Implementation with Statistical Evaluation

> **Course:** Advanced Statistics  
> **Reference Implementation:** [Code Ocean Capsule #1398208](https://codeocean.com/capsule/1398208/tree)  
> **Reference Dataset:** [Scientific Data (2024)](https://www.nature.com/articles/s41597-024-03950-5)  
> **Reference Paper:** A ballistocardiogram dataset with reference sensor signals in long-term natural sleep environments

---

## Table of Contents

- [Project Overview](#project-overview)
- [Repository Structure](#repository-structure)
- [Dependencies](#dependencies)
- [Installation](#installation)
- [Dataset](#dataset)
- [Usage Instructions](#usage-instructions)
- [Pipeline Description](#pipeline-description)
- [Evaluation Metrics](#evaluation-metrics)
- [Results](#results)
- [Known Limitations](#known-limitations)
- [References](#references)

---

## Project Overview

This project implements a **signal processing pipeline in Julia** to estimate Heart Rate (HR) from Ballistocardiogram (BCG) signals recorded during overnight sleep sessions. The implementation is a manual conversion of the original Python codebase published on Code Ocean, with several bugs identified and corrected during translation.

The pipeline is evaluated against ECG-derived reference RR interval signals across **22 patients** using seven statistical metrics.

### Source & Reference

| | Link |
|---|---|
| **CodeOcean Capsule** (original Python implementation) | https://codeocean.com/capsule/1398208/tree |
| **Reference Publication** (dataset paper) | https://doi.org/10.1038/s41597-024-03950-5 |

---

## Repository Structure

```
├── Code/
│   ├── mainEntry.jl                          # Main entry point — runs full pipeline
│   ├── data_preprocessing.jl                 # Signal processing pipeline functions
│   ├── metrics_computation.jl                # MAE, RMSE, MAPE computation
│   ├── pearson_relation_coeff.jl             # Pearson correlation coefficient
│   ├── p_value_calculations.jl               # p-value via t-distribution (manual)
│   └── regression_plot_and_bland_altman_plot.jl  # Visualisation functions
│
├── Deliverables/
│   ├── IEEE_Paper.docx                       # Two-page IEEE format report
│   ├── Presentation.pptx                     # Project slides
│   └── Demo_Recording.mp4                    # VMware demo walkthrough
│
├── Plots/
│   └── regression_bland_altman.png           # Regression & Bland–Altman plots
│
└── README.md                                 # This file
```

---

## Dependencies

### Julia Version
```
Julia ≥ 1.9
```

### Required Julia Packages

| Package | Purpose |
|---------|---------|
| `FFTW` | Fast Fourier Transform for wavelet decomposition |
| `DSP` | Digital signal processing (Chebyshev filters, filtfilt) |
| `Statistics` | Mean, standard deviation |
| `LinearAlgebra` | Matrix operations for filter initial conditions |
| `CSV` | Reading patient data files |
| `DataFrames` | Tabular data handling |
| `Dates` | Parsing RR reference file timestamps |
| `Plots` | Regression and Bland–Altman visualisation |

---

## Installation

### 1. Install Julia

Download and install Julia from [https://julialang.org/downloads/](https://julialang.org/downloads/)

### 2. Install Required Packages

Open a Julia REPL and run:

```julia
using Pkg

Pkg.add([
    "FFTW",
    "DSP",
    "Statistics",
    "LinearAlgebra",
    "CSV",
    "DataFrames",
    "Dates",
    "Plots"
])
```

---

## Dataset

This project uses the dataset described in:

> Chen et al., *"A ballistocardiogram dataset with reference sensor signals in long-term natural sleep environments"*, Scientific Data, 2024.  
> DOI: [10.1038/s41597-024-03950-5](https://doi.org/10.1038/s41597-024-03950-5)

### Expected Folder Structure

```
modified_signals/
├── 01/
│   ├── BCG/
│   │   └── 01_20231104_BCG.csv       ← BCG signal file
│   └── RR/
│       └── 01_20231104_RR.csv        ← Reference RR intervals
├── 02/
│   ├── BCG/
│   └── RR/
└── ...
```

### BCG File Format

```
BCG,Timestamp,fs
-86,1699022112866,140        ← Row 2: first value, start timestamp (ms), fs (Hz)
-90                          ← Rows 3+: BCG values only
-91
...
```

### RR Reference File Format

```
Timestamp,Heart Rate,RR Interval in seconds
2023/11/3 22:39:15,85,0.708
2023/11/3 22:39:15,85,0.657
...
```

---

## Usage Instructions

### 1. Set the Dataset Path

Open `Code/mainEntry.jl` and update the base directory path:

```julia
base_dir = "/path/to/your/modified_signals"
```

### 2. Run the Full Pipeline

From the terminal:

```bash
julia Code/mainEntry.jl
```

Or from a Jupyter notebook, run each cell in order starting from the includes:

```julia
include("Code/data_preprocessing.jl")
include("Code/metrics_computation.jl")
include("Code/pearson_relation_coeff.jl")
include("Code/p_value_calculations.jl")
include("Code/regression_plot_and_bland_altman_plot.jl")
```

Then call the main pipeline:

```julia
base_dir = raw"C:\path\to\modified_signals"
all_hr_rates, all_br_rates, all_hr_ref = process_patients_data(base_dir)
```

Results, metrics, and plots are printed and saved automatically by `mainEntry.jl`.

---

## Pipeline Description

### Functions in `data_preprocessing.jl`

| Function | Description |
|----------|-------------|
| `get_bior39_filters()` | Returns bior3.9 reconstruction filters (Lo, Hi) |
| `modwt(x, wname, J)` | Maximal Overlap Discrete Wavelet Transform |
| `modwtmra(w, wname)` | Multi-Resolution Analysis from MODWT coefficients |
| `detect_patterns(pt1, pt2, win_size, data, time)` | Segments signal into sleep/movement/empty windows |
| `scipy_filtfilt(filter_coefs, x)` | Zero-phase filter matching SciPy's filtfilt exactly |
| `band_pass_filtering(data, fs, filter_type)` | Chebyshev Type I bandpass for BCG or breathing |
| `detect_peaks(x; mpd, mph, ...)` | Peak detection with minimum distance filtering |
| `compute_rate(beats, fs, mpd)` | BPM from peak spacing in samples |
| `vitals(t1, t2, win_size, limit, sig, fs, mpd)` | Per-window BPM across full signal |
| `process_patient(file)` | Full pipeline for one BCG file |
| `process_patient_reference_data(rr_file, window_sec)` | Windowed HR from RR reference file |
| `find_rr_file(patient_dir, bcg_file)` | Locates matching RR file by date |
| `process_patients_data(base_dir)` | Iterates all patients, returns HR/BR dicts |
| `display_summary_for_data_analysis(...)` | Prints BCG vs reference comparison table |

### Functions in other modules

| File | Function | Description |
|------|----------|-------------|
| `metrics_computation.jl` | `compute_metrics(...)` | MAE, RMSE, MAPE across all patients |
| `pearson_relation_coeff.jl` | `compute_pearson(...)` | Pearson r (manual implementation) |
| `p_value_calculations.jl` | `compute_pvalue(r, n; tail)` | p-value via incomplete Beta function |
| `p_value_calculations.jl` | `regularised_beta(x, a, b)` | Regularised incomplete Beta |
| `p_value_calculations.jl` | `log_gamma(x)` | ln Γ(x) via Lanczos approximation |
| `regression_plot_and_bland_altman_plot.jl` | `plot_regression_and_bland_altman(...)` | Regression + Bland–Altman plots |

---

## Evaluation Metrics

All metrics are computed **manually** without external statistical libraries.

| Metric | Formula | Interpretation |
|--------|---------|----------------|
| **MAE** | `(1/n) Σ\|pred − ref\|` | Average absolute error in bpm |
| **RMSE** | `√((1/n) Σ(pred − ref)²)` | Penalises large errors more than MAE |
| **MAPE** | `(100/n) Σ\|pred − ref\| / \|ref\|` | Scale-independent percentage error |
| **Pearson r** | `Σ(p−p̄)(r−r̄) / √(Σ(p−p̄)²·Σ(r−r̄)²)` | Linear correlation in [−1, +1] |
| **p-value** | Right-tailed t-test, H₁: r > 0 | Statistical significance of correlation |
| **Regression** | Linear fit: BCG = slope × Ref + intercept | Systematic relationship |
| **Bland–Altman** | Bias ± 1.96 SD limits of agreement | Method agreement analysis |

---

## Results

Evaluated across **22 patients** (overnight BCG recordings, fs = 140 Hz):

| Metric | Value |
|--------|-------|
| MAE | 24.642 bpm |
| RMSE | 26.348 bpm |
| MAPE | 42.339 % |
| Pearson r | −0.5694 |
| t-statistic | −3.0977 (df = 20) |
| p-value | 0.002837 ✓ (significant) |
| Bias (Bland–Altman) | +24.64 bpm |
| Limits of Agreement | [5.93, 43.35] bpm |

![Regression and Bland-Altman Plots](Plots/regression_bland_altman.png)

---

## Known Limitations

- BCG estimates converge to ~91 bpm across patients regardless of true HR, indicating the wavelet level and peak detection minimum distance require further tuning for fs = 140 Hz.
- Breathing rate estimates (avg ~19 bpm) are physiologically plausible and appear more reliable than HR estimates.
- MODWTMRA reconstruction is imperfect (inherited from original Python implementation).
- RR reference file timestamps may contain gaps, affecting per-window reference HR accuracy.

---

## References

1. Chen et al., "A ballistocardiogram dataset with reference sensor signals in long-term natural sleep environments," *Scientific Data*, 2024. https://doi.org/10.1038/s41597-024-03950-5

2. Code Ocean Capsule #1398208 — BCG Signal Processing Reference Implementation. https://codeocean.com/capsule/1398208

3. M. Duarte, "detect_peaks.py," BMC Notebook, 2016. https://github.com/demotu/BMC
