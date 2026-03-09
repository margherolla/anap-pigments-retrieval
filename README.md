# aNAP and phytoplankton pigment retrieval from hyperspectral absorption

MATLAB workflow for estimating **non-algal particle absorption spectra (aNAP)** from cp(660), deriving **phytoplankton absorption spectra (aPH)**, performing **Gaussian decomposition**, and retrieving **phytoplankton pigment concentrations**.

The repository provides two main workflows:

1. **Standalone aNAP estimation**
2. **Full phytoplankton pigment retrieval workflow**

---

# Workflow overview

Full pipeline:

cp(660)  
↓  
aNAP(λ)  
↓  
aPH(λ) = aP(λ) − aNAP(λ)  
↓  
Gaussian decomposition of aPH(λ)  
↓  
Pigment concentration estimation  

If particulate absorption spectra are not available, the workflow can also be used only to estimate:

cp(660) → aNAP(λ)

---

# Repository structure

```
anap-pigment-retrieval/

main/
    Example scripts to run the workflow

functions/
    Core MATLAB functions

coefficients/
    Regression coefficients used by the algorithms

example_data/
    Example datasets

docs/
    Additional documentation
```

---

# Installation

Clone the repository:

```
git clone https://github.com/USERNAME/anap-pigment-retrieval.git
```

Add functions to the MATLAB path:

```
addpath(genpath('functions'))
```

---

# Quick start

## Example 1 — derive only aNAP

Run:

```
run main/run_example_anap_only.m
```

Required inputs:

- wavelength vector (`lambda`)
- cp(660)
- optional cp(660) uncertainty

Outputs:

- aNAP(λ)
- aNAP(400)
- S_NAP
- uncertainty statistics

---

## Example 2 — full workflow

Run:

```
run main/run_example_full_workflow.m
```

Required inputs:

- wavelength vector (`lambda`)
- cp(660)
- cp(660) uncertainty (optional)
- particulate absorption spectra `aP(λ)`
- uncertainty on `aP(λ)`

Outputs:

- aNAP(λ)
- aPH(λ)
- Gaussian peak amplitudes
- pigment concentrations

---

# Main functions

## function_anap_estimation_from_cp660

Estimate aNAP spectra from cp(660).

Supports:

- deterministic estimation
- coefficient uncertainty
- coefficient + cp(660) uncertainty

---

## spectral_decomp_aph

Performs Gaussian decomposition of phytoplankton absorption spectra.

Outputs:

- Gaussian peak amplitudes
- component spectra
- reconstructed spectra

---

## derive_pigm_gaus

Derives pigment concentrations from Gaussian peak amplitudes using regression coefficients.

---

# Example datasets

The repository includes small example datasets:

| File | Description |
|-----|-------------|
| example_dataset_anap_only.mat | minimal dataset for aNAP estimation |
| example_dataset_full.mat | dataset for the full workflow |

---

# Citation

If you use this code please cite:

Costanzo et al. (2026)  
*A novel approach for estimating non-algal particle absorption from AC-S and retrieving phytoplankton pigments.*

---

# Author

Margherita Costanzo