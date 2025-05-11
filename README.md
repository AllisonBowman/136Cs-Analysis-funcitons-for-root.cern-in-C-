# 136Cs ROOT Analysis Toolkit

This repository contains a set of ROOT-based C++ functions developed for the analysis of neutron time-of-flight (TOF) spectra from experiments conducted at the Triangle Universities Nuclear Laboratory (TUNL) investigating decay sequences of **Cesium-136 (136Cs)**.

## Project Overview

The primary goal of this project was to analyze β-decay transitions and neutron emission following Gamow-Teller transitions in 136Cs. The experimental data was collected using a gaseous xenon target bombarded with a 6.8 MeV proton beam, resulting in the formation of excited 136Cs nuclei. The TOF spectra from neutron and γ-ray detectors were analyzed to identify specific excited nuclear states and to search for potential missing transitions.

## Key Features

- **TOF Histogram Construction**: Code organizes raw ROOT files into `TH1F` objects by detector cluster (top, mid, bottom) and angle (0° to 12°).
- **Detector Synchronization**: Time alignment of signals using γ-ray peaks as a timing reference.
- **Polynomial Background Subtraction**: Allows subtraction of background using customizable polynomial models.
- **Peak Fitting**:
  - Single Gaussian and triple Gaussian peak fitting functions.
  - Fit centers anchored to theoretical TOF positions for known 136Cs states from literature.
  - Convolution with gamma-ray response shape using a parametrized triple-Gaussian model.
- **Compound Nuclear Modeling**: Includes simulated neutron smearing from TRIM modeling and non-selective decay models to explain angular independence.
- **Angular Comparison**: Enables comparison of transition intensities across angles to test for preferential state selection.

## Repository Contents

- `tunl_pn_2dplot.C`: Primary script for histogram generation and detector synchronization.
- `SMEAR_tunl_pn_2dplot_APBgammaDISfit_compoundLIT.C`: Compound nuclear model fitting with smearing and background modeling.
- Fitting and analysis functions tied to theoretical state energies extracted from:
  - TUNL 136Cs state evaluations
  - Rebeiro Ph.D. thesis (University of the Western Cape)
  - NNDC and AME2016 atomic mass evaluation data

## How to Use

1. **Install ROOT** (if not already):
   ```bash
   sudo apt install root-system
   ```
2. **Compile and run the main script**:
   ```bash
   root -l
   .x tunl_pn_2dplot.C
   ```
3. **Generate `rainbow#.root` histograms** for each angle.
4. **Fit histograms** using peak models in the `SMEAR_tunl_*.C` script.

## Notable Parameters and Models

- **Gaussian Fit Model**:
  ```cpp
  f(x) = A * exp(-0.5 * ((x - (μ + δμ)) / σ)^2)
  ```
  - Amplitude (A), width (σ), and offset (δμ) are optimized by the fitting routine.
- **Triple Gaussian Response Model**: Matches the γ-ray signal structure to model neutron peaks more accurately.
- **Smearing Function from TRIM**:
  ```math
  σ = a*μ³ - b*μ² + c*μ - d
  ```

## Findings

- **No angular dependence** was observed in the neutron spectrum, supporting the hypothesis of compound nuclear decay over selective transition rules.
- **Several previously unreported or unresolved states** were indicated in the neutron TOF fits, suggesting further exploration in the >2.6 MeV energy region may be fruitful.

## Citation

Please cite this code and project using:

> Bowman, A. Daniels, T PhD. (2023). *Analyzing Nuclear Decay Sequence of Cesium*. Triangle Universities Nuclear Laboratory.

## Acknowledgements

This project was developed under the mentorship of Dr. Daniels at TUNL. Support and data courtesy of the TUNL collaboration.
