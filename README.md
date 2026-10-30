# Endogenous Hawkes-Merton Degradation Framework for Rolling Element Bearings

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Paper: RESS](https://img.shields.io/badge/Journal-RESS%20Under%20Review-orange.svg)](#citation)

Reference implementation of the stochastic degradation and maintenance-risk model presented in:
> **An Endogenous Hawkes-Merton Degradation Framework for Rolling Element Bearings: Modeling, Calibration, and Maintenance-Risk Assessment**  
> *Ognjen Raketić, Phuc Do, Enrico Zio* (Submitted to *Reliability Engineering & System Safety*)

---

## Overview

Standard Prognostics and Health Management (PHM) models frequently rely on memoryless Markovian or homogeneous renewal assumptions that fail to reflect the path-dependent shock cascades characteristic of Rolling Contact Fatigue (RCF). 

This repository implements a coupled **Hawkes-Merton jump-diffusion Stochastic Differential Equation (SDE)** where:
1. Structural micro-impacts are modeled as a non-homogeneous point process modulated by an endogenous, state-dependent intensity $\lambda(t)$.
2. The excitation kick $\alpha(H(t))$ scales quadratically with cumulative fatigue damage $H(t)$, derived from Linear Elastic Fracture Mechanics (LEFM).
3. The Markovian exponential kernel enables exact $\mathcal{O}(1)$ discrete recursive hazard updates ($1.50\,\mu\text{s}$ execution footprint).
4. Run-to-failure telemetry from the **NASA Intelligent Maintenance Systems (IMS)** bearing benchmark is used to evaluate prognostic lead times ($\Delta\tau_{\text{LIH}}$) across a 40,000-path Monte Carlo battery.

---

## Mathematical Formulation

### 1. Coupled Jump-Diffusion State Equation
The latent structural damage index $H(t) \in [0.10, H_{\text{crit}}]$ evolves via:

$$dH(t) = a(H(t))dt + \sigma(H(t))dW(t) + d\left(\sum_{i=1}^{N(t)} Y_i\right)$$

- **Continuous Wear & Diffusion:** Deterministic drift $a(H(t))$ and Brownian driver $\sigma(H(t))$ represent background tribological wear and thermal vibration.
- **Counting Process $N(t)$:** Point process modulated by instantaneous conditional hazard rate $\lambda(t)$.
- **Jump Severity $Y_i$:** Lognormal damage severity per fracture event: $\ln(Y_i) \sim \mathcal{N}(\mu_J, \sigma_J^2), \quad Y_i > 0$.

### 2. State-Dependent Hawkes Intensity
The hazard rate $\lambda(t)$ incorporates past damage history via an exponential response kernel:

$$\lambda(t) = \mu_0 + \sum_{t_i < t} \alpha(H(t_i)) e^{-\beta(t - t_i)}$$

$$\alpha(H(t)) = \alpha_0 \left[1 + \gamma \left(\frac{\min(H(t), H_{\text{crit}})}{H_{\text{crit}}}\right)^2\right]$$

- $\mu_0$: Baseline hazard in nominal elastohydrodynamic regime.
- $\beta$: Elastodynamic stress-relaxation and shear dissipation rate ($h^{-1}$).
- $\alpha(H)$: Quadratic state-feedback amplification scaling with normalized crack area ($K \propto \sqrt{A} \implies G \propto K^2 \propto H$).
- Subcritical stability is governed by the branching ratio $n(H) = \alpha(H)/\beta < 1$.

### 3. Exact $\mathcal{O}(1)$ Recursive Hazard Updates
Over inspection epochs $\Delta t = t_k - t_{k-1}$ ($10\text{ min} = 0.1667\text{ h}$):

$$\lambda(t_k) = \mu_0 + (\lambda(t_{k-1}) - \mu_0)e^{-\beta \Delta t} + \sum_{t_{k-1} \le t_j < t_k} \alpha(H(t_j))$$

---

| Metric / Statistical Test | Empirical Value | Physical / Operational Significance |
| :--- | :---: | :--- |
| **Bearing 1 Mean Lead Time (E[Δτ_LIH])** | **61.99 h** | Average tactical decision margin before catastrophic seizure |
| **Bearing 1 Median Lead Time (Med(Δτ_LIH))** | **59.95 h** | Robust central tendency (95% CI: [57.10 h, 62.40 h]) |
| **24-Hour Buffer Compliance (P(Δτ_LIH ≥ 24 h))** | **96.92%** | Probability of meeting mobilization window τ_req = 24 h |
| **Bearing 3 False Positive Rate (Nominal)** | **0.00%** | Zero spurious alarms during undamaged operational life |
| **Bearings 2 & 4 False Positive Rate (Nominal)** | **0.00%** | Rejection of shaft crosstalk during nominal phase |
| **Mann-Whitney U Test (H₁: Δτ > 24 h)** | **U = 9.69 × 10⁷ (p < 10⁻¹⁵)** | Significant superiority over minimum operational notice |
| **Kolmogorov-Smirnov Test (Real Telemetry)** | **D = 0.4003 (p = 1.22 × 10⁻¹⁴²)** | Decisive rejection of Gaussian degradation increments |
| **99% Value-at-Risk (Reactive Replacement)** | **€ 2,660,627.42** | Uncontained catastrophic downtime tail loss |
| **99% Value-at-Risk (Prescriptive LIH Policy)** | **€ 288,637.01** | Truncated planned downtime cost |
| **Extremal Tail-Risk Mitigation (ΔVaR 99%)** | **89.15%** | Net financial exposure reduction |

---

## Repository Structure

```text
├── Hawkes_Merton_Bearing_Prognostics.ipynb  # Vectorized 40k Monte Carlo simulation & plotting
├── requirements.txt                         # Dependency specifications
├── figures/                                 # Generated publication-quality vector outputs (PDF)
│   ├── fig1_degradation_ensemble.pdf        # Ensemble damage trajectories H(t) across bearings
│   ├── fig2_intensity_evolution.pdf         # Temporal KDE snapshots of hazard density lambda(t)
│   ├── fig3_ecdf_lih.pdf                    # Empirical CDF of actionable lead time margin
│   ├── fig4_qq_plot.pdf                     # Exponential Q-Q plot proving jump heavy tails
│   ├── fig5_stability_surface.pdf           # Phase diagram of branching ratio n_0 = alpha / beta
│   ├── fig6_financial_tails.pdf             # VaR99% loss distributions: reactive vs. prescriptive
│   ├── fig7_roc_curves.pdf                  # Multi-cohort ROC/AUC curves for lead horizons
│   └── fig8_volatility_bifurcation.pdf      # Intensity volatility bifurcation (faulty vs. nominal)
└── README.md
