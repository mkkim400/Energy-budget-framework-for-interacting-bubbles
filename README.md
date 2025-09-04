# Bubble dynamics simulations with weak bubble-bubble interactions

This repository contains code for simulating and analyzing the dynamics of **interacting gas bubbles** in a liquid using the **Keller–Miksis equations**. The project includes:

- ✅ MATLAB driver for parameter sweeps
- ✅ MATLAB solver for coupled bubble dynamics (`KM_Ida_solver`)
- ✅ R script for post-processing, regression analysis, and 3D visualization


---

## Project Objective

This project models the **mutual interaction of two gas bubbles** in a compressible fluid. It focuses on the effect of:
- **Initial separation (D)**
- **Radius ratio (RR)**
- **Ambient-to-internal pressure ratio (PR)**

Outputs include:
- Bubble radius evolution
- Wall velocities and pressures
- Dimensional and non-dimensional metrics
- Power-law fits to max dynamics variables

---

## Methodology

### 1. MATLAB Driver (`main_driver.m`)
Sweeps through combinations of:
- Normalized separation distances `D`
- Radius ratios `RR`
- Pressure ratios `PR`

Calls the solver for each configuration, extracts results, and writes:
- `R_Q.txt` – Dimensional results
- `R_Q_non.txt` – Non-dimensional results

### 2. MATLAB Solver (`KM_Ida_solver.m`)
Solves the **coupled Keller–Miksis equations** for two interacting bubbles with:
- Compressibility
- Bubble–bubble coupling
- Adiabatic gas dynamics

Uses MATLAB’s `ode15s` stiff solver and nondimensionalizes the equations for numerical stability.

### 3. R Script (`regression_model.R`)
Performs:
- **Nonlinear power-law regression**:  
  \[
  \text{value} = \beta_0 + \beta_1 \cdot RR^{p_1} \cdot D^{p_2}
  \]
- **Grid search** to find optimal exponents `p1`, `p2`
- **Fit evaluation** via RMSE, AIC, BIC
- **3D plots** of predicted vs. actual values
- **Exports results** with confidence intervals

---

## Dependencies

### R
Install dependencies via:

```r
install.packages(c("ggplot2", "data.table", "dplyr", "readxl", "mgcv", "MASS", "plotly", "reshape2"))

