# Repository guidance

## Purpose

This repository contains interactive teaching modules for a Solid State Ionics
course. Favor physical transparency, projector readability, and short
self-contained notebooks over framework-building.

## Scientific invariants

- State canonical units at every model boundary.
- In Module 01, use K, eV per defect, k_B per defect, and dimensionless defect
  fraction. Preserve the exact finite-N binomial multiplicity, distinguish it
  from the Stirling limit, and introduce the dilute exponential only as a
  low-defect approximation.
- In Module 02, use K, bar, and cm^-3. Preserve the documented SrTiO3
  mass-action laws and exact electroneutrality equation.
- Calculate Module 02 concentration curves from the full equilibrium model.
  Limiting Brouwer balances and slopes may be annotations or tests, never
  solver inputs.
- Keep all solved concentrations strictly positive. Prefer logarithmic variables
  for equations spanning many decades.
- In Module 03, use the documented per-neighbor convention
  w0 = nu0 exp(-H_mig/kBT) and D = a^2 w0. Reconcile it explicitly with the
  lecture convention Gamma = 2 w0 and D = a^2 Gamma / 2.
- Keep particle-scale (k_B, e) and molar (R, F) transport equations visibly
  separated. For the ideal monovalent ambipolar model, preserve local
  electroneutrality, zero current, and D_chem = 2 Di De/(Di + De) Theta, while
  stating that the factor of two is model-specific.
- If a species set or thermodynamic constant changes, explain the scientific
  reason and update the displayed assumptions and validation checks together.

## Notebook structure

- Keep physics constants and solver functions in dedicated cells, independent of
  marimo controls and plotting.
- Keep plotting and prose downstream of the solved data.
- Use only marimo, NumPy, SciPy, and matplotlib unless a new dependency is
  clearly justified and compatible with browser/WASM execution.
- Do not create a separate physics package until multiple modules need shared
  kernels.

## Verification

Before handing off a change:

1. Compile all notebooks as Python.
2. Run strict marimo checks and execute or export all notebooks in a clean
   environment.
3. For Module 01, check the free-energy minimum, chemical-potential zero,
   large-N Stirling convergence, finite-N composition spacing, and dilute limit.
4. For Module 02, check mass-action residuals, scaled electroneutrality
   residuals, positivity, and limiting slopes where the selected window samples
   them.
5. For Module 03, check hopping-rate conventions, MSD, detailed balance,
   low-field drift, number conservation, Fick flux, electrochemical-potential
   cancellation, ambipolar zero current, thermodynamic-factor multiplication,
   positivity, and finiteness.
6. Run git diff --check and review all rendered figures at projector scale.

Do not weaken numerical tolerances simply to silence a failed physics check;
identify whether the model, regime test, or implementation is responsible.
