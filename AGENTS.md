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
- In Module 03, follow the lecture convention
  Gamma = nu exp(-Delta H_mig/kBT), where Gamma is the total one-dimensional
  hop frequency, and preserve D = a^2 Gamma / 2.
- Keep every spatial model in Module 03 one-dimensional. Keep particle-scale
  (k_B, e) and molar (R, F) transport equations visibly separated.
- For the dilute Li <-> Li+ + e- example, preserve local electroneutrality,
  local equilibrium, equal steady fluxes, zero current, and
  D_Li^delta = 2 D_Li+ D_e-/(D_Li+ + D_e-). Do not make a non-ideal
  thermodynamic factor a central control unless the teaching scope changes.
- In Module 04, keep the geometry planar and one-dimensional. Use K, V, cm^-3,
  nm, and microF/cm^2 at the teaching boundary, with SI conversions explicit in
  the physics functions.
- Preserve the exact planar Gouy-Chapman solution, the lecture's frozen-dopant
  Mott-Schottky depletion approximation, phi_infinity = 0, and the distinction
  between chemical, electrical, and electrochemical potentials.
- In the Gouy-Chapman-Stern and Frumkin sections, enforce one common interfacial
  charge, phi_0 = (phi_0 - phi_1) + phi_1, series differential capacitance, and
  reaction-plane concentrations evaluated at x = x_1. Use a signed reactant
  charge number in the Boltzmann factor.
- In Module 05, keep a one-dimensional ideal H <-> H+ + e- pair with local bulk
  electroneutrality c_i = c_e = c, two ion-blocking/electron-reversible
  electrodes, and the lecture voltage convention
  U = (mu_tilde_e(L) - mu_tilde_e(0))/F.
- Preserve mu = mu_i + mu_e = 2 RT ln(c/c_bar),
  D_delta = 2 D_i D_e/(D_i + D_e), total ion conservation, and zero ionic flux
  at both boundaries. Support constant current and constant potential without
  assuming the transient concentration or potential profiles.

- In Module 06, retain the Module 05 ideal H <-> H+ + e- pair,
  local electroneutrality, mu_H = mu_i + mu_e = 2 RT ln(c/c0), and
  D_delta = 2 D_i D_e/(D_i + D_e), but use complementary selective contacts:
  J_i(0,t) = 0 at the electronic current collector and J_e(L,t) = 0 at the
  ion-conducting electrolyte.
- Define positive current and voltage as extraction and state that convention.
  Generate PITT, GITT, and OCV concentration and potential profiles from the
  full finite one-dimensional diffusion model. At OCV set terminal current to
  zero, preserve the end-of-pulse profile, and conserve the mean composition.
- Present the Cottrell, square-root, first-mode, and late linear formulas only
  as classical one-sided, small-signal limits. Never use them as solver inputs,
  and keep theta = D_delta t/L^2 distinct from t/tau_delta.
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
5. For Module 03, check the Gamma convention, MSD, detailed balance,
   low-field drift, one-dimensional number conservation, Fick flux,
   electrochemical-potential cancellation, equal Li-ion/electron flux,
   zero current, conductivity-form equivalence, positivity, and finiteness.
6. For Module 04, check the Boltzmann and electrochemical-potential identities,
   exact Gouy-Chapman profile and Gauss law, Mott-Schottky boundaries and
   Poisson curvature, GCS charge and voltage matching, series capacitance,
   Frumkin reaction-plane consistency, positivity, and finiteness.
7. For Module 05, check the uniform initial state, total-ion conservation,
   positivity, the selected electrical drive, the initial total-conductivity
   response, zero ionic boundary flux, chemical-diffusivity identity,
   electrochemical-potential decomposition, measured voltage, and late-time
   steady state.
8. For Module 06, check the uniform initial state, positivity, PITT voltage and
   GITT current control, zero terminal current and constant mean composition
   during OCV, current-composition balance, the chemical-diffusivity identity,
   voltage reconstruction from chemical and electrical potentials, classical
   PITT/GITT short- and long-time limits, and first-mode OCV decay.
9. Run git diff --check and review all rendered figures at projector scale.

Do not weaken numerical tolerances simply to silence a failed physics check;
identify whether the model, regime test, or implementation is responsible.
