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
- Preserve mu = mu_i + mu_e = 2 RT ln(c/c0),
  D_delta = 2 D_i D_e/(D_i + D_e), total ion conservation, and zero ionic flux
  at both boundaries. Support constant current and constant potential without
  assuming the transient concentration or potential profiles.
- In Module 06, retain the Module 05 ideal neutral M <-> M+ + e- pair, with a
  Li/H label choice, local electroneutrality,
  mu_M = mu_i + mu_e = 2 RT ln(c/c0), and
  D_delta = 2 D_i D_e/(D_i + D_e). In the student-facing article coordinate,
  J_e(0,t) = 0 at the ion-conducting electrolyte and J_i(L,t) = 0 at the
  electronic current collector. A mirrored internal coordinate is acceptable
  only when stated and all face-flux/profile mappings are checked.
- Define positive current and voltage as extraction and state that convention.
  Generate PITT, GITT, and OCV concentration and potential profiles from the
  full finite one-dimensional diffusion model. At OCV set terminal current to
  zero, preserve the end-of-pulse profile, and conserve the mean composition.
  Enforce the two selective-contact face fluxes directly with a conservative
  spatial representation; do not reconstruct a nonzero boundary gradient from
  a basis whose derivative is forced to zero at the endpoints.
- Present the Cottrell, square-root, first-mode, and late linear formulas only
  as classical one-sided, small-signal limits. Never use them as solver inputs,
  and keep theta = D_delta t/L^2 distinct from t/tau_delta.
- Keep the finite-kinetics PITT extension separate from the core ideal boundary.
  Use Bi = k_delta L/D_delta, tau_d = L^2/D_delta, tau_ct = L/k_delta,
  lambda tan(lambda) = Bi, and show reaction-, mixed-, and diffusion-controlled
  limits without patching curves. Check fitting bias rather than claiming that
  a linear long-time plot proves diffusion control.
- In Module 07, use the e^(i omega t) convention, define Z = Z' + i Z'', and
  plot Z' against -Z'' on Nyquist axes. Preserve the distinction between
  frequency in Hz and angular frequency in rad/s.
- Derive Warburg curves from the one-dimensional chemical-diffusion equation.
  Keep the semi-infinite, fixed-composition finite-length, and zero-flux
  finite-length boundary conditions explicit; never insert a 45-degree line as
  a solver input.
  Pair every boundary name with its equation and state that an open diffusion
  boundary is a fixed-composition reservoir, not electrical open circuit.
- For the Module 07 transmission line, preserve the continuous two-rail MIEC
  equations, voltage-equivalent electrochemical potentials, total-current
  conservation, SI circuit units, and explicit terminal boundary conditions.
  Keep the chemical time scale distinct from a universally assigned peak
  frequency.
  Use u_e = -mu_tilde_e/F and u_i = +mu_tilde_i/F for monovalent carriers, so
  u_e - u_i = -mu_neutral/F. Distinguish distributed r_e, r_i, c_chem from
  total R_e, R_i, C_chem and reserve S for area.
- If a species set or thermodynamic constant changes, explain the scientific
  reason and update the displayed assumptions and validation checks together.

## Notebook structure

- Keep physics constants and solver functions in dedicated cells, independent of
  marimo controls and plotting.
- Keep plotting and prose downstream of the solved data.
- Use only marimo, NumPy, SciPy, and matplotlib unless a new dependency is
  clearly justified and compatible with browser/WASM execution.
- Begin every module with a guiding question, two or three learning goals, a
  notation/model-scope box, and a prediction before the main controls.
- Organize the core classroom path as prediction -> controls -> figure ->
  takeaway. Put extended derivations, implementation detail, advanced
  interpretation, and detailed checks in collapsed sections where practical.
- End every module with exactly three primary messages and cross-link the next
  relevant module or the shared notation guide.

- Do not create a separate physics package until multiple modules need shared
  kernels.
- Update `PROJECT_CONTEXT.md` only when scientific scope, repository
  architecture, shared notation, or the module roadmap changes. Routine visual
  polish and numerical bug fixes do not require a context rewrite.

## Visualization

- Use projector-readable typography, low-saturation colors, and approximately
  1.5--2.0 pt data lines. Pair color with line style, markers, or direct labels
  where curves must remain distinguishable without color.
- Capitalize plot labels and titles consistently. Typeset scientific quantities
  with math text, and state units or an explicit dimensionless normalization on
  every axis.
- Use equal horizontal and vertical scale for every Nyquist diagram.

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
5. For Module 03, check the Gamma convention, the fitted MSD slope against
   <x^2> = 2 D t and the analytical diffusivity, detailed balance, low-field
   drift, one-dimensional number conservation, Fick flux,
   electrochemical-potential cancellation, equal Li-ion/electron flux,
   zero current, conductivity-form equivalence, positivity, and finiteness.
6. For Module 04, check the Boltzmann and electrochemical-potential identities,
   exact Gouy-Chapman profile and Gauss law, the Mott-Schottky boundaries and
   Poisson curvature, GCS charge and voltage matching, series capacitance,
   Frumkin reaction-plane consistency, positivity, and finiteness.
7. For Module 05, check the uniform initial state, total-ion conservation,
   positivity, the selected electrical drive, the initial total-conductivity
   response, zero ionic boundary flux, chemical-diffusivity identity,
   electrochemical-potential decomposition, measured voltage, and late-time
   steady state.
8. For Module 06, check the uniform initial state, positivity, PITT voltage and
   GITT current control, zero terminal current and constant mean composition
   during OCV, selective-contact face fluxes, current-composition balance,
   spatial-grid convergence, the chemical-diffusivity identity,
   voltage reconstruction from chemical and electrical potentials, classical
   PITT/GITT short- and long-time limits, long-time OCV decay, finite-kinetics
   eigenvalue limits, positivity and conservation, and the direction and
   fitting-window dependence of diffusion-only bias versus Biot number.
9. For Module 07, check the phasor convention, resistor, capacitor, and ideal
   parallel-RC limits, the RC semicircle and apex, equal Nyquist scaling,
   DC/AC diffusion-length scaling, general-to-dilute Warburg resistance mapping,
   finite-length Warburg limits and passivity, TLM boundary residuals,
   voltage-equivalent-potential signs, distributed/total conversions,
   total-current conservation, reversible-contact limit, passivity, and
   finiteness.
10. Run git diff --check and review all rendered figures at projector scale.

Do not weaken numerical tolerances simply to silence a failed physics check;
identify whether the model, regime test, or implementation is responsible.
