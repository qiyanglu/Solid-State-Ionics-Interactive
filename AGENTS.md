# Repository guidance

## Purpose

This repository contains interactive teaching modules for a Solid State Ionics
course. Favor physical transparency, projector readability, and short
self-contained notebooks over framework-building.

## Scientific invariants

- State canonical units at every model boundary: K, bar, and cm^-3 in module 01.
- Preserve the exact mass-action laws and electroneutrality equation documented
  in the notebook.
- Calculate concentration curves from the full equilibrium model. Limiting
  Brouwer balances and slopes may be annotations or tests, never solver inputs.
- Keep all solved concentrations strictly positive. Prefer logarithmic variables
  for equations spanning many decades.
- If the species set or thermodynamic constants change, explain the scientific
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

1. Compile the notebook as Python.
2. Run the notebook or marimo's notebook tests in a clean environment.
3. Check mass-action residuals, scaled electroneutrality residuals, positivity,
   and limiting slopes where the requested parameter window samples them.
4. Run `git diff --check` and review the rendered figure at projector scale.

Do not weaken numerical tolerances simply to silence a failed physics check;
identify whether the model, regime test, or implementation is responsible.
