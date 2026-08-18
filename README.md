# Solid State Ionics Interactive

[![Deploy marimo apps to GitHub Pages](https://github.com/qiyanglu/Solid-State-Ionics-Interactive/actions/workflows/pages.yml/badge.svg)](https://github.com/qiyanglu/Solid-State-Ionics-Interactive/actions/workflows/pages.yml)

**[Browse the interactive course modules &rarr;](https://qiyanglu.github.io/Solid-State-Ionics-Interactive/)**

This repository is a growing series of interactive [marimo](https://marimo.io/)
notebooks for teaching the **Solid State Ionics** course at Westlake University.
The notebooks complement the lectures and course materials collected on the
[Solid State Ionics Lab teaching page](https://ssi-westlake.com/teaching/).

The series helps students explore how defect formation, defect chemistry, ionic
and electronic transport, electrochemical polarization, space charge, and
impedance emerge from governing equations. Each module favors transparent
physics, interactive controls, and numerical checks over preassembled textbook
curves.

## Notebooks

| Module | Topic | Interactive app |
|---|---|---|
| 01 | Defect Formation Thermodynamics | [Launch](https://qiyanglu.github.io/Solid-State-Ionics-Interactive/01-defect-formation/) |
| 02 | Brouwer Diagram for Acceptor-Doped SrTiO3 | [Launch](https://qiyanglu.github.io/Solid-State-Ionics-Interactive/02-brouwer-sto/) |

### Module 01: Defect Formation Thermodynamics

[01_defect_formation.py](01_defect_formation.py) follows one neutral defect
species on N equivalent lattice sites. It connects the exact finite-system
multiplicity

~~~text
Omega = binomial(N, n)
S_config = k_B ln(Omega)
~~~

to the Stirling-limit entropy, free-energy minimum, chemical-potential zero, and
equilibrium logistic occupancy. Interactive controls expose temperature,
formation enthalpy, non-configurational formation entropy, and finite lattice
size. The notebook compares exact finite-N combinatorics, the Stirling limit,
and the dilute Boltzmann approximation without introducing charged defects,
oxygen pressure, or electroneutrality.

### Module 02: Brouwer Diagram for Acceptor-Doped SrTiO3

[02_brouwer_sto.py](02_brouwer_sto.py) solves, independently at every oxygen
partial pressure,

~~~text
K_red = V n^2 pO2^(1/2)
K_eh  = n p
2 V + p = A + n
~~~

for oxygen vacancies V, electrons n, and holes p, with a fixed fully ionized
acceptor concentration A. The unique positive solution is found with a
bracketed log-space root solve. Brouwer regimes and slopes are measured only
after solving the full equations; they are never used to construct the curves.

The default state is 973 K, A = 10^18 cm^-3, and pO2 from 10^-25 to 1 bar.
Controls span 700-1500 K and log10(A/cm^-3) from 13 to 21.

## Run locally

With Python 3.11 or newer and [uv](https://docs.astral.sh/uv/) installed:

~~~console
uv sync
uv run marimo edit 01_defect_formation.py
uv run marimo edit 02_brouwer_sto.py
~~~

To present either notebook as an app, replace edit with run.

## Validation

~~~console
uv run marimo check --strict 01_defect_formation.py 02_brouwer_sto.py
uv run marimo export html 01_defect_formation.py -o defect-formation.html --no-include-code -f
uv run marimo export html 02_brouwer_sto.py -o brouwer.html --no-include-code -f
~~~

Module 01 displays checks for its free-energy minimum, chemical-potential zero,
finite-N rounding, large-N Stirling convergence, and dilute limit. Module 02
displays mass-action, electroneutrality, positivity, regime-coverage, and
limiting-slope checks. Its mass-action rows report logarithmic residuals with an
exact target of zero (equivalently, an unlogged equilibrium ratio of one); values
below `1e-12` are displayed consistently as numerical zero rather than giving
exact cancellation and floating-point roundoff different visual weight.

## Browser deployment

Both notebooks use only marimo, NumPy, SciPy, and matplotlib and perform no
network or filesystem access at runtime. The workflow in
[pages.yml](.github/workflows/pages.yml) exports each notebook as a
browser-hosted WASM app and deploys the generated static site to GitHub Pages on
every push to main. The project root is a stable module index; each notebook is
published under its own numbered path. After deployment, the workflow retries
both live module routes and verifies their browser entry bundles plus the
dynamically imported marimo run-page modules.

The former /01-brouwer-sto/ URL redirects to /02-brouwer-sto/.
