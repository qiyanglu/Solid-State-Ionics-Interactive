# Solid State Ionics Interactive

[![Deploy marimo app to GitHub Pages](https://github.com/qiyanglu/Solid-State-Ionics-Interactive/actions/workflows/pages.yml/badge.svg)](https://github.com/qiyanglu/Solid-State-Ionics-Interactive/actions/workflows/pages.yml)

**[Open the interactive notebook →](https://qiyanglu.github.io/Solid-State-Ionics-Interactive/)**

This repository is a growing series of interactive [marimo](https://marimo.io/)
notebooks for teaching the **Solid State Ionics** course at Westlake University.
The notebooks complement the lectures and course materials collected on the
[Solid State Ionics Lab teaching page](https://ssi-westlake.com/teaching/).

The series is designed to help students explore how defect chemistry, ionic and
electronic transport, electrochemical polarization, space charge, and impedance
emerge from governing equations. Each module favors transparent physics,
interactive controls, and numerical checks over preassembled textbook curves.

## Notebooks

| Module | Topic | Interactive app |
|---|---|---|
| 01 | Brouwer diagram for weakly acceptor-doped SrTiO3 | [Launch](https://qiyanglu.github.io/Solid-State-Ionics-Interactive/) |

### Module 01: Brouwer diagram for SrTiO3

[`01_brouwer_sto.py`](01_brouwer_sto.py) solves, independently at every oxygen
partial pressure,

```text
K_red = V n^2 pO2^(1/2)
K_eh  = n p
2 V + p = A + n
```

for oxygen vacancies `V`, electrons `n`, and holes `p`, with a fixed fully
ionized acceptor concentration `A`. The unique positive solution is found with
a bracketed log-space root solve. Brouwer regimes and slopes are measured only
after solving the full equations; they are never used to construct the curves.

The default state is 973 K, A = 10^18 cm^-3, and pO2 from 10^-25 to 1 bar.
Controls span 700-1500 K and log10(A/cm^-3) from 13 to 21. The plot distinguishes
majority-carrier character from charge-compensation regimes so the electron-hole
crossover is not mistaken for a change in electroneutrality.

## Run locally

With Python 3.11 or newer and [uv](https://docs.astral.sh/uv/) installed:

```console
uv sync
uv run marimo edit 01_brouwer_sto.py
```

To present it as an app:

```console
uv run marimo run 01_brouwer_sto.py
```

## Validation

```console
uv run marimo check --strict 01_brouwer_sto.py
uv run marimo export html 01_brouwer_sto.py -o brouwer.html --no-include-code -f
```

The notebook displays mass-action, electroneutrality, positivity, finiteness,
regime-coverage, and limiting-slope checks below the plot.

## Browser deployment

The notebook uses only marimo, NumPy, SciPy, and matplotlib and performs no
network or filesystem access at runtime. The workflow in
[`pages.yml`](.github/workflows/pages.yml) exports it as a browser-hosted WASM
app and deploys the generated static site to GitHub Pages on every push to
`main`.

## Scientific scope

This first model assumes dilute ideal defects and carriers, fixed fully ionized
acceptors, and the quantitative equilibrium constants written in the notebook.
It intentionally excludes cation defects, complexes, non-ideal activities,
finite site-density constraints, and an external physics package.
