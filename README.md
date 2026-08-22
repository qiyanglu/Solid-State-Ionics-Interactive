# Solid State Ionics Interactive

[![Deploy marimo apps to GitHub Pages](https://github.com/qiyanglu/Solid-State-Ionics-Interactive/actions/workflows/pages.yml/badge.svg)](https://github.com/qiyanglu/Solid-State-Ionics-Interactive/actions/workflows/pages.yml)

**[Open the interactive course modules →](https://qiyanglu.github.io/Solid-State-Ionics-Interactive/)**

This repository contains interactive [marimo](https://marimo.io/) notebooks for
teaching **Solid State Ionics** at Westlake University. They complement the
[course teaching materials](https://ssi-westlake.com/teaching/) and
[tutorial collection](https://ssi-westlake.com/tutorial/).

The modules are designed as short, self-contained readers. Students can change
physical parameters, inspect the resulting plots, and see how familiar
solid-state ionics relationships emerge from thermodynamics and transport equations.
No software installation is needed when using the hosted versions.

## Course modules

| Module | Guiding question | Interactive app |
|---:|---|---|
| 01 — Defect Formation Thermodynamics | Why can defects be present at equilibrium even when forming one costs energy? | [Launch Module 01](https://qiyanglu.github.io/Solid-State-Ionics-Interactive/01-defect-formation/) |
| 02 — Brouwer Diagram for Acceptor-Doped SrTiO3 | How do defect-chemistry regimes and slopes emerge from equilibrium and charge neutrality? | [Launch Module 02](https://qiyanglu.github.io/Solid-State-Ionics-Interactive/02-brouwer-sto/) |
| 03 — Defect Transport | How do atomic hops become diffusion, conductivity, and coupled chemical transport? | [Launch Module 03](https://qiyanglu.github.io/Solid-State-Ionics-Interactive/03-defect-transport/) |
| 04 — Space-Charge Layers and the Frumkin Effect | How does an interface redistribute charged defects and alter capacitance and reaction rate? | [Launch Module 04](https://qiyanglu.github.io/Solid-State-Ionics-Interactive/04-space-charge-frumkin/) |
| 05 — Stoichiometry Polarization | How do blocking electrodes create time-dependent composition and potential profiles? | [Launch Module 05](https://qiyanglu.github.io/Solid-State-Ionics-Interactive/05-stoichiometry-polarization/) |
| 06 — PITT and GITT | What can voltage steps, current steps, and open-circuit relaxation reveal about thermodynamics and diffusion? | [Launch Module 06](https://qiyanglu.github.io/Solid-State-Ionics-Interactive/06-pitt-gitt/) |
| 07 — Impedance Spectroscopy and Transmission Lines | How do time scales, diffusion boundaries, and mixed conduction shape an impedance spectrum? | [Launch Module 07](https://qiyanglu.github.io/Solid-State-Ionics-Interactive/07-impedance-tlm/) |

The sequence follows three broad stages:

1. **Defects and transport:** Modules 01–03 build from equilibrium defect
   populations to microscopic hopping and macroscopic diffusion.
2. **Interfaces and polarization:** Modules 04–05 introduce space charge,
   capacitance, interfacial kinetics, and blocking-electrode polarization.
3. **Electrochemical methods:** Modules 06–07 connect time-domain titration
   experiments with frequency-domain impedance and transmission lines.

## Using the notebooks

Each module begins with a short core lesson and a small set of controls. More
detailed derivations and parameter explorations are collapsed under **Explore
further**, while **Physical consistency checks** explain why the calculation is
trustworthy without interrupting the main narrative.

The apps are intended to work on phones and laptops and remain readable on a
classroom projector. For the clearest learning path, follow the modules in
numbered order and make a prediction before moving each control.

## Run locally

With Python 3.11 or newer and [uv](https://docs.astral.sh/uv/) installed:

~~~console
uv sync
uv run marimo edit 01_defect_formation.py
~~~

Replace the filename with any of the seven numbered notebooks. Use `marimo run`
instead of `marimo edit` to present a notebook as an app.

## Repository guide

- `01_defect_formation.py` through `07_impedance_tlm.py` are the self-contained
  teaching notebooks.
- `pages/index.html` is the course landing page deployed by GitHub Pages.
- [`docs/NOTATION.md`](docs/NOTATION.md) connects symbols and sign conventions
  across the course.
- [`docs/PROJECT_CONTEXT.md`](docs/PROJECT_CONTEXT.md) and
  [`docs/STYLE_GUIDE.md`](docs/STYLE_GUIDE.md) record the shared teaching and
  visual conventions.
- [`docs/STUDENT_FIRST_AUDIT.md`](docs/STUDENT_FIRST_AUDIT.md),
  [`docs/VISUAL_QA.md`](docs/VISUAL_QA.md), and
  [`docs/IMPLEMENTATION_REPORT.md`](docs/IMPLEMENTATION_REPORT.md) contain the
  detailed editorial and validation record.

All browser apps are rebuilt and deployed automatically when `main` is updated.
