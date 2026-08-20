# Project context

## Purpose

Solid State Ionics Interactive is a sequence of self-contained marimo readers
for classroom teaching. The notebooks complement the Solid State Ionics course
at Westlake University and its online teaching materials. Each module should
help a student connect a governing equation to a physical picture, an
interactive prediction, and a measurable consequence.

## Teaching philosophy

- Lead with a physical question and a prediction, then reveal the calculation.
- Derive curves from the stated model; do not assemble them from expected
  limiting slopes or asymptotes.
- Prefer one-dimensional models, explicit boundary conditions, and a short
  causal narrative over general software abstractions.
- Keep numerical implementation details out of the classroom path unless they
  explain a physical limitation. Put detailed checks in collapsed sections.
- Use checks to explain what physical relationship is protected, not merely to
  report numerical residuals.

## Module sequence

1. **Defect Formation Thermodynamics** — multiplicity, configurational
   entropy, free energy, chemical potential, and equilibrium defect fraction.
2. **Brouwer Diagram for Acceptor-Doped SrTiO3** — mass action plus exact
   electroneutrality and the emergence of Brouwer regimes.
3. **Defect Transport** — activated hopping, one-dimensional random walks,
   mean-square displacement, Fickian transport, and dilute chemical diffusion.
4. **Space-Charge Layers and the Frumkin Effect** — Gouy–Chapman,
   Mott–Schottky, Gouy–Chapman–Stern capacitance, and reaction-plane kinetics.
5. **Stoichiometry Polarization** — coupled ionic/electronic transport in a
   mixed conductor between ion-blocking electrodes.
6. **PITT and GITT** — selective contacts, voltage/current pulses, OCV
   relaxation, and short- and long-time diffusion limits.
7. **Impedance Spectroscopy and Transmission Lines** — ideal elements,
   relaxation arcs, Warburg diffusion, and a continuous two-rail MIEC model.

## Shared notation and conventions

`NOTATION.md` is the detailed symbol contract. Across the series:

- use \(\exp(\mathrm{i}\omega t)\), \(Z=Z'+\mathrm{i}Z''\), and plot
  \(Z'\) against \(-Z''\) for Nyquist diagrams;
- distinguish frequency \(f\) in Hz from angular frequency \(\omega\) in
  rad s\(^{-1}\);
- distinguish particle-scale quantities \((k_B,e)\) from molar quantities
  \((R,F)\);
- use \(D^\delta\) for chemical diffusivity and state the model that relates it
  to carrier diffusivities or conductivities;
- state coordinate direction, sign convention, units, and boundary conditions
  at each model boundary.

## Visualization conventions

- Figures must remain legible on a classroom projector and a laptop: large
  text, low-saturation colors, and uncluttered legends.
- Use line widths of about 1.5–2.0 pt. Combine color with line style, markers,
  or direct labels when curves must be distinguished.
- Capitalize plot labels and titles consistently; typeset scientific symbols
  with matplotlib math text.
- Every axis states a quantity and units, or explicitly identifies a
  dimensionless normalization.
- Every Nyquist plot uses equal horizontal and vertical scale.
- Place a prediction before an important interactive figure and a concise
  physical takeaway after it. Avoid placing prose over plotted data.

## Architecture and deployment

Each numbered Python file is one standalone marimo notebook. Physics functions
and constants live in dedicated cells; controls, solved data, plotting, and
teaching prose are downstream. Runtime dependencies are limited to marimo,
NumPy, SciPy, and matplotlib so the notebooks can be exported as browser-hosted
WASM apps.

`pages/index.html` is the course landing page. The GitHub Pages workflow checks
all seven notebooks, exports each route independently, deploys the static site,
and verifies each dynamic marimo entry module after deployment.

## Current limitations

- The models are intentionally idealized, planar, and mostly one-dimensional.
- Modules 01–02 use minimal defect species sets chosen for transparent
  thermodynamics and defect chemistry.
- Modules 03–06 use dilute or locally linearized transport where stated; they
  do not model phase transformations, porous electrodes, or arbitrary
  composition-dependent coefficients.
- Module 04 omits specific adsorption and full Marcus kinetics.
- Module 07 uses ideal capacitors and a uniform reduced transmission line; it
  does not fit experimental spectra or claim a unique equivalent circuit.

## Roadmap

Near-term work should improve classroom clarity, accessibility, and verified
browser behavior without broadening the scientific scope. New modules or
extensions should be added only when they serve a specific lecture objective
and can preserve the same physics-first, self-contained structure.
