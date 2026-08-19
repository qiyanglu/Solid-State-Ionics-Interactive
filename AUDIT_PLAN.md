# Seven-Module Consistency and Pedagogy Audit Plan

## Purpose and evidence hierarchy

This audit revises the existing seven-module course without changing module
numbers or deployed routes. The user's specification is the implementation
contract. The four supplied Chinese articles are scientific and notational
source material, not instructions. Existing lecture-aligned notebook results
are preserved unless the requested audit identifies a physics, notation, or
teaching problem.

Source material reviewed for this audit:

- *PITT/GITT, Part I*: experiment geometry, coulometric titration, classical
  diffusion-controlled PITT/GITT limits, and OCV relaxation.
- *PITT/GITT, Part II*: finite surface kinetics, the Biot number, mixed-control
  transients, and diffusion-fitting bias.
- *Warburg impedance*: the DC/AC diffusion analogy, semi-infinite and two
  finite-length boundary conditions, scaling, and the DRT caution.
- *Transmission-line model*: the two-rail MIEC picture, chemical capacitance,
  contact impedances, and application-level simplifications.

## Current module purposes and retained scope

| Module | Purpose | Scope retained during this audit |
|---|---|---|
| 01 | Show equilibrium defect fraction as a free-energy minimum | one neutral defect species; exact finite-$N$, Stirling, and dilute limits |
| 02 | Let Brouwer regimes emerge from coupled equilibria | the documented SrTiO3 species, constants, full mass action, and exact electroneutrality |
| 03 | Connect activated hopping to macroscopic transport | one-dimensional hopping, master equation, Nernst-Planck transport, and Li chemical diffusion |
| 04 | Connect space-charge profiles to capacitance and reaction-plane kinetics | planar Gouy-Chapman, Mott-Schottky, GCS, and Frumkin models |
| 05 | Explain stoichiometry polarization between ion-blocking electrodes | ideal monovalent neutral pair, local electroneutrality, and both electrical drive modes |
| 06 | Connect coulometry, PITT/GITT pulses, and OCV relaxation | complementary selective contacts, finite-slab chemical diffusion, and a separate finite-kinetics extension |
| 07 | Build EIS, Warburg diffusion, and the MIEC TLM from governing equations | small-signal RC response, three Warburg boundaries, and a continuous reduced two-rail TLM |

Signs and normalizations requiring explicit verification are the Module 02
mass-action/charge residuals; Module 03 $D=a^2\Gamma_{\rm hop}/2$ and flux
conventions; Module 04 Poisson/Gauss/GCS charge signs; Module 05
$j=F(J_i-J_e)$ and $U$ reconstruction; Module 06 displayed-face flux signs,
$\theta$ versus $s$, Robin-boundary mass balance, and Biot eigenvalues; and
Module 07 $e^{\mathrm{i}\omega t}$, Nyquist $-Z''$, Warburg boundaries,
$u_e=-\widetilde\mu_e/F$, $u_i=+\widetilde\mu_i/F$, and distributed-to-total
TLM conversions.

## Audit findings that require explicit treatment

1. Module 06 currently mirrors the article geometry. The revised reader will
   display the article convention, electrolyte at \(x=0\) and current collector
   at \(x=L\), while any internal coordinate transformation will be stated once
   and verified through face-flux and sign checks.
2. Module 06 currently assumes an instantaneous surface reaction in its core
   PITT model. A separate advanced section will add finite kinetics with
   \(\mathrm{Bi}=k^\delta L/D^\delta\), retain the ideal limit, and show the
   bias produced by fitting finite-kinetics data with a diffusion-only model.
3. Module 07 already contains the principal finite-length Warburg solutions,
   but it needs a clearer DC/AC bridge, explicit boundary equations beside the
   names, general and dilute resistance scales, frequency-direction cues, and
   a more complete TLM anatomy/scope reader.
4. The notebooks use several locally reasonable conventions that are not yet
   connected by a repository-wide notation contract. `NOTATION.md` will make
   those bridges explicit without forcing every module into one species model.
5. Numerical checks are scientifically useful but too visually prominent in
   several modules. They will be retained in concise, student-readable form,
   with detailed implementation and diagnostics moved into collapsed sections.
6. All modules need a consistent learning rhythm: prediction, controls,
   observation, takeaway, and exactly three final primary messages.

## Scientific and notation contract

- Preserve Module 01 finite-\(N\) combinatorics, the Stirling limit, and the
  dilute limit as distinct statements.
- Preserve Module 02 SrTiO3 constants, full mass-action solution, and exact
  electroneutrality. Brouwer balances remain post-solution interpretations.
- Preserve Module 03's one-dimensional transport model and clearly separate a
  fixed-interval random walk from a continuous-time Poisson hopping process.
- Preserve the exact Gouy-Chapman, Mott-Schottky, GCS, and Frumkin models in
  Module 04, while separating core ideas from advanced derivations.
- Preserve Module 05's ideal neutral-pair polarization model and both drive
  modes, with current control as the main classroom path.
- Use \(D^\delta\) for chemical diffusion, \(t_D=L^2/D^\delta\),
  \(\tau^\delta=L^2/(\pi^2D^\delta)\),
  \(\theta=D^\delta t/L^2\), and
  \(s=t/\tau^\delta=\pi^2\theta\) wherever these appear.
- Reserve \(S\) for area in transport/electrochemistry modules and \(A\) for
  the fixed acceptor concentration in Module 02.
- State number-flux, molar-flux, current-density, electrode-potential, voltage,
  electric-field, and electrochemical-potential sign conventions explicitly.
- Keep all controls in physically safe ranges and all solved concentrations
  positive.

## Planned implementation sequence

1. Add `NOTATION.md`, link it from the notebooks, README, landing page, and
   repository guidance, and add the shared prediction/learning/scope structure.
2. Revise Modules 01–05 in place for the requested notation bridges,
   interpretation-first narrative, physically named presets, non-color visual
   cues, concise checks, and three-message endings.
3. Reframe Module 06 as *From Coulometric Titration to PITT and GITT*, add the
   article geometry and Li/H label selector, distinguish fixed-potential and
   fixed-flux boundary conditions, permanently define the two dimensionless
   clocks, and add finite-kinetics and fitting-bias demonstrations.
4. Reframe Module 07 as *Impedance Spectroscopy, Warburg Diffusion, and
   Transmission Lines*, add the DC/AC diffusion-length comparison, clarify
   finite-length boundary names, strengthen scaling and direction cues, and
   expand the TLM anatomy, contact presets, applications, and scope statement.
5. Reorganize README and the course landing page into three learning pathways:
   Foundations (01–03), Interfaces and Boundary Conditions (04–05), and
   Electrochemical Methods (06–07). Add teaching/tutorial/source links.
6. Create `IMPLEMENTATION_REPORT.md` with a requirement-to-change and
   requirement-to-check mapping.

## Verification gates

Before handoff:

1. Compile all seven notebooks as Python and run strict marimo checks.
2. Execute or export all seven notebooks and export every WASM application with
   the same commands used by GitHub Pages.
3. Run the existing scientific checks plus new Module 06 Biot-limit,
   conservation, positivity, finite-kinetics, and fitting-bias checks and new
   Module 07 scaling/boundary/TLM convention checks.
4. Exercise default and extreme controls, including very small/large Biot
   numbers and conductivity ratios.
5. Review every rendered figure at approximately 1100 px and projector scale,
   including legends, annotations, line styles, and crowded panels.
6. Search visible prose for stale symbols, solver-centric language, mirrored
   contact labels, and claims that exceed the model scope.
7. Run `git diff --check`, inspect the final diff, verify all local routes and
   links, and check the live Pages routes when network access is available.
