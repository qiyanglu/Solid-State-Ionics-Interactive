# Student-first redesign implementation report

## Outcome

The seven-module course has been reorganized as a set of short, self-contained
student readers. The governing physics, numerical solvers, species sets, sign
conventions, units, and public routes were preserved. The redesign changes the
visible teaching sequence: one physical question, no more than three core
controls, and no more than two panels in a core figure.

Every notebook now has three layers:

1. a core reader intended for roughly 10--20 minutes;
2. collapsed **Explore further** material for derivations, parameter sweeps,
   alternate views, and advanced models; and
3. collapsed numerical and physical checks with a short explanation of why
   each check matters.

The editorial design is documented in
[STUDENT_FIRST_AUDIT.md](STUDENT_FIRST_AUDIT.md). Browser evidence and findings
are recorded in [VISUAL_QA.md](VISUAL_QA.md).

## Module-by-module result

| Module | Core reader after redesign | Material moved out of the core |
|---|---|---|
| 01 | Lattice picture, multiplicity/entropy, free-energy minimum, chemical-potential zero, and approximation limits | Finite-size controls, state tables, exact-point toggles, and detailed diagnostics |
| 02 | Four exact equilibrium curves from mass action plus \(2V+p=A+n\) | Equilibrium-constant table, regime spans, slope tables, and optional limiting guides |
| 03 | One activated hop, fixed-realization random walk, MSD, Fick flux, field bias, and the chemical-diffusion bottleneck | Periodic solver details, full electrochemical cancellation derivation, conductivity mapping, and slab relaxation |
| 04 | One GC/MS/compare selector with potential and concentration profiles | Electrochemical-potential decomposition, screening-length sweeps, complete GCS capacitance, and Frumkin kinetics |
| 05 | Blocking-electrode geometry, one selected concentration profile, and its measured response | Raw material controls, constant-potential extension, ratio sweeps, histories, and potential decomposition |
| 06 | Coulometric overview, selective contacts, pulse/rest timeline, one selected PITT or GITT profile, and one measured response | Potential decomposition, classical asymptotes, fitting windows, and finite-kinetics/Biot analysis |
| 07 | Fixed-time waveform, series/parallel RC, two-panel Warburg views, clean two-rail TLM, one spectrum, and one selected internal state | Overlapping relaxations, Bode views, physical scaling controls, boundary comparison, distributed TLM parameters, and application cards |

## Important classroom fixes

- Module 02 opens without decorative regime tints, plateau bands, corner labels,
  or six simultaneous slope annotations. The full-equilibrium curves remain the
  only concentration data.
- Module 03 keeps the random realization fixed while physical controls change,
  so a visible trend is not confused with sampling noise.
- Module 04 compares Gouy--Chapman and Mott--Schottky like-for-like in one
  two-panel figure. The Butler--Volmer transfer coefficient remains confined to
  the advanced kinetic section and does not alter equilibrium GCS capacitance.
- Module 05 begins with the actual blocking-electrode geometry and shows one
  time state rather than a dashboard.
- Module 06 presents PITT or GITT one at a time and makes current interruption
  plus continued OCV redistribution visually explicit.
- Module 07 uses a fixed two-second waveform window, so frequency visibly
  changes cycle count. Its finite Warburg and TLM readers no longer use
  three simultaneous viewports or a core twin axis.
- Browser review found that the default cross-selective TLM blocking tail
  compressed a single equal-scale Nyquist view. The default was changed to
  electron-reversible, ion-blocked contacts, producing a readable first arc;
  cross-selective and fully reversible limits remain selectable.

## Scientific scope preserved

No limiting Brouwer law, asymptotic PITT/GITT expression, Warburg reference
line, or TLM feature is used as solver input. All concentration, potential,
current, and impedance curves continue to come from their full documented
models.

The redesign did not add new species, a shared physics package, charged defects
to Module 01, extra dimensions, a porous-electrode model, a fitting framework,
or experimental parameter claims. Module-specific invariants in
[AGENTS.md](AGENTS.md) and notation in [NOTATION.md](NOTATION.md) remain the
release contract.

## Functional audit

The repository contains 81 marimo controls after the redesign:

| Module | Controls |
|---:|---:|
| 01 | 6 |
| 02 | 3 |
| 03 | 14 |
| 04 | 10 |
| 05 | 9 |
| 06 | 14 |
| 07 | 25 |

The higher counts in Modules 03, 06, and 07 are primarily collapsed advanced
controls. Strict marimo checks pass, and core sliders were exercised at default,
low/contrasting, and high/contrasting states in real browser-hosted WASM apps.

## Release validation

| Gate | Result |
|---|---|
| Python compilation | All seven notebooks compile. |
| Strict marimo check | All seven pass with no findings. |
| Clean static execution | All seven export successfully in isolated environments. |
| Embedded physics results | 52 PASS, 0 CHECK, 0 cell errors. |
| WASM export | All seven export in run mode with the Pages options. |
| WASM route assets | Every route has an index bundle and its dynamically imported `run-page` module. |
| Automated figure-layout scan | Every rendered figure checks title, label, tick, legend, colorbar, and annotation bounding boxes; conservative warnings are adjudicated against browser screenshots. |
| Browser visual QA | Default, low/contrasting, and high/contrasting states reviewed at mobile, laptop, and projector sizes; no release-blocking overlap or clipping remains. |
| Patch integrity | `git diff --check` passes. |

Per-module clean-execution results:

| Module | PASS | CHECK | Cell errors |
|---:|---:|---:|---:|
| 01 | 6 | 0 | 0 |
| 02 | 3 | 0 | 0 |
| 03 | 7 | 0 | 0 |
| 04 | 6 | 0 | 0 |
| 05 | 6 | 0 | 0 |
| 06 | 12 | 0 | 0 |
| 07 | 12 | 0 | 0 |

Generated validation exports and screenshots remain under the ignored
`dist/` directory. The evidence filenames and qualitative findings—not the
generated binaries—are versioned in [VISUAL_QA.md](VISUAL_QA.md).

## Commit and validation record

| Commit | Scope | Validation status |
|---|---|---|
| `473d6fc` | Student-first rules, audit map, style guide, and context | Documentation and repository-diff review passed. |
| `6a82aaf` | Reusable browser visual-QA helper | Python compilation and real-browser capture workflow passed. |
| `dcfd031` | Module 01 redesign | Strict check, static/WASM export, 6/0 embedded checks, and three-size browser review passed. |
| `d4a2818` | Module 03 redesign | Strict check, static/WASM export, 7/0 embedded checks, and three-size browser review passed. |
| `29a884a` | Module 05 redesign | Strict check, static/WASM export, 6/0 embedded checks, and three-size browser review passed. |
| `673e927` | Module 06 redesign | Strict check, static/WASM export, 12/0 embedded checks, and three-size browser review passed. |
| `4704b6d` | Module 07 redesign | Strict check, static/WASM export, 12/0 embedded checks, and three-size browser review passed. |
| `96c2082` | Module 02 redesign | Strict check, static/WASM export, 3/0 embedded checks, and three-size browser review passed. |
| `40c7024` | Module 04 redesign | Strict check, static/WASM export, 6/0 embedded checks, and three-size browser review passed. |
| This commit | Release-time matplotlib warning tool, final documentation, and release gate | Full seven-module compile, strict check, static/WASM export, 52/0 embedded checks, asset verification, and browser review passed. |
## Preserved limitations

- The models remain idealized and predominantly one-dimensional.
- Module 04 omits specific adsorption and full Marcus kinetics.
- Module 06 treats the classical short- and long-time formulas only as
  one-sided, small-signal limits.
- Module 07 uses ideal capacitors and a uniform reduced TLM, not a unique fit to
  an experimental spectrum.
- Public GitHub Pages content changes only after the committed main branch is
  pushed and the Pages workflow completes.