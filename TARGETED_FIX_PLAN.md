# Targeted notebook correction plan

## Scope and invariants

This pass corrects the specific scientific, notation, plotting, and boundary-
condition issues identified after the student-first redesign. It preserves the
seven module filenames and routes, the existing scientific model scope, the
browser/WASM deployment, and the marimo + NumPy/SciPy + matplotlib dependency
set. Limiting formulas remain explanations or checks, never solver inputs.

## Module-by-module changes

### Module 01 — Defect Formation Thermodynamics

- Rename the visible entropy control to **Non-configurational formation
  entropy** and display its unit with $k_B$.
- Leave the physics, figures, and reader structure unchanged.

### Module 02 — Brouwer Diagram for Acceptor-Doped SrTiO3

- Repair all mixed raw/LaTeX species and unit notation in prose, controls,
  axes, legends, and checks.
- Use a plain-language acceptor-concentration control with the mathematical
  definition beside it.
- Replace the all-at-once limiting-guide checkbox with a default-off selector
  for one optional regime guide at a time.
- Update marimo generation metadata while preserving the exact equilibrium
  solver and its numerical concentration curves.

### Module 03 — Diffusion and Transport

- Separate the random-walk lesson into hop-space trajectories and a physical-
  time mean-square-displacement plot.
- Generate a fixed random ensemble independent of temperature, barrier, hop
  distance, and attempt frequency; map hop number to time only through the
  selected total hop frequency.
- Compare simulated and analytical diffusivities and retain the $a^2$ and
  $\Gamma$ scaling callouts.
- Preserve the existing one-dimensional Fick and mixed-conductor sections,
  while cleaning student-facing notation.

### Module 04 — Space Charge and the Frumkin Effect

- Give Gouy–Chapman and Mott–Schottky separate physically valid potential
  controls; Mott–Schottky remains on the positive depletion branch.
- Expand the Gouy–Chapman range to $-0.50$ to $+0.50$ V and add an explicit
  ideal-dilute-model warning at large dimensionless surface potential.
- Rename diffuse/space-charge capacitance to $C_{\mathrm{sc}}$ and distinguish
  it from total Gouy–Chapman–Stern capacitance.
- Rework the advanced capacitance view so the diffuse-only and full series
  GCS interpretations are not conflated, and extend GCS/Frumkin sweeps where
  numerically stable.
- Add direct checks of $C_{\mathrm{sc,pZC}}=\varepsilon/\lambda_D$, series
  capacitance, positivity, and finiteness.

### Module 05 — Stoichiometry Polarization

- Replace the single selected transient with a representative family at
  $t/\tau^\delta=0.02,0.10,0.50,2.0$ plus the steady profile.
- Use sequential low-saturation styling and highlight the profile nearest any
  retained time selector without drawing an extra curve.
- Label $x=0$, $L/2$, and $L$, mark both ion-blocking faces outside the data,
  and keep the response in its own figure.
- Verify the displayed profiles against the analytical boundary gradients.

### Module 06 — PITT and GITT

- Shorten and resize the titration figure; move explanation into Markdown.
- Correct the selective-contact schematic colors, width, and displayed
  coordinate labels.
- Replace the single before/after profile view with a pulse/OCV stage selector:
  four pulse profiles or three dashed OCV-rest profiles for the selected PITT
  or GITT experiment.
- Remove the erroneous profile reversal so displayed $x=0$ is the ion
  electrolyte and $x=L$ is the current collector.
- State only boundary properties imposed by the actual finite-volume model,
  and add coordinate-aware face-flux tests.
- Clean student-facing notation while keeping advanced asymptotic analysis
  collapsed.

### Module 07 — Impedance, Warburg Diffusion, and Transmission Lines

- Add one shared `set_equal_nyquist_limits` helper and apply it to every
  Nyquist plot, always retaining equal units and a visible negative real-axis
  margin when data approach zero.
- Insert a two-panel pure-resistor/pure-capacitor foundation before the RC
  comparison; keep the parallel-RC semicircle central and series RC secondary.
- Clean Warburg controls, labels, and phase notation.
- Keep only the two meaningful visible TLM contact presets; retain the fully
  reversible $R_e\parallel R_i$ limit only as a collapsed regression check.
- Keep the original TLM tool's two-rail/capacitor schematic anatomy and terminal
  mapping, with non-crossing capacitor plates and explicit distributed/total
  units.
- Add regression checks for a pure-capacitor line and the capacitive TLM case
  remaining inside equal-scale Nyquist axes.

## Cross-module notation and visual audit

- Search all visible prose, control labels, axes, titles, legends, and check
  tables for raw notation and unit strings identified in the correction brief.
- Use plain-language widget labels where math rendering is unreliable, with
  adjacent Markdown definitions.
- Preserve the current projector typography and low-saturation palette.

## Verification and evidence

1. Compile all seven notebooks.
2. Run strict marimo checks and execute/export all notebooks.
3. Export all seven WASM applications.
4. Run the existing and newly added scientific regression checks.
5. Run `git diff --check`.
6. Browser-test every changed core interaction at laptop and projector sizes,
   capturing default, low/contrasting, and high/contrasting states.
7. Record screenshots and findings in `VISUAL_QA.md`, then summarize the
   correction and validation results in `IMPLEMENTATION_REPORT.md`.

## Commit sequence

Use focused commits for: planning and Module 02; Module 03; Module 04;
Modules 05–06; Module 07; and final notation/browser QA documentation. Work
directly on `main`; do not create a feature branch.


## Completion status — 2026-08-21

All scoped corrections are implemented directly on `main`. The final gate
compiled and strictly checked all seven notebooks, executed all seven static
exports, rebuilt all seven WASM routes, verified every route's dynamic
`run-page` asset, rendered 52 PASS and 0 CHECK validation rows, and passed
`git diff --check`. Fresh browser evidence covers mobile, laptop, and projector
views plus default and contrasting control states; findings are recorded in
`VISUAL_QA.md`.
