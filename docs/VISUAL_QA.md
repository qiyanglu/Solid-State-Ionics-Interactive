# Browser visual QA record

## Protocol

Every notebook is exported as the same WASM app used by GitHub Pages and served
over HTTP. Each core interactive figure is inspected in a real browser at:

- mobile: 390 × 844;
- laptop: 1280 × 800;
- projector: 1920 × 1080.

For each interaction, screenshots cover the default state and low or contrasting
and high or contrasting settings. Screenshot paths are relative to the
generated, ignored dist/visual-qa/ directory. The screenshots are evidence
rather than repository source files.

The review checks responsive control wrapping, clipped text, title/label/tick
collisions, legends and annotations over data, meaningful parameter response,
and equal-scale Nyquist geometry. Release exports enable the lightweight
`scripts/matplotlib_layout_qa/sitecustomize.py` draw-time scan through
`PYTHONPATH`. It inspects titles, labels, tick labels, legends, colorbars, and
annotations without adding QA code to the notebooks or their WASM runtime. The
scan stores possible intersections or clipping on each figure and emits warnings
only. Endpoint tick boxes can generate conservative warnings, so every warning
is adjudicated in the browser; screenshots remain the acceptance evidence.

## Status legend

- **Pending redesign** — baseline finding recorded; implementation not ready.
- **Pass** — all three browser sizes and contrasting states inspected.
- **Blocked** — a visible defect remains and the module must not be committed.

## Module records

The initial entries below define the screenshot matrix. Paths and findings are
filled during each focused module pass.

| Module | Core figure and physical question | Core controls and expected change | Screenshot evidence | Findings | Status |
|---|---|---|---|---|---|
| 01 | Finite lattice — what does equilibrium occupancy look like? | Thermodynamic controls change defect count; finite \(N\) changes spacing in the advanced view. | `module-01/default-laptop.png`; `module-01/free-energy-low-mobile.jpg`; `module-01/free-energy-high-projector.jpg` | Concise title and caption keep \(N,n,x\) visible. Randomized sites remain deterministic for comparison. A mobile menu-button warning was visually inspected and is a false positive; the button clears the title. | **Pass** |
| 01 | Multiplicity and entropy — why do mixed states dominate? | \(N\) makes exact points denser toward Stirling. | `module-01/entropy-laptop.png`; `module-01/free-energy-low-mobile.jpg`; `module-01/free-energy-high-projector.jpg` | Exact points and Stirling lines remain distinguishable; labels, legends, and axes do not overlap at the three sizes. | **Pass** |
| 01 | Free energy and chemical potential — where is equilibrium? | \(T,\Delta h_f,\Delta s_f^0\) move the minimum and zero together. | `module-01/free-energy-default-laptop.jpg`; `module-01/free-energy-low-mobile.jpg`; `module-01/free-energy-high-projector.jpg` | The minimum and zero crossing move together from dilute to high-occupancy states. Projector labels are clear; the mobile view is compact but unclipped and zoomable. | **Pass** |
| 01 | Approximation limits — when is dilute valid? | The selected driving force and advanced \(N\) reveal convergence and discreteness. | `module-01/free-energy-low-mobile.jpg`; `module-01/free-energy-high-projector.jpg` | The long in-axes paragraph is gone; the dilute and finite-\(N\) comparison remains readable without obscuring data. | **Pass** |
| 02 | Exact Brouwer curves — how do regimes emerge? | (T,A) shift all full-equilibrium curves; optional guides compare limiting slopes afterward. | `module-02/brouwer-default-laptop.jpg`; `module-02/brouwer-low-mobile.jpg`; `module-02/brouwer-high-projector.jpg` | Decorative regime tints, corner labels, and the heavy default guide layer are removed. Four low-saturation curves, one derived (n=p) marker, and a compact legend remain clear across all three sizes with no clipping or overlap warnings. | **Pass** |
| 03 | One activated hop | \(T,\Delta H_{\rm mig}\) change the activated rate. | `module-03/hopping-default-laptop.jpg`; `module-03/hopping-low-mobile.jpg`; `module-03/hopping-high-projector.jpg` | The barrier remains clear while its height changes; no clipped labels or overlaps at all three sizes. | **Pass** |
| 03 | Many hops give diffusivity | A fixed realization plus \(a\) and \(\Gamma\) gives explicit expected \(D\) ratios. | `module-03/hopping-default-laptop.jpg`; `module-03/hopping-low-mobile.jpg`; `module-03/hopping-high-projector.jpg` | The random realization remains fixed. Physical time and the analytical/fit ratio respond without a redundant bar panel. | **Pass** |
| 03 | Net exchange gives Fick flux | Population contrast and time change a finite profile and interface flux. | `module-03/fick-default-laptop.jpg`; `module-03/fick-low-mobile.jpg`; `module-03/fick-high-projector.jpg` | Reflecting ends remove the wrap-around peak. The profile smooths with time and bond exchange has physical units; controls and axes are unclipped. | **Pass** |
| 03 | Field-biased hopping | Field sign and magnitude swap directional preference. | `module-03/field-default-laptop.jpg`; `module-03/field-low-mobile.jpg`; `module-03/field-high-projector.jpg` | A single barrier plot makes the magnitude response legible. Direction reversal was also exercised and reverses the site-energy tilt. | **Pass** |
| 03 | Driving-force cancellation | The advanced gradient control changes equal-and-opposite chemical/electrical terms. | `module-03/field-default-laptop.jpg`; `module-03/field-low-mobile.jpg`; `module-03/field-high-projector.jpg` | Moved into a collapsed two-panel exploration, so it no longer competes with the core barrier lesson. The zero-flux identity remains checked. | **Pass** |
| 03 | Chemical-diffusion bottleneck | The diffusivity ratio changes \(D^\delta/D_i\) and matched carrier flux. | `module-03/chemical-default-laptop.jpg`; `module-03/chemical-low-mobile.jpg`; `module-03/chemical-high-projector.jpg` | The bottleneck question precedes the plot. Electron-slow and electron-fast limits are distinct; the mobile plot is compact but unclipped. | **Pass** |
| 04 | GC/MS profiles — how does interface charge redistribute defects? | Model selector, core potential, and bulk concentration change potential, concentrations, and screening length. | `module-04/profiles-default-laptop.jpg`; `module-04/profiles-low-mobile.jpg`; `module-04/profiles-high-projector.jpg` | One two-panel comparison replaces two simultaneous three-panel dashboards. GC and MS use distinct low-saturation styles; high concentration visibly narrows both layers and high core potential increases redistribution. Detailed electrochemical-potential, screening-length, GCS, and Frumkin views remain collapsed. No clipping or overlap warnings occur. | **Pass** |
| 05 | Polarization geometry | Boundary roles remain visually obvious. | `module-05/geometry-default-laptop.jpg`; `module-05/polarization-early-mobile.jpg`; `module-05/polarization-late-projector.jpg` | The slab and two ion-blocking/electron-passing contacts appear before equations; current direction and coordinate are unclipped. | **Pass** |
| 05 | Selected concentration profile | Strength preset, conductivity-ratio preset, and time change the gradient. | `module-05/polarization-default-laptop.jpg`; `module-05/polarization-early-mobile.jpg`; `module-05/polarization-late-projector.jpg` | The early state is nearly uniform and the late state approaches the blocked linear profile. One selected curve replaces the former profile family and heat map. | **Pass** |
| 05 | Constant-current voltage response | Strength and ratio change magnitude; time selects the matching state. | `module-05/polarization-default-laptop.jpg`; `module-05/polarization-early-mobile.jpg`; `module-05/polarization-late-projector.jpg` | The single voltage curve rises above its initial ohmic value and shares a selected-time marker with the profile. No twin axis, clipping, or overlap remains. | **Pass** |
| 06 | Coulometric titration | Composition points change equilibrium potential. | `module-06/titration-default-laptop.jpg` | Equilibrium composition and potential are introduced before the pulse methods; the illustrative curve and charge-to-composition narrative remain unclipped. | **Pass** |
| 06 | Pulse/rest timeline | Mode and duration define pulse, interruption, and rest. | `module-06/timeline-default-laptop.jpg` | A single low-ink timeline makes current interruption and continued OCV relaxation explicit without numerical detail. | **Pass** |
| 06 | Selected PITT or GITT profile | Mode, pulse size, duration, and (D^\delta) change the inherited profile. | `module-06/pitt-default-laptop.jpg`; `module-06/pitt-low-mobile.jpg`; `module-06/pitt-high-projector.jpg` | One experiment is shown at a time. End-of-pulse and rested profiles remain distinct from low to high drive; no clipping or overlap warnings occur. | **Pass** |
| 06 | Selected measured response | PITT current or GITT voltage responds to pulse and rest. | `module-06/pitt-default-laptop.jpg`; `module-06/pitt-low-mobile.jpg`; `module-06/pitt-high-projector.jpg` | The core response uses one y-axis: current decay and exact interruption for PITT, or potential evolution and OCV relaxation for GITT. Mobile text is readable and projector geometry is clear. | **Pass** |
| 07 | Waveform and phasor | Frequency changes visible cycle count in a fixed two-second window; phase changes offset. | `module-07/waveform-default-laptop.jpg`; `module-07/waveform-low-mobile.jpg`; `module-07/waveform-high-projector.jpg` | The low state shows about one cycle and the high state twenty cycles on the same physical time axis. Phasor angle stays independent of frequency; no overlap or clipping warnings occur. | **Pass** |
| 07 | Series and parallel RC | The same ideal (R,C) values produce a vertical series response or a parallel semicircle. | `module-07/rc-default-laptop.jpg` | Both Nyquist panels use equal scale, the (omega RC=1) markers are clear, and the pure-element limits are stated immediately above. Overlapping arcs and Bode views are optional. | **Pass** |
| 07 | Semi-infinite Warburg | The concentration-wave picture produces the (45^\circ) Nyquist line. | `module-07/warburg-default-laptop.jpg`; `module-07/warburg-low-mobile.jpg`; `module-07/warburg-high-projector.jpg` | The former three-panel/twin-axis dashboard is now a two-panel causal link from spatial decay to impedance. Labels remain clear at all three sizes. | **Pass** |
| 07 | Finite Warburg boundary | Boundary, frequency, and (D^\delta) change penetration and the low-frequency termination. | `module-07/warburg-default-laptop.jpg`; `module-07/warburg-low-mobile.jpg`; `module-07/warburg-high-projector.jpg` | Three core controls drive one profile and one equal-scale Nyquist plot. Optional phase, geometry, comparison, scale, and limiting formulas stay collapsed. | **Pass** |
| 07 | TLM schematic and spectrum | Contact, conductivity ratio, and frequency change one equal-scale Nyquist spectrum. | `module-07/tlm-schematic-default-laptop.jpg`; `module-07/tlm-spectrum-default-laptop.jpg`; `module-07/tlm-low-mobile.jpg`; `module-07/tlm-high-projector.jpg` | The schematic preserves the original two-rail teaching anatomy without vertical strokes through the (c_{\rm chem}) label. Browser review caught and fixed a compressed default: electron-reversible/ion-blocked contacts now open with a legible arc, while cross-selective blocking remains selectable. | **Pass** |
| 07 | Selected TLM internal state | Frequency and view selector change composition, potentials, or currents. | `module-07/tlm-internal-default-laptop.jpg`; `module-07/tlm-internal-low-mobile.jpg`; `module-07/tlm-internal-high-projector.jpg` | One selected internal observable replaces the former simultaneous three-panel view. The composition envelope remains the default bridge to Warburg; no clipping or overlap warnings occur. | **Pass** |


## Targeted correction pass — 2026-08-21

This section supersedes the changed figure descriptions in the historical
student-first table above. The final apps were rebuilt from the current source
and inspected through headless Chrome over HTTP. The in-app browser controller
could not attach on this Windows host because of a local ACL restriction, so
the repository's Chrome DevTools capture helper provided equivalent real-browser
screenshots. Every reported layout check returned no clipping. On mobile, the
lightweight overlap scan sometimes associated marimo's top-right menu button
with a nearby heading; direct inspection confirmed that the button clears the
text in every case.

| Module | Final targeted evidence | Finding | Status |
|---:|---|---|---:|
| 01 | `targeted-fixes/module-01-default-mobile.jpg`; `targeted-fixes/module-01-default-laptop.jpg` | The entropy widget now has a plain-language label and the adjacent text defines $\Delta s_f^0/k_B$, including the meaning of a value of 3. The randomized lattice and controls remain clear. | **Pass** |
| 02 | `targeted-fixes/module-02-default-mobile.jpg`; `targeted-fixes/module-02-default-laptop.jpg`; `targeted-fixes/module-02-low-laptop.jpg`; `targeted-fixes/module-02-high-projector.jpg`; `targeted-fixes/module-02-oxidizing-guide-projector.jpg` | Species and units render consistently; the acceptor control is defined beside the widget; one default-off limiting guide can be shown without crowding or altering exact solved curves. | **Pass** |
| 03 | `targeted-fixes/module-03-msd-default-mobile.jpg`; `targeted-fixes/module-03-default-laptop.jpg`; `targeted-fixes/module-03-low-laptop.jpg`; `targeted-fixes/module-03-high-projector.jpg`; `targeted-fixes/module-03-msd-default-laptop.jpg` | Fixed hop-space trajectories no longer change with physical parameters. A separate physical-time MSD figure follows $\langle x^2\rangle=2Dt$ and retains readable labels at all sizes. | **Pass** |
| 04 | `targeted-fixes/module-04-default-mobile.jpg`; `targeted-fixes/module-04-default-laptop.jpg`; `targeted-fixes/module-04-low-laptop.jpg`; `targeted-fixes/module-04-high-projector.jpg` | Gouy--Chapman and Mott--Schottky controls stay in their valid potential domains. The large-field ideal-model warning is visible, and $C_{\rm sc}$ is distinguished from the series GCS capacitance. | **Pass** |
| 05 | `targeted-fixes/module-05-default-mobile.jpg`; `targeted-fixes/module-05-default-laptop.jpg`; `targeted-fixes/module-05-early-laptop.jpg`; `targeted-fixes/module-05-late-projector.jpg` | The representative transient family makes profile evolution visible, the selector highlights rather than adds a curve, both blocking faces are labelled, and the voltage response remains separate. | **Pass** |
| 06 | `targeted-fixes/module-06-pulse-default-mobile.jpg`; `targeted-fixes/module-06-contacts-projector.jpg`; `targeted-fixes/module-06-pulse-default-laptop.jpg`; `targeted-fixes/module-06-pulse-low-laptop.jpg`; `targeted-fixes/module-06-pulse-high-projector.jpg`; `targeted-fixes/module-06-ocv-default-laptop.jpg`; `targeted-fixes/module-06-gitt-pulse-laptop.jpg`; `targeted-fixes/module-06-gitt-ocv-projector.jpg` | The public coordinate is consistent from ion electrolyte to current collector. Pulse and OCV stages show truthful profile families for both PITT and GITT, with no misleading numerical-method narrative. | **Pass** |
| 07 | `targeted-fixes/module-07-wave-default-laptop.jpg`; `targeted-fixes/module-07-wave-low-laptop.jpg`; `targeted-fixes/module-07-wave-high-projector.jpg`; `targeted-fixes/module-07-ideal-elements-mobile.jpg`; `targeted-fixes/module-07-ideal-elements-laptop.jpg`; `targeted-fixes/module-07-parallel-rc-projector.jpg`; `targeted-fixes/module-07-warburg-default-laptop.jpg`; `targeted-fixes/module-07-warburg-blocked-projector.jpg`; `targeted-fixes/module-07-tlm-cross-low-laptop.jpg`; `targeted-fixes/module-07-tlm-cross-high-projector.jpg`; `targeted-fixes/module-07-tlm-schematic-projector.jpg` | Frequency changes cycle count on a seconds axis; ideal elements precede the parallel-RC semicircle; every Nyquist view preserves equal scale and visible padding; finite Warburg and the two visible TLM contact cases remain readable; capacitor strokes do not cross the $c_{\rm chem}$ label. | **Pass** |

The complete targeted screenshot set remains under the ignored
`dist/visual-qa/targeted-fixes/` directory. Low/default/high states were
exercised for each changed core interaction; the table lists the most useful
acceptance views rather than every generated file.


## Module 04 Gouy–Chapman range follow-up — 2026-08-21

The core reader now opens with Gouy–Chapman over
$-1.00\leq\phi_0\leq1.00\ \mathrm{V}$ and overlays the exact planar solution
with the Debye–Hückel low-potential exponential by default. The Mott–Schottky
control remains restricted to its positive depletion branch.

| State | Evidence | Finding | Status |
|---|---|---|---:|
| Default, $\phi_0=0.20$ V | `module-04-gc-range/default-laptop.jpg` | The two potential curves have begun to separate, while controls, legends, and both panels remain unclipped. | **Pass** |
| Low potential, $\phi_0=0.01$ V | `module-04-gc-range/low-potential-mobile.jpg`; `module-04-gc-range/low-potential-laptop.jpg` | Exact and Debye–Hückel potential profiles are visually coincident, as required by $|ze\phi_0/(k_BT)|\ll1$. The mobile menu-button overlap warning is a visually adjudicated false positive; the button clears the heading. | **Pass** |
| High potential, $\phi_0=1.00$ V | `module-04-gc-range/high-potential-projector.jpg`; `module-04-gc-range/high-potential-caution-projector.jpg` | The exact profile screens much more sharply than the exponential guide, counter-ion enrichment spans many decades, and the ideal-dilute-model caution is fully visible. | **Pass** |

All captures used the exported WASM app in a real Chromium browser. Layout
checks reported no clipping, and the executed notebook rendered 7 PASS and
0 CHECK physical-consistency rows, including explicit low-potential convergence
and high-potential separation tests.


## Full interactive-control and readability audit — 2026-08-21

All seven final WASM apps were rebuilt, served over HTTP, and executed one at a
time in Chromium. Sequential loading avoids mistaking a resource-contended
Pyodide hourglass state for a finished notebook. Acceptance captures were taken
only after the actual matplotlib outputs were present.

A source-level reactive audit found **83 controls and zero controls without a
value consumer**. Browser checks then exercised the conditional and relocated
controls most likely to appear as ghost widgets:

- Module 03 shows charge sign, field direction, and field magnitude beside the
  hopping figure; the concentration-gradient slider is beside the advanced
  cancellation figure it changes.
- Module 04 shows the Debye–Hückel guide only for Gouy–Chapman states. In a pure
  Mott–Schottky state, both the checkbox and the words “linear guide” disappear.
- Module 05 shows only the active current or potential strength in the core
  controls and updates the section narrative with the drive mode.
- Module 06 places stage progress directly above the composition and potential
  profiles it changes.
- Module 07 hides the second-relaxation ratio and separation sliders until the
  second parallel-RC branch is enabled; the TLM snapshot phase sits beside its
  internal-view selector and figure.

| Audit view | Final evidence | Finding | Status |
|---|---|---|---:|
| Rebuilt default state, all modules | `full-audit-2026-08-21/01-final-default-laptop.jpg` through `07-final-default-laptop.jpg` | Every app finished executing with its controls and figures present. No confirmed clipping or overlap remains. | **Pass** |
| Module 03 field and cancellation controls | `full-audit-2026-08-21/03-field-default-laptop.jpg`; `03-field-low-mobile.jpg`; `03-field-high-projector.jpg`; `03-gradient-laptop.jpg` | Field magnitude changes the hopping bias, and the advanced gradient control is co-located with the equal-and-opposite driving-force plot. | **Pass** |
| Module 04 model-dependent controls and GCS | `full-audit-2026-08-21/04-gcs-laptop.jpg`; `04-ms-projector.jpg`; `module-04-gc-range/low-potential-mobile.jpg`; `module-04-gc-range/high-potential-projector.jpg` | The exact and low-potential GC curves coincide at low field and separate at high field. Pure MS has no inapplicable linear-guide widget. GCS controls, equations, and figure form one continuous advanced section. | **Pass** |
| Module 05 drive modes | `full-audit-2026-08-21/05-current-laptop.jpg`; `05-potential-projector.jpg` | Switching drive mode replaces the active strength widget and narrative; no inactive current or voltage widget remains in the core row. | **Pass** |
| Module 06 OCV potentials | `full-audit-2026-08-21/06-ocv-potentials-projector.jpg` | Stage progress is visibly attached to the three advanced profiles; the OCV state retains zero terminal current while internal chemical relaxation remains visible. | **Pass** |
| Module 07 optional RC and TLM states | `full-audit-2026-08-21/07-rc-second-on-laptop.jpg`; `07-rc-second-off-mobile.jpg`; `07-tlm-phase0-mobile.jpg`; `07-tlm-phase300-projector.jpg`; `07-final-tlm-schematic-projector.jpg` | Disabled second-arc controls disappear, the phase slider changes the selected spatial state, and the final two-rail schematic has clean capacitor connections and no line through the $c_{\rm chem}$ label. | **Pass** |

The lightweight mobile overlap scan flagged two slider value bubbles near their
section headings. Direct inspection of `03-field-low-mobile.jpg` and
`07-tlm-phase0-mobile.jpg` confirmed clear vertical separation, so both were
false positives. All other final layout scans reported no clipping or overlap.

The final browser-executed physical checks report **53 PASS and 0 CHECK**:
Modules 01–07 contribute 6, 3, 7, 7, 6, 12, and 12 passing rows. Every exported
route also contains its entry bundle and dynamically imported `run-page` asset.
