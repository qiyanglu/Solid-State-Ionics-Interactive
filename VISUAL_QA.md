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
and equal-scale Nyquist geometry. A matplotlib bounding-box scan supplies
warnings only; the browser screenshot remains the acceptance evidence.

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
| 04 | GC/MS profiles — how does interface charge redistribute defects? | Model selector, core potential, and bulk-state preset change profile and width. | Pending | Two current 1×3 figures repeat concepts. | Pending redesign |
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
