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
| 01 | Finite lattice — what does equilibrium occupancy look like? | Thermodynamic controls change defect count; finite \(N\) changes spacing in the advanced view. | Pending | Baseline title is too long and the state table delays the picture. | Pending redesign |
| 01 | Multiplicity and entropy — why do mixed states dominate? | \(N\) makes exact points denser toward Stirling. | Pending | The two-panel concept is suitable after control simplification. | Pending redesign |
| 01 | Free energy and chemical potential — where is equilibrium? | \(T,\Delta h_f,\Delta s_f^0\) move the minimum and zero together. | Pending | The central figure is strong; initial controls and toggles need layering. | Pending redesign |
| 01 | Approximation limits — when is dilute valid? | The selected driving force and advanced \(N\) reveal convergence and discreteness. | Pending | Remove prose from the axes. | Pending redesign |
| 02 | Exact Brouwer curves — how do regimes emerge? | \(T,A\) shift curves; optional interpretation reveals one guide layer. | Pending | Baseline is annotation-heavy and uses 2.8 pt curves. | Pending redesign |
| 03 | One activated hop | \(T,\Delta H_{\rm mig}\) change the activated rate. | Pending | Split from the existing 2×2 pipeline. | Pending redesign |
| 03 | Many hops give diffusivity | A fixed realization plus \(a\) and \(\Gamma\) gives explicit expected \(D\) ratios. | Pending | Baseline seed and autoscaling obscure comparison. | Pending redesign |
| 03 | Net exchange gives Fick flux | Population contrast and time change a finite profile and interface flux. | Pending | Periodic wrap and normalized flux are confusing. | Pending redesign |
| 03 | Field-biased hopping | Field sign and magnitude swap directional preference. | Pending | Separate from electrochemical cancellation. | Pending redesign |
| 03 | Driving-force cancellation | Chemical/electrical imbalance changes total flux; balance gives zero. | Pending | Reduce controls and panels. | Pending redesign |
| 03 | Chemical-diffusion bottleneck | Conductivity-ratio preset changes \(D^\delta/D_i\) and matched flux. | Pending | Lead with the physical question. | Pending redesign |
| 04 | GC/MS profiles — how does interface charge redistribute defects? | Model selector, core potential, and bulk-state preset change profile and width. | Pending | Two current 1×3 figures repeat concepts. | Pending redesign |
| 05 | Polarization geometry | Boundary roles remain visually obvious. | Pending | A new core schematic is required. | Pending redesign |
| 05 | Selected concentration profile | Strength preset, ratio preset, and time change the gradient. | Pending | Baseline dashboard combines profiles, a heat map, and response. | Pending redesign |
| 05 | Constant-current voltage response | Strength and ratio change magnitude and relaxation. | Pending | Separate from current/constant-potential twin-axis mode. | Pending redesign |
| 06 | Coulometric titration | Composition points change equilibrium potential. | Pending | Baseline concept is usable. | Pending redesign |
| 06 | Pulse/rest timeline | Mode and duration change pulse, interruption, and rest spans. | Pending | A new simple schematic is required. | Pending redesign |
| 06 | Selected PITT or GITT profile | Pulse size, duration, and speed preset change penetration. | Pending | Baseline shows both methods simultaneously. | Pending redesign |
| 06 | Selected measured response | PITT current or GITT voltage responds to pulse and rest. | Pending | Remove core twin axes. | Pending redesign |
| 07 | Waveform and phasor | Frequency changes visible cycle count in fixed seconds; phase changes offset. | Pending | Baseline frequency appears inert. | Pending redesign |
| 07 | Ideal R, C, and one RC arc | Element selector or topology changes faithful Nyquist geometry. | Pending | Pure R and C need visible figures. | Pending redesign |
| 07 | Semi-infinite Warburg | Frequency changes penetration depth and the selected Nyquist location. | Pending | Baseline 1×3 plus twin axis is too dense. | Pending redesign |
| 07 | Finite Warburg boundary | Boundary, frequency, and speed preset change profile and low-frequency termination. | Pending | Baseline has nine controls and three panels. | Pending redesign |
| 07 | TLM schematic and spectrum | Contact and ratio presets change one equal-scale Nyquist spectrum. | Pending | Remove three simultaneous viewports and colorbar. | Pending redesign |
| 07 | Selected TLM internal state | Frequency and view selector change composition, potentials, or currents. | Pending | Baseline displays all three observables at once. | Pending redesign |
