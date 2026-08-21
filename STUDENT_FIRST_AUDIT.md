# Student-first editorial and visual audit

## Audit contract

This document maps the existing seven notebooks to a student-first layered
reader. The governing physics, signs, units, and numerical checks remain
protected by AGENTS.md and NOTATION.md. The redesign does not add scientific
scope. It changes what students see first, how many ideas share one figure, and
where derivations and diagnostics live.

The intended core audience is an undergraduate or beginning graduate student
in Solid State Ionics who has not read the source code, lecture slides, or
tutorial articles. Each core reader should take roughly 10--20 minutes.

## Module 01 — Defect formation thermodynamics

**Main question.** Why can a crystal contain defects at equilibrium even when
forming each defect costs energy?

**Minimum core concepts.** A lattice picture; multiplicity; configurational
entropy; formation and entropy contributions to \(G\); the free-energy minimum;
the matching zero of \(\mu_D\); and the dilute limit.

**Move to Explore further.** The long thermodynamic-state table, finite-\(N\)
mode/ensemble-mean nuance, lattice-size control, exact-point display toggles,
and detailed approximation diagnostics.

### Existing visible figures and redesign

| Figure | Question and controls | Expected visible change | Current problem | Redesign |
|---|---|---|---|---|
| Randomized finite lattice | What does a defect fraction look like? \(T,\Delta h_f,\Delta s_f^0,N\) indirectly set \(n\). | Thermodynamic controls change the number of open defect sites; \(N\) changes grid resolution. | The long title embeds \(N,n,x\); the state table delays the physical picture. | Put this first. Title it “A finite lattice at its most probable composition”; move \(N,n,x\) to a caption and state that it is a configurational schematic, not an atomistic simulation. |
| Multiplicity and entropy, two panels | Why are mixed states numerous? \(N\). | Exact points become denser and approach the Stirling curves as \(N\) grows. | Appropriate conceptually, but finite-size control is mixed with the main thermodynamic controls. | Keep two panels in the core after the lattice; place \(N\) in a finite-size accordion. |
| Free energy and chemical potential, two panels | Where is equilibrium, and why do \(G\) minimum and \(\mu_D=0\) agree? \(T,\Delta h_f,\Delta s_f^0\); optional components/exact points. | The marked minimum and zero crossing move together; temperature and formation entropy favor larger \(x\). | This is the strongest figure, but toggles crowd the initial controls and the same derivation is repeated in prose. | Make it the central core figure with the three thermodynamic controls. Keep total \(G\) primary; put contribution/exact-point toggles in Explore further. |
| Exact/Stirling/dilute comparison | When is each approximation valid? \(N\), selected state. | Finite steps shrink with \(N\); the dilute curve converges at large \(\Delta g_f^0/k_BT\). | A useful final comparison, but its long annotation and repeated explanation compete with the core ending. | Keep as the final core or first advanced figure with a short explanation and no paragraph inside the axes. |

**Core narrative.** Physical question → lattice → multiplicity/entropy → total
free energy and minimum → derivative/chemical potential → approximation limits.

## Module 02 — Brouwer diagram for acceptor-doped SrTiO3

**Main question.** Can defect-chemistry regimes emerge from mass action and
exact electroneutrality rather than being drawn as assumed power laws?

**Minimum core concepts.** Four species, the two mass-action laws, exact
\(2V+p=A+n\), exact calculated concentrations, and one worked dominant balance.

**Move to Explore further.** Numerical equilibrium-constant values, activity
conventions, the thermodynamic-state table, detailed regime spans and slope
tables, and simultaneous display of all limiting annotations.

### Existing visible figure and redesign

| Figure | Question and controls | Expected visible change | Current problem | Redesign |
|---|---|---|---|---|
| Exact Brouwer explorer | How do \(V,n,p,A\) respond to oxygen pressure? \(T\), \(\log_{10}A\), slope-guide toggle. | Temperature shifts equilibria; acceptor level moves the vacancy plateau and crossovers; guides reveal limiting slopes. | Four curves are combined with two background tints, a plateau band, corner labels, a crossover label, external legend, and up to six slope annotations. Local line width is 2.8 pt. | Default to the four exact curves only, with direct labels or a compact external legend and 1.7 pt lines. Keep the \(n=p\) marker subtle. One Show limiting interpretation control reveals either regime shading or slope guides, not both. Explain the three balances below as concise cards. |

**Core narrative.** Species and exact neutrality → exact solution → identify one
plateau from the plotted concentrations → optionally reveal other limiting
guides → explain that bends arise because the full balance changes smoothly.

## Module 03 — Defect transport

**Main question.** How does a thermally activated one-dimensional hop become a
macroscopic diffusivity, and why must ions and electrons move together during
chemical diffusion?

**Minimum core concepts.** Activated hop, \(D=a^2\Gamma/2\), MSD, Fick's law,
field-biased hopping, electrochemical-potential cancellation, and the
slower-carrier bottleneck in \(D_{\rm Li}^{\delta}\).

**Move to Explore further.** Attempt frequency and jump distance controls,
periodic Fourier/master-equation implementation, detailed barrier/rate tables,
Haven-ratio detail, the conductivity-form derivation, and macroscopic slab
relaxation sweeps.

### Existing visible figures and redesign

| Figure | Question and controls | Expected visible change | Current problem | Redesign |
|---|---|---|---|---|
| Microscopic 2×2 pipeline | How do barrier, paths, MSD, and fitted \(D\) connect? \(T,\Delta H_{\rm mig},\nu,a\). | Higher \(T\) or lower barrier raises \(\Gamma\); \(D\propto a^2\Gamma\). | Four concepts share one figure. The random seed depends on physical controls, and automatic axes make causal changes hard to see. | Split into Figure A (one barrier plus a short 1D lattice/path) and Figure B (trajectories/distribution plus MSD). Use a fixed realization. Show analytical \(D\), fitted \(D\), and their ratio below; use baseline/fixed references so \(a\to2a\) and \(\Gamma\to10\Gamma\) are unmistakable. |
| Periodic profile and normalized bond flux | Why does random exchange give Fick's law? Time and step contrast. | The step smooths and the bond flux follows \(-Ddc/dx\). | The periodic wrap creates a second peak; self-normalized flux has no physical unit and confuses the first Fick-law explanation. | Replace the core with two adjacent populations and the difference between expected left-to-right and right-to-left crossings, then a finite non-periodic smoothing profile with stated reflecting boundaries. Retain the periodic solver only in Explore further. |
| Field/electrochemical 1×3 figure | How do field-biased barriers and chemical/electrical forces combine? Charge, field sign/magnitude, concentration gradient, balance toggle. | Field swaps directional rates; at selected balance the chemical and electrical terms cancel and total flux vanishes. | Barrier tilt, gradient decomposition, and flux bars are three lessons in one figure with five controls. | Split into two lessons: barrier tilt with field sign/magnitude; then cancellation of chemical and electrical driving forces. Put charge sign and gradient details in advanced controls. |
| Ambipolar two-panel figure | Why is chemical diffusion limited by the slower carrier? \(D_i\), \(D_e/D_i\). | \(D^\delta/D_i\) approaches \(2r/(1+r)\); the internal field makes ion and electron fluxes equal. | The physical bottleneck appears after substantial derivation and a bar comparison mixes uncoupled/coupled fluxes. | Lead with the bottleneck question. Keep one ratio plot and one simple matched-flux comparison; move the full conductivity derivation into Explore further. |
| Relaxation two-panel figure | How does \(D^\delta\) set a macroscopic time? \(L\), Fourier number. | Time grows as \(L^2\); the profile decays with Fourier time. | Correct but extends beyond the shortest core path and adds another parameter sweep. | Keep as Explore further — from diffusivity to sample response, preserving \(t_D\) versus \(\tau^\delta\). |

**Core narrative.** One hop → many fixed-realization walks → MSD and \(D\) →
net exchange and Fick flux → field bias → cancellation at equilibrium →
ion/electron bottleneck.

## Module 04 — Space-charge layers and the Frumkin effect

**Main question.** How does an interfacial charge redistribute mobile defects,
and how does that redistribution affect capacitance and reaction rate?

**Minimum core concepts.** Planar potential and concentration profiles,
electrochemical equilibrium, Gouy–Chapman versus Mott–Schottky, and the
screening/depletion length.

**Move to Explore further.** Chemical/electrical/electrochemical cancellation
plot, screening-length sweeps, full nonlinear derivations, the complete GCS
model, and the complete Frumkin section.

### Existing visible figures and redesign

| Figure | Question and controls | Expected visible change | Current problem | Redesign |
|---|---|---|---|---|
| Gouy–Chapman 1×3 | How do potential, concentrations, and potential components vary? \(T,c_\infty,\epsilon_r,z,\phi_0\), linear-overlay toggle. | Concentration/potential decay length and nonlinear surface enrichment change. | Three panels and six controls lead the core; cancellation is mixed into the profile lesson. | Use one two-panel core figure: potential and mobile-defect concentration. Add a model selector for GC/MS/compare. Move cancellation to a collapsed subsection. |
| Mott–Schottky 1×3 | How does frozen dopant charge create depletion? Same core controls. | Width and parabolic potential change with concentration, permittivity, and surface potential. | Repeats another dense three-panel block before direct comparison. | Use the same two-panel figure/mode selector so students compare like with like. |
| Screening/depletion lengths, two panels | How do bulk concentration and permittivity set width? \(c_\infty,\epsilon_r\). | Both length scales shrink with concentration and grow with permittivity. | Useful sweep but interrupts the core profile story. | Move to Explore further — what sets the length scale. |
| GCS potential/capacitance, two panels | How is voltage split between Stern and diffuse layers, and how do capacitances combine? \(C_s\) plus core state. | Near pZC the total follows the series combination and is Stern-limited when \(C_d\gg C_s\). | Only the introduction is accordion-wrapped; downstream controls, figure, and prose can remain visible. | Put controls, calculation, figure, and prose inside one collapsed GCS section. Preserve common charge and differential series capacitance. |
| Frumkin 2×2 | How do reaction-plane concentration and potential modify rate? \(\alpha,z_R\) plus GCS controls. | Signed reactant charge changes enrichment/depletion; \(\alpha\) changes kinetic weighting, not equilibrium profiles or capacitance. | Four panels are too dense and the transfer coefficient appears disconnected when students are looking at equilibrium GCS. | Put the complete section in one collapsed accordion. Split into at most two panels per view and state explicitly that \(\alpha\) is kinetic and does not alter GCS equilibrium or capacitance. |

**Core narrative.** Charged interface → linked potential/concentration → choose
GC or MS → compare length scales qualitatively → optional advanced GCS/Frumkin.

## Module 05 — Stoichiometry polarization

**Main question.** What happens to composition and measured voltage when ions
are blocked but electrons can pass through a mixed conductor?

**Minimum core concepts.** Geometry, immediate current partition, a selected
concentration profile, measured response over time, and the late steady state.

**Move to Explore further.** Factor-of-two derivation, constant-potential mode,
raw material controls, conductivity-ratio sweep, heat map and multi-time family,
and ion/electron potential decomposition.

### Existing visible figures and redesign

| Figure | Question and controls | Expected visible change | Current problem | Redesign |
|---|---|---|---|---|
| Conductivity-ratio two-panel sweep with twin axis | Which carrier controls current partition and \(D^\delta\)? Conductivity ratio plus material state. | Transference shifts; \(D^\delta\) and \(\tau^\delta\) reveal the slower carrier. | Appears before the sample geometry and uses a twin axis for diffusivity/time. | Move to Explore further; separate diffusivity and time if retained. |
| Transient 1×3 with twin response axis | How do profiles, history, and electrical response evolve? Nine controls including drive mode and time. | Selected time changes profile/heat map; current or voltage relaxes. | Six profiles, heat map, and voltage/current share one figure; twin axis and nine controls make it a dashboard. | Add a geometry schematic first. Default constant current. Core figure: selected \(c(x,t)/c_0\) profile plus \(U(t)\), at most three controls (strength preset, conductivity-ratio preset, time). Move heat map/profile family to advanced. |
| Potential decomposition 1×3 | How do chemical, electrical, and electrochemical potentials reconstruct voltage? Drive/time/material controls. | The components rearrange while their sums obey boundary conditions. | Three technical panels interrupt the core and require substantial notation. | Move to advanced. Use a selector for ion, electron, or comparison and show at most two panels. |

**Core narrative.** Geometry → immediate response → selected composition
profile → voltage relaxation at constant current → steady polarized state.

## Module 06 — PITT and GITT

**Main question.** How do controlled potential/current pulses and subsequent
OCV relaxation reveal chemical diffusion in a mixed conductor?

**Minimum core concepts.** Coulometric titration, selective-contact geometry,
one experiment at a time, pulse/rest timeline, one concentration profile, and
one measured response.

**Move to Explore further.** PITT/GITT comparison dashboard, potential
decomposition, asymptotic fitting tools, raw normalized controls, and all
finite-kinetics/Biot analysis.

### Existing visible figures and redesign

| Figure | Question and controls | Expected visible change | Current problem | Redesign |
|---|---|---|---|---|
| Coulometric titration curve | What equilibrium state does a composition step establish? Species label. | Li/H labels change while the same ideal neutral-pair relation remains. | Good opening concept, but the experiment transition can be shorter. | Keep as the overview figure with a concise bridge to transient methods. |
| Selective-contact schematic | Which carrier passes at each boundary? Species label. | Labels change consistently; boundary roles stay fixed. | Correct and useful, but separated from the pulse timeline. | Keep, followed immediately by a simple pulse/interruption/OCV timeline. |
| PITT/GITT 2×2 dashboard with twin axes | How do both experiments evolve? Eleven controls. | Pulse size, time, and diffusivity alter profiles and current/voltage transients. | Both techniques are shown simultaneously; four panels and twin axes make a reference dashboard. | Add a top-level PITT/GITT selector. Default to PITT. Show one concentration figure and one response figure with pulse size, duration, and material-speed preset. Compare only in advanced mode. |
| Composition/ion/electron potential 1×3 | How do potential components evolve in a selected case? Experiment, stage, and time controls. | Selected pulse/rest state changes composition and potential gradients. | Too technical for the core and three panels at once. | Move to advanced; provide a species/view selector and at most two panels. |
| PITT/GITT asymptotic comparison, two panels | Where do short- and long-time approximations work? Pulse parameters and \(D^\delta\). | Exact series approaches square-root and first-mode/linear limits in their domains. | Useful analysis appears before each experiment is individually understood. | Move to Explore further — classical limits, retaining notation consistent with the course articles. |
| Surface-kinetics 1×3 | When does finite exchange bias a diffusion-only fit? Biot/profile-time controls. | Small Bi slows current and changes profiles and fitted \(D^\delta\). | Advanced expert dashboard. | Keep all controls, prose, and figures inside Advanced analysis — when surface kinetics matters. |

**Core narrative.** Coulometric equilibrium → selective contacts → choose PITT
or GITT → pulse/rest timeline → profile → measured response → OCV relaxation.

## Module 07 — Impedance spectroscopy, Warburg diffusion, and TLM

**Main question.** How do phase, diffusion length, and spatial current transfer
appear in an impedance spectrum?

**Minimum core concepts.** \(e^{i\omega t}\), ideal \(R\) and \(C\), one
relaxation arc, semi-infinite Warburg penetration and 45-degree response, one
finite boundary at a time, and a clean two-rail TLM with one selected internal
view.

**Move to Explore further.** Series RC, two relaxations, Bode views, physical
Warburg scaling controls, boundary comparison, distributed TLM parameters,
zoomed spectra, current/potential/composition alternatives, and detailed TLM
scope equations.

The original qiyanglu/TLM-teaching-tool repository was reviewed as the
benchmark. Its most useful teaching features are the clear two-rail geometry,
explicit terminal conditions, linked spectrum/profile frequency markers, and
separation of a selected internal state from the dense spectrum. Its three
simultaneous Nyquist viewports, detailed four-terminal editor, dielectric
shunt, and implementation discussion are intentionally not all appropriate for
this course notebook's core reader.

### Existing visible figures and redesign

| Figure | Question and controls | Expected visible change | Current problem | Redesign |
|---|---|---|---|---|
| Waveform and phasor, two panels | How do frequency and phase appear in time and phasor views? Frequency and phase lead. | In a fixed time window, frequency changes cycle count; phase changes relative offset. | The current implementation uses a fixed two-cycle window, so frequency appears visually inert despite being numerically connected. | Use a fixed physical time window. Keep frequency and phase controls; optionally move normalized phase to advanced. |
| Series/parallel RC comparison, two Nyquist panels | How does topology change impedance? Fixed example values. | Series RC is a vertical line; parallel RC is a semicircle. | Pure resistor and capacitor are not taught as visible concepts first. | Core sequence: ideal resistor, ideal capacitor, then \(R\parallel C\). Keep series RC optional or advanced. |
| One/two-relaxation Nyquist and Bode with twin axis | What does each time constant add? Several raw RC controls. | A second time constant creates another arc or feature. | Too many raw controls and twin-axis Bode clutter. | Keep one \(R\parallel C\) core semicircle. Move two relaxations and Bode to advanced; use stacked Bode axes or a selector. |
| Semi-infinite Warburg 1×3 with twin Bode axis | Why is the Warburg line 45 degrees? Diffusion state. | Penetration shrinks with frequency; real and imaginary impedance remain equal. | Concentration, Nyquist, magnitude, and phase are shown simultaneously after dense formulas. | Start with the physical penetration picture. Use two panels: concentration-wave penetration and Nyquist. Put Bode in an optional view and explain 45 degrees after the figure. |
| Finite-length Warburg 1×3 | How does the far boundary change the low-frequency response? Boundary, frequency, phase, comparison, \(D^\delta,L,T,c,S\). | Fixed composition terminates resistively; zero flux turns upward capacitively; frequency changes penetration. | Nine controls and three panels compare too many states at once. | Use boundary, frequency, and a diffusion-speed preset. Show selected concentration profile and selected Nyquist response. Put the boundary equation beside the selector; move Bode and physical scale controls to advanced. |
| TLM schematic | Why are two rails coupled by distributed chemical storage? Contact preset. | Terminal symbols change with contact geometry. | The current schematic is useful, but it precedes a dense parameter/control block. | Retain a two-rail schematic closely matching the original tool, with vertical connections ending cleanly on capacitor plates and visually obvious terminal conditions. |
| TLM three-viewport Nyquist plus Bode | How do contacts and bulk ratio shape the spectrum? Contact case, conductivity ratio, distributed \(r,c,L\), frequency. | Contact preset and ratio alter arc/low-frequency behavior; selected frequency maps to profiles. | Full/high/low Nyquist viewports, dense colored scatter, colorbar, and separate Bode figure recreate a technical app rather than a lesson. | Use one equal-scale Nyquist plot with a Full/High/Low zoom selector, a few markers, and visible contact/ratio presets. Move distributed parameters and Bode to advanced. |
| TLM composition/potential/current 1×3 | What happens inside at the selected frequency? Frequency and phase plus TLM controls. | Frequency moves chemical penetration/current transfer; phase changes snapshots. | Three different internal observables and a long suptitle are shown at once. | Add a view selector: composition response, rail potentials, or rail currents. Show one or two panels and connect the selected marker directly to the Nyquist spectrum. |

**Core narrative.** Time-domain phase → ideal \(R\), \(C\), and one RC arc →
diffusion penetration and Warburg → finite boundary choice → why a lumped RC
cannot show carrier transfer → TLM schematic → one spectrum → one selected
internal state.

## Repository-wide acceptance map

- Core figures: no more than two panels and no twin y-axes.
- Core interactions: no more than three visible controls per figure.
- Random demonstrations: fixed realization under physical parameter changes.
- Checks: collapsed at the end with one-line physical explanations.
- Advanced content: calculations, controls, prose, and figures all contained in
  the same collapsed section rather than only collapsing an introduction.
- Browser QA: default, low/contrasting, and high/contrasting states at mobile,
  laptop, and projector sizes, recorded in VISUAL_QA.md.
- Physics: preserve all invariants and regression checks in AGENTS.md; never
  use limiting laws as solver inputs.
