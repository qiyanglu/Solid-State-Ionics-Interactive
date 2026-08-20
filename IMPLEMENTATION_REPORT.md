# Second-Stage Scientific, Functional, and Visual QA Report

## Scope and outcome

This pass audited all seven existing marimo teaching notebooks without adding a
new scientific model or changing any public route. The work focused on
scientific consistency, control responsiveness, classroom narrative, visual
legibility, and browser/WASM deployment.

The original physical scope was preserved. The substantial teaching revisions
are confined to Module 03 and Module 07; the remaining modules received shared
typography, line-weight, label, color, figure-order, and overlap corrections.

## Stable project context

PROJECT_CONTEXT.md now records the project purpose, teaching philosophy, module
sequence, shared notation, visualization rules, architecture, limitations, and
roadmap. AGENTS.md now states that this context file is updated only when
scientific scope, repository architecture, shared notation, or the roadmap
changes.

## Repository-wide figure QA

- Increased the shared matplotlib typography to improve laptop and projector
  readability.
- Standardized data lines to approximately 1.5--2.0 pt.
- Capitalized titles and axis labels consistently and retained units or an
  explicit dimensionless normalization.
- Used restrained colors together with dashed, dotted, marker, or hatch cues
  where curves need a non-color distinction.
- Removed annotations that collided with data or legends and shortened or
  wrapped crowded titles.
- Enforced equal horizontal and vertical data scale for every Nyquist diagram.
- Reviewed all 33 executed figures in module contact sheets and inspected the
  revised Module 03 and Module 07 figures at full resolution.

## Module 03: activated hopping to diffusivity

The microscopic sequence is now explicit:

\[
\text{activated hop}
\rightarrow \text{random walk}
\rightarrow \langle x^2\rangle
\rightarrow \text{MSD slope}
\rightarrow D.
\]

The trajectory panel now uses physical position and physical time. A
deterministic seed derived from \(T\), \(\Delta H_{\rm mig}\), \(\nu\), and
\(a\) regenerates the displayed realization whenever any microscopic control
changes. Six thousand walkers reduce sampling noise while remaining suitable
for browser execution.

The MSD is fitted through the physical origin and compared with
\(\langle x^2\rangle=2Dt\). A separate panel compares the fitted diffusivity
with the analytical one-dimensional result \(D=a^2\Gamma/2\).

The prose now reads as one standalone article rather than a commentary on
lecture slides. Lecture-dependent headings and the disconnected advanced aside were
removed. The Li chemical diffusivity is derived in one sequence
from equal carrier fluxes, through the neutral chemical-potential gradient, to
the ideal dilute harmonic-mean result.

Independent low, default, and high parameter checkpoints changed the displayed
time span from \(3.77\times10^{-1}\) s to \(8.81\times10^{-11}\) s and the
largest sampled displacement from 17.5 to 62.3 nm. The fitted diffusivities
agreed with the analytical values within 0.31--2.67%, with
\(R^2\ge 0.99944\).

## Module 07: ideal elements, Warburg diffusion, and TLM

The teaching order is now:

1. time-domain sinusoid and phasor convention;
2. series- and parallel-RC Nyquist signatures;
3. one or two parallel-RC relaxations;
4. semi-infinite Warburg concentration, Nyquist, and Bode views;
5. interactive finite fixed-composition and zero-flux boundaries;
6. transport-equation-to-circuit mapping for the continuous two-rail MIEC TLM.

The resistor and capacitor limits are calculated explicitly from
\(Z_R=R\) and \(Z_C=1/(\mathrm{i}\omega C)\). The semi-infinite Warburg figure
shows the concentration wave, the derived 45-degree Nyquist response, the
\(\widetilde\omega^{-1/2}\) magnitude, and the \(-45^\circ\) phase without
using a reference line as solver input.

The finite-length reader keeps
\(\widehat{\Delta c}(L)=0\) and
\(d\widehat{\Delta c}/dx|_L=0\) separate. Independent checkpoints verified the
fixed-composition endpoint, zero-gradient blocked endpoint, high-frequency
semi-infinite limit, and low-frequency real parts \(1\) and \(1/3\).

The TLM introduction now explains why a lumped circuit loses spatial current
transfer information and maps each transport equation to its rail resistance,
chemical storage element, or terminal boundary condition. The reduced model's
included and excluded physics remain explicit. Its schematic, distributed
controls, frequency-colored Nyquist/Bode views, and chemical-storage/current
profiles now closely follow the original TLM teaching tool, while the notebook
retains three simple ideal contact presets for classroom use.

## Other module refinements

- Module 01 retains exact finite-\(N\) combinatorics, Stirling thermodynamics,
  and the dilute limit; only shared figure styling and label consistency
  changed.
- Module 02 retains the full mass-action and exact electroneutrality solver;
  Brouwer slopes remain post-solution annotations.
- Module 04 retains the nonlinear Gouy--Chapman, Mott--Schottky, GCS, and
  Frumkin equations. Core takeaways now follow their figures, while
  \(\alpha\) remains a kinetic parameter that does not change equilibrium GCS
  profiles or capacitance.
- Module 05 retains its ideal one-dimensional pair model and both drive modes.
  Crowded multi-panel titles were shortened or wrapped.
- Module 06 retains its conservative finite-slab PITT/GITT and OCV model plus
  the collapsed finite-kinetics extension. Numerical implementation detail
  remains outside the core classroom path.

## Functional audit

A read-only AST audit found 99 marimo controls across the seven notebooks and
no disconnected widget variables. Default executed exports completed for every
module. The revised Module 03 and Module 07 controls were additionally exercised
at multiple parameter values and limiting cases.

## Validation results

| Gate | Result |
|---|---|
| Python compilation | All seven notebooks compile. |
| Strict marimo checks | All seven pass marimo check --strict in an isolated marimo 0.24.0 environment. |
| Executed static exports | All seven complete without cell exceptions; 33 figures were rendered and reviewed. |
| Module 03 checks | The \(\Gamma\) convention, \(\langle x^2\rangle=2Dt\), fitted-versus-analytical \(D\), and existing transport checks pass. |
| Module 07 checks | Waveform period, series/parallel RC, finite Warburg profiles, equal-Nyquist-axis, and TLM boundary, passivity, unit, and conservation checks pass. |
| WASM exports | All seven routes export in run mode with the Pages workflow options. |
| Route assets | Every local route contains its index bundle and dynamically imported run-page module. |
| Widget audit | 99 controls found; no disconnected controls. |
| Repository hygiene | Temporary QA helpers and generated files were removed; git diff --check passes. |

## Preserved limitations

- The models remain intentionally idealized and mostly one-dimensional.
- No new species, phase transformations, porous-electrode model, fitting
  framework, or shared physics package was introduced.
- Module 04 still omits specific adsorption and full Marcus kinetics.
- Module 07 still uses ideal capacitors and a uniform reduced TLM rather than a
  unique experimental equivalent-circuit fit.
- The public Pages content updates only after this commit is deployed.
