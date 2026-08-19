# Seven-Module Audit Implementation Report

## Delivery summary

The seven existing module filenames, numbering, conceptual order, and public
routes were preserved. Each notebook remains self-contained and limited to
marimo, NumPy, SciPy, and matplotlib. No shared runtime physics package or new
notebook was introduced.

The four supplied Chinese articles were treated as scientific and notational
references rather than instructions. Their source hierarchy and the staged
implementation plan are recorded in [AUDIT_PLAN.md](AUDIT_PLAN.md), while the
cross-module student notation is recorded in [NOTATION.md](NOTATION.md).

## Files changed

| File | Main purpose of the revision |
|---|---|
| `01_defect_formation.py` | Added the prediction/learning/scope rhythm, per-defect notation and reaction-reference explanation, non-color cues, concise figure takeaways, and collapsed checks. |
| `02_brouwer_sto.py` | Added Kröger-Vink/shorthand and equilibrium-constant mappings, activity assumptions, oxygen chemical-potential bridge, regime logic, prediction prompt, and collapsed checks without changing the exact solver. |
| `03_defect_transport.py` | Clarified the embedded random walk, flux/current bridge, Haven ratio, periodic two-front geometry, general chemical-diffusion bridge, and the distinction between diffusion and first-mode clocks. |
| `04_space_charge_frumkin.py` | Reorganized the core reader around electrochemical equilibrium, Gouy-Chapman, and Mott-Schottky; placed GCS and Frumkin material in advanced collapsed readers; clarified reaction-plane signs. |
| `05_stoichiometry_polarization.py` | Made current control the core path, added physical presets and current-density mapping, measurement decomposition, transference-number language, and the total/distributed chemical-capacitance bridge. |
| `06_pitt_gitt.py` | Reframed the module from coulometric titration through PITT/GITT and OCV relaxation; aligned the displayed geometry; added Li/H labels, permanent time normalizations, finite surface kinetics, Biot regimes, and fitting-bias demonstrations. |
| `07_impedance_tlm.py` | Added the DC/AC diffusion bridge, explicit finite-length Warburg boundaries and resistance scales, frequency direction, DRT caution, general TLM anatomy, voltage-equivalent potentials, contact mapping, applications, and scope limits. |
| `README.md` | Organized the collection into three conceptual pathways and updated each module description, validation summary, source framing, and notation links. |
| `pages/index.html` | Updated course navigation, pathways, module summaries, tutorial/teaching links, and notation-guide access without changing any route. |
| `AGENTS.md` | Extended the scientific invariants and verification gates to all seven modules and made `NOTATION.md` the student-facing notation contract. |
| `AUDIT_PLAN.md` | Recorded the pre-edit scope, additions, sign/normalization risks, implementation stages, and verification gates. |
| `NOTATION.md` | Defined concentration, stoichiometry, area, flux/current, potential, diffusion, kinetic, time-scale, and EIS conventions and their context-dependent mappings. |

The existing Pages workflow already exported all seven required routes with the
correct filenames and deployment settings, so it did not need a behavioral
change.

## Notation decisions

- Module 02 retains full Kröger-Vink species in explanations and defines
  (V,n,p,A) only as shorthand. (A) remains reserved for acceptor
  concentration; electrochemical modules use (S) for area.
- Particle-level equations use (k_B,e); molar equations use (R,F). Number
  concentration (c_N), molar concentration (c), and dimensionless
  stoichiometry (delta) are distinguished and mapped where they meet.
- (J_N) is number flux, (J) is molar species flux, (j) is conventional
  current density, and (I) is total current. For the monovalent pair,
  (j=F(J_i-J_e)), with positive directions stated locally.
- (phi) is local electrostatic potential, (mathcal E=-dphi/dx) is electric
  field, (E) is an electrode potential versus reference in PITT/GITT, and
  (U) is a two-terminal voltage in Module 05.
- (D^*), (D^q), and (D^delta) denote tracer, conductivity-derived, and
  chemical diffusivities; (Gamma_{\rm hop}) remains a hopping frequency and
  is not reused as a thermodynamic factor.
- (t_D=L^2/D^delta),
  (	au^delta=L^2/(\pi^2D^delta)),
  (	heta=D^delta t/L^2), and
  (s=t/	au^delta=\pi^2\theta) are kept distinct.
- Module 07 uses (e^{\mathrm{i}\omega t}),
  (Z=Z'+\mathrm{i}Z''), and Nyquist axes (Z') versus (-Z'').

## Physics retained and clarified

- Module 01 still separates exact finite-(N) binomial combinatorics, the
  Stirling limit, and the dilute limit. Its most-probable macrostate remains
  distinct from the ensemble mean.
- Module 02 still calculates every concentration from the documented SrTiO3
  mass-action equations and exact electroneutrality. Brouwer balances and
  slopes remain post-solution interpretations, never inputs.
- Module 03 retains the microscopic hopping, one-dimensional master-equation,
  Nernst-Planck, and ideal ambipolar calculations. The general non-ideal bridge
  is additional context, not a replacement for the lecture's dilute result.
- Module 04 retains the exact planar nonlinear Gouy-Chapman, Mott-Schottky
  depletion, GCS, and Frumkin calculations. The low-potential exponential is
  identified explicitly as the Debye-Hückel limit.
- Module 05 retains the ideal monovalent neutral-pair model, local
  electroneutrality, two ion-blocking electrodes, and both drive modes. The
  constant-potential case is still available as an extension.
- Module 06 retains its finite-slab PITT/GITT and OCV-relaxation solutions. The
  finite-kinetics reader is a separate verified extension rather than a change
  to the introductory diffusion-controlled limit.
- Module 07 retains the RC, three Warburg boundary models, and continuous
  dual-rail TLM solution. The general four-terminal anatomy is presented as a
  map around the deliberately reduced interactive model.

## New PITT/GITT content

The visible sequence now begins with what potentiostatic, galvanostatic,
intermittent, and titration mean experimentally. Charge integration connects
(Q=\int I\,dt) to a small stoichiometry change, an illustrative (E(\delta))
curve, and the thermodynamic slope used in data analysis. The displayed
coordinate follows the tutorial geometry: electrolyte/MIEC at (x=0) and
metal current collector/MIEC at (x=L).

PITT is identified as a prescribed surface chemical-potential/composition
boundary, while GITT is a prescribed surface flux/current boundary. The pulse
and OCV steps show concentration, chemical, electrical, and electrochemical
potential profiles for the neutral Li or H pair. Classical short- and
long-time expressions remain comparisons to the full finite-slab result.

The collapsed finite-kinetics extension introduces
(mathrm{Bi}=k^delta L/D^delta), its reaction-, mixed-, and
diffusion-controlled regimes, and representative values
(infty,100,1,0.01). A verified Robin-boundary eigenfunction model generates
the synthetic profiles and currents. Diffusion-only Cottrell and long-time fits
are then applied to those full finite-kinetics transients to show the inferred
(D^delta/D^delta_{\rm true}) bias. The GITT reader also separates
equilibrium surface potential, charge-transfer overpotential, and Ohmic drop.

## New Warburg and TLM content

Module 07 now links (sqrt{D^delta t}) in a DC step to
(sqrt{D^delta/\omega}) in a small-signal experiment. The semi-infinite,
open-boundary finite-length, and blocked-boundary finite-length Warburg models
are paired with their actual far-face equations. General
(R_D=L|\partial E/\partial c|/(zFS D^delta)) and ideal-dilute resistance
scales are shown separately, with (L), (D^delta), and (S) scaling made
visible. Frequency direction is annotated, and the reader explicitly states
that a 45-degree segment is an asymptote rather than the definition of a
Warburg response.

Before the reduced TLM solver, a schematic identifies electronic and ionic
rails, distributed (r_{\rm eon}), (r_{\rm ion}), and (c_{\rm chem}), four
generic terminal impedances, and the optional dielectric capacitance not
implemented by the solver. The notebook derives
(u_e=-\widetilde\mu_e/F), (u_i=+\widetilde\mu_i/F), maps all three contact
presets, distinguishes distributed and total quantities, and labels broader
applications as conceptual reductions rather than computed features.

## Validation results

| Gate | Result |
|---|---|
| Python compilation | All seven notebooks compiled successfully. |
| Strict marimo checks | All seven passed `marimo check --strict`. |
| Executed static HTML exports | All seven completed without cell exceptions or failed physical-check output. |
| WASM exports | All seven completed using the same route names and run-mode options as the Pages workflow. |
| Module 06 additions | Mass balance, positivity, Biot limiting behavior, eigenvalue behavior, and fitting-bias direction passed. |
| Module 07 additions | DC/AC diffusion-length matching, general/dilute resistance mapping, voltage-equivalent potential signs, and distributed/total conversions passed. |
| Route assets | Every local exported route returned its index and dynamically imported run-page asset successfully. |
| Live route smoke test | All seven existing GitHub Pages routes and their dynamic run-page assets returned HTTP 200 at audit time. |
| Figure review | Thirty-one rendered figures were inspected in module contact sheets at approximately 1100-pixel width and projector-scale typography; legends, labels, annotations, line styles, and boundary profiles were checked. |
| Repository hygiene | Visible prose searches and `git diff --check` were included in the final audit. |

## Browser export results

The local WASM routes verified were:

- `/01-defect-formation/`
- `/02-brouwer-sto/`
- `/03-defect-transport/`
- `/04-space-charge-frumkin/`
- `/05-stoichiometry-polarization/`
- `/06-pitt-gitt/`
- `/07-impedance-tlm/`

No redirects or route changes were added. The existing dynamically imported
marimo run-page bundle was fetched successfully for every local and live route.

## Remaining limitations

- The source DOCX files were read from extracted document text because the
  local visual DOCX renderer was unavailable. Equations and notation were
  independently checked against the implemented models rather than copied
  mechanically.
- Modules 03, 05, 06, and 07 deliberately use one-dimensional idealized models.
  Composition-dependent transport, stress coupling, multi-species reactions,
  and nonlinear large-signal effects remain outside their stated scope.
- Module 06's finite-kinetics extension is linearized with constant
  (k^delta). It does not claim to be a full concentration-dependent
  Butler-Volmer or phase-transforming electrode model.
- Module 07's interactive TLM implements a continuous two-rail interior with
  three ideal contact presets. It does not implement arbitrary four-terminal
  impedances, general reaction/interface RC elements, dielectric/stray
  capacitance, or a DRT inversion.
- Public routes currently verify deployment integrity. The revised notebook
  content will appear on the public site only after these repository changes
  are committed and deployed.
