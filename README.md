# Solid State Ionics Interactive

[![Deploy marimo apps to GitHub Pages](https://github.com/qiyanglu/Solid-State-Ionics-Interactive/actions/workflows/pages.yml/badge.svg)](https://github.com/qiyanglu/Solid-State-Ionics-Interactive/actions/workflows/pages.yml)

**[Browse the interactive course modules &rarr;](https://qiyanglu.github.io/Solid-State-Ionics-Interactive/)**

This repository is a growing series of interactive [marimo](https://marimo.io/)
notebooks for teaching the **Solid State Ionics** course at Westlake University.
The notebooks complement the lectures and course materials collected on the
[Solid State Ionics Lab teaching page](https://ssi-westlake.com/teaching/) and
the growing [tutorial collection](https://ssi-westlake.com/tutorial/).

The series helps students explore how defect formation, defect chemistry, ionic
and electronic transport, electrochemical polarization, space charge, titration
transients, and impedance spectra emerge from governing equations. Each module favors transparent
physics, interactive controls, and numerical checks over preassembled textbook
curves.

This is a permanent conceptual collection assembled from multiple course years,
lecture decks, and independently audited tutorial articles—not a literal
transcription of one year's slides. Symbols and sign conventions are connected
in [NOTATION.md](NOTATION.md).

## Learning pathways

### Foundations (Modules 01–03)

Build equilibrium defects, defect-chemistry regimes, and transport from atomic
hopping to coupled chemical diffusion.

### Interfaces and boundary conditions (Modules 04–05)

See how interfacial charge and blocking electrodes reshape potential,
capacitance, reaction rate, and stoichiometry.

### Electrochemical methods (Modules 06–07)

Connect time-domain titration to frequency-domain impedance, with boundary
conditions and model assumptions kept visible.

| Module | Topic | Interactive app |
|---|---|---|
| 01 | Defect Formation Thermodynamics | [Launch](https://qiyanglu.github.io/Solid-State-Ionics-Interactive/01-defect-formation/) |
| 02 | Brouwer Diagram for Acceptor-Doped SrTiO3 | [Launch](https://qiyanglu.github.io/Solid-State-Ionics-Interactive/02-brouwer-sto/) |
| 03 | Defect Transport: From Atomic Hopping to Chemical Diffusion | [Launch](https://qiyanglu.github.io/Solid-State-Ionics-Interactive/03-defect-transport/) |
| 04 | Space-Charge Layers and the Frumkin Effect | [Launch](https://qiyanglu.github.io/Solid-State-Ionics-Interactive/04-space-charge-frumkin/) |
| 05 | Stoichiometry Polarization in a Mixed Conductor | [Launch](https://qiyanglu.github.io/Solid-State-Ionics-Interactive/05-stoichiometry-polarization/) |
| 06 | From Coulometric Titration to PITT and GITT | [Launch](https://qiyanglu.github.io/Solid-State-Ionics-Interactive/06-pitt-gitt/) |
| 07 | Impedance Spectroscopy, Warburg Diffusion, and Transmission Lines | [Launch](https://qiyanglu.github.io/Solid-State-Ionics-Interactive/07-impedance-tlm/) |

### Module 01: Defect Formation Thermodynamics

[01_defect_formation.py](01_defect_formation.py) follows one neutral defect
species on \(N\) equivalent lattice sites. It connects the exact finite-system
multiplicity

~~~text
Omega = binomial(N, n)
S_config = k_B ln(Omega)
~~~

to the Stirling-limit entropy, free-energy minimum, chemical-potential zero, and
thermodynamic logistic occupancy. For a finite lattice, the discrete minimum of
\(G(n)\) is identified precisely as the **most probable finite-\(N\)
macrostate**, because

~~~text
P(n) is proportional to binomial(N, n) exp[-n Delta g_f^0/(k_B T)].
~~~

For independent equivalent sites, the ensemble mean remains
\(\langle n\rangle/N=x_{\rm eq}\), even when the most probable small-\(N\)
macrostate is \(n=0\). Interactive controls expose temperature, formation
enthalpy, non-configurational formation entropy (for example, vibrational
entropy), and finite lattice size. The notebook keeps exact finite-\(N\)
combinatorics, the Stirling limit, and the dilute Boltzmann approximation
conceptually distinct.

The default state is 1000 K, Delta h_f = 0.45 eV/defect,
Delta s_f^0 = 3 k_B/defect, and N = 200. It gives x_eq approximately 0.098 and
a visibly populated randomized finite lattice, making the curvature around the
free-energy minimum easy to see. This intentionally large fraction is a
teaching choice, not a typical dilute defect concentration in an oxide.

### Module 02: Brouwer Diagram for Acceptor-Doped SrTiO3

[02_brouwer_sto.py](02_brouwer_sto.py) solves, independently at every oxygen
partial pressure,

~~~text
K_red = V n^2 pO2^(1/2)
K_eh  = n p
2 V + p = A + n
~~~

for oxygen vacancies V, electrons n, and holes p, with a fixed fully ionized
acceptor concentration A. The full mass-action and charge-neutrality equations
are satisfied together. Brouwer regimes and slopes are identified only after
the equilibrium curves are calculated; they are never used to construct them.

The default state is 973 K, A = 10^18 cm^-3, and pO2 from 10^-25 to 1 bar.
Controls span 700-1500 K and log10(A/cm^-3) from 13 to 21.

### Module 03: Defect Transport

[03_defect_transport.py](03_defect_transport.py) follows the one-dimensional
sequence used in Lectures 3–5, from a thermally activated jump to macroscopic
composition relaxation:

~~~text
Gamma = nu exp[-Delta H_mig/(k_B T)]
D = a^2 Gamma / 2
<x^2> = 2 D t
J = -D dc/dx
~~~

Here Gamma is the total hop frequency and 1/Gamma is the mean time between
hops. A discrete master equation recovers Fick's law using the same Gamma
throughout. A symmetric electric-field bias then produces
detailed balance, low-field Nernst-Einstein drift, and
mu_tilde = mu + z F phi. An interactive equilibrium example shows nonzero
chemical and electrical contributions cancelling to give
d(mu_tilde)/dx = 0 and zero flux.

The final sections use the lecture notation for tracer diffusivity \(D^*\),
conductivity-derived self-diffusivity \(D^q\), and Li chemical diffusivity
\(D_{\rm Li}^{\delta}\). Correlated motion is introduced through the Haven
ratio \(D^*=H D^q\). A collapsed non-ideal bridge gives the dimensionally
complete conductivity-thermodynamic-susceptibility expression before
specializing to the dilute reaction
\(Li \rightleftharpoons Li^+ + e^-\), local equilibrium and local charge
neutrality require \(J_{\rm Li}=J_{\rm Li^+}=J_{e^-}\), giving

~~~text
D_Li^delta = 2 D_Li+ D_e- / (D_Li+ + D_e-)
t_D        = L^2 / D_Li^delta
tau^delta  = t_D / pi^2
~~~

The internal field is solved from the equal-flux condition. Controls expose the
Li-ion/electron mobility contrast and sample length. The introductory module
stays with the dilute Li derivation developed in the slides.

### Module 04: Space-Charge Layers and the Frumkin Effect

[04_space_charge_frumkin.py](04_space_charge_frumkin.py) uses the planar,
one-dimensional interface convention from the space-charge lectures. The bulk
reference is phi_infinity = 0 and the core is positively charged. Students first
see a flat electrochemical potential generate the Boltzmann concentration
profile, then use Poisson's equation to compare the two lecture limits:

~~~text
Gouy-Chapman:  both +ze and -ze defects are mobile
Mott-Schottky: negative majority dopants are frozen
~~~

The Gouy-Chapman profiles use the exact nonlinear planar solution; the familiar
exponential is shown only as its small-potential approximation. The
Mott-Schottky section keeps the frozen-dopant depletion approximation and its
parabolic potential over the finite width lambda. Both cases display the
chemical and electrical contributions cancelling to keep the mobile-defect
electrochemical potential flat.

The advanced, initially collapsed sections turn the Gouy-Chapman
charge-potential relation into the differential capacitance C_d, introduce a
constant-capacitance Stern layer, and
solve the Gouy-Chapman-Stern voltage split

~~~text
Q_core = C_s (phi_0 - phi_1) = Q_GC(phi_1)
1/C_tot = 1/C_s + 1/C_d.
~~~

The Frumkin section then evaluates both concentration and potential at the
reaction plane x = x_1. A signed reactant charge number z_R makes the Boltzmann
factor unambiguous and shows why the two Frumkin contributions can either
reinforce or oppose one another. The default z_R = +1 represents the
proton-like reactant depleted by a positive core in the lecture example. For a
single transparent voltage coordinate, the teaching plot takes the formal
potential to coincide with the point of zero charge; the notebook states that
this is a reference choice, not a general identity. The other defaults are
T = 800 K, c_i,infinity = 10^18 cm^-3, epsilon_r = 100, z = 1,
phi_0 = 0.16 V, and C_s = 20 microF/cm^2. Marcus theory and specific adsorption
are deliberately left outside this module.

### Module 05: Stoichiometry Polarization

[05_stoichiometry_polarization.py](05_stoichiometry_polarization.py) follows
the lecture's one-dimensional slab from x = 0 to L. It uses an explicit ideal
pair reaction

~~~text
H <-> H+ + e-
c_i = c_e = c(x,t)
~~~

with two electrodes that pass electrons but block ions. The notebook derives
the coupled transport model from the two electrochemical-potential gradients,
giving

~~~text
J_i = t_i j/F - D_delta dc/dx
D_delta = 2 D_i D_e / (D_i + D_e)
J_i(0,t) = J_i(L,t) = 0.
~~~

This explicit species model gives mu = mu_i + mu_e = 2 RT ln(c/c0). The
notebook therefore retains the factor 1/2 in the small-polarization boundary
coefficient instead of importing a coefficient that belongs to a different
chemical-potential model.

Constant current is the core classroom path, with weak, moderate, strong,
and reverse physical presets that report both beta and the corresponding
current density. Constant potential remains available in a collapsed extension.
The conductivity ratio sigma_e/sigma_i controls the carrier bottleneck, while
the total conductivity sets the initial Ohmic scale. Interactive figures show
the full concentration history, voltage/current relaxation, and separate
chemical, electrical, and electrochemical potentials. A measurement panel
separates imposed, measured, immediate Ohmic, evolving chemical-polarization,
and late Nernstian contributions, then bridges total C_chem to the distributed
c_chem used in Module 07.

The defaults are T = 800 K, c0 = 10^20 cm^-3, L = 100 micrometers,
sigma_i + sigma_e = 10^-3 S/cm, sigma_e/sigma_i = 100, and beta = 0.8.
Electron-blocking, reversible, and Hebb-Wagner electrode cases are deliberately
reserved for later extensions with their own boundary conditions.

### Module 06: From Coulometric Titration to PITT and GITT

[06_pitt_gitt.py](06_pitt_gitt.py) begins with coulometric composition steps and
an equilibrium $E(\delta)$ titration curve, then extends the ideal
one-dimensional pair model from Module 05 to complementary selective contacts.
The student-facing article coordinate is

~~~text
ion electrolyte at x = 0:    J_e = 0
current collector at x = L:  J_i = 0
~~~

Students may label the neutral pair as Li or H. PITT fixes a small electrode-
potential step and solves for the current;
GITT fixes a current step and calculates the voltage. In both cases, the full
finite-slab chemical-diffusion equation generates the transient concentration,
chemical-potential, electrical-potential, and electrochemical-potential
profiles. The selective contact fluxes are enforced directly at both faces, so
the steep electrolyte-side concentration response remains smooth and physically
consistent. The final pulse profile becomes the initial condition for the OCV
relaxation, where terminal current is zero but equal internal ion and electron
fluxes can continue.

The notebook then derives the classical one-sided, small-signal PITT and GITT
series and their short- and long-time limits. These formulas are displayed as
asymptotic comparisons, not used to construct the full solutions. Controls
expose temperature, concentration, thickness, chemical diffusivity,
electronic-to-ionic conductivity ratio, pulse size, pulse duration, and rest
duration. The defaults are T = 800 K, c0 = 10^20 cm^-3, L = 100 micrometers,
D_delta = 10^-8 cm^2/s, and sigma_e/sigma_i = 100.

A collapsed advanced reader replaces the instantaneous PITT boundary with a
finite surface-exchange condition. It compares Bi = infinity, 100, 1, and 0.01,
shows current and concentration profiles, and quantifies the bias produced when
finite-kinetics data are fitted with the ideal diffusion-only long-time slope.
The same section separates the GITT measurement into equilibrium surface
potential, charge-transfer overpotential, and Ohmic drop.

### Module 07: Impedance Spectroscopy, Warburg Diffusion, and Transmission Lines

[07_impedance_tlm.py](07_impedance_tlm.py) begins with the lecture convention

~~~text
exp(i omega t),   Z = V_hat / I_hat = Z' + i Z'',
Nyquist axes: Z' versus -Z''.
~~~

Students first connect voltage/current phase to complex impedance and explore
how one or two ideal parallel-RC relaxations produce separated or overlapping
arcs. The diffusion section then solves the one-dimensional frequency-domain
chemical-diffusion equation. Semi-infinite, fixed-composition finite-length,
and zero-flux finite-length Warburg responses are kept distinct, and the
controls show how L and D_delta set the laboratory diffusion frequency.
A direct DC/AC comparison connects $\sqrt{D^\delta t}$ with
$\sqrt{D^\delta/\omega}$, and every finite-length name is paired with its
boundary equation. General and dilute resistance scales are stated separately.

The final section incorporates the continuous dual-rail model from the
[TLM teaching tool](https://qiyanglu.github.io/TLM-teaching-tool/):

~~~text
du_e/dx = -r_e I_e                    dI_e/dx = -i omega c_chem (u_e - u_i)
du_i/dx = -r_i I_i                    dI_i/dx = +i omega c_chem (u_e - u_i)
~~~

Here u_e and u_i are voltage-equivalent electrochemical potentials. Three
transparent ideal contact cases reveal how the same interior MIEC can look
conducting, chemically polarized, or blocking. Interactive internal profiles
show current transferring between rails while I_e + I_i remains conserved.
The general anatomy labels all four contact impedances and distinguishes the
interactive model from conceptual extensions involving dielectric capacitance,
surface reaction resistance, and more general contacts.

## Run locally

With Python 3.11 or newer and [uv](https://docs.astral.sh/uv/) installed:

~~~console
uv sync
uv run marimo edit 01_defect_formation.py
uv run marimo edit 02_brouwer_sto.py
uv run marimo edit 03_defect_transport.py
uv run marimo edit 04_space_charge_frumkin.py
uv run marimo edit 05_stoichiometry_polarization.py
uv run marimo edit 06_pitt_gitt.py
uv run marimo edit 07_impedance_tlm.py
~~~

To present a notebook as an app, replace edit with run.

## Validation

Generated validation exports belong in the ignored `dist/` directory; they are
not course source files.

~~~console
uv run marimo check --strict 01_defect_formation.py 02_brouwer_sto.py 03_defect_transport.py 04_space_charge_frumkin.py 05_stoichiometry_polarization.py 06_pitt_gitt.py 07_impedance_tlm.py
uv run marimo export html 01_defect_formation.py -o dist/defect-formation.html --no-include-code -f
uv run marimo export html 02_brouwer_sto.py -o dist/brouwer.html --no-include-code -f
uv run marimo export html 03_defect_transport.py -o dist/defect-transport.html --no-include-code -f
uv run marimo export html 04_space_charge_frumkin.py -o dist/space-charge-frumkin.html --no-include-code -f
uv run marimo export html 05_stoichiometry_polarization.py -o dist/stoichiometry-polarization.html --no-include-code -f
uv run marimo export html 06_pitt_gitt.py -o dist/pitt-gitt.html --no-include-code -f
uv run marimo export html 07_impedance_tlm.py -o dist/impedance-tlm.html --no-include-code -f
~~~

Module 01 checks its thermodynamic free-energy minimum, chemical-potential zero,
finite-\(N\) macrostate spacing, large-\(N\) Stirling convergence, and dilute
limit. Module 02 checks mass action, electroneutrality, positivity, regime
coverage, and limiting slopes. Student-facing tables report the physical result
and why it matters, while detailed numerical tolerances remain internal.

Module 03 checks the lecture identity \(D=a^2\Gamma/2\), the stochastic MSD
fit, detailed balance, low-field Nernst-Einstein drift, one-dimensional
master-equation conservation, Fick flux, electrochemical cancellation,
equal Li-ion/electron flux, zero current, and the agreement between the
conductivity and diffusivity forms of \(D_{\rm Li}^{\delta}\).

Module 04 checks the Boltzmann and electrochemical-potential identities, the
exact Gouy-Chapman profile and Gauss law, the Mott-Schottky boundary conditions
and Poisson curvature, GCS charge and voltage matching, series capacitance,
Frumkin reaction-plane consistency, positivity, and finiteness. Each displayed
check includes a short explanation of the physical link it protects.

Module 05 checks the uniform initial state, total-ion conservation, positivity,
constant-current or constant-potential control, the initial total-conductivity
response, zero ionic boundary flux, the ambipolar diffusivity identity, both
electrochemical-potential decompositions, the measured voltage, and the
late-time linear profile and flat ionic electrochemical potential.


Module 06 checks the uniform initial state, positivity, fixed-voltage PITT and
fixed-current GITT control, zero terminal current and conserved mean composition
during OCV, selective-contact fluxes and pulse mass balance, the Module 05 chemical-diffusivity identity,
potential reconstruction of the measured voltage, the PITT/GITT limiting
series, long-rest decay, and spatial-grid convergence.
It additionally checks the finite-kinetics reaction- and diffusion-controlled
eigenvalue limits, positive concentrations, and the direction of the
diffusion-only fitting bias across Biot number.


Module 07 checks the phasor sign convention, the ideal parallel-RC semicircle
and apex, DC/AC diffusion-length scaling, the general-to-dilute Warburg
resistance reduction, finite-length Warburg limits and passivity, TLM boundary
residuals, voltage-equivalent potential signs, distributed/total conversions,
total-current conservation, the reversible-contact limit, passivity, and
finiteness.


## Repository layout

The seven numbered Python files are the self-contained marimo modules. The
`pages/` directory is also source: `pages/index.html` is the course landing page.
The Pages workflow copies it before exporting the seven notebooks.
`NOTATION.md` is the cross-module symbol contract; `AUDIT_PLAN.md` and
`IMPLEMENTATION_REPORT.md` record the current comprehensive audit.


Local `.venv/`, `__marimo__/`, `__pycache__/`, and `dist/` directories are
generated working files. They are ignored by Git and may be removed whenever a
clean checkout is desired.

## Browser deployment

All seven notebooks use only marimo, NumPy, SciPy, and matplotlib and perform no
network or filesystem access at runtime. The workflow in
[pages.yml](.github/workflows/pages.yml) exports each notebook as a
browser-hosted WASM app and deploys the generated static site to GitHub Pages on
every push to main. The project root is a stable module index; each notebook is
published under its own numbered path. After deployment, the workflow retries
all live module routes and verifies their browser entry bundles plus the
dynamically imported marimo run-page modules.
