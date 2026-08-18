# Solid State Ionics Interactive

[![Deploy marimo apps to GitHub Pages](https://github.com/qiyanglu/Solid-State-Ionics-Interactive/actions/workflows/pages.yml/badge.svg)](https://github.com/qiyanglu/Solid-State-Ionics-Interactive/actions/workflows/pages.yml)

**[Browse the interactive course modules &rarr;](https://qiyanglu.github.io/Solid-State-Ionics-Interactive/)**

This repository is a growing series of interactive [marimo](https://marimo.io/)
notebooks for teaching the **Solid State Ionics** course at Westlake University.
The notebooks complement the lectures and course materials collected on the
[Solid State Ionics Lab teaching page](https://ssi-westlake.com/teaching/).

The series helps students explore how defect formation, defect chemistry, ionic
and electronic transport, electrochemical polarization, space charge, and
impedance emerge from governing equations. Each module favors transparent
physics, interactive controls, and numerical checks over preassembled textbook
curves.

## Notebooks

| Module | Topic | Interactive app |
|---|---|---|
| 01 | Defect Formation Thermodynamics | [Launch](https://qiyanglu.github.io/Solid-State-Ionics-Interactive/01-defect-formation/) |
| 02 | Brouwer Diagram for Acceptor-Doped SrTiO3 | [Launch](https://qiyanglu.github.io/Solid-State-Ionics-Interactive/02-brouwer-sto/) |
| 03 | Defect Transport: From Atomic Hopping to Chemical Diffusion | [Launch](https://qiyanglu.github.io/Solid-State-Ionics-Interactive/03-defect-transport/) |

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
acceptor concentration A. The unique positive solution is found with a
bracketed log-space root solve. Brouwer regimes and slopes are measured only
after solving the full equations; they are never used to construct the curves.

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
\(D_{\rm Li}^{\delta}\). For the dilute reaction
\(Li \rightleftharpoons Li^+ + e^-\), local equilibrium and local charge
neutrality require \(J_{\rm Li}=J_{\rm Li^+}=J_{e^-}\), giving

~~~text
D_Li^delta = 2 D_Li+ D_e- / (D_Li+ + D_e-)
tau^delta  ~ L^2 / D_Li^delta
~~~

The internal field is solved from the equal-flux condition. Controls expose the
Li-ion/electron mobility contrast and sample length. The introductory module
stays with the dilute Li derivation developed in the slides.

## Run locally

With Python 3.11 or newer and [uv](https://docs.astral.sh/uv/) installed:

~~~console
uv sync
uv run marimo edit 01_defect_formation.py
uv run marimo edit 02_brouwer_sto.py
uv run marimo edit 03_defect_transport.py
~~~

To present a notebook as an app, replace edit with run.

## Validation

~~~console
uv run marimo check --strict 01_defect_formation.py 02_brouwer_sto.py 03_defect_transport.py
uv run marimo export html 01_defect_formation.py -o defect-formation.html --no-include-code -f
uv run marimo export html 02_brouwer_sto.py -o brouwer.html --no-include-code -f
uv run marimo export html 03_defect_transport.py -o defect-transport.html --no-include-code -f
~~~

Module 01 checks its thermodynamic free-energy minimum, chemical-potential zero,
finite-\(N\) macrostate spacing, large-\(N\) Stirling convergence, and dilute
limit. Module 02 checks mass action, electroneutrality, positivity, regime
coverage, and limiting slopes. Its mass-action rows report logarithmic residuals
with an exact target of zero (equivalently, an unlogged equilibrium ratio of
one); values below 1e-12 are displayed consistently as numerical zero.

Module 03 checks the lecture identity \(D=a^2\Gamma/2\), the stochastic MSD
fit, detailed balance, low-field Nernst-Einstein drift, one-dimensional
master-equation conservation, Fick flux, electrochemical cancellation,
equal Li-ion/electron flux, zero current, and the agreement between the
conductivity and diffusivity forms of \(D_{\rm Li}^{\delta}\).

## Browser deployment

All three notebooks use only marimo, NumPy, SciPy, and matplotlib and perform no
network or filesystem access at runtime. The workflow in
[pages.yml](.github/workflows/pages.yml) exports each notebook as a
browser-hosted WASM app and deploys the generated static site to GitHub Pages on
every push to main. The project root is a stable module index; each notebook is
published under its own numbered path. After deployment, the workflow retries
all live module routes and verifies their browser entry bundles plus the
dynamically imported marimo run-page modules.

The former /01-brouwer-sto/ URL redirects to /02-brouwer-sto/.
