# Notation and sign conventions

This page is the bridge between the eight notebooks. Each notebook remains
self-contained and introduces the symbols it needs; this page explains how the
same physical ideas are named across modules.

## Defects and concentrations

- A defect is first written in full Kröger–Vink notation and may then be given
  a short label. In Module 02,
  \(V=[V_{\mathrm O}^{\bullet\bullet}]\),
  \(n=[e']\), \(p=[h^\bullet]\), and
  \(A=[A_{\mathrm{Ti}}']\).
- \(\delta\) is a dimensionless stoichiometry or defect fraction.
- \(c_N\) is a number concentration (for example, cm\(^{-3}\)); \(c\) is a
  molar concentration (mol m\(^{-3}\)) in transport and electrochemical
  formulas. When a host molar volume \(V_m\) is used, the teaching convention is
  \(c=\delta/V_m\).
- \(c_0\) is an initial or bulk reference concentration. A subscript \(s\)
  denotes a surface value.
- Chemical diffusion changes an electrically neutral composition species, such
  as Li, H, or O. It does not describe an isolated charged ion moving alone.

## Flux and current

- \(\Phi_{j+1/2}\): net particle exchange rate across a discrete bond in
  Module 03, s\(^{-1}\). It is not the continuum flux \(J_N\).
- \(J_N\): number flux, particles m\(^{-2}\) s\(^{-1}\).
- \(J\): molar flux, mol m\(^{-2}\) s\(^{-1}\).
- \(I\): total current, A. \(j=I/S\): conventional current density, A m\(^{-2}\).
- \(S\): specimen or electrode area. The letter \(A\) is reserved for the
  acceptor concentration in Module 02.
- For monovalent positive ions and electrons, positive particle/molar flux is
  in the +\(x\) direction and conventional charge current is

  \[
  j_i=FJ_i,\qquad j_e=-FJ_e,\qquad j=F(J_i-J_e).
  \]

  A notebook states separately whether positive terminal current means
  insertion or extraction.

## Potentials and electric field

- \(\mu_k\): chemical potential; \(\widetilde\mu_k\): electrochemical potential;
  \(\phi\): electrostatic potential.
- The electric field is \(\mathcal E=-\partial\phi/\partial x\). The letter
  \(E\) is not used for electric field in the electrochemical-method modules.
- In Module 06, \(E\) is an electrode potential relative to a reference
  electrode. In Module 05, \(U=V_{\mathrm{left}}-V_{\mathrm{right}}\) and
  \(U=[\widetilde\mu_e(L)-\widetilde\mu_e(0)]/F\). In Module 08, \(V(t)\)
  and \(\widehat V\) are the time-domain and phasor voltages.
- For monovalent carriers, the voltage-equivalent electrochemical potentials in
  the transmission-line model are

  \[
  u_e=-\widetilde\mu_e/F,\qquad
  u_i=+\widetilde\mu_i/F,\qquad
  u_e-u_i=-\mu_{\mathrm{neutral}}/F.
  \]

## Transport coefficients and time scales

- \(D^*\): tracer diffusivity measured by following labeled particles.
- \(D^q\): conductivity-derived self-diffusivity. Correlated hopping gives
  \(D^*=H D^q\), with Haven ratio \(H\).
- \(D^\delta\): chemical diffusivity of the neutral composition variable.
- \(k^\delta\): linearized surface exchange coefficient for that composition,
  in m s\(^{-1}\).
- \(\Gamma_{\mathrm{hop}}\): total one-dimensional hop frequency in Module 03.
- \(f_{\mathrm{th}}=[c/(RT)](\partial\mu/\partial c)\): dimensionless
  thermodynamic factor used in Module 07. It is introduced only after the
  chemical-potential derivative that defines chemical capacitance.
- \(t_i=\sigma_i/(\sigma_i+\sigma_e)\) and
  \(t_e=\sigma_e/(\sigma_i+\sigma_e)\) are ionic and electronic transference
  numbers.

Two diffusion times appear for complementary purposes:

\[
t_D=\frac{L^2}{D^\delta},\qquad
\tau^\delta=\frac{L^2}{\pi^2D^\delta}.
\]

The first is the direct scaling time; the second is the slowest-mode time for
the common finite-slab relaxation. Correspondingly,

\[
\theta=\frac{D^\delta t}{L^2},\qquad
s=\frac{t}{\tau^\delta}=\pi^2\theta.
\]

For finite surface kinetics,
\(\tau_d=L^2/D^\delta\), \(\tau_{ct}=L/k^\delta\), and
\(\mathrm{Bi}=k^\delta L/D^\delta=\tau_d/\tau_{ct}\).

## Chemical capacitance

Module 07 defines total chemical capacitance from the local chemical-potential
slope:

\[
C_{\mathrm{chem}}=z^2F^2V_{\mathrm{act}}
\left(\frac{\partial c}{\partial\mu}\right)_T .
\]

Its volumetric form is
\(c_{\mathrm{chem}}^V=C_{\mathrm{chem}}/V_{\mathrm{act}}\), in
F m\(^{-3}\). For a uniform slab,
\(c_{\mathrm{chem}}=c_{\mathrm{chem}}^V S=C_{\mathrm{chem}}/L\), in
F m\(^{-1}\), is the distributed capacitance used by the Module 08
transmission line. The chemical diffusivity bridge is

\[
D^\delta=\frac{\sigma_{\mathrm{amb}}}{c_{\mathrm{chem}}^V}.
\]

## Impedance convention

Module 08 uses \(e^{\mathrm{i}\omega t}\), where \(\omega=2\pi f\), and

\[
Z=\widehat V/\widehat I=Z'+\mathrm{i}Z''.
\]

A passive capacitor has \(Z''<0\). Nyquist plots therefore show \(Z'\) on the
horizontal axis and \(-Z''\) on the vertical axis.
