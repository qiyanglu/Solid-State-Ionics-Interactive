# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "marimo>=0.23.14",
#     "matplotlib>=3.8",
#     "numpy>=1.26",
#     "scipy>=1.12",
# ]
# ///
import marimo

__generated_with = "0.24.0"
app = marimo.App(width="full")


@app.cell
def _():
    import marimo as mo
    import matplotlib.pyplot as plt
    import numpy as np
    from scipy.integrate import cumulative_trapezoid, solve_ivp

    plt.rcParams.update(
        {
            "font.size": 13,
            "axes.titlesize": 15,
            "axes.labelsize": 13,
            "xtick.labelsize": 11,
            "ytick.labelsize": 11,
            "legend.fontsize": 10.5,
            "axes.facecolor": "#FCFCFA",
            "figure.facecolor": "white",
            "grid.color": "#C7CCD1",
            "grid.alpha": 0.28,
            "axes.titlepad": 10,
            "axes.labelpad": 6,
            "legend.frameon": False,
            "axes.prop_cycle": plt.cycler(
                color=[
                    "#4C7C86",
                    "#B8734A",
                    "#7C6A91",
                    "#6B86A5",
                    "#B77A82",
                    "#5F8A6B",
                    "#C49345",
                ]
            ),
            "lines.solid_capstyle": "round",
            "figure.dpi": 115,
        }
    )

    return cumulative_trapezoid, mo, np, plt, solve_ivp


@app.cell
def _(mo):
    mo.Html(r"""
    <style>
      .markdown.prose { font-size: 1.12rem !important; line-height: 1.70 !important; }
      .markdown.prose table { font-size: 1.02rem !important; }
      .marimo-cell-output label,
      .marimo-cell-output button,
      .marimo-cell-output input,
      .marimo-cell-output select { font-size: 1rem !important; }
    </style>
    """)
    return


@app.cell
def _(mo):
    mo.md(r"""
    # From coulometric titration to PITT and GITT

    **How can a short electrical pulse reveal both equilibrium composition and
    chemical-diffusion speed?**

    A titration step has two jobs. Integrated charge tells how much neutral
    material entered or left the MIEC. After a long open-circuit rest, the
    electrode potential \(E\) gives one equilibrium point. The transient
    between those states reveals how quickly composition redistributes.

    - **PITT** applies a potential step and measures the current decay.
    - **GITT** applies a current step and measures the potential response.
    - **OCV rest** sets terminal current to zero while an inherited internal
      concentration gradient may continue to relax.

    The core reader follows one experiment at a time. All profiles come from
    the finite one-dimensional diffusion model; the familiar short- and
    long-time formulas appear only in advanced analysis.
    """)
    return

@app.cell
def _(mo):
    species_label_06 = mo.ui.dropdown(
        options={
            "Lithium host: Li ⇌ Li⁺ + e⁻": "Li",
            "Hydrogen host: H ⇌ H⁺ + e⁻": "H",
        },
        value="Lithium host: Li ⇌ Li⁺ + e⁻",
        label="neutral composition label",
    )
    species_label_06
    return (species_label_06,)


@app.cell
def _(mo, np, plt, species_label_06):
    _delta = np.linspace(0.03, 0.97, 400)
    _reduced_potential = -np.log(_delta / (1.0 - _delta))
    _selected = np.array([0.18, 0.36, 0.58, 0.78])
    _selected_potential = -np.log(_selected / (1.0 - _selected))
    _figure, _axis = plt.subplots(figsize=(11.8, 4.1), constrained_layout=True)
    _axis.plot(_delta, _reduced_potential, color="#4F7881", lw=1.9)
    _axis.scatter(
        _selected, _selected_potential, marker="D", s=65, color="#B8734A",
        edgecolor="white", zorder=4, label="successive long-rest states",
    )
    for _left, _right in zip(_selected[:-1], _selected[1:]):
        _axis.annotate(
            "", xy=(_right, -np.log(_right / (1.0 - _right))),
            xytext=(_left, -np.log(_left / (1.0 - _left))),
            arrowprops={"arrowstyle": "->", "color": "#7C6A91", "lw": 1.8},
        )
    _axis.set(
        xlabel=rf"neutral {species_label_06.value} stoichiometry, $\delta$",
        ylabel=r"Illustrative $F(E-E^0)/(RT)$",
        title="Coulometry locates composition; long-rest potential locates equilibrium",
    )
    _axis.grid(alpha=0.22)
    _axis.legend(loc="best")
    mo.vstack([
        mo.md(rf"""
        ## 1. First build an equilibrium titration curve

        A pulse moves neutral **{species_label_06.value}**. For a monovalent
        insertion reaction,

        $$Q=\int I\,dt=F n_{{\rm host}}(\delta_2-\delta_1),\qquad
        c=\frac{{\delta}}{{V_m}}.$$

        Integrated charge gives the composition step. After a long OCV rest,
        the uniform specimen gives one point on $E(\delta)$. Repeating small
        pulse–rest steps maps equilibrium; each transient contains transport
        information. The curve below is illustrative rather than a material fit.
        """),
        _figure,
        mo.md("Coulometry answers *where equilibrium is*; PITT/GITT asks *how the solid approaches it*."),
    ])
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## 2. The selective-contact experiment

    Write the selectable neutral species as
    $$M \rightleftharpoons M^+ + e^-, \qquad c_i=c_e=c(x,t).$$

    Following the articles, the **ion electrolyte is at $x=0$** and the
    **electronic current collector is at $x=L$**.

    - At $x=0$, the electrolyte passes $M^+$ and blocks electrons: $J_e(0,t)=0$.
    - At $x=L$, the collector passes electrons and blocks $M^+$: $J_i(L,t)=0$.

    The carriers cross opposite faces, while bulk electroneutrality couples them
    into a neutral composition change. We define $j=F(J_i-J_e)$. With the sign
    convention used here, positive drive extracts $M$, so the mean concentration
    falls. Every displayed profile uses the left-to-right coordinate shown above.
    """)
    return


@app.cell
def _(np, solve_ivp):
    GAS_CONSTANT_J_PER_MOL_K = 8.314462618
    FARADAY_C_PER_MOL = 96485.33212
    AVOGADRO_PER_MOL = 6.02214076e23

    def _positive(name, value):
        number = float(value)
        if not np.isfinite(number) or number <= 0.0:
            raise ValueError(f"{name} must be positive and finite")
        return number

    def pitt_gitt_parameters(
        temperature_k,
        reference_concentration_cm3,
        length_m,
        chemical_diffusivity_cm2_per_s,
        electronic_to_ionic_ratio,
    ):
        """Return a self-consistent ideal H+/e- MIEC parameter set."""
        temperature = _positive("temperature_k", temperature_k)
        concentration_cm3 = _positive("reference_concentration_cm3", reference_concentration_cm3)
        length = _positive("length_m", length_m)
        diffusivity_cm2_per_s = _positive(
            "chemical_diffusivity_cm2_per_s", chemical_diffusivity_cm2_per_s
        )
        ratio = _positive("electronic_to_ionic_ratio", electronic_to_ionic_ratio)

        concentration_mol_per_m3 = concentration_cm3 * 1.0e6 / AVOGADRO_PER_MOL
        diffusivity_m2_per_s = diffusivity_cm2_per_s * 1.0e-4
        ionic_fraction = 1.0 / (1.0 + ratio)
        electronic_fraction = ratio / (1.0 + ratio)
        tau_delta_s = length**2 / (np.pi**2 * diffusivity_m2_per_s)

        # The ideal-pair identity used in Module 05 is
        # D_delta = 2 RT sigma_i sigma_e / [F^2 c (sigma_i + sigma_e)].
        conductivity_total_s_per_m = (
            FARADAY_C_PER_MOL**2
            * concentration_mol_per_m3
            * diffusivity_m2_per_s
            / (2.0 * GAS_CONSTANT_J_PER_MOL_K * temperature * ionic_fraction * electronic_fraction)
        )
        return {
            "temperature_k": temperature,
            "concentration_cm3": concentration_cm3,
            "concentration_mol_per_m3": concentration_mol_per_m3,
            "length_m": length,
            "chemical_diffusivity_cm2_per_s": diffusivity_cm2_per_s,
            "chemical_diffusivity_m2_per_s": diffusivity_m2_per_s,
            "conductivity_ratio": ratio,
            "ionic_fraction": ionic_fraction,
            "electronic_fraction": electronic_fraction,
            "conductivity_total_s_per_m": conductivity_total_s_per_m,
            "conductivity_i_s_per_m": ionic_fraction * conductivity_total_s_per_m,
            "conductivity_e_s_per_m": electronic_fraction * conductivity_total_s_per_m,
            "tau_delta_s": tau_delta_s,
        }

    def _transport_operator(positions, params):
        """Conservative diffusion operator with selective-contact face fluxes."""
        xi = np.asarray(positions, dtype=float)
        if xi.ndim != 1 or xi.size < 9:
            raise ValueError("positions must contain at least nine points")
        spacing = np.diff(xi)
        if not np.allclose(spacing, spacing[0], rtol=1.0e-11, atol=1.0e-13):
            raise ValueError("positions must be uniformly spaced")
        dx = float(spacing[0])
        if abs(xi[0]) > 1.0e-13 or abs(xi[-1] - 1.0) > 1.0e-13:
            raise ValueError("positions must span xi = 0 to 1")

        point_count = xi.size
        operator = np.zeros((point_count, point_count), dtype=float)
        scale = 1.0 / (np.pi**2 * dx**2)
        interior = np.arange(1, point_count - 1)
        operator[interior, interior - 1] = scale
        operator[interior, interior] = -2.0 * scale
        operator[interior, interior + 1] = scale
        operator[0, 0] = -2.0 * scale
        operator[0, 1] = 2.0 * scale
        operator[-1, -2] = 2.0 * scale
        operator[-1, -1] = -2.0 * scale

        contact_forcing = np.zeros(point_count, dtype=float)
        contact_forcing[0] = -2.0 * params["ionic_fraction"] / (np.pi**2 * dx)
        contact_forcing[-1] = -2.0 * params["electronic_fraction"] / (np.pi**2 * dx)
        quadrature_weights = np.full(point_count, dx, dtype=float)
        quadrature_weights[[0, -1]] *= 0.5
        return operator, contact_forcing, quadrature_weights

    def _voltage_step_current(profile, target_voltage, positions, params):
        safe_profile = np.maximum(np.asarray(profile, dtype=float), 1.0e-10)
        resistance_integral = np.trapezoid(1.0 / safe_profile, positions)
        t_i = params["ionic_fraction"]
        t_e = params["electronic_fraction"]
        return (
            target_voltage + t_i * np.log(safe_profile[0]) + t_e * np.log(safe_profile[-1])
        ) / (t_i * t_e * resistance_integral)

    def _voltage_from_profile(profile, q_value, positions, params):
        profile = np.asarray(profile, dtype=float)
        if np.min(profile) <= 0.0:
            return np.nan
        t_i = params["ionic_fraction"]
        t_e = params["electronic_fraction"]
        resistance_integral = np.trapezoid(1.0 / profile, positions)
        return (
            -t_i * np.log(profile[0])
            - t_e * np.log(profile[-1])
            + q_value * t_i * t_e * resistance_integral
        )

    def _linear_profile_evolution(initial_profile, reduced_times, operator, forcing):
        """Advance a constant-forcing diffusion problem by matrix modes."""
        times = np.asarray(reduced_times, dtype=float)
        point_count = initial_profile.size
        metric = np.ones(point_count, dtype=float)
        metric[[0, -1]] = 0.5
        sqrt_metric = np.sqrt(metric)
        symmetric_operator = (
            sqrt_metric[:, None] * operator / sqrt_metric[None, :]
        )
        eigenvalues, eigenvectors = np.linalg.eigh(symmetric_operator)
        initial_coefficients = eigenvectors.T @ (sqrt_metric * initial_profile)
        forcing_coefficients = eigenvectors.T @ (sqrt_metric * forcing)
        exponentials = np.exp(np.outer(times, eigenvalues))
        response = np.empty_like(exponentials)
        stationary = np.abs(eigenvalues) < 1.0e-12
        response[:, stationary] = times[:, None]
        response[:, ~stationary] = np.expm1(
            np.outer(times, eigenvalues[~stationary])
        ) / eigenvalues[None, ~stationary]
        coefficients = (
            exponentials * initial_coefficients[None, :]
            + response * forcing_coefficients[None, :]
        )
        return (coefficients @ eigenvectors.T) / sqrt_metric[None, :]

    def simulate_pitt(
        reduced_pulse_times,
        reduced_rest_times,
        reduced_voltage,
        positions,
        params,
    ):
        """Finite-slab PITT pulse followed by an open-circuit relaxation."""
        pulse_times = np.asarray(reduced_pulse_times, dtype=float)
        rest_times = np.asarray(reduced_rest_times, dtype=float)
        xi = np.asarray(positions, dtype=float)
        operator, contact_forcing, quadrature_weights = _transport_operator(xi, params)
        initial_profile = np.ones(xi.size, dtype=float)
        t_i = params["ionic_fraction"]
        t_e = params["electronic_fraction"]

        def current_and_gradient(profile):
            safe_profile = np.maximum(profile, 1.0e-10)
            resistance_integral = np.sum(quadrature_weights / safe_profile)
            numerator = (
                reduced_voltage
                + t_i * np.log(safe_profile[0])
                + t_e * np.log(safe_profile[-1])
            )
            q_value = numerator / (t_i * t_e * resistance_integral)
            numerator_gradient = np.zeros_like(profile)
            numerator_gradient[0] = t_i / safe_profile[0]
            numerator_gradient[-1] = t_e / safe_profile[-1]
            integral_gradient = -quadrature_weights / safe_profile**2
            q_gradient = (
                numerator_gradient / (t_i * t_e * resistance_integral)
                - q_value * integral_gradient / resistance_integral
            )
            return q_value, q_gradient

        def voltage_step_rhs(_reduced_time, profile_now):
            q_now, _ = current_and_gradient(profile_now)
            return operator @ profile_now + contact_forcing * q_now

        def voltage_step_jacobian(_reduced_time, profile_now):
            _, q_gradient = current_and_gradient(profile_now)
            return operator + np.outer(contact_forcing, q_gradient)

        pulse_solution = solve_ivp(
            voltage_step_rhs,
            (float(pulse_times[0]), float(pulse_times[-1])),
            initial_profile,
            t_eval=pulse_times,
            method="BDF",
            jac=voltage_step_jacobian,
            rtol=2.0e-8,
            atol=2.0e-10,
        )
        if not pulse_solution.success:
            raise RuntimeError(str(pulse_solution.message))
        pulse_profiles = pulse_solution.y.T
        pulse_current = np.array(
            [_voltage_step_current(profile, reduced_voltage, xi, params) for profile in pulse_profiles]
        )
        pulse_voltage = np.array(
            [_voltage_from_profile(profile, q, xi, params) for profile, q in zip(pulse_profiles, pulse_current)]
        )

        rest_profiles = _linear_profile_evolution(
            pulse_profiles[-1], rest_times, operator, np.zeros_like(contact_forcing)
        )
        rest_current = np.zeros(rest_times.size)
        rest_voltage = np.array(
            [_voltage_from_profile(profile, 0.0, xi, params) for profile in rest_profiles]
        )
        return {
            "pulse_profiles": pulse_profiles,
            "rest_profiles": rest_profiles,
            "pulse_q": pulse_current,
            "rest_q": rest_current,
            "pulse_u": pulse_voltage,
            "rest_u": rest_voltage,
        }

    def simulate_gitt(
        reduced_pulse_times,
        reduced_rest_times,
        q_step,
        positions,
        params,
    ):
        """Finite-slab GITT pulse followed by an open-circuit relaxation."""
        pulse_times = np.asarray(reduced_pulse_times, dtype=float)
        rest_times = np.asarray(reduced_rest_times, dtype=float)
        xi = np.asarray(positions, dtype=float)
        operator, contact_forcing, _ = _transport_operator(xi, params)
        initial_profile = np.ones(xi.size, dtype=float)
        pulse_profiles = _linear_profile_evolution(
            initial_profile, pulse_times, operator, contact_forcing * q_step
        )
        pulse_q = np.full(pulse_times.size, float(q_step))
        pulse_u = np.array(
            [_voltage_from_profile(profile, q_step, xi, params) for profile in pulse_profiles]
        )

        rest_profiles = _linear_profile_evolution(
            pulse_profiles[-1], rest_times, operator, np.zeros_like(contact_forcing)
        )
        rest_q = np.zeros(rest_times.size)
        rest_u = np.array(
            [_voltage_from_profile(profile, 0.0, xi, params) for profile in rest_profiles]
        )
        return {
            "pulse_profiles": pulse_profiles,
            "rest_profiles": rest_profiles,
            "pulse_q": pulse_q,
            "rest_q": rest_q,
            "pulse_u": pulse_u,
            "rest_u": rest_u,
        }

    def potential_decomposition(profile, q_value, positions, params):
        """Chemical, electrical, and electrochemical profiles in RT units."""
        xi = np.asarray(positions, dtype=float)
        profile = np.asarray(profile, dtype=float)
        t_i = params["ionic_fraction"]
        t_e = params["electronic_fraction"]
        mu_i = np.log(profile)
        mu_e = np.log(profile)
        resistance_integral = np.zeros_like(xi)
        resistance_integral[1:] = np.cumsum(
            0.5 * (1.0 / profile[1:] + 1.0 / profile[:-1]) * np.diff(xi)
        )
        # Integrated electrochemical-potential relation for the selective contacts.
        tilde_mu_i = (
            mu_i[0]
            + 2.0 * t_e * (np.log(profile) - np.log(profile[0]))
            - 2.0 * q_value * t_i * t_e * resistance_integral
        )
        tilde_mu_e = 2.0 * np.log(profile) - tilde_mu_i
        electrical_i = tilde_mu_i - mu_i
        electrical_e = -electrical_i
        reconstructed_u = -0.5 * (tilde_mu_i[-1] + tilde_mu_e[0])
        return {
            "profile": profile,
            "mu_i": mu_i,
            "mu_e": mu_e,
            "electrical_i": electrical_i,
            "electrical_e": electrical_e,
            "tilde_mu_i": tilde_mu_i,
            "tilde_mu_e": tilde_mu_e,
            "reconstructed_u": reconstructed_u,
        }

    def classical_pitt_series(reduced_times, concentration_step, terms=500):
        """Neumann/Dirichlet finite-slab PITT current in the classical limit."""
        s_values = np.asarray(reduced_times, dtype=float)
        half_modes = np.arange(int(terms), dtype=float) + 0.5
        return 2.0 * concentration_step * np.sum(np.exp(-np.outer(s_values, half_modes**2)), axis=1)

    def classical_gitt_series(reduced_times, q_step, terms=500):
        """Small-signal surface voltage during a finite-slab current pulse."""
        s_values = np.asarray(reduced_times, dtype=float)
        modes = np.arange(1, int(terms) + 1, dtype=float)
        theta = s_values / np.pi**2
        transient_sum = np.sum(
            np.exp(-np.outer(s_values, modes**2)) / modes[None, :] ** 2,
            axis=1,
        )
        return q_step * (theta + 1.0 / 3.0 - 2.0 * transient_sum / np.pi**2)

    return (
        AVOGADRO_PER_MOL,
        FARADAY_C_PER_MOL,
        GAS_CONSTANT_J_PER_MOL_K,
        classical_gitt_series,
        classical_pitt_series,
        pitt_gitt_parameters,
        potential_decomposition,
        simulate_gitt,
        simulate_pitt,
    )



@app.cell
def _(np):
    from scipy.optimize import brentq

    def finite_pitt_modes_06(biot_number, mode_count=100):
        biot = float(biot_number)
        if biot <= 0.0 and not np.isinf(biot):
            raise ValueError("Bi must be positive")
        modes = np.arange(int(mode_count), dtype=float)
        if np.isinf(biot):
            roots = (modes + 0.5) * np.pi
        else:
            roots = np.array([
                brentq(
                    lambda root: root * np.tan(root) - biot,
                    mode * np.pi + 1.0e-12,
                    mode * np.pi + 0.5 * np.pi - 1.0e-12,
                ) for mode in modes
            ])
        coefficients = 4.0 * np.sin(roots) / (2.0 * roots + np.sin(2.0 * roots))
        return roots, coefficients

    def finite_pitt_kinetics_06(theta, position, biot_number, mode_count=100):
        theta_values = np.atleast_1d(np.asarray(theta, dtype=float))
        positions = np.atleast_1d(np.asarray(position, dtype=float))
        if np.any(theta_values <= 0.0):
            raise ValueError("theta must be positive")
        if np.any((positions < 0.0) | (positions > 1.0)):
            raise ValueError("position must lie between zero and one")
        roots, coefficients = finite_pitt_modes_06(biot_number, mode_count)
        decays = np.exp(-np.outer(theta_values, roots**2))
        shapes = np.cos(roots[:, None] * (1.0 - positions[None, :]))
        return {
            "theta": theta_values,
            "position": positions,
            "Bi": float(biot_number),
            "roots": roots,
            "profile": (decays * coefficients[None, :]) @ shapes,
            "current": decays @ (coefficients * roots * np.sin(roots)),
            "mean_profile": decays @ (coefficients * np.sin(roots) / roots),
        }

    def finite_pitt_fit_bias_06(biot_values, fitting_windows, mode_count=70):
        biot_grid = np.asarray(biot_values, dtype=float)
        results = {name: np.empty_like(biot_grid) for name in fitting_windows}
        theta_grid = np.geomspace(8.0e-3, 1.6, 420)
        for index, biot in enumerate(biot_grid):
            response = finite_pitt_kinetics_06(
                theta_grid, np.array([0.0, 1.0]), biot, mode_count
            )
            log_current = np.log(response["current"])
            for name, (lower, upper) in fitting_windows.items():
                selected = (theta_grid >= lower) & (theta_grid <= upper)
                slope = np.polyfit(theta_grid[selected], log_current[selected], 1)[0]
                results[name][index] = -4.0 * slope / np.pi**2
        return results

    return finite_pitt_fit_bias_06, finite_pitt_kinetics_06, finite_pitt_modes_06

@app.cell
def _(mo, plt, species_label_06):
    _fig, _ax = plt.subplots(figsize=(12.0, 2.7), constrained_layout=True)
    _ax.set_xlim(-0.12, 1.12)
    _ax.set_ylim(-0.42, 0.48)
    _ax.axvspan(-0.1, 0.0, color="#697386", alpha=0.95)
    _ax.axvspan(0.0, 1.0, color="#DDEBDD", alpha=0.95)
    _ax.axvspan(1.0, 1.1, color="#B8DDE3", alpha=0.95)
    _ax.text(
        -0.05, 0.12, "ion\nelectrolyte", ha="center", va="center", color="white", weight="bold"
    )
    _ax.text(
        0.5,
        0.15,
        rf"MIEC:  ${species_label_06.value} \rightleftharpoons {species_label_06.value}^+ + e^-$",
        ha="center",
        va="center",
        fontsize=15,
        weight="bold",
    )
    _ax.text(
        1.05, 0.12, "current\ncollector", ha="center", va="center", color="#405E66", weight="bold"
    )
    _ax.annotate(
        rf"${species_label_06.value}^+$ passes",
        xy=(-0.02, -0.12),
        xytext=(0.20, -0.12),
        arrowprops={"arrowstyle": "->", "lw": 2},
        ha="center",
    )
    _ax.annotate(
        "$e^-$ passes",
        xy=(1.02, -0.12),
        xytext=(0.80, -0.12),
        arrowprops={"arrowstyle": "->", "lw": 2},
        ha="center",
    )
    _ax.text(0.02, -0.31, "$J_e(0,t)=0$", ha="left", color="#A65E5E", weight="bold")
    _ax.text(0.98, -0.31, "$J_i(L,t)=0$", ha="right", color="#A65E5E", weight="bold")
    _ax.text(0.0, 0.37, "$x=0$", ha="center")
    _ax.text(1.0, 0.37, "$x=L$", ha="center")
    _ax.axis("off")
    mo.vstack(
        [
            _fig,
            mo.md(
                "The contacts are complementary: neither carrier can cross the whole "
                "cell by itself, yet the neutral pair can be inserted or extracted through "
                "the external circuit."
            ),
        ]
    )
    return


@app.cell
def _(mo):
    derivation = mo.md(r"""
    ## 3. From carrier fluxes to one diffusion equation

    Let $c_0=c(x,0)$ be the initially uniform pair concentration, consistent
    with Module 05. With ideal dilute chemical potentials

    $$
    \mu_i=\mu_i^0+RT\ln(c/c_0),\qquad
    \mu_e=\mu_e^0+RT\ln(c/c_0),
    $$

    the neutral-pair chemical potential is
    $\mu_M=\mu_i+\mu_e=\mu_M^0+2RT\ln(c/c_0)$. Eliminating the
    internal electric field gives exactly the same chemical diffusivity as in
    Module 05,

    $$
    D^\delta=\frac{2D_iD_e}{D_i+D_e},\qquad
    t_D=\frac{L^2}{D^\delta},\qquad
    \tau^\delta=\frac{L^2}{\pi^2D^\delta}.
    $$

    We use both clocks: $\theta=D^\delta t/L^2$ and
    $s=t/\tau^\delta=\pi^2\theta$. The first exposes scaling; the second makes
    the slowest finite-slab mode decay as $e^{-s}$.

    The coupled fluxes become

    $$
    J_i=t_i\frac{j}{F}-D^\delta\frac{\partial c}{\partial x},\qquad
    J_e=-t_e\frac{j}{F}-D^\delta\frac{\partial c}{\partial x},
    $$

    where $t_i=\sigma_i/(\sigma_i+\sigma_e)$ and $t_e=1-t_i$.
    Conservation gives the full finite-slab equation

    $$\frac{\partial c}{\partial t}=D^\delta\frac{\partial^2c}{\partial x^2}.$$

    For compact notation, introduce $y=c/c_0$, $\xi=x/L$,
    $s=t/\tau^\delta$, and

    $$\hat j=\frac{jL}{F D^\delta c_0}.$$

    The selective contacts impose

    $$y_\xi(0)=t_i\hat j,\qquad y_\xi(1)=-t_e\hat j,$$

    and the mean composition obeys $d\bar y/ds=-\hat j/\pi^2$. Thus a
    measured current is also a direct composition balance.

    The voltage measured between the selective terminals is the neutral-pair
    electrochemical drive. Relative to the initial state, define

    $$
    \hat E=\frac{F\Delta E}{2RT}
    =-t_i\ln y(0)-t_e\ln y(1)
    +\hat j t_i t_e\int_0^1\frac{d\xi}{y}.
    $$

    The first two terms are concentration polarization at the two selective
    faces; the integral is the internal Ohmic contribution. PITT holds $\hat E$
    fixed and solves this relation for $\hat j(t)$. GITT holds $\hat j$ fixed and
    evaluates $\hat E(t)$. At OCV, $\hat j=0$ while the inherited concentration profile relaxes.
    """)
    mo.accordion({"Model details - selective-contact diffusion equation": derivation})
    return


@app.cell
def _(mo):
    experiment_mode_06 = mo.ui.dropdown(
        options=["PITT", "GITT"],
        value="PITT",
        label="Experiment",
    )
    temperature_06 = mo.ui.slider(
        600, 1200, value=800, step=25, label="Temperature T (K)", show_value=True
    )
    log_concentration_06 = mo.ui.slider(
        17.0, 22.0, value=20.0, step=0.25, label="Initial concentration, log10(c0 / cm^-3)", show_value=True
    )
    length_06 = mo.ui.dropdown(
        options={
            "25 micrometers": 25.0,
            "50 micrometers": 50.0,
            "100 micrometers": 100.0,
            "250 micrometers": 250.0,
        },
        value="100 micrometers",
        label="MIEC thickness, L",
    )
    log_diffusivity_06 = mo.ui.slider(
        -12.0, -6.0, value=-8.0, step=0.25, label="Chemical diffusivity, log10(D-delta / cm^2 s^-1)", show_value=True
    )
    log_ratio_06 = mo.ui.slider(
        -1.0, 4.0, value=2.0, step=0.25, label="Conductivity ratio, log10(sigma_e / sigma_i)", show_value=True
    )
    pitt_voltage_06 = mo.ui.slider(
        5.0, 35.0, value=15.0, step=1.0, label="PITT potential step (mV, extraction)", show_value=True
    )
    gitt_current_06 = mo.ui.slider(
        0.05, 0.60, value=0.30, step=0.025, label="GITT current-step strength", show_value=True
    )
    pulse_duration_06 = mo.ui.slider(
        0.20, 3.0, value=1.20, step=0.10, label="Pulse duration / tau-delta", show_value=True
    )
    rest_duration_06 = mo.ui.slider(
        0.25, 5.0, value=2.50, step=0.25, label="OCV rest duration / tau-delta", show_value=True
    )
    potential_case_06 = mo.ui.dropdown(
        options=[
            "PITT: voltage step",
            "PITT: OCV relaxation",
            "GITT: current step",
            "GITT: OCV relaxation",
        ],
        value="GITT: current step",
        label="Advanced profile case",
    )
    potential_time_06 = mo.ui.slider(
        0, 100, value=65, step=5, label="Stage progress (%)", show_value=True
    )
    return (
        experiment_mode_06,
        gitt_current_06,
        length_06,
        log_concentration_06,
        log_diffusivity_06,
        log_ratio_06,
        pitt_voltage_06,
        potential_case_06,
        potential_time_06,
        pulse_duration_06,
        rest_duration_06,
        temperature_06,
    )


@app.cell
def _(
    experiment_mode_06,
    gitt_current_06,
    length_06,
    log_concentration_06,
    log_diffusivity_06,
    log_ratio_06,
    mo,
    pitt_voltage_06,
    potential_case_06,
    potential_time_06,
    pulse_duration_06,
    rest_duration_06,
    temperature_06,
):
    pulse_size_control = (
        pitt_voltage_06 if experiment_mode_06.value == "PITT" else gitt_current_06
    )
    core_controls_06 = mo.hstack(
        [
            experiment_mode_06,
            pulse_size_control,
            pulse_duration_06,
            log_diffusivity_06,
        ],
        justify="start",
        align="center",
        wrap=True,
        gap=1.2,
    )
    advanced_controls_06 = mo.vstack([
        mo.hstack(
            [temperature_06, log_concentration_06, length_06],
            justify="start", align="center", wrap=True, gap=1.2,
        ),
        mo.hstack(
            [log_ratio_06, rest_duration_06],
            justify="start", align="center", wrap=True, gap=1.2,
        ),
        mo.hstack(
            [potential_case_06, potential_time_06],
            justify="start", align="center", wrap=True, gap=1.2,
        ),
    ])
    mo.vstack([
        mo.md(r"""
        ## 3. Choose one pulse experiment

        Select PITT or GITT, then change the pulse size, duration, or chemical
        diffusivity. A larger \(D^\delta\) shortens the physical response time;
        it does not change the dimensionless boundary condition.
        """),
        core_controls_06,
        mo.accordion({"Explore further - material and profile controls": advanced_controls_06}),
    ])
    return (core_controls_06,)

@app.cell
def _(
    FARADAY_C_PER_MOL,
    GAS_CONSTANT_J_PER_MOL_K,
    gitt_current_06,
    length_06,
    log_concentration_06,
    log_diffusivity_06,
    log_ratio_06,
    np,
    pitt_gitt_parameters,
    pitt_voltage_06,
    pulse_duration_06,
    rest_duration_06,
    simulate_gitt,
    simulate_pitt,
    temperature_06,
):
    parameters_06 = pitt_gitt_parameters(
        temperature_k=temperature_06.value,
        reference_concentration_cm3=10.0**log_concentration_06.value,
        length_m=length_06.value * 1.0e-6,
        chemical_diffusivity_cm2_per_s=10.0**log_diffusivity_06.value,
        electronic_to_ionic_ratio=10.0**log_ratio_06.value,
    )
    positions_06 = np.linspace(0.0, 1.0, 81)
    pulse_times_06 = np.concatenate(([0.0], np.geomspace(1.0e-5, pulse_duration_06.value, 239)))
    rest_times_06 = np.concatenate(([0.0], np.geomspace(1.0e-5, rest_duration_06.value, 179)))
    pitt_reduced_voltage_06 = (
        FARADAY_C_PER_MOL
        * pitt_voltage_06.value
        * 1.0e-3
        / (2.0 * GAS_CONSTANT_J_PER_MOL_K * parameters_06["temperature_k"])
    )
    pitt_result_06 = simulate_pitt(
        pulse_times_06,
        rest_times_06,
        pitt_reduced_voltage_06,
        positions_06,
        parameters_06,
    )
    gitt_result_06 = simulate_gitt(
        pulse_times_06,
        rest_times_06,
        gitt_current_06.value,
        positions_06,
        parameters_06,
    )

    voltage_scale_v_06 = (
        2.0 * GAS_CONSTANT_J_PER_MOL_K * parameters_06["temperature_k"] / FARADAY_C_PER_MOL
    )
    current_scale_a_per_m2_06 = (
        FARADAY_C_PER_MOL
        * parameters_06["chemical_diffusivity_m2_per_s"]
        * parameters_06["concentration_mol_per_m3"]
        / parameters_06["length_m"]
    )
    return (
        current_scale_a_per_m2_06,
        gitt_result_06,
        parameters_06,
        pitt_reduced_voltage_06,
        pitt_result_06,
        positions_06,
        pulse_times_06,
        rest_times_06,
        voltage_scale_v_06,
    )


@app.cell
def _(current_scale_a_per_m2_06, mo, parameters_06):
    _tau = parameters_06["tau_delta_s"]
    _time_text = f"{_tau:.2e} s"
    if 1.0 <= _tau < 3600.0:
        _time_text = f"{_tau:.1f} s"
    elif _tau >= 3600.0:
        _time_text = f"{_tau / 3600.0:.2f} h"
    state_details = mo.md(
        f"""
        ### What the controls imply

        - $t_i={parameters_06["ionic_fraction"]:.4f}$ and
          $t_e={parameters_06["electronic_fraction"]:.4f}$.
        - $D^\\delta={parameters_06["chemical_diffusivity_cm2_per_s"]:.2e}$
          cm² s⁻¹ and $\\tau^\\delta=L^2/(\\pi^2D^\\delta)={_time_text}$.
        - The ideal-pair identity implies
          $\\sigma_i+\\sigma_e={parameters_06["conductivity_total_s_per_m"] / 100.0:.2e}$
          S cm⁻¹ for this $c_0$, $T$, $D^\\delta$, and conductivity ratio.
        - One unit of normalized pulse strength
          $\\hat j=j_0L/(FD^\\delta c_0)$ corresponds to
          ${0.1 * current_scale_a_per_m2_06:.3e}$ mA cm⁻².

        Changing $L$ or $D^\\delta$ stretches the physical time axis through
        $\\tau^\\delta$. Changing $\\sigma_e/\\sigma_i$ changes how the applied
        voltage is divided between concentration polarization and internal
        Ohmic loss.
        """
    )
    mo.accordion({"Model details - selected time and transport scales": state_details})
    return


@app.cell
def _(experiment_mode_06, mo, plt):
    timeline_figure, timeline_axis = plt.subplots(figsize=(11.8, 2.6), dpi=120)
    timeline_axis.set_xlim(0.0, 3.0)
    timeline_axis.set_ylim(-0.15, 1.0)
    timeline_axis.fill_between([0.0, 1.35], 0.15, 0.65, color="#D9E6E7")
    timeline_axis.fill_between([1.35, 3.0], 0.15, 0.65, color="#EEE4D5")
    timeline_axis.axvline(1.35, color="#6D747A", lw=1.3, ls="--")
    timeline_axis.text(
        0.675, 0.40,
        "Potential step" if experiment_mode_06.value == "PITT" else "Current step",
        ha="center", va="center", fontsize=14,
    )
    timeline_axis.text(1.35, 0.82, "Current interruption", ha="center")
    timeline_axis.text(2.175, 0.40, "OCV rest", ha="center", va="center", fontsize=14)
    timeline_axis.annotate(
        "", xy=(2.95, 0.05), xytext=(0.05, 0.05),
        arrowprops={"arrowstyle": "->", "lw": 1.4, "color": "#4C7C86"},
    )
    timeline_axis.text(1.50, -0.04, "Time", ha="center", va="top")
    timeline_axis.axis("off")
    timeline_axis.set_title("Pulse, interruption, and internal relaxation")
    timeline_figure.tight_layout()
    plt.close(timeline_figure)
    mo.vstack([
        mo.md(r"""
        ## 4. Follow the pulse into OCV rest

        Interruption sets terminal current to zero immediately. It does not
        erase the end-of-pulse concentration profile. During OCV, equal ion and
        electron fluxes can continue to flatten that profile without carrying
        electrical current.
        """),
        timeline_figure,
    ])
    return (timeline_figure,)

@app.cell
def _(
    current_scale_a_per_m2_06,
    experiment_mode_06,
    gitt_result_06,
    mo,
    np,
    parameters_06,
    pitt_result_06,
    plt,
    positions_06,
    pulse_times_06,
    rest_times_06,
    voltage_scale_v_06,
):
    selected_result = (
        pitt_result_06 if experiment_mode_06.value == "PITT" else gitt_result_06
    )
    tau_s = parameters_06["tau_delta_s"]
    x_um = positions_06 * parameters_06["length_m"] * 1.0e6
    pulse_seconds = pulse_times_06 * tau_s
    rest_seconds = (pulse_times_06[-1] + rest_times_06) * tau_s

    profile_figure, profile_axis = plt.subplots(figsize=(11.8, 4.2), dpi=120)
    profile_axis.axhline(
        1.0, color="#8D949A", lw=1.1, ls=":", label="Initial state"
    )
    profile_axis.plot(
        x_um, selected_result["pulse_profiles"][-1][::-1],
        color="#4C7C86", lw=1.7, label="At current interruption",
    )
    profile_axis.plot(
        x_um, selected_result["rest_profiles"][-1][::-1],
        color="#B8734A", lw=1.4, ls="--", label="After OCV rest",
    )
    profile_axis.set(
        xlabel=r"Position, $x$ (micrometers)",
        ylabel=r"Composition, $c/c_0$",
        title=f"{experiment_mode_06.value}: the pulse creates a profile",
    )
    profile_axis.grid(alpha=0.22)
    profile_axis.legend(frameon=False)
    profile_figure.tight_layout()
    plt.close(profile_figure)

    response_figure, response_axis = plt.subplots(figsize=(11.8, 4.2), dpi=120)
    if experiment_mode_06.value == "PITT":
        pulse_response = (
            0.1 * current_scale_a_per_m2_06 * selected_result["pulse_q"]
        )
        rest_response = np.zeros_like(rest_seconds)
        response_axis.plot(
            pulse_seconds, pulse_response, color="#4C7C86", lw=1.7,
            label="Measured current",
        )
        response_axis.plot(
            rest_seconds, rest_response, color="#4C7C86", lw=1.4,
        )
        response_axis.set_ylabel(r"Current density, $j$ (mA cm$^{-2}$)")
        response_axis.set_title("PITT current decays, then is interrupted")
        explanation = (
            "The imposed potential is held during the pulse, so the measured "
            "current decays. At interruption the terminal current is exactly zero."
        )
    else:
        pulse_response = 1.0e3 * voltage_scale_v_06 * selected_result["pulse_u"]
        rest_response = 1.0e3 * voltage_scale_v_06 * selected_result["rest_u"]
        response_axis.plot(
            pulse_seconds, pulse_response, color="#4C7C86", lw=1.7,
            label="During current pulse",
        )
        response_axis.plot(
            rest_seconds, rest_response, color="#B8734A", lw=1.5, ls="--",
            label="OCV relaxation",
        )
        response_axis.set_ylabel(r"Electrode-potential change, $\Delta E$ (mV)")
        response_axis.set_title("GITT potential evolves and relaxes")
        response_axis.legend(frameon=False)
        explanation = (
            "The current is fixed during the pulse. After interruption, the "
            "potential relaxes toward the equilibrium value of the new mean composition."
        )
    response_axis.axvline(
        pulse_seconds[-1], color="#6D747A", lw=1.2, ls=":",
    )
    response_axis.set_xlabel("Time (s)")
    response_axis.grid(alpha=0.22)
    response_figure.tight_layout()
    plt.close(response_figure)

    mo.vstack([
        mo.md(f"### {experiment_mode_06.value}: concentration first, response second"),
        profile_figure,
        mo.md(
            "The solid line is the inherited end-of-pulse profile; the dashed "
            "line shows how OCV rest flattens it while preserving its mean."
        ),
        response_figure,
        mo.md(explanation),
    ])
    return profile_figure, response_figure

@app.cell
def _(mo):
    mo.md(r"""
    ### Read the interruption carefully

    The voltage generally jumps when current is switched off because the Ohmic
    contribution disappears immediately. The remaining voltage then changes
    gradually as the concentration gradient relaxes. At OCV,

    $$j=F(J_i-J_e)=0,$$

    but this does **not** require $J_i=J_e=0$ everywhere. Instead,
    $J_i=J_e=-D^\delta\partial c/\partial x$ can be nonzero: neutral pairs
    redistribute internally without carrying net electrical current. Only after
    a sufficiently long rest are the composition and both chemical potentials
    spatially uniform.
    """)
    return


@app.cell
def _(
    gitt_result_06,
    np,
    pitt_result_06,
    potential_case_06,
    potential_decomposition,
    potential_time_06,
    positions_06,
    pulse_times_06,
    rest_times_06,
    parameters_06,
):
    _progress = potential_time_06.value / 100.0
    if "PITT" in potential_case_06.value:
        _selected_result = pitt_result_06
    else:
        _selected_result = gitt_result_06
    if "OCV" in potential_case_06.value:
        _selected_times = rest_times_06
        _selected_profiles = _selected_result["rest_profiles"]
        _selected_q = _selected_result["rest_q"]
        _stage_name = "OCV relaxation"
    else:
        _selected_times = pulse_times_06
        _selected_profiles = _selected_result["pulse_profiles"]
        _selected_q = _selected_result["pulse_q"]
        _stage_name = "driven pulse"
    potential_index_06 = int(np.argmin(np.abs(_selected_times - _progress * _selected_times[-1])))
    selected_profile_06 = _selected_profiles[potential_index_06]
    selected_q_06 = float(_selected_q[potential_index_06])
    selected_reduced_time_06 = float(_selected_times[potential_index_06])
    selected_stage_06 = _stage_name
    selected_potentials_06 = potential_decomposition(
        selected_profile_06,
        selected_q_06,
        positions_06,
        parameters_06,
    )
    return (
        potential_index_06,
        selected_potentials_06,
        selected_q_06,
        selected_reduced_time_06,
        selected_stage_06,
        selected_profile_06,
    )


@app.cell
def _(
    GAS_CONSTANT_J_PER_MOL_K,
    mo,
    parameters_06,
    plt,
    positions_06,
    potential_case_06,
    selected_potentials_06,
    selected_q_06,
    selected_reduced_time_06,
    selected_stage_06,
):
    _energy_scale = GAS_CONSTANT_J_PER_MOL_K * parameters_06["temperature_k"] / 1000.0
    _x_um = positions_06 * parameters_06["length_m"] * 1.0e6
    _fig, _axes = plt.subplots(1, 3, figsize=(15.0, 4.5), constrained_layout=True)
    _axes[0].plot(_x_um, selected_potentials_06["profile"][::-1], color="#4F7881", lw=1.9)
    _axes[0].set_title("Composition")
    _axes[0].set_ylabel("c / c0")

    _axes[1].plot(
        _x_um,
        _energy_scale * selected_potentials_06["mu_i"][::-1],
        color="#5F8F8D",
        lw=1.7,
        label="$\\mu_i$",
    )
    _axes[1].plot(
        _x_um,
        _energy_scale * selected_potentials_06["electrical_i"][::-1],
        color="#C49345",
        lw=1.7,
        label="$+F\\phi$",
    )
    _axes[1].plot(
        _x_um,
        _energy_scale * selected_potentials_06["tilde_mu_i"][::-1],
        color="#A65E5E",
        lw=1.9,
        label="$\\widetilde\\mu_i$",
    )
    _axes[1].set_title("Ion: chemical + electrical")
    _axes[1].set_ylabel("Change (kJ mol$^{-1}$)")
    _axes[1].legend(fontsize=10)

    _axes[2].plot(
        _x_um,
        _energy_scale * selected_potentials_06["mu_e"][::-1],
        color="#5F8F8D",
        lw=1.7,
        label="$\\mu_e$",
    )
    _axes[2].plot(
        _x_um,
        _energy_scale * selected_potentials_06["electrical_e"][::-1],
        color="#C49345",
        lw=1.7,
        label="$-F\\phi$",
    )
    _axes[2].plot(
        _x_um,
        _energy_scale * selected_potentials_06["tilde_mu_e"][::-1],
        color="#A65E5E",
        lw=1.9,
        label="$\\widetilde\\mu_e$",
    )
    _axes[2].set_title("Electron: chemical + electrical")
    _axes[2].set_ylabel("Change (kJ mol$^{-1}$)")
    _axes[2].legend(fontsize=10)
    for _axis in _axes:
        _axis.set_xlabel("Position x (micrometers)")
        _axis.grid(alpha=0.22)
    _fig.suptitle(
        f"{potential_case_06.value} | {selected_stage_06}, "
        f"t/tau = {selected_reduced_time_06:.3f}, "
        f"normalized current = {selected_q_06:.3f}",
        fontsize=14,
        weight="bold",
    )
    mo.accordion({
        "Explore further - composition and electrochemical-potential profiles":
        mo.vstack(
            [
                _fig,
                mo.md(
                    r"The $+F\phi$ and $-F\phi$ lines show the same electrostatic "
                    r"potential in molar-energy form. During OCV, sloped "
                    r"electrochemical potentials can still drive equal ion and "
                    r"electron fluxes even though electrical current is zero."
                ),
            ]
        )
    })
    return


@app.cell
def _(mo):
    approximation_reader = mo.md(r"""
    ## 5. Where the textbook approximations come from

    To isolate the standard formulas, now take the **classical one-sided
    limit**: electronic transport in the MIEC is much faster than ionic
    transport ($t_e\rightarrow1$), the perturbation is small, the slab is
    planar, and interfacial and series resistances are negligible. Let $A$ be
    the active area, $I=Aj$ the measured current, and $\Delta c$ the imposed
    change of surface concentration. We use $D^\delta$ for chemical
    diffusivity, as in Module 05.

    ### PITT: current after a surface-concentration step

    The finite-slab solution is

    $$
    I(t)=\frac{2FAD^\delta\Delta c}{L}
    \sum_{m=0}^{\infty}
    \exp\!\left[-(m+\tfrac12)^2\pi^2
    \frac{D^\delta t}{L^2}\right].
    $$

    Its short- and long-time limits are

    $$
    \underbrace{I(t)\simeq FA\Delta c
    \sqrt{\frac{D^\delta}{\pi t}}}_{D^\delta t/L^2\ll1\;\text{(Cottrell)}},
    \qquad
    \underbrace{I(t)\simeq\frac{2FAD^\delta\Delta c}{L}
    e^{-\pi^2D^\delta t/(4L^2)}}_{D^\delta t/L^2\gg1\;\text{(first mode)}}.
    $$

    The current sign depends on insertion versus extraction; the comparison
    plot shows its magnitude. The factor $1/4$ in the long-time exponent follows
    from a blocking boundary at one face and fixed composition at the other.

    ### GITT: voltage during a current pulse

    Let $I_0=Aj_0$ be the applied current. In this ideal planar model, the
    finite-slab diffusion-voltage response can be written

    $$
    \Delta E(t)=\frac{2RT}{F}\frac{j_0L}{FD^\delta c_0}
    \left[\frac{D^\delta t}{L^2}+\frac13-
    \frac{2}{\pi^2}\sum_{n=1}^{\infty}
    \frac{e^{-n^2\pi^2D^\delta t/L^2}}{n^2}\right].
    $$

    At short time,

    $$
    \Delta E(t)\simeq\frac{4RTj_0}{F^2c_0}
    \sqrt{\frac{t}{\pi D^\delta}},
    $$

    whereas the long-time expression is

    $$
    \Delta E(t)\simeq\frac{2RT}{F}\frac{j_0L}{FD^\delta c_0}
    \left(\frac{D^\delta t}{L^2}+\frac13\right).
    $$

    These are limiting descriptions of the same finite-slab diffusion problem;
    the full transient curves are not assembled from them.
    """)
    mo.accordion({"Advanced analysis - classical PITT and GITT limits": approximation_reader})
    return


@app.cell
def _(
    classical_gitt_series,
    classical_pitt_series,
    gitt_current_06,
    np,
    pitt_reduced_voltage_06,
    plt,
):
    approximation_times_06 = np.logspace(-3.0, 2.0, 360)
    approximation_fourier_time_06 = approximation_times_06 / np.pi**2
    classical_concentration_step_06 = 1.0 - np.exp(-pitt_reduced_voltage_06)
    pitt_series_06 = classical_pitt_series(approximation_times_06, classical_concentration_step_06)
    pitt_short_06 = classical_concentration_step_06 / np.sqrt(np.pi * approximation_fourier_time_06)
    pitt_long_06 = 2.0 * classical_concentration_step_06 * np.exp(-approximation_times_06 / 4.0)
    gitt_series_06 = classical_gitt_series(approximation_times_06, gitt_current_06.value)
    gitt_short_06 = 2.0 * gitt_current_06.value * np.sqrt(approximation_fourier_time_06 / np.pi)
    gitt_long_06 = gitt_current_06.value * (approximation_fourier_time_06 + 1.0 / 3.0)

    _fig, _axes = plt.subplots(1, 2, figsize=(13.8, 4.8), constrained_layout=True)
    _axes[0].loglog(
        approximation_fourier_time_06,
        pitt_series_06 / classical_concentration_step_06,
        color="#4F7881",
        lw=1.9,
        label="finite-slab solution",
    )
    _axes[0].loglog(
        approximation_fourier_time_06,
        pitt_short_06 / classical_concentration_step_06,
        "--",
        color="#C49345",
        lw=1.7,
        label="short-time expression",
    )
    _axes[0].loglog(
        approximation_fourier_time_06,
        pitt_long_06 / classical_concentration_step_06,
        ":",
        color="#A65E5E",
        lw=1.8,
        label="long-time expression",
    )
    _axes[0].set_title("PITT current")
    _axes[0].set_ylabel(r"$|I|L/(FAD^\delta|\Delta c|)$")

    _axes[1].loglog(
        approximation_fourier_time_06,
        gitt_series_06 / gitt_current_06.value,
        color="#4F7881",
        lw=1.9,
        label="finite-slab solution",
    )
    _axes[1].loglog(
        approximation_fourier_time_06,
        gitt_short_06 / gitt_current_06.value,
        "--",
        color="#C49345",
        lw=1.7,
        label="short-time expression",
    )
    _axes[1].loglog(
        approximation_fourier_time_06,
        gitt_long_06 / gitt_current_06.value,
        ":",
        color="#A65E5E",
        lw=1.8,
        label="long-time expression",
    )
    _axes[1].set_title("GITT voltage")
    _axes[1].set_ylabel(r"$F^2D^\delta c_0\Delta E/(2RTj_0L)$")
    for _axis in _axes:
        _axis.set_xlabel(r"Fourier time, $D^\delta t/L^2$")
        _axis.grid(which="both", alpha=0.22)
        _axis.legend()
    _fig.suptitle(
        "Short and long formulas are asymptotes, not universal fits", fontsize=15, weight="bold"
    )
    mo.accordion({"Advanced analysis - asymptotic response curves": _fig})
    return (
        approximation_fourier_time_06,
        approximation_times_06,
        classical_concentration_step_06,
        gitt_long_06,
        gitt_series_06,
        gitt_short_06,
        pitt_long_06,
        pitt_series_06,
        pitt_short_06,
    )


@app.cell
def _(mo):
    fitting_reader = mo.md(r"""
    ### Reading the usual experimental plots

    The short- and long-time formulas touch the same exact
    finite-slab response only in their own asymptotic windows; they are not
    interchangeable global fits.

    The asymptotes suggest simple diagnostics, provided the assumptions above
    have been tested. For PITT, the short-time plot of $I$ versus $t^{-1/2}$ is
    linear and

    $$D^\delta=\pi\left[\frac{I\sqrt{t}}
    {FA\Delta c}\right]^2.$$

    At long PITT times, a plot of $\ln|I|$ versus $t$ has slope
    $-\pi^2D^\delta/(4L^2)$.

    For a short GITT pulse of duration $\tau$, let $\Delta E_t$ be the gradual
    diffusion voltage during the pulse after removing the instantaneous Ohmic
    jump, and let $\Delta E_s$ be the equilibrium OCV change caused by that
    pulse. The planar small-signal result is

    $$
    D^\delta\simeq\frac{4L^2}{\pi\tau}
    \left(\frac{\Delta E_s}{\Delta E_t}\right)^2,
    \qquad \frac{D^\delta\tau}{L^2}\ll1.
    $$

    This familiar GITT formula is a **short-pulse result**, not a definition of
    diffusivity. The surface area, diffusion length, OCV slope, and removal of
    non-diffusive voltage drops must all be consistent with the specimen.

    ### OCV after a pulse remembers how the pulse was made

    Immediately after interruption, the concentration profile is the final
    profile of the voltage or current pulse. The OCV stage must therefore use
    that profile as its initial condition.

    For a **short GITT pulse** in a semi-infinite planar solid, the surface
    relaxation has the time dependence

    $$
    \Delta c_s(t_r)\propto
    \sqrt{t_r+t_p}-\sqrt{t_r},
    $$

    where $t_p$ is the pulse duration and $t_r$ starts at current interruption.
    This is not simply a universal $1/\sqrt{t_r}$ law. For a finite slab at long
    rest times, the slowest remaining spatial variation gives

    $$c(x,t_r)-c_{\rm eq}\propto\exp(-t_r/\tau^\delta),$$

    where $c_{\rm eq}=\langle c\rangle$ is the uniform concentration reached
    after that pulse; it is generally not the initial value $c_0$. The notebook
    carries the complete end-of-pulse profile into the OCV stage,
    so it does not need either approximation to generate the relaxation curves.
    """)
    mo.accordion({"Advanced analysis - fitting windows and OCV limits": fitting_reader})
    return




@app.cell
def _(mo):
    kinetics_profile_time_06 = mo.ui.slider(
        -3.0, 0.0, value=-1.3, step=0.1,
        label=r"profile time, $\log_{10}\theta$", show_value=True,
    )
    return (kinetics_profile_time_06,)


@app.cell
def _(
    cumulative_trapezoid,
    finite_pitt_fit_bias_06,
    finite_pitt_kinetics_06,
    finite_pitt_modes_06,
    np,
):
    finite_biot_cases_06 = (np.inf, 100.0, 1.0, 0.01)
    finite_theta_06 = np.geomspace(1.0e-4, 2.0, 340)
    finite_position_06 = np.linspace(0.0, 1.0, 260)
    finite_responses_06 = {
        biot: finite_pitt_kinetics_06(
            finite_theta_06, finite_position_06, biot, mode_count=120
        ) for biot in finite_biot_cases_06
    }
    finite_mass_balance_error_06 = max(
        np.max(np.abs(
            finite_responses_06[biot]["mean_profile"][0]
            - finite_responses_06[biot]["mean_profile"]
            - cumulative_trapezoid(
                finite_responses_06[biot]["current"], finite_theta_06, initial=0.0
            )
        ))
        for biot in finite_biot_cases_06
    )
    fit_biot_grid_06 = np.geomspace(1.0e-3, 1.0e3, 76)
    fitting_windows_06 = {
        "transition window, 0.03–0.12": (0.03, 0.12),
        "late window, 0.30–1.00": (0.30, 1.00),
    }
    fit_bias_06 = finite_pitt_fit_bias_06(
        fit_biot_grid_06, fitting_windows_06, mode_count=70
    )
    _roots_fast, _ = finite_pitt_modes_06(1.0e6, 1)
    _roots_slow, _ = finite_pitt_modes_06(1.0e-6, 1)
    finite_fast_limit_error_06 = abs(_roots_fast[0] / (0.5 * np.pi) - 1.0)
    finite_slow_limit_error_06 = abs(_roots_slow[0] ** 2 / 1.0e-6 - 1.0)
    finite_bias_direction_pass_06 = bool(
        fit_bias_06["late window, 0.30–1.00"][0] < 0.02
        and fit_bias_06["late window, 0.30–1.00"][-1] > 0.94
    )
    return (
        fit_bias_06,
        fit_biot_grid_06,
        finite_bias_direction_pass_06,
        finite_biot_cases_06,
        finite_mass_balance_error_06,
        finite_fast_limit_error_06,
        finite_position_06,
        finite_responses_06,
        finite_slow_limit_error_06,
        finite_theta_06,
        fitting_windows_06,
    )


@app.cell
def _(
    finite_biot_cases_06,
    finite_pitt_kinetics_06,
    finite_position_06,
    kinetics_profile_time_06,
    np,
):
    selected_finite_theta_06 = 10.0 ** kinetics_profile_time_06.value
    selected_finite_profiles_06 = {
        biot: finite_pitt_kinetics_06(
            np.array([selected_finite_theta_06]),
            finite_position_06,
            biot,
            mode_count=140,
        )["profile"][0]
        for biot in finite_biot_cases_06
    }
    finite_minimum_concentration_06 = min(
        np.min(0.85 + 0.15 * profile)
        for profile in selected_finite_profiles_06.values()
    )
    return (
        finite_minimum_concentration_06,
        selected_finite_profiles_06,
        selected_finite_theta_06,
    )


@app.cell
def _(
    fit_bias_06,
    fit_biot_grid_06,
    finite_bias_direction_pass_06,
    finite_biot_cases_06,
    finite_mass_balance_error_06,
    finite_fast_limit_error_06,
    finite_minimum_concentration_06,
    finite_position_06,
    finite_responses_06,
    finite_slow_limit_error_06,
    finite_theta_06,
    fitting_windows_06,
    kinetics_profile_time_06,
    mo,
    np,
    plt,
    selected_finite_profiles_06,
    selected_finite_theta_06,
):
    _colors = ("#4F7881", "#6B86A5", "#B8734A", "#7C6A91")
    _styles = ("-", "--", "-.", ":")
    _labels = (
        r"$\mathrm{Bi}=\infty$", r"$\mathrm{Bi}=100$",
        r"$\mathrm{Bi}=1$", r"$\mathrm{Bi}=0.01$",
    )
    _figure, _axes = plt.subplots(1, 3, figsize=(15.0, 4.8), constrained_layout=True)
    for _biot, _color, _style, _label in zip(
        finite_biot_cases_06, _colors, _styles, _labels
    ):
        _axes[0].loglog(
            finite_theta_06, finite_responses_06[_biot]["current"],
            color=_color, ls=_style, lw=1.9, label=_label,
        )
        _axes[1].plot(
            finite_position_06,
            0.85 + 0.15 * selected_finite_profiles_06[_biot],
            color=_color, ls=_style, lw=1.9, label=_label,
        )
    _axes[0].set(
        xlabel=r"Diffusion time, $\theta=D^\delta t/L^2$",
        ylabel="Normalized PITT current",
        title="Current separates diffusion and reaction control",
    )
    _axes[1].set(
        xlabel=r"Article coordinate, $x/L$",
        ylabel=r"Illustrative $c/c_0$",
        title=rf"Profiles at $\theta={selected_finite_theta_06:.3g}$",
    )
    for _name, _color, _style in zip(
        fitting_windows_06, ("#B8734A", "#4F7881"), ("--", "-")
    ):
        _axes[2].semilogx(
            fit_biot_grid_06, fit_bias_06[_name], color=_color,
            ls=_style, lw=1.9, label=_name,
        )
    _axes[2].axhline(1.0, color="#73808C", lw=1.2, ls=":")
    _axes[2].set(
        xlabel=r"Surface Biot number, $\mathrm{Bi}$",
        ylabel=r"$D^\delta_{\rm inferred}/D^\delta_{\rm true}$",
        title="A diffusion-only fit can underestimate diffusivity",
        ylim=(0.0, 1.35),
    )
    for _axis in _axes:
        _axis.grid(alpha=0.22, which="both")
        _axis.legend(fontsize=10.5, loc="best")
    _checks_pass = (
        finite_fast_limit_error_06 < 2.0e-6
        and finite_slow_limit_error_06 < 2.0e-6
        and finite_mass_balance_error_06 < 2.0e-4
        and finite_minimum_concentration_06 > 0.0
        and finite_bias_direction_pass_06
    )
    _content = mo.vstack([
        mo.md(r"""
        ### Advanced extension: finite surface kinetics

        Ideal PITT fixes the surface chemical potential instantly. A finite
        linearized surface reaction instead imposes

        $$-D^\delta\left.\frac{\partial c}{\partial x}\right|_{0}
        =k^\delta[c_s-c(0,t)],\qquad
        \left.\frac{\partial c}{\partial x}\right|_{L}=0.$$

        The competition is summarized by

        $$\mathrm{Bi}=\frac{k^\delta L}{D^\delta}
        =\frac{\tau_d}{\tau_{ct}},\qquad
        \tau_d=\frac{L^2}{D^\delta},\quad
        \tau_{ct}=\frac{L}{k^\delta}.$$

        - $\mathrm{Bi}<0.1$: nearly uniform, reaction-limited response.
        - $0.1\lesssim\mathrm{Bi}\lesssim10$: mixed control.
        - $\mathrm{Bi}>10$: the diffusion-controlled limit is approached.

        The finite-slab eigenvalues obey $\lambda\tan\lambda=\mathrm{Bi}$.
        The curves are full transients, not patched limiting expressions.
        """),
        kinetics_profile_time_06,
        _figure,
        mo.md(rf"""
        A straight line in
        $\ln|I|$ versus $t$ does not prove diffusion control. Finite exchange
        makes the slowest eigenvalue smaller than $\pi/2$; using the ideal PITT
        slope then returns a diffusivity that is too small. The inferred value
        also depends on the fitting window.

        **GITT measurement decomposition.** Under a current step,

        $$E_{{\rm meas}}(t)=E_{{\rm eq}}[c(0,t)]+\eta_{{ct}}(t)+I R_\Omega.$$

        The instantaneous $IR_\Omega$ jump, gradual equilibrium-composition
        term, and surface overpotential should not be assigned to one diffusion
        formula. Constant $\eta_{{ct}}$ may leave a short-time $\sqrt t$ slope
        unchanged; composition-dependent kinetics need not.

        **Finite-kinetics checks: {'PASS' if _checks_pass else 'CHECK'}.** The
        reaction- and diffusion-controlled eigenvalue limits agree; integrated
        surface flux matches the composition change (maximum normalized error
        {finite_mass_balance_error_06:.2e}); concentrations remain positive;
        and the late-window fitting bias moves from severe underestimation at
        small Bi toward unity at large Bi.
        """),
    ])
    mo.accordion({"Advanced reader — surface kinetics and fitting bias": _content})
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## 6. What can PITT and GITT determine?

    Under the assumptions of this notebook, the measurements separate two
    kinds of information:

    - **Thermodynamics:** the long-rest OCV gives the equilibrium potential at
      the new mean composition. Repeating steps traces an equilibrium titration
      curve and, through its slope, a differential chemical capacity.
    - **Transport:** the pulse and relaxation time scales give the chemical
      diffusivity $D^\delta$. The relevant length must be the actual diffusion
      length of the specimen.
    - **Carrier balance:** the immediate and gradual parts of the voltage
      response reveal how electronic and ionic transport share the drive.

    **Assumptions checklist before extracting $D^\delta$**

    - Is the step small enough that $D^\delta$ and $dE/d\delta$ are locally
      constant?
    - Is the active geometry one-dimensional, with known $L$ and $S$?
    - Is the fitting window inside the short- or long-time regime used by the
      formula?
    - Are $IR_\Omega$, surface overpotential, and early capacitive response
      separated from diffusion?
    - Does the inferred value remain stable when pulse size, duration, and fit
      window are changed?

    A fitted number is not automatically a material constant. The standard
    formulas assume a small step, a single phase, constant $D^\delta$, planar
    one-dimensional transport, known active area and diffusion length, and
    negligible charge-transfer, contact, and uncompensated series resistances.
    Porous composite electrodes, phase transformations, strongly
    composition-dependent diffusivity, and early-time capacitive transients
    require a more complete model. A useful experimental habit is to vary pulse
    size and duration: a real diffusion coefficient should not depend strongly
    on either within the valid regime.
    """)
    return


@app.cell
def _(
    FARADAY_C_PER_MOL,
    GAS_CONSTANT_J_PER_MOL_K,
    approximation_times_06,
    classical_concentration_step_06,
    gitt_result_06,
    gitt_series_06,
    gitt_short_06,
    parameters_06,
    pitt_reduced_voltage_06,
    pitt_result_06,
    pitt_series_06,
    pitt_short_06,
    positions_06,
    potential_decomposition,
    pulse_times_06,
    rest_times_06,
    np,
):
    initial_uniform_error_06 = max(
        np.max(np.abs(pitt_result_06["pulse_profiles"][0] - 1.0)),
        np.max(np.abs(gitt_result_06["pulse_profiles"][0] - 1.0)),
    )
    minimum_concentration_06 = min(
        np.min(pitt_result_06["pulse_profiles"]),
        np.min(pitt_result_06["rest_profiles"]),
        np.min(gitt_result_06["pulse_profiles"]),
        np.min(gitt_result_06["rest_profiles"]),
    )
    pitt_voltage_control_error_06 = np.max(
        np.abs(pitt_result_06["pulse_u"] - pitt_reduced_voltage_06)
    )
    gitt_current_control_error_06 = np.max(
        np.abs(gitt_result_06["pulse_q"] - gitt_result_06["pulse_q"][0])
    )
    pitt_rest_current_error_06 = np.max(np.abs(pitt_result_06["rest_q"]))
    gitt_rest_current_error_06 = np.max(np.abs(gitt_result_06["rest_q"]))

    pitt_spatial_mean_06 = np.trapezoid(pitt_result_06["pulse_profiles"], positions_06, axis=1)
    gitt_spatial_mean_06 = np.trapezoid(gitt_result_06["pulse_profiles"], positions_06, axis=1)
    pitt_expected_mean_06 = 1.0 - cumulative_trapezoid(
        pitt_result_06["pulse_q"], pulse_times_06, initial=0.0
    ) / np.pi**2
    pitt_mass_balance_error_06 = np.max(
        np.abs(pitt_spatial_mean_06 - pitt_expected_mean_06)
    )
    gitt_mass_balance_error_06 = np.max(
        np.abs(
            gitt_spatial_mean_06
            - (1.0 - gitt_result_06["pulse_q"][0] * pulse_times_06 / np.pi**2)
        )
    )
    ocv_mean_drift_06 = max(
        np.ptp(np.trapezoid(pitt_result_06["rest_profiles"], positions_06, axis=1)),
        np.ptp(np.trapezoid(gitt_result_06["rest_profiles"], positions_06, axis=1)),
    )
    composition_balance_error_06 = max(
        pitt_mass_balance_error_06,
        gitt_mass_balance_error_06,
        ocv_mean_drift_06,
    )

    identity_diffusivity_06 = (
        2.0
        * GAS_CONSTANT_J_PER_MOL_K
        * parameters_06["temperature_k"]
        * parameters_06["conductivity_i_s_per_m"]
        * parameters_06["conductivity_e_s_per_m"]
        / (
            FARADAY_C_PER_MOL**2
            * parameters_06["concentration_mol_per_m3"]
            * parameters_06["conductivity_total_s_per_m"]
        )
    )
    diffusivity_identity_error_06 = abs(
        identity_diffusivity_06 / parameters_06["chemical_diffusivity_m2_per_s"] - 1.0
    )

    _sample_indices = np.unique(np.linspace(0, pulse_times_06.size - 1, 9, dtype=int))
    _reconstruction_errors = []
    for _sample_index in _sample_indices:
        _potential_data = potential_decomposition(
            pitt_result_06["pulse_profiles"][_sample_index],
            pitt_result_06["pulse_q"][_sample_index],
            positions_06,
            parameters_06,
        )
        _reconstruction_errors.append(
            abs(_potential_data["reconstructed_u"] - pitt_result_06["pulse_u"][_sample_index])
        )
    voltage_reconstruction_error_06 = max(_reconstruction_errors)

    _short_index = int(np.argmin(np.abs(approximation_times_06 - 1.0e-3)))
    _long_index = int(np.argmin(np.abs(approximation_times_06 - 40.0)))
    pitt_short_error_06 = abs(pitt_series_06[_short_index] / pitt_short_06[_short_index] - 1.0)
    pitt_long_reference_06 = (
        2.0 * classical_concentration_step_06 * np.exp(-approximation_times_06[_long_index] / 4.0)
    )
    pitt_long_error_06 = abs(pitt_series_06[_long_index] / pitt_long_reference_06 - 1.0)
    gitt_short_error_06 = abs(gitt_series_06[_short_index] / gitt_short_06[_short_index] - 1.0)
    gitt_long_reference_06 = gitt_result_06["pulse_q"][0] * (
        approximation_times_06[_long_index] / np.pi**2 + 1.0 / 3.0
    )
    gitt_long_error_06 = abs(gitt_series_06[_long_index] / gitt_long_reference_06 - 1.0)
    _first_shape = np.cos(np.pi * positions_06)
    _rest_means = np.trapezoid(gitt_result_06["rest_profiles"], positions_06, axis=1)
    _first_amplitudes = np.trapezoid(
        (gitt_result_06["rest_profiles"] - _rest_means[:, None]) * _first_shape[None, :],
        positions_06,
        axis=1,
    )
    first_mode_relaxation_error_06 = abs(
        _first_amplitudes[-1] / _first_amplitudes[0] - np.exp(-rest_times_06[-1])
    )
    return (
        composition_balance_error_06,
        diffusivity_identity_error_06,
        first_mode_relaxation_error_06,
        gitt_current_control_error_06,
        gitt_long_error_06,
        gitt_mass_balance_error_06,
        gitt_rest_current_error_06,
        gitt_short_error_06,
        initial_uniform_error_06,
        minimum_concentration_06,
        ocv_mean_drift_06,
        pitt_long_error_06,
        pitt_rest_current_error_06,
        pitt_short_error_06,
        pitt_voltage_control_error_06,
        voltage_reconstruction_error_06,
    )


@app.cell
def _(
    composition_balance_error_06,
    diffusivity_identity_error_06,
    first_mode_relaxation_error_06,
    gitt_current_control_error_06,
    gitt_long_error_06,
    gitt_rest_current_error_06,
    gitt_short_error_06,
    initial_uniform_error_06,
    minimum_concentration_06,
    mo,
    pitt_long_error_06,
    pitt_rest_current_error_06,
    pitt_short_error_06,
    pitt_voltage_control_error_06,
    voltage_reconstruction_error_06,
):
    _checks = [
        (
            "Uniform initial specimen",
            initial_uniform_error_06 < 1.0e-12,
            "Both experiments must begin from the same equilibrated composition.",
        ),
        (
            "Positive concentrations",
            minimum_concentration_06 > 0.0,
            "Logarithmic chemical potentials require c > 0 everywhere.",
        ),
        (
            "PITT voltage held fixed",
            pitt_voltage_control_error_06 < 2.0e-9,
            "The current must adjust so that the requested voltage step is maintained.",
        ),
        (
            "GITT current held fixed",
            gitt_current_control_error_06 < 1.0e-14,
            "The voltage, not the imposed current, is allowed to evolve during GITT.",
        ),
        (
            "Zero current during both OCV rests",
            max(pitt_rest_current_error_06, gitt_rest_current_error_06) < 1.0e-14,
            "Open circuit stops charge transfer through the terminals while "
            "internal diffusion continues.",
        ),
        (
            "Composition balance",
            composition_balance_error_06 < 2.0e-4,
            "Integrated current changes the mean during a pulse; the mean is "
            "conserved during rest.",
        ),
        (
            "Chemical-diffusivity identity",
            diffusivity_identity_error_06 < 1.0e-12,
            "The conductivity and diffusivity parameters must describe the same ideal H+/e- pair.",
        ),
        (
            "Potential decomposition reconstructs voltage",
            voltage_reconstruction_error_06 < 3.0e-4,
            "Chemical and electrical contributions must add back to the measured terminal drive.",
        ),
        (
            "PITT short- and long-time limits",
            max(pitt_short_error_06, pitt_long_error_06) < 2.0e-3,
            "The short- and long-time expressions must approach the finite-slab "
            "response in their own regimes.",
        ),
        (
            "GITT short- and long-time limits",
            max(gitt_short_error_06, gitt_long_error_06) < 2.0e-3,
            "The square-root and late linear forms must approach the same finite-slab solution.",
        ),
        (
            "OCV long-time relaxation",
            first_mode_relaxation_error_06 < 2.0e-4,
            "At long rest, the slowest spatial relaxation sets the time constant tau_delta.",
        ),
    ]
    _rows = []
    for _name, _passed, _why in _checks:
        _status = "PASS" if _passed else "CHECK"
        _rows.append(f"| {_status} | {_name} | {_why} |")
    _table_lines = [
        "## 7. Physical consistency checks",
        "",
        "Each check protects a physical link between the imposed pulse, the "
        "composition profile, and the measured response.",
        "",
        "| status | check | why it matters |",
        "|---:|---|---|",
    ]
    _checks_table = mo.md("\n".join(_table_lines + _rows))
    mo.accordion({"Physical consistency checks": _checks_table})
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## What to carry forward

    $$
    \boxed{
    \text{controlled pulse}
    \rightarrow \text{selective carrier fluxes}
    \rightarrow c(x,t)
    \rightarrow \mu_i,\mu_e,\phi
    \rightarrow E(t)\text{ or }j(t)
    }
    $$

    1. **Titration separates state from rate.** Integrated charge gives
       $\Delta\delta$, long-rest $E(\delta)$ gives thermodynamics, and the
       transient carries kinetic information.
    2. **PITT and GITT impose different boundaries.** PITT fixes a surface
       chemical potential, GITT fixes a flux, and OCV sets terminal current to
       zero while internal neutral diffusion may continue.
    3. **$D^\delta$ is inferred through assumptions.** Finite surface kinetics,
       ohmic drop, geometry, and fitting window can bias it; consistency tests
       are part of the measurement.

    ### Sources and further reading

    - Q. Lu, [Solid State Ionics tutorials](https://ssi-westlake.com/tutorial/),
      including [PITT/GITT Part I](https://mp.weixin.qq.com/s/ktu9MiGhfYrE6l563pE38g)
      and [Part II](https://mp.weixin.qq.com/s/AzyRL3cZv6heEcB40wZplg), which
      motivated this English notebook.
    - W. Weppner and R. A. Huggins, “Determination of the Kinetic Parameters of
      Mixed-Conducting Electrodes and Application to the System Li3Sb,”
      *Journal of The Electrochemical Society* **124** (1977),
      [doi:10.1149/1.2133112](https://doi.org/10.1149/1.2133112).
    - C. J. Wen, B. A. Boukamp, R. A. Huggins, and W. Weppner,
      “Thermodynamic and Mass Transport Properties of LiAl,”
      *Journal of The Electrochemical Society* **126** (1979),
      [doi:10.1149/1.2128939](https://doi.org/10.1149/1.2128939).
    - S. D. Kang and W. C. Chueh, “Galvanostatic Intermittent Titration
      Technique Reinvented: Part I. A Critical Review,” *Journal of The
      Electrochemical Society* **168** (2021),
      [open article and record](https://www.osti.gov/biblio/1838037).
    - BioLogic, [Application Note 70: diffusion coefficients by EIS, PITT, and GITT](https://www.biologic.net/documents/determination-diffusion-coefficient-inserted-species-host-electrode-eis-pitt-gitt-techniques-an70/).

    The online sources were used as evidence, not copied as instructions. All
    equations here were rederived for the selective-contact geometry and sign
    convention stated in this notebook.
    """)
    return


@app.cell
def _(plt):
    plt.rcParams.update(
        {
            "font.size": 15,
            "axes.titlesize": 17,
            "axes.labelsize": 15,
            "legend.fontsize": 12,
            "lines.solid_capstyle": "round",
            "figure.dpi": 115,
        }
    )
    return


@app.cell
def _(mo):
    mo.md(r"""
    **Continue:** [Module 07 — Impedance, Warburg Diffusion, and Transmission Lines](https://qiyanglu.github.io/Solid-State-Ionics-Interactive/07-impedance-tlm/)
    """)
    return



if __name__ == "__main__":
    app.run()
