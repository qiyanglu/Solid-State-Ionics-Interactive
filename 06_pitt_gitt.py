import marimo

__generated_with = "0.24.0"
app = marimo.App(width="full")


@app.cell
def _():
    import marimo as mo
    import matplotlib.pyplot as plt
    import numpy as np
    from scipy.integrate import cumulative_trapezoid, solve_ivp

    return cumulative_trapezoid, mo, np, plt, solve_ivp


@app.cell
def _(mo):
    mo.md(r"""
    # PITT and GITT: Watching Chemical Diffusion

    **Guiding question.** How can a voltage or current pulse reveal both the
    equilibrium composition and the rate of chemical diffusion in a mixed
    ionic-electronic conductor (MIEC)?

    PITT and GITT use the same pulse--relaxation idea with different controls.
    PITT is the **potentiostatic intermittent titration technique**; GITT is the
    corresponding **galvanostatic intermittent titration technique**.

    | technique | controlled during the pulse | measured response | especially useful for |
    |---|---|---|---|
    | **PITT** | voltage step | current decay | $D^\delta$, differential capacity |
    | **GITT** | current step | voltage transient | OCV curve, $D^\delta$ |

    After each pulse, the circuit is opened. The current becomes zero, but a
    nonuniform composition can still relax internally. Repeating small steps
    maps an equilibrium titration curve while the transients contain kinetic
    information.

    This notebook does more than draw the standard straight-line fits. It solves
    the finite, one-dimensional diffusion problem first, then asks when the
    familiar short- and long-time approximations are actually justified.
    """)
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## 1. The selective-contact experiment

    We keep the ideal pair model and notation of Module 05,

    $$H \rightleftharpoons H^+ + e^-, \qquad c_i=c_e=c(x,t),$$

    but change the contacts. The MIEC occupies $0\le x\le L$.

    - At $x=0$, the **current collector passes electrons and blocks ions**:
      $J_i(0,t)=0$.
    - At $x=L$, the **electrolyte passes ions and blocks electrons**:
      $J_e(L,t)=0$.

    The two carriers therefore enter or leave through opposite faces. Local
    electroneutrality couples them into a neutral composition change in the
    MIEC. Thin interfacial layers and charge-transfer kinetics are not resolved
    in this introductory bulk-diffusion model.

    We define

    $$j=F(J_i-J_e),$$

    and choose **positive $j$ and positive $\Delta U$ to extract $H$**, so the
    mean concentration decreases. This sign choice is stated here because
    battery literature often uses the opposite charge/discharge convention.
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

    def _mode_numbers(mode_count):
        count = int(mode_count)
        if count < 8:
            raise ValueError("mode_count must be at least eight")
        return np.arange(1, count + 1, dtype=float)

    def _state_profiles(states, positions, modes):
        state_array = np.asarray(states, dtype=float)
        cosine_matrix = np.cos(np.pi * np.outer(positions, modes))
        if state_array.ndim == 1:
            return state_array[0] + cosine_matrix @ state_array[1:]
        return state_array[:, :1] + state_array[:, 1:] @ cosine_matrix.T

    def _state_voltage_and_current(state, target_voltage, positions, modes, params):
        profile = _state_profiles(state, positions, modes)
        safe_profile = np.maximum(profile, 1.0e-10)
        resistance_integral = np.trapezoid(1.0 / safe_profile, positions)
        t_i = params["ionic_fraction"]
        t_e = params["electronic_fraction"]
        if target_voltage is None:
            raise ValueError("target_voltage is required")
        q_value = (
            target_voltage + t_i * np.log(safe_profile[0]) + t_e * np.log(safe_profile[-1])
        ) / (t_i * t_e * resistance_integral)
        return q_value

    def _voltage_from_state(state, q_value, positions, modes, params):
        profile = _state_profiles(state, positions, modes)
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

    def _modal_rhs_for_q(state, q_value, modes, params):
        t_i = params["ionic_fraction"]
        t_e = params["electronic_fraction"]
        forcing = -(2.0 * q_value / np.pi**2) * (t_i + t_e * (-1.0) ** modes)
        return np.concatenate(([-q_value / np.pi**2], -(modes**2) * state[1:] + forcing))

    def simulate_pitt(
        reduced_pulse_times,
        reduced_rest_times,
        reduced_voltage,
        positions,
        params,
        mode_count=36,
    ):
        """Finite-slab PITT pulse followed by an open-circuit relaxation."""
        pulse_times = np.asarray(reduced_pulse_times, dtype=float)
        rest_times = np.asarray(reduced_rest_times, dtype=float)
        xi = np.asarray(positions, dtype=float)
        modes = _mode_numbers(mode_count)
        initial_state = np.concatenate(([1.0], np.zeros(modes.size)))

        def voltage_step_rhs(_reduced_time, state_now):
            q_now = _state_voltage_and_current(state_now, reduced_voltage, xi, modes, params)
            return _modal_rhs_for_q(state_now, q_now, modes, params)

        pulse_solution = solve_ivp(
            voltage_step_rhs,
            (float(pulse_times[0]), float(pulse_times[-1])),
            initial_state,
            t_eval=pulse_times,
            method="BDF",
            rtol=2.0e-8,
            atol=2.0e-10,
        )
        if not pulse_solution.success:
            raise RuntimeError(str(pulse_solution.message))
        pulse_states = pulse_solution.y.T
        pulse_current = np.array(
            [
                _state_voltage_and_current(state, reduced_voltage, xi, modes, params)
                for state in pulse_states
            ]
        )
        pulse_voltage = np.array(
            [
                _voltage_from_state(state, q, xi, modes, params)
                for state, q in zip(pulse_states, pulse_current)
            ]
        )

        final_state = pulse_states[-1]
        rest_states = np.empty((rest_times.size, final_state.size))
        rest_states[:, 0] = final_state[0]
        rest_states[:, 1:] = final_state[1:] * np.exp(-np.outer(rest_times, modes**2))
        rest_current = np.zeros(rest_times.size)
        rest_voltage = np.array(
            [_voltage_from_state(state, 0.0, xi, modes, params) for state in rest_states]
        )
        return {
            "modes": modes,
            "pulse_states": pulse_states,
            "rest_states": rest_states,
            "pulse_profiles": _state_profiles(pulse_states, xi, modes),
            "rest_profiles": _state_profiles(rest_states, xi, modes),
            "pulse_q": pulse_current,
            "rest_q": rest_current,
            "pulse_u": pulse_voltage,
            "rest_u": rest_voltage,
            "solver_message": str(pulse_solution.message),
        }

    def simulate_gitt(
        reduced_pulse_times,
        reduced_rest_times,
        q_step,
        positions,
        params,
        mode_count=36,
    ):
        """Analytical finite-slab GITT pulse followed by open circuit."""
        pulse_times = np.asarray(reduced_pulse_times, dtype=float)
        rest_times = np.asarray(reduced_rest_times, dtype=float)
        xi = np.asarray(positions, dtype=float)
        modes = _mode_numbers(mode_count)
        forcing = -(2.0 * q_step / np.pi**2) * (
            params["ionic_fraction"] + params["electronic_fraction"] * (-1.0) ** modes
        )
        pulse_states = np.empty((pulse_times.size, modes.size + 1))
        pulse_states[:, 0] = 1.0 - q_step * pulse_times / np.pi**2
        pulse_states[:, 1:] = (
            forcing[None, :]
            / modes[None, :] ** 2
            * (1.0 - np.exp(-np.outer(pulse_times, modes**2)))
        )
        pulse_q = np.full(pulse_times.size, float(q_step))
        pulse_u = np.array(
            [_voltage_from_state(state, q_step, xi, modes, params) for state in pulse_states]
        )

        final_state = pulse_states[-1]
        rest_states = np.empty((rest_times.size, final_state.size))
        rest_states[:, 0] = final_state[0]
        rest_states[:, 1:] = final_state[1:] * np.exp(-np.outer(rest_times, modes**2))
        rest_q = np.zeros(rest_times.size)
        rest_u = np.array(
            [_voltage_from_state(state, 0.0, xi, modes, params) for state in rest_states]
        )
        return {
            "modes": modes,
            "pulse_states": pulse_states,
            "rest_states": rest_states,
            "pulse_profiles": _state_profiles(pulse_states, xi, modes),
            "rest_profiles": _state_profiles(rest_states, xi, modes),
            "pulse_q": pulse_q,
            "rest_q": rest_q,
            "pulse_u": pulse_u,
            "rest_u": rest_u,
            "solver_message": "analytical Fourier solution",
        }

    def potential_decomposition(state, q_value, positions, modes, params):
        """Chemical, electrical, and electrochemical profiles in RT units."""
        xi = np.asarray(positions, dtype=float)
        profile = _state_profiles(state, xi, modes)
        t_i = params["ionic_fraction"]
        t_e = params["electronic_fraction"]
        mu_i = np.log(profile)
        mu_e = np.log(profile)
        resistance_integral = np.zeros_like(xi)
        resistance_integral[1:] = np.cumsum(
            0.5 * (1.0 / profile[1:] + 1.0 / profile[:-1]) * np.diff(xi)
        )
        # Integrated form of
        # d(tilde_mu_i/RT)/dxi = 2 t_e d(ln y)/dxi - 2 q t_i t_e/y.
        # It remains accurate at selective boundaries without differentiating a
        # truncated Fourier series.
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
def _(mo, plt):
    _fig, _ax = plt.subplots(figsize=(12.0, 2.7), constrained_layout=True)
    _ax.set_xlim(-0.12, 1.12)
    _ax.set_ylim(-0.42, 0.48)
    _ax.axvspan(-0.1, 0.0, color="#4f5d75", alpha=0.95)
    _ax.axvspan(0.0, 1.0, color="#d8f3dc", alpha=0.95)
    _ax.axvspan(1.0, 1.1, color="#90e0ef", alpha=0.95)
    _ax.text(
        -0.05, 0.12, "current\ncollector", ha="center", va="center", color="white", weight="bold"
    )
    _ax.text(
        0.5,
        0.15,
        "MIEC:  $H \u21cc H^+ + e^-$",
        ha="center",
        va="center",
        fontsize=15,
        weight="bold",
    )
    _ax.text(
        1.05, 0.12, "ion\nelectrolyte", ha="center", va="center", color="#073b4c", weight="bold"
    )
    _ax.annotate(
        "$e^-$ passes",
        xy=(-0.02, -0.12),
        xytext=(0.20, -0.12),
        arrowprops={"arrowstyle": "->", "lw": 2},
        ha="center",
    )
    _ax.annotate(
        "$H^+$ passes",
        xy=(1.02, -0.12),
        xytext=(0.80, -0.12),
        arrowprops={"arrowstyle": "->", "lw": 2},
        ha="center",
    )
    _ax.text(0.02, -0.31, "$J_i(0,t)=0$", ha="left", color="#9b2226", weight="bold")
    _ax.text(0.98, -0.31, "$J_e(L,t)=0$", ha="right", color="#9b2226", weight="bold")
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
    mo.md(r"""
    ## 2. From carrier fluxes to one diffusion equation

    With ideal dilute chemical potentials

    $$
    \mu_i=\mu_i^0+RT\ln(c/\bar c),\qquad
    \mu_e=\mu_e^0+RT\ln(c/\bar c),
    $$

    the neutral-pair chemical potential is
    $\mu_H=\mu_i+\mu_e=\mu_H^0+2RT\ln(c/\bar c)$. Eliminating the
    internal electric field gives exactly the same chemical diffusivity as in
    Module 05,

    $$
    D^\delta=\frac{2D_iD_e}{D_i+D_e},\qquad
    \tau^\delta=\frac{L^2}{\pi^2D^\delta}.
    $$

    The coupled fluxes become

    $$
    J_i=t_i\frac{j}{F}-D^\delta\frac{\partial c}{\partial x},\qquad
    J_e=-t_e\frac{j}{F}-D^\delta\frac{\partial c}{\partial x},
    $$

    where $t_i=\sigma_i/(\sigma_i+\sigma_e)$ and $t_e=1-t_i$.
    Conservation gives the full finite-slab equation

    $$\frac{\partial c}{\partial t}=D^\delta\frac{\partial^2c}{\partial x^2}.$$

    Introduce $y=c/\bar c_0$, $\xi=x/L$, $s=t/\tau^\delta$, and

    $$q=\frac{jL}{F D^\delta\bar c_0}.$$

    The selective contacts impose

    $$y_\xi(0)=t_iq,\qquad y_\xi(1)=-t_eq,$$

    and the mean composition obeys $d\bar y/ds=-q/\pi^2$. Thus a
    measured current is also a direct composition balance.

    The voltage measured between the selective terminals is the neutral-pair
    electrochemical drive. Relative to the initial state, define

    $$
    u=\frac{F\Delta U}{2RT}
    =-t_i\ln y(0)-t_e\ln y(1)
    +q t_i t_e\int_0^1\frac{d\xi}{y}.
    $$

    The first two terms are concentration polarization at the two selective
    faces; the integral is the internal Ohmic contribution. PITT holds $u$
    fixed and solves this relation for $q(t)$. GITT holds $q$ fixed and evaluates
    $u(t)$. At OCV, $q=0$ while the inherited concentration profile relaxes.
    """)
    return


@app.cell
def _(mo):
    temperature_06 = mo.ui.slider(
        600, 1200, value=800, step=25, label="Temperature T (K)", show_value=True
    )
    log_concentration_06 = mo.ui.slider(
        17.0, 22.0, value=20.0, step=0.25, label="log10(c0 / cm^-3)", show_value=True
    )
    length_06 = mo.ui.dropdown(
        options={
            "25 micrometers": 25.0,
            "50 micrometers": 50.0,
            "100 micrometers": 100.0,
            "250 micrometers": 250.0,
        },
        value="100 micrometers",
        label="MIEC thickness L",
    )
    log_diffusivity_06 = mo.ui.slider(
        -12.0, -6.0, value=-8.0, step=0.25, label="log10(D_delta / cm^2 s^-1)", show_value=True
    )
    log_ratio_06 = mo.ui.slider(
        -1.0, 4.0, value=2.0, step=0.25, label="log10(sigma_e / sigma_i)", show_value=True
    )
    pitt_voltage_06 = mo.ui.slider(
        5.0, 35.0, value=15.0, step=1.0, label="PITT voltage step (mV, extraction)", show_value=True
    )
    gitt_current_06 = mo.ui.slider(
        0.05, 0.60, value=0.30, step=0.025, label="GITT reduced current q", show_value=True
    )
    pulse_duration_06 = mo.ui.slider(
        0.20, 3.0, value=1.20, step=0.10, label="Pulse duration / tau_delta", show_value=True
    )
    rest_duration_06 = mo.ui.slider(
        0.25, 5.0, value=2.50, step=0.25, label="OCV duration / tau_delta", show_value=True
    )
    potential_case_06 = mo.ui.dropdown(
        options=[
            "PITT: voltage step",
            "PITT: OCV relaxation",
            "GITT: current step",
            "GITT: OCV relaxation",
        ],
        value="GITT: current step",
        label="Potential-profile case",
    )
    potential_time_06 = mo.ui.slider(
        0, 100, value=65, step=5, label="Progress through selected stage (%)", show_value=True
    )
    return (
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
    control_panel_06 = mo.vstack(
        [
            mo.md("### Classroom controls"),
            mo.hstack([temperature_06, log_concentration_06, length_06], justify="start", gap=1.2),
            mo.hstack([log_diffusivity_06, log_ratio_06], justify="start", gap=1.2),
            mo.hstack([pitt_voltage_06, gitt_current_06], justify="start", gap=1.2),
            mo.hstack([pulse_duration_06, rest_duration_06], justify="start", gap=1.2),
            mo.hstack([potential_case_06, potential_time_06], justify="start", gap=1.2),
        ]
    )
    control_panel_06
    return


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
    positions_06 = np.linspace(0.0, 1.0, 181)
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
    mo.md(
        f"""
        ### What the controls imply

        - $t_i={parameters_06["ionic_fraction"]:.4f}$ and
          $t_e={parameters_06["electronic_fraction"]:.4f}$.
        - $D^\\delta={parameters_06["chemical_diffusivity_cm2_per_s"]:.2e}$
          cm² s⁻¹ and $\\tau^\\delta=L^2/(\\pi^2D^\\delta)={_time_text}$.
        - The ideal-pair identity implies
          $\\sigma_i+\\sigma_e={parameters_06["conductivity_total_s_per_m"] / 100.0:.2e}$
          S cm⁻¹ for this $c_0$, $T$, $D^\\delta$, and conductivity ratio.
        - One unit of reduced current $q$ corresponds to
          ${0.1 * current_scale_a_per_m2_06:.3e}$ mA cm⁻².

        Changing $L$ or $D^\\delta$ stretches the physical time axis through
        $\\tau^\\delta$. Changing $\\sigma_e/\\sigma_i$ changes how the applied
        voltage is divided between concentration polarization and internal
        Ohmic loss.
        """
    )
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## 3. The three stages: pulse, interruption, relaxation

    The plots below use the **full finite-slab solution**.

    1. During a PITT voltage step, $U$ is fixed and $j(t)$ adjusts.
    2. During a GITT current step, $j$ is fixed and $U(t)$ adjusts.
    3. At interruption, $j$ becomes exactly zero. The mean composition is then
       conserved, while internal ion and electron fluxes flatten $c(x,t)$.

    Solid concentration curves are during the pulse; dashed curves are during
    the following OCV relaxation. The vertical line in each response plot marks
    current interruption.
    """)
    return


@app.cell
def _(
    current_scale_a_per_m2_06,
    gitt_result_06,
    np,
    parameters_06,
    pitt_result_06,
    plt,
    positions_06,
    pulse_times_06,
    rest_times_06,
    voltage_scale_v_06,
):
    _tau = parameters_06["tau_delta_s"]
    _x_um = positions_06 * parameters_06["length_m"] * 1.0e6
    _pulse_seconds = pulse_times_06 * _tau
    _rest_seconds = (pulse_times_06[-1] + rest_times_06) * _tau
    _profile_fractions = [0.08, 0.30, 1.0]
    _rest_fractions = [0.08, 0.35, 1.0]
    _pulse_indices = [
        int(np.argmin(np.abs(pulse_times_06 - _fraction * pulse_times_06[-1])))
        for _fraction in _profile_fractions
    ]
    _rest_indices = [
        int(np.argmin(np.abs(rest_times_06 - _fraction * rest_times_06[-1])))
        for _fraction in _rest_fractions
    ]
    _colors = ["#0a9396", "#0077b6", "#023e8a"]
    _rest_colors = ["#ee9b00", "#ca6702", "#9b2226"]

    _fig, _axes = plt.subplots(2, 2, figsize=(14.2, 9.0), constrained_layout=True)
    for _index, _color in zip(_pulse_indices, _colors):
        _axes[0, 0].plot(
            _x_um,
            pitt_result_06["pulse_profiles"][_index],
            color=_color,
            lw=2.3,
            label=f"pulse: {pulse_times_06[_index]:.2f} tau",
        )
    for _index, _color in zip(_rest_indices, _rest_colors):
        _axes[0, 0].plot(
            _x_um,
            pitt_result_06["rest_profiles"][_index],
            color=_color,
            lw=2.2,
            ls="--",
            label=f"OCV: {rest_times_06[_index]:.2f} tau after",
        )
    _axes[0, 0].set_title("PITT: concentration profiles")
    _axes[0, 0].set_xlabel("position x (micrometers)")
    _axes[0, 0].set_ylabel("c / c0")
    _axes[0, 0].legend(fontsize=8.5, ncol=2)

    _pitt_voltage_mv = 1000.0 * voltage_scale_v_06 * pitt_result_06["pulse_u"]
    _pitt_rest_voltage_mv = 1000.0 * voltage_scale_v_06 * pitt_result_06["rest_u"]
    _pitt_current_milliamp_cm2 = 0.1 * current_scale_a_per_m2_06 * pitt_result_06["pulse_q"]
    _pitt_current_rest = 0.1 * current_scale_a_per_m2_06 * pitt_result_06["rest_q"]
    _axes[0, 1].plot(_pulse_seconds, _pitt_voltage_mv, color="#9b2226", lw=2.5)
    _axes[0, 1].plot(_rest_seconds, _pitt_rest_voltage_mv, color="#9b2226", lw=2.5)
    _axes[0, 1].axvline(_pulse_seconds[-1], color="0.45", ls=":", lw=1.7)
    _axes[0, 1].set_title("PITT response: fixed voltage, relaxing current")
    _axes[0, 1].set_xlabel("time (s)")
    _axes[0, 1].set_ylabel("voltage change (mV)", color="#9b2226")
    _axes[0, 1].tick_params(axis="y", labelcolor="#9b2226")
    _current_axis = _axes[0, 1].twinx()
    _current_axis.plot(_pulse_seconds, _pitt_current_milliamp_cm2, color="#005f73", lw=2.3)
    _current_axis.plot(_rest_seconds, _pitt_current_rest, color="#005f73", lw=2.3)
    _current_axis.set_ylabel("current density (mA cm$^{-2}$)", color="#005f73")
    _current_axis.tick_params(axis="y", labelcolor="#005f73")

    for _index, _color in zip(_pulse_indices, _colors):
        _axes[1, 0].plot(
            _x_um,
            gitt_result_06["pulse_profiles"][_index],
            color=_color,
            lw=2.3,
            label=f"pulse: {pulse_times_06[_index]:.2f} tau",
        )
    for _index, _color in zip(_rest_indices, _rest_colors):
        _axes[1, 0].plot(
            _x_um,
            gitt_result_06["rest_profiles"][_index],
            color=_color,
            lw=2.2,
            ls="--",
            label=f"OCV: {rest_times_06[_index]:.2f} tau after",
        )
    _axes[1, 0].set_title("GITT: concentration profiles")
    _axes[1, 0].set_xlabel("position x (micrometers)")
    _axes[1, 0].set_ylabel("c / c0")
    _axes[1, 0].legend(fontsize=8.5, ncol=2)

    _gitt_voltage_mv = 1000.0 * voltage_scale_v_06 * gitt_result_06["pulse_u"]
    _gitt_rest_voltage_mv = 1000.0 * voltage_scale_v_06 * gitt_result_06["rest_u"]
    _gitt_current_milliamp_cm2 = 0.1 * current_scale_a_per_m2_06 * gitt_result_06["pulse_q"]
    _axes[1, 1].plot(_pulse_seconds, _gitt_voltage_mv, color="#9b2226", lw=2.5)
    _axes[1, 1].plot(_rest_seconds, _gitt_rest_voltage_mv, color="#9b2226", lw=2.5)
    _axes[1, 1].axvline(_pulse_seconds[-1], color="0.45", ls=":", lw=1.7)
    _axes[1, 1].set_title("GITT response: fixed current, relaxing voltage")
    _axes[1, 1].set_xlabel("time (s)")
    _axes[1, 1].set_ylabel("voltage change (mV)", color="#9b2226")
    _axes[1, 1].tick_params(axis="y", labelcolor="#9b2226")
    _gitt_current_axis = _axes[1, 1].twinx()
    _gitt_current_axis.plot(_pulse_seconds, _gitt_current_milliamp_cm2, color="#005f73", lw=2.3)
    _gitt_current_axis.plot(_rest_seconds, np.zeros_like(_rest_seconds), color="#005f73", lw=2.3)
    _gitt_current_axis.set_ylabel("current density (mA cm$^{-2}$)", color="#005f73")
    _gitt_current_axis.tick_params(axis="y", labelcolor="#005f73")

    for _axis in _axes.flat:
        _axis.grid(alpha=0.22)
    _fig.suptitle(
        "The same diffusion equation under voltage, current, and OCV control",
        fontsize=16,
        weight="bold",
    )
    _fig
    return


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
        _selected_states = _selected_result["rest_states"]
        _selected_q = _selected_result["rest_q"]
        _stage_name = "OCV relaxation"
    else:
        _selected_times = pulse_times_06
        _selected_states = _selected_result["pulse_states"]
        _selected_q = _selected_result["pulse_q"]
        _stage_name = "driven pulse"
    potential_index_06 = int(np.argmin(np.abs(_selected_times - _progress * _selected_times[-1])))
    selected_state_06 = _selected_states[potential_index_06]
    selected_q_06 = float(_selected_q[potential_index_06])
    selected_reduced_time_06 = float(_selected_times[potential_index_06])
    selected_stage_06 = _stage_name
    selected_potentials_06 = potential_decomposition(
        selected_state_06,
        selected_q_06,
        positions_06,
        _selected_result["modes"],
        parameters_06,
    )
    return (
        potential_index_06,
        selected_potentials_06,
        selected_q_06,
        selected_reduced_time_06,
        selected_stage_06,
        selected_state_06,
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
    _axes[0].plot(_x_um, selected_potentials_06["profile"], color="#005f73", lw=2.8)
    _axes[0].set_title("Composition")
    _axes[0].set_ylabel("c / c0")

    _axes[1].plot(
        _x_um,
        _energy_scale * selected_potentials_06["mu_i"],
        color="#0a9396",
        lw=2.2,
        label="$\\mu_i$",
    )
    _axes[1].plot(
        _x_um,
        _energy_scale * selected_potentials_06["electrical_i"],
        color="#ee9b00",
        lw=2.2,
        label="$+F\\phi$",
    )
    _axes[1].plot(
        _x_um,
        _energy_scale * selected_potentials_06["tilde_mu_i"],
        color="#9b2226",
        lw=2.8,
        label="$\\widetilde\\mu_i$",
    )
    _axes[1].set_title("Ion: chemical + electrical")
    _axes[1].set_ylabel("change (kJ mol$^{-1}$)")
    _axes[1].legend(fontsize=9)

    _axes[2].plot(
        _x_um,
        _energy_scale * selected_potentials_06["mu_e"],
        color="#0a9396",
        lw=2.2,
        label="$\\mu_e$",
    )
    _axes[2].plot(
        _x_um,
        _energy_scale * selected_potentials_06["electrical_e"],
        color="#ee9b00",
        lw=2.2,
        label="$-F\\phi$",
    )
    _axes[2].plot(
        _x_um,
        _energy_scale * selected_potentials_06["tilde_mu_e"],
        color="#9b2226",
        lw=2.8,
        label="$\\widetilde\\mu_e$",
    )
    _axes[2].set_title("Electron: chemical + electrical")
    _axes[2].set_ylabel("change (kJ mol$^{-1}$)")
    _axes[2].legend(fontsize=9)
    for _axis in _axes:
        _axis.set_xlabel("position x (micrometers)")
        _axis.grid(alpha=0.22)
    _fig.suptitle(
        f"{potential_case_06.value} | {selected_stage_06}, "
        f"t/tau = {selected_reduced_time_06:.3f}, q = {selected_q_06:.3f}",
        fontsize=14,
        weight="bold",
    )
    mo.vstack(
        [
            _fig,
            mo.md(
                r"Move the **progress** slider or switch cases. The orange curves show "
                r"the electrostatic potential in molar-energy form, $+F\phi$ for the ion "
                r"and $-F\phi$ for the electron; divide by $F$ to express $\phi$ in volts. "
                r"During OCV the electrical current is zero, but sloped electrochemical "
                r"potentials can still drive equal ion and electron fluxes."
            ),
        ]
    )
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## 4. Where the textbook approximations come from

    To isolate the standard formulas, now take the **classical one-sided
    limit**: electronic transport in the MIEC is much faster than ionic
    transport ($t_e\rightarrow1$), the perturbation is small, the slab is
    planar, and interfacial and series resistances are negligible. Define

    $$\theta=\frac{D^\delta t}{L^2}=\frac{s}{\pi^2}.$$

    ### PITT: current after a surface-composition step

    If the voltage step changes the electrolyte-side surface composition by
    $\Delta y$, the finite-slab series is

    $$
    q_{\rm PITT}=2\Delta y\sum_{m=0}^{\infty}
    \exp[-(m+\tfrac12)^2\pi^2\theta].
    $$

    Its two limits are

    $$
    \underbrace{q\simeq\frac{\Delta y}{\sqrt{\pi\theta}}}_{\theta\ll1\;\text{(Cottrell)}},
    \qquad
    \underbrace{q\simeq2\Delta y\exp(-\pi^2\theta/4)}_{\theta\gg1\;\text{(first mode)}}.
    $$

    The factor $1/4$ in the long-time exponent is not optional: it follows from
    a blocking boundary at one face and a fixed surface composition at the
    other.

    ### GITT: voltage during a current pulse

    For a small signal, voltage follows the surface concentration. The exact
    one-sided finite-slab result is

    $$
    u_{\rm GITT}=q\left[\theta+\frac13-
    \frac{2}{\pi^2}\sum_{n=1}^{\infty}\frac{e^{-n^2\pi^2\theta}}{n^2}\right],
    $$

    with

    $$
    \underbrace{u\simeq2q\sqrt{\theta/\pi}}_{\theta\ll1},
    \qquad
    \underbrace{u\simeq q(\theta+1/3)}_{\theta\gg1}.
    $$

    These approximations are comparison curves, never inputs to the full
    numerical solution above.
    """)
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
    approximation_theta_06 = approximation_times_06 / np.pi**2
    classical_concentration_step_06 = 1.0 - np.exp(-pitt_reduced_voltage_06)
    pitt_series_06 = classical_pitt_series(approximation_times_06, classical_concentration_step_06)
    pitt_short_06 = classical_concentration_step_06 / np.sqrt(np.pi * approximation_theta_06)
    pitt_long_06 = 2.0 * classical_concentration_step_06 * np.exp(-approximation_times_06 / 4.0)
    gitt_series_06 = classical_gitt_series(approximation_times_06, gitt_current_06.value)
    gitt_short_06 = 2.0 * gitt_current_06.value * np.sqrt(approximation_theta_06 / np.pi)
    gitt_long_06 = gitt_current_06.value * (approximation_theta_06 + 1.0 / 3.0)

    _fig, _axes = plt.subplots(1, 2, figsize=(13.8, 4.8), constrained_layout=True)
    _axes[0].loglog(
        approximation_theta_06,
        pitt_series_06 / classical_concentration_step_06,
        color="#005f73",
        lw=3.0,
        label="finite-slab series",
    )
    _axes[0].loglog(
        approximation_theta_06,
        pitt_short_06 / classical_concentration_step_06,
        "--",
        color="#ee9b00",
        lw=2.2,
        label="short time",
    )
    _axes[0].loglog(
        approximation_theta_06,
        pitt_long_06 / classical_concentration_step_06,
        ":",
        color="#9b2226",
        lw=2.5,
        label="long time",
    )
    _axes[0].set_title("PITT current")
    _axes[0].set_ylabel(r"$q / \Delta y$")

    _axes[1].loglog(
        approximation_theta_06,
        gitt_series_06 / gitt_current_06.value,
        color="#005f73",
        lw=3.0,
        label="finite-slab series",
    )
    _axes[1].loglog(
        approximation_theta_06,
        gitt_short_06 / gitt_current_06.value,
        "--",
        color="#ee9b00",
        lw=2.2,
        label="short time",
    )
    _axes[1].loglog(
        approximation_theta_06,
        gitt_long_06 / gitt_current_06.value,
        ":",
        color="#9b2226",
        lw=2.5,
        label="long time",
    )
    _axes[1].set_title("GITT voltage")
    _axes[1].set_ylabel("u / q")
    for _axis in _axes:
        _axis.set_xlabel(r"dimensionless time $\theta=D^\delta t/L^2$")
        _axis.grid(which="both", alpha=0.22)
        _axis.legend()
    _fig.suptitle(
        "Short and long formulas are asymptotes, not universal fits", fontsize=15, weight="bold"
    )
    _fig
    return (
        approximation_theta_06,
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
    mo.md(r"""
    ### Reading the usual experimental plots

    The asymptotes suggest simple diagnostics, provided the assumptions above
    have been tested. In dimensional form, the PITT short-time current density
    is

    $$
    j(t)\simeq F\bar c_0\Delta y
    \sqrt{\frac{D^\delta}{\pi t}},
    $$

    so $j$ versus $t^{-1/2}$ is linear and

    $$D^\delta=\pi\left[\frac{j\sqrt{t}}
    {F\bar c_0\Delta y}\right]^2.$$

    At long PITT times, a plot of $\ln|j|$ versus $t$ has slope
    $-\pi^2D^\delta/(4L^2)$.

    For a short GITT pulse of duration $t_p$, let $\Delta U_t$ be the gradual
    diffusion voltage during the pulse after removing the instantaneous Ohmic
    jump, and let $\Delta U_s$ be the equilibrium OCV change caused by that
    pulse. The planar small-signal result is

    $$
    D^\delta\simeq\frac{4L^2}{\pi t_p}
    \left(\frac{\Delta U_s}{\Delta U_t}\right)^2,
    \qquad \frac{D^\delta t_p}{L^2}\ll1.
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
    rest times, every nonuniform Fourier mode decays and the slowest one gives

    $$c(x,t_r)-\bar c\propto\exp(-t_r/\tau^\delta).$$

    The notebook carries the complete end-of-pulse profile into the OCV stage,
    so it does not need either approximation to generate the relaxation curves.
    """)
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## 5. What can PITT and GITT determine?

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
    modal_mean_error_06 = max(
        np.max(np.abs(pitt_spatial_mean_06 - pitt_result_06["pulse_states"][:, 0])),
        np.max(np.abs(gitt_spatial_mean_06 - gitt_result_06["pulse_states"][:, 0])),
    )
    gitt_mass_balance_error_06 = np.max(
        np.abs(
            gitt_result_06["pulse_states"][:, 0]
            - (1.0 - gitt_result_06["pulse_q"][0] * pulse_times_06 / np.pi**2)
        )
    )
    ocv_mean_drift_06 = max(
        np.ptp(pitt_result_06["rest_states"][:, 0]),
        np.ptp(gitt_result_06["rest_states"][:, 0]),
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
            pitt_result_06["pulse_states"][_sample_index],
            pitt_result_06["pulse_q"][_sample_index],
            positions_06,
            pitt_result_06["modes"],
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
    first_mode_relaxation_error_06 = abs(
        gitt_result_06["rest_states"][-1, 1] / gitt_result_06["rest_states"][0, 1]
        - np.exp(-rest_times_06[-1])
    )
    return (
        diffusivity_identity_error_06,
        first_mode_relaxation_error_06,
        gitt_current_control_error_06,
        gitt_long_error_06,
        gitt_mass_balance_error_06,
        gitt_rest_current_error_06,
        gitt_short_error_06,
        initial_uniform_error_06,
        minimum_concentration_06,
        modal_mean_error_06,
        ocv_mean_drift_06,
        pitt_long_error_06,
        pitt_rest_current_error_06,
        pitt_short_error_06,
        pitt_voltage_control_error_06,
        voltage_reconstruction_error_06,
    )


@app.cell
def _(
    diffusivity_identity_error_06,
    first_mode_relaxation_error_06,
    gitt_current_control_error_06,
    gitt_long_error_06,
    gitt_mass_balance_error_06,
    gitt_rest_current_error_06,
    gitt_short_error_06,
    initial_uniform_error_06,
    minimum_concentration_06,
    mo,
    modal_mean_error_06,
    ocv_mean_drift_06,
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
            initial_uniform_error_06,
            "Both experiments must begin from the same equilibrated composition.",
        ),
        (
            "Positive concentrations",
            minimum_concentration_06 > 0.0,
            minimum_concentration_06,
            "Logarithmic chemical potentials require c > 0 everywhere.",
        ),
        (
            "PITT voltage held fixed",
            pitt_voltage_control_error_06 < 2.0e-9,
            pitt_voltage_control_error_06,
            "The current must adjust so that the requested voltage step is maintained.",
        ),
        (
            "GITT current held fixed",
            gitt_current_control_error_06 < 1.0e-14,
            gitt_current_control_error_06,
            "The voltage, not the imposed current, is allowed to evolve during GITT.",
        ),
        (
            "Zero current during both OCV rests",
            max(pitt_rest_current_error_06, gitt_rest_current_error_06) < 1.0e-14,
            max(pitt_rest_current_error_06, gitt_rest_current_error_06),
            "Open circuit stops charge transfer through the terminals while "
            "internal diffusion continues.",
        ),
        (
            "Composition balance",
            max(modal_mean_error_06, gitt_mass_balance_error_06, ocv_mean_drift_06) < 2.0e-5,
            max(modal_mean_error_06, gitt_mass_balance_error_06, ocv_mean_drift_06),
            "Integrated current changes the mean during a pulse; the mean is "
            "conserved during rest.",
        ),
        (
            "Chemical-diffusivity identity",
            diffusivity_identity_error_06 < 1.0e-12,
            diffusivity_identity_error_06,
            "The conductivity and diffusivity parameters must describe the same ideal H+/e- pair.",
        ),
        (
            "Potential decomposition reconstructs voltage",
            voltage_reconstruction_error_06 < 3.0e-4,
            voltage_reconstruction_error_06,
            "Chemical and electrical contributions must add back to the measured terminal drive.",
        ),
        (
            "PITT short- and long-time limits",
            max(pitt_short_error_06, pitt_long_error_06) < 2.0e-3,
            max(pitt_short_error_06, pitt_long_error_06),
            "The Cottrell and first-mode formulas must approach the finite-slab "
            "series in their own regimes.",
        ),
        (
            "GITT short- and long-time limits",
            max(gitt_short_error_06, gitt_long_error_06) < 2.0e-3,
            max(gitt_short_error_06, gitt_long_error_06),
            "The square-root and late linear forms must approach the same finite-slab solution.",
        ),
        (
            "OCV slow mode",
            first_mode_relaxation_error_06 < 1.0e-12,
            first_mode_relaxation_error_06,
            "At long rest, the first Fourier mode sets the time constant tau_delta.",
        ),
    ]
    _rows = []
    for _name, _passed, _value, _why in _checks:
        _status = "PASS" if _passed else "CHECK"
        _rows.append(f"| {_status} | {_name} | {_value:.2e} | {_why} |")
    _table_lines = [
        "## 6. Numerical and physical sanity checks",
        "",
        "These checks are not extra decoration. Each one protects a link in the "
        "experimental interpretation.",
        "",
        "| status | check | diagnostic | why it matters |",
        "|---|---|---:|---|",
    ]
    mo.md("\n".join(_table_lines + _rows))
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## 7. Take-home picture

    $$
    \boxed{
    \text{controlled pulse}
    \rightarrow \text{selective carrier fluxes}
    \rightarrow c(x,t)
    \rightarrow \mu_i,\mu_e,\phi
    \rightarrow U(t)\text{ or }j(t)
    }
    $$

    PITT and GITT do not measure $D^\delta$ directly. They measure a transient,
    and $D^\delta$ is inferred through a model with stated geometry, boundary
    conditions, and thermodynamics. Short-time square-root behavior and
    long-time exponential or linear behavior emerge from that model; they
    should not be imposed outside their regimes.

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
            "font.size": 12,
            "axes.titlesize": 13,
            "axes.labelsize": 12,
            "legend.fontsize": 10,
            "lines.solid_capstyle": "round",
            "figure.dpi": 115,
        }
    )
    return


if __name__ == "__main__":
    app.run()
