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
    from scipy.integrate import solve_ivp

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
        }
    )

    return mo, np, plt, solve_ivp


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
def _(mo, plt):
    introduction = mo.md(r"""
    # Stoichiometry polarization in a mixed conductor

    **What happens inside a mixed conductor when both electrodes block ions but
    still pass electrons?**

    Immediately after a current is applied, ions and electrons both respond.
    The ions then meet blocking faces: they cannot leave, so the stoichiometry
    redistributes until an internal chemical-potential gradient balances their
    electrical drive.

    The model uses one ideal reaction,
    \(\mathrm{H}\rightleftharpoons\mathrm{H^+}+e^-\), with local bulk
    electroneutrality \(c_i=c_e=c(x,t)\). The initially uniform concentration is
    \(c_0\), and both blocking faces keep its spatial average constant.
    """)

    geometry_figure, geometry_axis = plt.subplots(figsize=(11.8, 3.1), dpi=120)
    geometry_axis.set_xlim(-0.08, 1.08)
    geometry_axis.set_ylim(-0.45, 0.65)
    geometry_axis.fill_between(
        [0.0, 1.0], [-0.22, -0.22], [0.30, 0.30],
        color="#DDE7E8", edgecolor="#4C7C86", linewidth=1.5,
    )
    geometry_axis.plot([0.0, 0.0], [-0.28, 0.36], color="#6D747A", lw=7)
    geometry_axis.plot([1.0, 1.0], [-0.28, 0.36], color="#6D747A", lw=7)
    geometry_axis.annotate(
        "", xy=(0.82, 0.46), xytext=(0.18, 0.46),
        arrowprops={"arrowstyle": "->", "lw": 1.6, "color": "#B8734A"},
    )
    geometry_axis.text(0.50, 0.51, "Positive conventional current", ha="center")
    geometry_axis.text(0.50, 0.02, r"MIEC: mobile H$^+$ and e$^-$", ha="center", fontsize=14)
    geometry_axis.text(0.0, -0.37, "Ion-blocking\nmetal", ha="center", va="top")
    geometry_axis.text(1.0, -0.37, "Ion-blocking\nmetal", ha="center", va="top")
    geometry_axis.text(0.50, -0.29, r"$0\leq x\leq L$", ha="center", va="top")
    geometry_axis.axis("off")
    geometry_axis.set_title("Ions remain in the slab; electrons cross both contacts")
    geometry_figure.tight_layout()
    plt.close(geometry_figure)

    mo.vstack([
        introduction,
        geometry_figure,
        mo.md(r"""
        Positive current points toward \(+x\). We define the measured voltage as
        \(U=[\widetilde{\mu}_e(L)-\widetilde{\mu}_e(0)]/F\).
        """),
    ])
    return (geometry_figure,)

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

    def mixed_conductor_parameters(
        temperature_k,
        mean_concentration_cm3,
        total_conductivity_s_per_cm,
        electronic_to_ionic_ratio,
        length_m,
    ):
        """Ideal H+/e- pair with both conductivities proportional to c."""
        temperature = _positive("temperature_k", temperature_k)
        concentration_cm3 = _positive(
            "mean_concentration_cm3", mean_concentration_cm3
        )
        sigma_total_s_per_cm = _positive(
            "total_conductivity_s_per_cm", total_conductivity_s_per_cm
        )
        ratio = _positive("electronic_to_ionic_ratio", electronic_to_ionic_ratio)
        length = _positive("length_m", length_m)

        concentration_mol_per_m3 = concentration_cm3 * 1.0e6 / AVOGADRO_PER_MOL
        sigma_total_s_per_m = sigma_total_s_per_cm * 100.0
        ionic_fraction = 1.0 / (1.0 + ratio)
        electronic_fraction = ratio / (1.0 + ratio)
        sigma_i_s_per_m = ionic_fraction * sigma_total_s_per_m
        sigma_e_s_per_m = electronic_fraction * sigma_total_s_per_m
        diffusivity_i_m2_per_s = (
            GAS_CONSTANT_J_PER_MOL_K
            * temperature
            * sigma_i_s_per_m
            / (FARADAY_C_PER_MOL**2 * concentration_mol_per_m3)
        )
        diffusivity_e_m2_per_s = (
            GAS_CONSTANT_J_PER_MOL_K
            * temperature
            * sigma_e_s_per_m
            / (FARADAY_C_PER_MOL**2 * concentration_mol_per_m3)
        )
        chemical_diffusivity_m2_per_s = (
            2.0
            * diffusivity_i_m2_per_s
            * diffusivity_e_m2_per_s
            / (diffusivity_i_m2_per_s + diffusivity_e_m2_per_s)
        )
        tau_delta_s = length**2 / (
            np.pi**2 * chemical_diffusivity_m2_per_s
        )
        return {
            "temperature_k": temperature,
            "concentration_cm3": concentration_cm3,
            "concentration_mol_per_m3": concentration_mol_per_m3,
            "length_m": length,
            "sigma_total_s_per_m": sigma_total_s_per_m,
            "sigma_i_s_per_m": sigma_i_s_per_m,
            "sigma_e_s_per_m": sigma_e_s_per_m,
            "ionic_fraction": ionic_fraction,
            "electronic_fraction": electronic_fraction,
            "diffusivity_i_m2_per_s": diffusivity_i_m2_per_s,
            "diffusivity_e_m2_per_s": diffusivity_e_m2_per_s,
            "chemical_diffusivity_m2_per_s": chemical_diffusivity_m2_per_s,
            "tau_delta_s": tau_delta_s,
        }

    def _odd_modes(mode_count):
        count = int(mode_count)
        if count < 4:
            raise ValueError("mode_count must be at least four")
        return np.arange(1, 2 * count, 2, dtype=float)

    def _profiles_from_modes(coefficients, positions, modes):
        cosine_matrix = np.cos(np.pi * np.outer(positions, modes))
        coefficient_array = np.asarray(coefficients, dtype=float)
        if coefficient_array.ndim == 1:
            return 1.0 + cosine_matrix @ coefficient_array
        return 1.0 + coefficient_array @ cosine_matrix.T

    def simulate_stoichiometry_polarization(
        drive_mode,
        current_strength_beta,
        applied_voltage_v,
        parameters,
        reduced_times,
        positions,
        mode_count=50,
    ):
        """Solve dc/dt=D_delta*d2c/dx2 with two ion-blocking electrodes."""
        time_ratios = np.asarray(reduced_times, dtype=float)
        xi = np.asarray(positions, dtype=float)
        if time_ratios[0] != 0.0 or np.any(np.diff(time_ratios) <= 0.0):
            raise ValueError("reduced_times must start at zero and increase")
        modes = _odd_modes(mode_count)
        cosine_matrix = np.cos(np.pi * np.outer(xi, modes))
        ionic_fraction = parameters["ionic_fraction"]
        electronic_fraction = parameters["electronic_fraction"]
        temperature = parameters["temperature_k"]

        if drive_mode == "current":
            beta_value = float(current_strength_beta)
            decay = np.exp(-np.outer(time_ratios, modes**2))
            mode_coefficients = (
                -4.0
                * beta_value
                / (np.pi**2 * modes**2)
                * (1.0 - decay)
            )
            beta_history = np.full(time_ratios.shape, beta_value)
            solver_success = True
            solver_message = "analytical Fourier solution"
            target_reduced_voltage = np.nan
        elif drive_mode == "voltage":
            target_reduced_voltage = (
                FARADAY_C_PER_MOL
                * float(applied_voltage_v)
                / (2.0 * GAS_CONSTANT_J_PER_MOL_K * temperature)
            )

            def beta_from_coefficients(coefficients_now):
                ratio_profile = 1.0 + cosine_matrix @ coefficients_now
                # This floor is used only during the implicit solver's internal
                # trial steps. Returned physical profiles are never clipped.
                safe_profile = np.maximum(ratio_profile, 1.0e-10)
                resistance_integral = np.trapezoid(1.0 / safe_profile, xi)
                chemical_jump = np.log(safe_profile[-1] / safe_profile[0])
                return (
                    target_reduced_voltage - ionic_fraction * chemical_jump
                ) / (electronic_fraction * resistance_integral)

            def modal_rhs(_time_ratio, coefficients_now):
                beta_now = beta_from_coefficients(coefficients_now)
                return (
                    -(modes**2) * coefficients_now
                    - 4.0 * beta_now / np.pi**2
                )

            modal_solution = solve_ivp(
                modal_rhs,
                (float(time_ratios[0]), float(time_ratios[-1])),
                np.zeros_like(modes),
                t_eval=time_ratios,
                method="BDF",
                rtol=2.0e-8,
                atol=2.0e-10,
            )
            solver_success = bool(modal_solution.success)
            solver_message = str(modal_solution.message)
            if not solver_success:
                raise RuntimeError(solver_message)
            mode_coefficients = modal_solution.y.T
            beta_history = np.array(
                [
                    beta_from_coefficients(coefficients_now)
                    for coefficients_now in mode_coefficients
                ]
            )
        else:
            raise ValueError("drive_mode must be 'current' or 'voltage'")

        concentration_ratio = _profiles_from_modes(
            mode_coefficients,
            xi,
            modes,
        )
        resistance_integral = np.trapezoid(
            1.0 / concentration_ratio,
            xi,
            axis=1,
        )
        chemical_jump = np.log(
            concentration_ratio[:, -1] / concentration_ratio[:, 0]
        )
        reduced_voltage = (
            beta_history * electronic_fraction * resistance_integral
            + ionic_fraction * chemical_jump
        )
        voltage_v = (
            2.0
            * GAS_CONSTANT_J_PER_MOL_K
            * temperature
            / FARADAY_C_PER_MOL
            * reduced_voltage
        )
        current_density_a_per_m2 = (
            beta_history
            * 2.0
            * GAS_CONSTANT_J_PER_MOL_K
            * temperature
            * parameters["sigma_e_s_per_m"]
            / (FARADAY_C_PER_MOL * parameters["length_m"])
        )
        return {
            "drive_mode": drive_mode,
            "times_over_tau_delta": time_ratios,
            "positions_over_length": xi,
            "modes": modes,
            "mode_coefficients": mode_coefficients,
            "concentration_ratio": concentration_ratio,
            "beta": beta_history,
            "resistance_integral": resistance_integral,
            "chemical_jump": chemical_jump,
            "reduced_voltage": reduced_voltage,
            "voltage_v": voltage_v,
            "current_density_a_per_m2": current_density_a_per_m2,
            "target_reduced_voltage": target_reduced_voltage,
            "solver_success": solver_success,
            "solver_message": solver_message,
        }

    def reconstruct_potential_profiles(
        coefficients,
        beta,
        current_density_a_per_m2,
        parameters,
        positions,
        modes,
    ):
        """Reconstruct chemical, electrical, and electrochemical potentials."""
        xi = np.asarray(positions, dtype=float)
        coefficients_now = np.asarray(coefficients, dtype=float)
        ratio_profile = _profiles_from_modes(coefficients_now, xi, modes)
        sine_matrix = np.sin(np.pi * np.outer(xi, modes))
        ratio_gradient_per_xi = sine_matrix @ (
            -np.pi * modes * coefficients_now
        )
        # The cosine representation converges to the non-zero boundary slope
        # only as a one-sided limit; impose the exact blocking condition there.
        ratio_gradient_per_xi[0] = float(beta)
        ratio_gradient_per_xi[-1] = float(beta)

        temperature = parameters["temperature_k"]
        length = parameters["length_m"]
        mu_i_j_per_mol = (
            GAS_CONSTANT_J_PER_MOL_K * temperature * np.log(ratio_profile)
        )
        mu_e_j_per_mol = mu_i_j_per_mol.copy()
        mu_pair_j_per_mol = mu_i_j_per_mol + mu_e_j_per_mol
        sigma_total_local = (
            parameters["sigma_total_s_per_m"] * ratio_profile
        )
        positions_m = xi * length
        resistance_segment_ohm_m2 = (
            0.5
            * (1.0 / sigma_total_local[1:] + 1.0 / sigma_total_local[:-1])
            * np.diff(positions_m)
        )
        cumulative_resistance_ohm_m2 = np.concatenate(
            ([0.0], np.cumsum(resistance_segment_ohm_m2))
        )
        tilde_mu_i_j_per_mol = (
            mu_i_j_per_mol[0]
            + parameters["electronic_fraction"]
            * (mu_pair_j_per_mol - mu_pair_j_per_mol[0])
            - FARADAY_C_PER_MOL
            * float(current_density_a_per_m2)
            * cumulative_resistance_ohm_m2
        )
        electrostatic_potential_v = (
            tilde_mu_i_j_per_mol - mu_i_j_per_mol
        ) / FARADAY_C_PER_MOL
        tilde_mu_e_j_per_mol = (
            mu_e_j_per_mol
            - FARADAY_C_PER_MOL * electrostatic_potential_v
        )
        reconstructed_voltage_v = (
            tilde_mu_e_j_per_mol[-1] - tilde_mu_e_j_per_mol[0]
        ) / FARADAY_C_PER_MOL
        return {
            "concentration_ratio": ratio_profile,
            "concentration_gradient_per_xi": ratio_gradient_per_xi,
            "mu_i_j_per_mol": mu_i_j_per_mol,
            "mu_e_j_per_mol": mu_e_j_per_mol,
            "mu_pair_j_per_mol": mu_pair_j_per_mol,
            "electrostatic_potential_v": electrostatic_potential_v,
            "electrical_i_j_per_mol": (
                FARADAY_C_PER_MOL * electrostatic_potential_v
            ),
            "electrical_e_j_per_mol": (
                -FARADAY_C_PER_MOL * electrostatic_potential_v
            ),
            "tilde_mu_i_j_per_mol": tilde_mu_i_j_per_mol,
            "tilde_mu_e_j_per_mol": tilde_mu_e_j_per_mol,
            "reconstructed_voltage_v": reconstructed_voltage_v,
        }

    return (
        AVOGADRO_PER_MOL,
        FARADAY_C_PER_MOL,
        GAS_CONSTANT_J_PER_MOL_K,
        mixed_conductor_parameters,
        reconstruct_potential_profiles,
        simulate_stoichiometry_polarization,
    )


@app.cell
def _(mo):
    derivation = mo.md(r"""
    ## 1. Derive the transient model from the two carrier fluxes

    For ideal dilute $H^+$ and $e^-$, measured relative to the uniform
    reference state,

    $$
    \mu_i=RT\ln\frac{c}{c_0},\qquad
    \mu_e=RT\ln\frac{c}{c_0},\qquad
    \mu=\mu_i+\mu_e=2RT\ln\frac{c}{c_0}.
    $$

    The electrochemical potentials and molar fluxes are

    $$
    \widetilde\mu_i=\mu_i+F\phi,\qquad
    \widetilde\mu_e=\mu_e-F\phi,
    $$

    $$
    J_i=-\frac{\sigma_i}{F^2}\frac{\partial\widetilde\mu_i}{\partial x},
    \qquad
    J_e=-\frac{\sigma_e}{F^2}\frac{\partial\widetilde\mu_e}{\partial x},
    \qquad
    j=F(J_i-J_e).
    $$

    Eliminating the internal electric field gives

    $$
    J_i=t_i\frac{j}{F}-D^\delta\frac{\partial c}{\partial x},
    \qquad
    D^\delta=\frac{2D_iD_e}{D_i+D_e}
    =\frac{2RT}{F^2c_0}
    \frac{\sigma_i\sigma_e}{\sigma_i+\sigma_e}.
    $$

    Both electrodes block ions:

    $$
    J_i(0,t)=J_i(L,t)=0.
    $$

    Conservation then produces the diffusion equation

    $$
    \frac{\partial c}{\partial t}
    =D^\delta\frac{\partial^2c}{\partial x^2},
    \qquad
    \left.\frac{\partial c}{\partial x}\right|_{0,L}
    =\frac{t_i j}{F D^\delta}.
    $$

    ### Why a factor of two matters here

    The general blocking condition can be written using
    $(\partial\mu/\partial c)^{-1}$. For the explicit ideal pair used here,
    $\mu=\mu_i+\mu_e=2RT\ln c+\text{constant}$, so the
    small-polarization boundary slope is

    $$
    \left.\frac{\partial c}{\partial x}\right|_{0,L}
    \approx \frac{Fjc_0}{2RT\sigma_e}.
    $$

    A coefficient without the $1/2$ corresponds to a different chemical-
    potential model. This is why the boundary equation must be re-derived from
    the chosen species rather than copied independently of the defect reaction.
    """)
    mo.accordion({"Model details - from carrier fluxes to diffusion": derivation})
    return


@app.cell
def _(mo):
    temperature_control = mo.ui.slider(
        start=400,
        stop=1200,
        step=25,
        value=800,
        label="Temperature, T (K)",
        show_value=True,
    )
    log_concentration_control = mo.ui.slider(
        start=17.0,
        stop=21.0,
        step=0.25,
        value=20.0,
        label="Initial-concentration exponent",
        show_value=True,
    )
    log_total_conductivity_control = mo.ui.slider(
        start=-6.0,
        stop=0.0,
        step=0.25,
        value=-3.0,
        label="Total-conductivity exponent",
        show_value=True,
    )
    log_conductivity_ratio_control = mo.ui.dropdown(
        options={
            "Ions faster (sigma_e / sigma_i = 0.01)": -2.0,
            "Balanced carriers (ratio = 1)": 0.0,
            "Electrons faster (ratio = 100)": 2.0,
            "Electrons much faster (ratio = 1000)": 3.0,
        },
        value="Electrons faster (ratio = 100)",
        label="Conductivity ratio",
    )
    length_control = mo.ui.dropdown(
        options={
            "50 micrometers": 50.0e-6,
            "100 micrometers": 100.0e-6,
            "250 micrometers": 250.0e-6,
            "500 micrometers": 500.0e-6,
            "1 millimeter": 1.0e-3,
        },
        value="100 micrometers",
        label="Sample thickness, L",
    )
    drive_mode_control = mo.ui.dropdown(
        options={
            "constant current": "current",
            "constant potential": "voltage",
        },
        value="constant current",
        label="Electrical drive",
    )
    current_strength_control = mo.ui.dropdown(
        options={
            "weak forward polarization (beta = 0.25)": 0.25,
            "moderate forward polarization (beta = 0.80)": 0.80,
            "strong forward polarization (beta = 1.40)": 1.40,
            "moderate reverse polarization (beta = -0.80)": -0.80,
        },
        value="moderate forward polarization (beta = 0.80)",
        label="Polarization strength",
    )
    applied_voltage_control = mo.ui.slider(
        start=-120,
        stop=120,
        step=5,
        value=60,
        label="Applied voltage, U (mV)",
        show_value=True,
    )
    log_time_control = mo.ui.slider(
        start=-3.0,
        stop=0.7,
        step=0.1,
        value=-0.3,
        label="Reduced-time exponent",
        show_value=True,
    )
    return (
        applied_voltage_control,
        current_strength_control,
        drive_mode_control,
        length_control,
        log_concentration_control,
        log_conductivity_ratio_control,
        log_time_control,
        log_total_conductivity_control,
        temperature_control,
    )


@app.cell
def _(
    applied_voltage_control,
    current_strength_control,
    drive_mode_control,
    length_control,
    log_concentration_control,
    log_conductivity_ratio_control,
    log_time_control,
    log_total_conductivity_control,
    mo,
    temperature_control,
):
    if drive_mode_control.value == "current":
        active_drive_control = current_strength_control
        drive_intro = mo.md(r"""
        ## 1. Apply a constant current

        The polarization preset fixes the current strength. The conductivity
        ratio decides which carrier is slower, and the time control lets the
        concentration profile grow from uniform toward its blocked steady
        state.

        \[
        \beta=\frac{FjL}{2RT\sigma_e},\qquad
        \frac{c(x,\infty)}{c_0}=1+\beta\left(\frac{x}{L}-\frac12\right).
        \]
        """)
    else:
        active_drive_control = applied_voltage_control
        drive_intro = mo.md(r"""
        ## 1. Apply a constant potential

        The selected terminal voltage is held fixed while the current relaxes.
        The conductivity ratio decides which carrier is slower, and the time
        control follows the growing stoichiometry polarization.
        """)
    core_controls = mo.hstack(
        [active_drive_control, log_conductivity_ratio_control, log_time_control],
        justify="start",
        align="center",
        wrap=True,
        gap=1.2,
    )
    advanced_controls = mo.vstack(
        [
            mo.hstack(
                [temperature_control, log_concentration_control, length_control],
                justify="start", align="center", wrap=True, gap=1.2,
            ),
            mo.hstack(
                [log_total_conductivity_control, drive_mode_control],
                justify="start", align="center", wrap=True, gap=1.2,
            ),
            mo.md(
                "Switch the electrical-drive selector to compare fixed-current "
                "and fixed-potential experiments. Only the active drive strength "
                "is shown in the core controls."
            ),
        ]
    )
    mo.vstack([
        drive_intro,
        core_controls,
        mo.accordion({"Explore further - material and drive controls": advanced_controls}),
    ])
    return (core_controls,)

@app.cell
def _(
    applied_voltage_control,
    current_strength_control,
    drive_mode_control,
    length_control,
    log_concentration_control,
    log_conductivity_ratio_control,
    log_time_control,
    log_total_conductivity_control,
    mixed_conductor_parameters,
    np,
    reconstruct_potential_profiles,
    simulate_stoichiometry_polarization,
    temperature_control,
):
    selected_temperature_k = float(temperature_control.value)
    selected_concentration_cm3 = 10.0 ** float(log_concentration_control.value)
    selected_total_conductivity_s_per_cm = 10.0 ** float(
        log_total_conductivity_control.value
    )
    selected_conductivity_ratio = 10.0 ** float(
        log_conductivity_ratio_control.value
    )
    selected_length_m = float(length_control.value)
    selected_drive_mode = str(drive_mode_control.value)
    selected_beta_setpoint = float(current_strength_control.value)
    selected_applied_voltage_v = float(applied_voltage_control.value) * 1.0e-3

    selected_parameters = mixed_conductor_parameters(
        selected_temperature_k,
        selected_concentration_cm3,
        selected_total_conductivity_s_per_cm,
        selected_conductivity_ratio,
        selected_length_m,
    )
    polarization_time_ratios = np.concatenate(
        ([0.0], np.geomspace(1.0e-4, 5.0, 180))
    )
    polarization_positions = np.linspace(0.0, 1.0, 401)
    polarization_solution = simulate_stoichiometry_polarization(
        selected_drive_mode,
        selected_beta_setpoint,
        selected_applied_voltage_v,
        selected_parameters,
        polarization_time_ratios,
        polarization_positions,
        mode_count=50,
    )
    profile_time_targets = np.array([0.02, 0.10, 0.50, 2.0])
    profile_time_indices = np.array([
        int(np.argmin(np.abs(polarization_time_ratios - target)))
        for target in profile_time_targets
    ])
    requested_time_ratio = 10.0 ** float(log_time_control.value)
    highlighted_profile_slot = int(
        np.argmin(np.abs(np.log(profile_time_targets / requested_time_ratio)))
    )
    selected_time_index = int(
        np.argmin(np.abs(polarization_time_ratios - requested_time_ratio))
    )
    selected_time_ratio = float(polarization_time_ratios[selected_time_index])
    selected_time_s = selected_time_ratio * selected_parameters["tau_delta_s"]
    selected_current_density_a_per_m2 = float(
        polarization_solution["current_density_a_per_m2"][selected_time_index]
    )
    selected_voltage_v = float(
        polarization_solution["voltage_v"][selected_time_index]
    )
    selected_potential_profiles = reconstruct_potential_profiles(
        polarization_solution["mode_coefficients"][selected_time_index],
        polarization_solution["beta"][selected_time_index],
        selected_current_density_a_per_m2,
        selected_parameters,
        polarization_positions,
        polarization_solution["modes"],
    )
    return (
        highlighted_profile_slot,
        polarization_positions,
        polarization_solution,
        polarization_time_ratios,
        profile_time_indices,
        profile_time_targets,
        requested_time_ratio,
        selected_applied_voltage_v,
        selected_beta_setpoint,
        selected_concentration_cm3,
        selected_conductivity_ratio,
        selected_current_density_a_per_m2,
        selected_drive_mode,
        selected_length_m,
        selected_parameters,
        selected_potential_profiles,
        selected_temperature_k,
        selected_time_index,
        selected_time_ratio,
        selected_time_s,
        selected_total_conductivity_s_per_cm,
        selected_voltage_v,
    )


@app.cell
def _(
    mo,
    selected_current_density_a_per_m2,
    selected_drive_mode,
    selected_parameters,
    selected_beta_setpoint,
    selected_time_ratio,
    selected_time_s,
    selected_voltage_v,
):
    selected_tau_minutes = selected_parameters["tau_delta_s"] / 60.0
    if selected_drive_mode == "current":
        _drive_details = (
            f"The selected preset has $\\beta={selected_beta_setpoint:.2f}$ and "
            f"therefore $j={selected_current_density_a_per_m2 * 0.1:.3g}$ mA/cm²."
        )
    else:
        _drive_details = (
            "The constant-potential extension is active; the classroom current "
            "preset is retained but is not imposed."
        )
    state_details = mo.md(
        rf"""
        ### What the selected material parameters imply

        $$
        t_i=\frac{{\sigma_i}}{{\sigma_i+\sigma_e}}
        =\mathbf{{{selected_parameters['ionic_fraction']:.3f}}},\qquad
        t_e=\mathbf{{{selected_parameters['electronic_fraction']:.3f}}}.
        $$

        $$
        D^\delta=\mathbf{{{selected_parameters['chemical_diffusivity_m2_per_s'] * 1.0e4:.2e}}}
        \ \mathrm{{cm^2\,s^{{-1}}}},\qquad
        \tau^\delta=\mathbf{{{selected_tau_minutes:.2f}}}\ \mathrm{{min}}.
        $$

        At the selected $t/\tau^\delta={selected_time_ratio:.3g}$
        ($t={selected_time_s:.3g}$ s), the model gives
        $U=\mathbf{{{selected_voltage_v * 1.0e3:.2f}}}$ mV and
        $j=\mathbf{{{selected_current_density_a_per_m2 * 0.1:.3g}}}$ mA/cm².
        The active drive is **{selected_drive_mode.replace('_', ' ')}**.
        {_drive_details}

        The slow carrier limits $D^\delta$: making one conductivity enormous
        cannot compensate for making the other one vanishingly small.
        """
    )
    mo.accordion({"Model details - selected material state": state_details})
    return


@app.cell
def _(
    mixed_conductor_parameters,
    np,
    plt,
    selected_conductivity_ratio,
    selected_concentration_cm3,
    selected_length_m,
    selected_temperature_k,
    selected_total_conductivity_s_per_cm,
):
    ratio_sweep = np.logspace(-3.0, 3.0, 300)
    sweep_parameters = [
        mixed_conductor_parameters(
            selected_temperature_k,
            selected_concentration_cm3,
            selected_total_conductivity_s_per_cm,
            ratio_value,
            selected_length_m,
        )
        for ratio_value in ratio_sweep
    ]
    ionic_fraction_sweep = np.array(
        [item["ionic_fraction"] for item in sweep_parameters]
    )
    electronic_fraction_sweep = np.array(
        [item["electronic_fraction"] for item in sweep_parameters]
    )
    chemical_diffusivity_sweep_cm2 = np.array(
        [item["chemical_diffusivity_m2_per_s"] for item in sweep_parameters]
    ) * 1.0e4
    tau_sweep_s = np.array([item["tau_delta_s"] for item in sweep_parameters])

    ratio_figure, ratio_axes = plt.subplots(1, 2, figsize=(13.2, 4.6), dpi=120)
    ratio_axes[0].semilogx(
        ratio_sweep,
        ionic_fraction_sweep,
        color="#B65C4A",
        lw=1.9,
        label=r"ionic $t_i$",
    )
    ratio_axes[0].semilogx(
        ratio_sweep,
        electronic_fraction_sweep,
        color="#4C7C86",
        lw=1.9,
        label=r"electronic $t_e$",
    )
    ratio_axes[0].axvline(
        selected_conductivity_ratio,
        color="#C49345",
        lw=1.5,
        ls="--",
        label="selected ratio",
    )
    ratio_axes[0].set(
        xlabel=r"Conductivity ratio, $\sigma_e/\sigma_i$",
        ylabel="Conductivity fraction",
        title="The conductivity ratio partitions the current",
        ylim=(-0.03, 1.03),
    )
    ratio_axes[0].grid(alpha=0.22)
    ratio_axes[0].legend(frameon=False)

    ratio_axes[1].loglog(
        ratio_sweep,
        chemical_diffusivity_sweep_cm2,
        color="#4C7C86",
        lw=1.9,
        label=r"$D^\delta$",
    )
    ratio_time_axis = ratio_axes[1].twinx()
    ratio_time_axis.loglog(
        ratio_sweep,
        tau_sweep_s,
        color="#B65C4A",
        lw=1.8,
        ls="--",
        label=r"$\tau^\delta$",
    )
    ratio_axes[1].axvline(
        selected_conductivity_ratio,
        color="#C49345",
        lw=1.5,
        ls="--",
    )
    ratio_axes[1].set(
        xlabel=r"Conductivity ratio, $\sigma_e/\sigma_i$",
        ylabel=r"$D^\delta$ (cm$^2$ s$^{-1}$)",
        title="The slower carrier sets the polarization time",
    )
    ratio_time_axis.set_ylabel(r"$\tau^\delta$ (s)", color="#B65C4A")
    ratio_axes[1].grid(which="both", alpha=0.22)
    ratio_lines = ratio_axes[1].get_lines()[:1] + ratio_time_axis.get_lines()[:1]
    ratio_axes[1].legend(
        ratio_lines,
        [line.get_label() for line in ratio_lines],
        frameon=False,
        loc="best",
    )
    ratio_figure.tight_layout()
    plt.close(ratio_figure)
    mo.accordion({
        "Explore further - conductivity partition and time scale": mo.vstack([
            ratio_figure,
            mo.md("Transference numbers partition the immediate current, while the slower carrier controls the chemical diffusivity and polarization time."),
        ])
    })
    return (ratio_figure,)


@app.cell
def _(
    highlighted_profile_slot,
    mo,
    np,
    plt,
    polarization_positions,
    polarization_solution,
    polarization_time_ratios,
    profile_time_indices,
    profile_time_targets,
    selected_drive_mode,
    selected_parameters,
    selected_time_index,
):
    profile_figure, profile_axis = plt.subplots(figsize=(9.2, 5.2), dpi=120)
    profile_colors = ("#B8CCD0", "#91AEB4", "#6D9098", "#4D737C")
    profile_styles = ("-", "--", "-.", ":")
    for slot, (target, index, color, style) in enumerate(
        zip(profile_time_targets, profile_time_indices, profile_colors, profile_styles)
    ):
        highlighted = slot == highlighted_profile_slot
        profile_axis.plot(
            polarization_positions,
            polarization_solution["concentration_ratio"][index],
            color=color,
            ls=style,
            lw=2.8 if highlighted else 1.7,
            label=rf"$t/\tau^\delta={target:.2f}$"
            + (" (nearest selection)" if highlighted else ""),
            zorder=4 if highlighted else 2,
        )
    profile_axis.plot(
        polarization_positions,
        polarization_solution["concentration_ratio"][-1],
        color="#A56F55",
        lw=2.0,
        label="Steady state",
        zorder=3,
    )
    profile_axis.axvline(0.0, color="#7F878C", lw=1.0, alpha=0.55)
    profile_axis.axvline(1.0, color="#7F878C", lw=1.0, alpha=0.55)
    profile_axis.set(
        xlabel=r"Position, $x$",
        ylabel=r"Stoichiometry, $c/c_0$",
        title="The concentration gradient grows toward the blocked steady state",
        xticks=[0.0, 0.5, 1.0],
        xticklabels=[r"$0$", r"$L/2$", r"$L$"],
    )
    profile_axis.text(
        0.0, -0.19, "Ion-blocking electrode", transform=profile_axis.get_xaxis_transform(),
        ha="left", va="top", color="#526173", clip_on=False,
    )
    profile_axis.text(
        1.0, -0.19, "Ion-blocking electrode", transform=profile_axis.get_xaxis_transform(),
        ha="right", va="top", color="#526173", clip_on=False,
    )
    profile_axis.grid(alpha=0.22)
    profile_axis.legend(frameon=False, ncols=2, loc="best")
    profile_figure.subplots_adjust(bottom=0.25)
    plt.close(profile_figure)

    physical_time_s = polarization_time_ratios * selected_parameters["tau_delta_s"]
    plotted_time_s = np.maximum(physical_time_s, physical_time_s[1] * 0.2)
    selected_plot_time_s = plotted_time_s[selected_time_index]
    response_figure, response_axis = plt.subplots(figsize=(9.2, 4.5), dpi=120)
    if selected_drive_mode == "current":
        response = 1.0e3 * polarization_solution["voltage_v"]
        response_axis.semilogx(plotted_time_s, response, color="#4C7C86", lw=1.9)
        response_axis.set_ylabel(r"Measured voltage, $U$ (mV)")
        response_axis.set_title("Voltage rises at fixed current")
        response_sentence = (
            "The imposed current is fixed. As the composition gradient grows, "
            "chemical polarization adds to the measured voltage."
        )
    else:
        response = 0.1 * polarization_solution["current_density_a_per_m2"]
        response_axis.semilogx(plotted_time_s, response, color="#B8734A", lw=1.9)
        response_axis.set_ylabel(r"Current density, $j$ (mA cm$^{-2}$)")
        response_axis.set_title("Current relaxes at fixed voltage")
        response_sentence = (
            "The imposed voltage is fixed. Chemical polarization takes over "
            "part of the drive, so the current relaxes."
        )
    response_axis.scatter(
        [selected_plot_time_s], [response[selected_time_index]], s=65,
        color="#C49345", edgecolor="#40464D", zorder=4,
    )
    response_axis.set_xlabel("Time (s)")
    response_axis.grid(alpha=0.22)
    response_figure.tight_layout()
    plt.close(response_figure)

    mo.vstack([
        mo.md(r"""
        ## 2. Watch the stoichiometry and electrical response evolve

        The total amount of ion is unchanged: depletion at one side equals
        accumulation at the other. The selected time control highlights the
        nearest representative profile; it does not add another curve.
        """),
        profile_figure,
        mo.md(
            "Every transient profile has the same imposed ion-blocking boundary "
            "condition. Their growing separation shows the approach to steady state."
        ),
        response_figure,
        mo.md(response_sentence),
    ])
    return profile_figure, response_figure

@app.cell
def _(
    FARADAY_C_PER_MOL,
    mo,
    plt,
    polarization_positions,
    selected_potential_profiles,
    selected_time_ratio,
):
    potential_figure, (ion_axis, electron_axis) = plt.subplots(
        1, 2, figsize=(13.2, 4.7), dpi=120
    )
    ion_axis.plot(
        polarization_positions,
        selected_potential_profiles["mu_i_j_per_mol"] / 1.0e3,
        color="#B8734A", lw=1.6, label=r"Chemical, $\mu_i$",
    )
    ion_axis.plot(
        polarization_positions,
        selected_potential_profiles["electrical_i_j_per_mol"] / 1.0e3,
        color="#4C7C86", lw=1.4, ls="--", label=r"Electrical, $+F\phi$",
    )
    ion_axis.plot(
        polarization_positions,
        selected_potential_profiles["tilde_mu_i_j_per_mol"] / 1.0e3,
        color="#40464D", lw=1.6, ls="-.", label=r"Sum, $\widetilde\mu_i$",
    )
    ion_axis.set(
        xlabel=r"Position, $x/L$",
        ylabel="Molar energy (kJ mol$^{-1}$)",
        title="Ion potential balance",
    )
    ion_axis.grid(alpha=0.22)
    ion_axis.legend(frameon=False)

    electron_axis.plot(
        polarization_positions,
        selected_potential_profiles["mu_e_j_per_mol"] / 1.0e3,
        color="#B8734A", lw=1.6, label=r"Chemical, $\mu_e$",
    )
    electron_axis.plot(
        polarization_positions,
        selected_potential_profiles["electrical_e_j_per_mol"] / 1.0e3,
        color="#4C7C86", lw=1.4, ls="--", label=r"Electrical, $-F\phi$",
    )
    electron_axis.plot(
        polarization_positions,
        selected_potential_profiles["tilde_mu_e_j_per_mol"] / 1.0e3,
        color="#40464D", lw=1.6, ls="-.", label=r"Sum, $\widetilde\mu_e$",
    )
    electron_axis.set(
        xlabel=r"Position, $x/L$",
        ylabel="Molar energy (kJ mol$^{-1}$)",
        title="Electron potential balance",
    )
    electron_axis.grid(alpha=0.22)
    electron_axis.legend(frameon=False)
    potential_figure.tight_layout()
    plt.close(potential_figure)

    reconstructed_voltage_mv = (
        selected_potential_profiles["reconstructed_voltage_v"] * 1.0e3
    )
    details = mo.vstack([
        mo.md(rf"""
        ### Why the measured voltage contains chemical polarization

        The selected profile at \(t/\tau^\delta={selected_time_ratio:.2g}\)
        fixes both carrier chemical potentials. The internal electrostatic
        potential follows from the same flux equations. Their sums reconstruct
        the measured voltage,
        \([\widetilde\mu_e(L)-\widetilde\mu_e(0)]/F
        ={reconstructed_voltage_mv:.2f}\) mV.
        """),
        potential_figure,
    ])
    mo.accordion({"Explore further - chemical and electrical potentials": details})
    return (potential_figure,)

@app.cell
def _(
    FARADAY_C_PER_MOL,
    GAS_CONSTANT_J_PER_MOL_K,
    mo,
    polarization_solution,
    profile_time_indices,
    selected_drive_mode,
    selected_parameters,
    selected_time_index,
):
    _scale_v = (
        2.0
        * GAS_CONSTANT_J_PER_MOL_K
        * selected_parameters["temperature_k"]
        / FARADAY_C_PER_MOL
    )
    _beta_now = float(polarization_solution["beta"][selected_time_index])
    _ohmic_now_v = (
        _scale_v
        * _beta_now
        * selected_parameters["electronic_fraction"]
        * float(polarization_solution["resistance_integral"][selected_time_index])
    )
    _chemical_now_v = (
        _scale_v
        * selected_parameters["ionic_fraction"]
        * float(polarization_solution["chemical_jump"][selected_time_index])
    )
    _chemical_final_v = (
        _scale_v
        * selected_parameters["ionic_fraction"]
        * float(polarization_solution["chemical_jump"][-1])
    )
    _initial_voltage_mv = 1.0e3 * float(polarization_solution["voltage_v"][0])
    _initial_current_ma_cm2 = 0.1 * float(
        polarization_solution["current_density_a_per_m2"][0]
    )
    _imposed = "current density $j$" if selected_drive_mode == "current" else "two-terminal voltage $U$"
    _measured = "voltage $U(t)$" if selected_drive_mode == "current" else "current density $j(t)$"
    decomposition = mo.md(
        rf"""
        ### Measurement decomposition

        | part of the experiment | quantity in this model |
        |---|---|
        | imposed | {_imposed} |
        | measured | {_measured} |
        | immediate Ohmic response | $U(0^+)={_initial_voltage_mv:.3g}$ mV, $j(0^+)={_initial_current_ma_cm2:.3g}$ mA cm$^{{-2}}$ |
        | selected-time transport contribution | ${1.0e3 * _ohmic_now_v:.3g}$ mV |
        | selected-time chemical-polarization contribution | ${1.0e3 * _chemical_now_v:.3g}$ mV |
        | late chemical / Nernstian contribution | ${1.0e3 * _chemical_final_v:.3g}$ mV |

        The displayed contributions add to the same $U(t)$ calculated from the
        carrier electrochemical potentials. The chemical term begins at zero
        and grows only after stoichiometry redistributes.

        **Connection to chemical capacitance.** For molar concentration $c$ and
        active volume $V_{{\rm act}}=SL$,

        $$
        C_{{\rm chem}}=z^2F^2V_{{\rm act}}
        \left(\frac{{\partial c}}{{\partial\mu_{{\rm neutral}}}}\right),
        \qquad c_{{\rm chem}}=\frac{{C_{{\rm chem}}}}{{L}}.
        $$

        Chemical capacitance is differential storage of **neutral composition**.
        It scales with active volume; Module 07 uses the distributed quantity
        $c_{{\rm chem}}$ in F m$^{{-1}}$. Here $z=1$, $c$ is mol m$^{{-3}}$, and
        $\mu_{{\rm neutral}}$ is molar (J mol$^{{-1}}$).
        """
    )
    mo.accordion({"Model details - measurement decomposition": decomposition})
    return


@app.cell
def _(
    FARADAY_C_PER_MOL,
    GAS_CONSTANT_J_PER_MOL_K,
    np,
    polarization_positions,
    polarization_solution,
    reconstruct_potential_profiles,
    selected_applied_voltage_v,
    selected_beta_setpoint,
    selected_drive_mode,
    selected_parameters,
    selected_potential_profiles,
    selected_temperature_k,
    selected_time_index,
    selected_voltage_v,
):
    concentration_history = polarization_solution["concentration_ratio"]
    mass_residual = float(
        np.max(
            np.abs(
                np.trapezoid(
                    concentration_history - 1.0,
                    polarization_positions,
                    axis=1,
                )
            )
        )
    )
    initial_uniform_residual = float(
        np.max(np.abs(concentration_history[0] - 1.0))
    )
    minimum_concentration_ratio = float(np.min(concentration_history))

    current_history = polarization_solution["current_density_a_per_m2"]
    voltage_history = polarization_solution["voltage_v"]
    if selected_drive_mode == "current":
        drive_residual = float(
            np.max(np.abs(current_history - current_history[0]))
            / max(abs(current_history[0]), 1.0e-300)
        )
        expected_steady_beta = selected_beta_setpoint
    else:
        drive_residual = float(
            np.max(np.abs(voltage_history - selected_applied_voltage_v))
            / max(abs(selected_applied_voltage_v), 1.0e-300)
        )
        expected_steady_beta = 2.0 * np.tanh(
            FARADAY_C_PER_MOL
            * selected_applied_voltage_v
            / (4.0 * GAS_CONSTANT_J_PER_MOL_K * selected_temperature_k)
        )

    steady_profile_expected = 1.0 + expected_steady_beta * (
        polarization_positions - 0.5
    )
    steady_profile_residual = float(
        np.max(np.abs(concentration_history[-1] - steady_profile_expected))
        / max(abs(expected_steady_beta), 1.0)
    )
    steady_beta_residual = float(
        abs(polarization_solution["beta"][-1] - expected_steady_beta)
        / max(abs(expected_steady_beta), 1.0)
    )

    initial_ohmic_voltage = (
        current_history[0]
        * selected_parameters["length_m"]
        / selected_parameters["sigma_total_s_per_m"]
    )
    initial_ohmic_residual = float(
        abs(voltage_history[0] - initial_ohmic_voltage)
        / max(abs(voltage_history[0]), 1.0e-300)
    )
    blocking_flux_scale = max(
        abs(
            selected_parameters["ionic_fraction"]
            * current_history[selected_time_index]
            / FARADAY_C_PER_MOL
        ),
        1.0e-300,
    )
    blocking_flux_residual = float(
        abs(
            selected_parameters["ionic_fraction"]
            * current_history[selected_time_index]
            / FARADAY_C_PER_MOL
            - selected_parameters["chemical_diffusivity_m2_per_s"]
            * selected_parameters["concentration_mol_per_m3"]
            / selected_parameters["length_m"]
            * polarization_solution["beta"][selected_time_index]
        )
        / blocking_flux_scale
    )
    displayed_boundary_residual = float(max(
        abs(
            selected_parameters["ionic_fraction"]
            * current_history[index]
            / FARADAY_C_PER_MOL
            - selected_parameters["chemical_diffusivity_m2_per_s"]
            * selected_parameters["concentration_mol_per_m3"]
            / selected_parameters["length_m"]
            * polarization_solution["beta"][index]
        )
        / max(
            abs(selected_parameters["ionic_fraction"] * current_history[index]
                / FARADAY_C_PER_MOL),
            1.0e-300,
        )
        for index in profile_time_indices
    ))

    potential_decomposition_residual = float(
        max(
            np.max(
                np.abs(
                    selected_potential_profiles["tilde_mu_i_j_per_mol"]
                    - selected_potential_profiles["mu_i_j_per_mol"]
                    - selected_potential_profiles["electrical_i_j_per_mol"]
                )
            ),
            np.max(
                np.abs(
                    selected_potential_profiles["tilde_mu_e_j_per_mol"]
                    - selected_potential_profiles["mu_e_j_per_mol"]
                    - selected_potential_profiles["electrical_e_j_per_mol"]
                )
            ),
        )
    )
    reconstructed_voltage_residual = float(
        abs(
            selected_potential_profiles["reconstructed_voltage_v"]
            - selected_voltage_v
        )
    )
    chemical_diffusivity_identity = (
        2.0
        * selected_parameters["diffusivity_i_m2_per_s"]
        * selected_parameters["diffusivity_e_m2_per_s"]
        / (
            selected_parameters["diffusivity_i_m2_per_s"]
            + selected_parameters["diffusivity_e_m2_per_s"]
        )
    )
    diffusivity_residual = float(
        abs(
            selected_parameters["chemical_diffusivity_m2_per_s"]
            - chemical_diffusivity_identity
        )
        / selected_parameters["chemical_diffusivity_m2_per_s"]
    )

    final_potential_profiles = reconstruct_potential_profiles(
        polarization_solution["mode_coefficients"][-1],
        polarization_solution["beta"][-1],
        current_history[-1],
        selected_parameters,
        polarization_positions,
        polarization_solution["modes"],
    )
    final_ion_tilde_range = float(
        np.ptp(final_potential_profiles["tilde_mu_i_j_per_mol"])
    )
    final_chemical_scale = max(
        float(np.ptp(final_potential_profiles["mu_pair_j_per_mol"])),
        1.0,
    )
    final_ion_flatness = final_ion_tilde_range / final_chemical_scale

    module05_validation = {
        "solver_pass": bool(polarization_solution["solver_success"]),
        "initial_uniform_residual": initial_uniform_residual,
        "initial_uniform_pass": initial_uniform_residual < 1.0e-14,
        "mass_residual": mass_residual,
        "mass_pass": mass_residual < 2.0e-13,
        "minimum_concentration_ratio": minimum_concentration_ratio,
        "positivity_pass": minimum_concentration_ratio > 0.0,
        "drive_residual": drive_residual,
        "drive_pass": drive_residual < 2.0e-9,
        "initial_ohmic_residual": initial_ohmic_residual,
        "initial_ohmic_pass": initial_ohmic_residual < 2.0e-13,
        "blocking_flux_residual": blocking_flux_residual,
        "blocking_flux_pass": blocking_flux_residual < 2.0e-13,
        "displayed_boundary_residual": displayed_boundary_residual,
        "displayed_boundary_pass": displayed_boundary_residual < 2.0e-13,
        "diffusivity_residual": diffusivity_residual,
        "diffusivity_pass": diffusivity_residual < 2.0e-14,
        "potential_decomposition_residual": potential_decomposition_residual,
        "potential_decomposition_pass": potential_decomposition_residual < 1.0e-10,
        "reconstructed_voltage_residual": reconstructed_voltage_residual,
        "reconstructed_voltage_pass": reconstructed_voltage_residual < 2.0e-5,
        "steady_profile_residual": steady_profile_residual,
        "steady_beta_residual": steady_beta_residual,
        "steady_limit_pass": max(
            steady_profile_residual,
            steady_beta_residual,
        ) < 8.0e-3,
        "final_ion_flatness": final_ion_flatness,
        "final_ion_flatness_pass": final_ion_flatness < 1.2e-2,
    }
    return (module05_validation,)


@app.cell
def _(mo, module05_validation):
    def _check_status(passed):
        return "PASS" if passed else "CHECK"

    _checks = mo.md(
        rf"""
        ## Physical consistency checks

        | status | physical statement | why it matters |
        |---:|---|---|
        | {_check_status(module05_validation['solver_pass'] and module05_validation['initial_uniform_pass'] and module05_validation['mass_pass'] and module05_validation['positivity_pass'])} | the sample begins uniform, remains positive, and conserves its total ion content | blocking electrodes redistribute stoichiometry without adding or removing ions |
        | {_check_status(module05_validation['drive_pass'])} | the selected current or voltage is held constant | the two experimental controls remain distinct |
        | {_check_status(module05_validation['initial_ohmic_pass'] and module05_validation['blocking_flux_pass'] and module05_validation['displayed_boundary_pass'])} | the initial response uses \(\sigma_i+\sigma_e\), and every displayed profile obeys zero ionic flux at both electrodes | the bulk and contact conditions match the stated experiment |
        | {_check_status(module05_validation['diffusivity_pass'])} | \(D^\delta=2D_iD_e/(D_i+D_e)\) | ion and electron motion combine into one chemical diffusivity |
        | {_check_status(module05_validation['potential_decomposition_pass'] and module05_validation['reconstructed_voltage_pass'])} | chemical and electrical contributions reproduce both electrochemical potentials and the terminal voltage | all plotted potentials belong to the same physical state |
        | {_check_status(module05_validation['steady_limit_pass'] and module05_validation['final_ion_flatness_pass'])} | the late profile reaches the expected steady state and \(\widetilde\mu_i\) becomes flat | the long-time limit agrees with ion blocking |

        These checks protect the physical interpretation without adding new
        equations to the model.
        """
    )
    mo.accordion({"Physical consistency checks": _checks})
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## What to carry forward

    1. **One neutral stoichiometry variable couples two carriers.** Local bulk
       electroneutrality gives $c_i=c_e=c$, while the slower carrier limits
       $D^\delta=2D_iD_e/(D_i+D_e)$ in this ideal pair model.
    2. **Blocking changes profiles, not total ion content.** Ion-blocking faces
       force a gradient to grow while conserving the spatial average. Under
       constant current the voltage rises; under constant potential the current
       relaxes.
    3. **The measured voltage has several contributions.** Chemical and
       electrical terms add to each carrier's electrochemical potential. At
       steady state $\widetilde\mu_i$ is flat and the remaining terminal current
       is electronic.

    **Model boundary.** This first stoichiometry-polarization module is planar
    and one-dimensional. It uses an ideal $H/H^+/e^-$ pair, local bulk
    electroneutrality, conductivities proportional to $c$, and two perfectly
    ion-blocking but electronically reversible electrodes. It does not resolve
    interfacial space charge, electrode kinetics, non-ideal concentration
    effects, concentration-dependent mobilities, electron-blocking electrodes,
    or Hebb–Wagner boundary conditions. Those require new boundary or
    constitutive equations and should be added as later cases rather than hidden
    inside this model.
    """)
    return


@app.cell
def _(mo):
    mo.md(r"""
    **Continue:** [Module 06 — From Coulometric Titration to PITT and GITT](https://qiyanglu.github.io/Solid-State-Ionics-Interactive/06-pitt-gitt/)
    """)
    return



if __name__ == "__main__":
    app.run()
