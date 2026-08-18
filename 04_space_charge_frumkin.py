import marimo

__generated_with = "0.24.0"
app = marimo.App(width="full")


@app.cell
def _():
    import marimo as mo
    import matplotlib.pyplot as plt
    import numpy as np

    return mo, np, plt


@app.cell
def _(mo):
    mo.md(r"""
    # Space-Charge Layers and the Frumkin Effect

    **Guiding question.** How does an interface redistribute charged defects,
    and why does that redistribution change both capacitance and reaction rate?

    A surface, grain boundary, or phase boundary can have a different defect
    formation energy from the bulk. Defects may therefore segregate into a thin
    **core**. If that core acquires net charge, nearby mobile defects rearrange
    to screen it: this surrounding region is the **space-charge layer**. We do
    not model the atomic structure of the core; its potential \(\phi_0\) is the
    boundary condition for the continuum region.

    This reader follows the lecture's one-dimensional path:

    \[
    \widetilde\mu_i=\mu_i+z_i e\phi=\text{constant}
    \rightarrow c_i(x)
    \rightarrow \rho(x)
    \rightarrow \phi(x)
    \rightarrow Q_{\rm core}
    \rightarrow C_d
    \rightarrow \text{GCS}
    \rightarrow \text{Frumkin effect}.
    \]

    We compare the two limits used in the space-charge lecture:

    - **Gouy–Chapman:** both positive and negative defects are mobile.
    - **Mott–Schottky:** the majority negative dopants are frozen, while the
      positive defects remain mobile.

    | symbol | meaning in this notebook |
    |---|---|
    | \(x\) | distance from the core; \(x=0\) at the interface |
    | \(c_{i,\infty}\) | bulk number concentration of each charge carrier |
    | \(\phi_0\), \(\phi_1\), \(\phi_\infty\) | core, reaction-plane, and bulk potentials |
    | \(Q_{\rm core}\) | core charge **per unit area** |
    | \(\lambda_D\), \(\lambda\) | Gouy–Chapman Debye length and Mott–Schottky depletion width |
    | \(C_d\), \(C_s\) | diffuse-layer and Stern-layer differential capacitances per area |

    The bulk is the reference, \(\phi_\infty=0\), and the main profile sections
    use a positively charged core, as in the slides. Inputs are shown in
    **K, V, cm⁻³, nm, and μF/cm²**; calculations use SI internally. We use
    \(k_B,e\) for energies per defect and \(R,F\) for the equivalent molar
    kinetic equations.
    """)
    return


@app.cell
def _(np):
    KB_J_PER_K = 1.380649e-23
    E_CHARGE_C = 1.602176634e-19
    EPSILON_0_F_PER_M = 8.8541878128e-12
    FARADAY_C_PER_MOL = 96485.33212
    GAS_CONSTANT_J_PER_MOL_K = 8.314462618

    def _positive(name, value):
        number = float(value)
        if not np.isfinite(number) or number <= 0.0:
            raise ValueError(f"{name} must be positive and finite")
        return number

    def debye_length_m(temperature_k, concentration_cm3, epsilon_r, charge_magnitude):
        """Debye length for two mobile species with charges +ze and -ze."""
        temperature = _positive("temperature_k", temperature_k)
        concentration_m3 = _positive("concentration_cm3", concentration_cm3) * 1.0e6
        relative_permittivity = _positive("epsilon_r", epsilon_r)
        z_value = _positive("charge_magnitude", charge_magnitude)
        permittivity = EPSILON_0_F_PER_M * relative_permittivity
        return np.sqrt(
            permittivity * KB_J_PER_K * temperature
            / (2.0 * z_value**2 * E_CHARGE_C**2 * concentration_m3)
        )

    def gc_surface_charge_c_per_m2(
        surface_potential_v,
        temperature_k,
        concentration_cm3,
        epsilon_r,
        charge_magnitude,
    ):
        """Charge per core area required by the planar Gouy-Chapman layer."""
        temperature = _positive("temperature_k", temperature_k)
        relative_permittivity = _positive("epsilon_r", epsilon_r)
        z_value = _positive("charge_magnitude", charge_magnitude)
        lambda_d = debye_length_m(
            temperature,
            concentration_cm3,
            relative_permittivity,
            z_value,
        )
        permittivity = EPSILON_0_F_PER_M * relative_permittivity
        thermal_voltage = KB_J_PER_K * temperature / (z_value * E_CHARGE_C)
        return (
            2.0
            * permittivity
            * thermal_voltage
            / lambda_d
            * np.sinh(float(surface_potential_v) / (2.0 * thermal_voltage))
        )

    def gc_differential_capacitance_f_per_m2(
        surface_potential_v,
        temperature_k,
        concentration_cm3,
        epsilon_r,
        charge_magnitude,
    ):
        """dQ_core/dphi_0 for the planar Gouy-Chapman solution."""
        temperature = _positive("temperature_k", temperature_k)
        relative_permittivity = _positive("epsilon_r", epsilon_r)
        z_value = _positive("charge_magnitude", charge_magnitude)
        lambda_d = debye_length_m(
            temperature,
            concentration_cm3,
            relative_permittivity,
            z_value,
        )
        permittivity = EPSILON_0_F_PER_M * relative_permittivity
        thermal_voltage = KB_J_PER_K * temperature / (z_value * E_CHARGE_C)
        return permittivity / lambda_d * np.cosh(
            np.asarray(surface_potential_v, dtype=float) / (2.0 * thermal_voltage)
        )

    def gouy_chapman_profile(
        distance_m,
        surface_potential_v,
        temperature_k,
        concentration_cm3,
        epsilon_r,
        charge_magnitude,
    ):
        """Exact planar Poisson-Boltzmann profile for a symmetric defect pair."""
        distance = np.asarray(distance_m, dtype=float)
        temperature = _positive("temperature_k", temperature_k)
        concentration = _positive("concentration_cm3", concentration_cm3)
        z_value = _positive("charge_magnitude", charge_magnitude)
        lambda_d = debye_length_m(
            temperature,
            concentration,
            epsilon_r,
            z_value,
        )
        thermal_voltage = KB_J_PER_K * temperature / (z_value * E_CHARGE_C)
        gamma = np.tanh(float(surface_potential_v) / (4.0 * thermal_voltage))
        atanh_argument = np.clip(
            gamma * np.exp(-distance / lambda_d),
            -0.999999999999,
            0.999999999999,
        )
        potential_v = 4.0 * thermal_voltage * np.arctanh(atanh_argument)
        reduced_potential = potential_v / thermal_voltage
        positive_ratio = np.exp(np.clip(-reduced_potential, -700.0, 700.0))
        negative_ratio = np.exp(np.clip(+reduced_potential, -700.0, 700.0))
        linear_potential_v = float(surface_potential_v) * np.exp(-distance / lambda_d)
        chemical_positive_ev = -z_value * potential_v
        electrical_positive_ev = +z_value * potential_v
        return {
            "lambda_d_m": lambda_d,
            "potential_v": potential_v,
            "linear_potential_v": linear_potential_v,
            "positive_ratio": positive_ratio,
            "negative_ratio": negative_ratio,
            "chemical_positive_ev": chemical_positive_ev,
            "electrical_positive_ev": electrical_positive_ev,
            "electrochemical_positive_ev": chemical_positive_ev + electrical_positive_ev,
        }

    def mott_schottky_width_m(
        surface_potential_v,
        concentration_cm3,
        epsilon_r,
        charge_magnitude,
    ):
        """Depletion width lambda for a positive core and frozen negative dopants."""
        potential = _positive("surface_potential_v", surface_potential_v)
        concentration_m3 = _positive("concentration_cm3", concentration_cm3) * 1.0e6
        relative_permittivity = _positive("epsilon_r", epsilon_r)
        z_value = _positive("charge_magnitude", charge_magnitude)
        permittivity = EPSILON_0_F_PER_M * relative_permittivity
        return np.sqrt(
            2.0 * permittivity * potential
            / (z_value * E_CHARGE_C * concentration_m3)
        )

    def mott_schottky_profile(
        distance_m,
        surface_potential_v,
        temperature_k,
        concentration_cm3,
        epsilon_r,
        charge_magnitude,
    ):
        """Lecture depletion profile with a Boltzmann mobile positive defect."""
        distance = np.asarray(distance_m, dtype=float)
        temperature = _positive("temperature_k", temperature_k)
        z_value = _positive("charge_magnitude", charge_magnitude)
        width = mott_schottky_width_m(
            surface_potential_v,
            concentration_cm3,
            epsilon_r,
            z_value,
        )
        inside = distance <= width
        potential_v = np.where(
            inside,
            float(surface_potential_v) * (distance / width - 1.0) ** 2,
            0.0,
        )
        reduced_potential = z_value * E_CHARGE_C * potential_v / (
            KB_J_PER_K * temperature
        )
        positive_ratio = np.exp(np.clip(-reduced_potential, -700.0, 700.0))
        negative_ratio = np.ones_like(distance)
        chemical_positive_ev = -z_value * potential_v
        electrical_positive_ev = +z_value * potential_v
        return {
            "width_m": width,
            "potential_v": potential_v,
            "positive_ratio": positive_ratio,
            "negative_ratio": negative_ratio,
            "chemical_positive_ev": chemical_positive_ev,
            "electrical_positive_ev": electrical_positive_ev,
            "electrochemical_positive_ev": chemical_positive_ev + electrical_positive_ev,
        }

    def ms_surface_charge_c_per_m2(
        surface_potential_v,
        concentration_cm3,
        epsilon_r,
        charge_magnitude,
    ):
        concentration_m3 = _positive("concentration_cm3", concentration_cm3) * 1.0e6
        z_value = _positive("charge_magnitude", charge_magnitude)
        width = mott_schottky_width_m(
            surface_potential_v,
            concentration_cm3,
            epsilon_r,
            z_value,
        )
        return z_value * E_CHARGE_C * concentration_m3 * width


    def gcs_diffuse_potential_v(
        core_potential_v,
        temperature_k,
        concentration_cm3,
        epsilon_r,
        charge_magnitude,
        stern_capacitance_f_per_m2,
    ):
        """Solve C_s(phi_0-phi_1)=Q_GC(phi_1) by monotone bisection."""
        phi0 = float(core_potential_v)
        capacitance = _positive(
            "stern_capacitance_f_per_m2",
            stern_capacitance_f_per_m2,
        )
        if phi0 == 0.0:
            return 0.0

        def mismatch(phi1):
            return capacitance * (phi0 - phi1) - gc_surface_charge_c_per_m2(
                phi1,
                temperature_k,
                concentration_cm3,
                epsilon_r,
                charge_magnitude,
            )

        lower = min(0.0, phi0)
        upper = max(0.0, phi0)
        f_lower = mismatch(lower)
        for _iteration in range(90):
            midpoint = 0.5 * (lower + upper)
            f_midpoint = mismatch(midpoint)
            if f_lower * f_midpoint <= 0.0:
                upper = midpoint
            else:
                lower = midpoint
                f_lower = f_midpoint
        return 0.5 * (lower + upper)

    def gcs_state(
        core_potential_v,
        temperature_k,
        concentration_cm3,
        epsilon_r,
        charge_magnitude,
        stern_capacitance_f_per_m2,
    ):
        phi1 = gcs_diffuse_potential_v(
            core_potential_v,
            temperature_k,
            concentration_cm3,
            epsilon_r,
            charge_magnitude,
            stern_capacitance_f_per_m2,
        )
        charge = gc_surface_charge_c_per_m2(
            phi1,
            temperature_k,
            concentration_cm3,
            epsilon_r,
            charge_magnitude,
        )
        diffuse_capacitance = gc_differential_capacitance_f_per_m2(
            phi1,
            temperature_k,
            concentration_cm3,
            epsilon_r,
            charge_magnitude,
        )
        total_capacitance = 1.0 / (
            1.0 / stern_capacitance_f_per_m2 + 1.0 / diffuse_capacitance
        )
        return {
            "phi1_v": phi1,
            "stern_drop_v": float(core_potential_v) - phi1,
            "charge_c_per_m2": charge,
            "diffuse_capacitance_f_per_m2": diffuse_capacitance,
            "total_capacitance_f_per_m2": total_capacitance,
        }

    def frumkin_log10_anodic_current(
        core_potentials_v,
        temperature_k,
        concentration_cm3,
        epsilon_r,
        charge_magnitude,
        stern_capacitance_f_per_m2,
        transfer_coefficient,
        reactant_charge_number,
    ):
        """Normalized anodic current using the reaction-plane concentration and phi_1."""
        potentials = np.asarray(core_potentials_v, dtype=float)
        alpha = float(transfer_coefficient)
        z_reactant = float(reactant_charge_number)
        phi1_values = np.array(
            [
                gcs_diffuse_potential_v(
                    phi0,
                    temperature_k,
                    concentration_cm3,
                    epsilon_r,
                    charge_magnitude,
                    stern_capacitance_f_per_m2,
                )
                for phi0 in potentials
            ]
        )
        molar_scale = FARADAY_C_PER_MOL / (
            GAS_CONSTANT_J_PER_MOL_K * float(temperature_k)
        )
        log10_local_ratio = -z_reactant * molar_scale * phi1_values / np.log(10.0)
        log10_naive_current = (
            (1.0 - alpha) * molar_scale * potentials / np.log(10.0)
        )
        log10_corrected_current = (
            (1.0 - alpha) * molar_scale * (potentials - phi1_values)
            - z_reactant * molar_scale * phi1_values
        ) / np.log(10.0)
        return {
            "phi1_v": phi1_values,
            "stern_drop_v": potentials - phi1_values,
            "log10_local_ratio": log10_local_ratio,
            "log10_naive_current": log10_naive_current,
            "log10_corrected_current": log10_corrected_current,
        }

    return (
        E_CHARGE_C,
        EPSILON_0_F_PER_M,
        FARADAY_C_PER_MOL,
        GAS_CONSTANT_J_PER_MOL_K,
        KB_J_PER_K,
        debye_length_m,
        frumkin_log10_anodic_current,
        gc_differential_capacitance_f_per_m2,
        gc_surface_charge_c_per_m2,
        gcs_diffuse_potential_v,
        gcs_state,
        gouy_chapman_profile,
        mott_schottky_profile,
        mott_schottky_width_m,
        ms_surface_charge_c_per_m2,
    )


@app.cell
def _(mo):
    temperature_control = mo.ui.slider(
        start=300,
        stop=1200,
        step=25,
        value=800,
        label="temperature, T (K)",
        show_value=True,
    )
    log_concentration_control = mo.ui.slider(
        start=16.0,
        stop=21.0,
        step=0.25,
        value=18.0,
        label="log10 bulk concentration, c_i,infinity (cm^-3)",
        show_value=True,
    )
    epsilon_r_control = mo.ui.slider(
        start=5,
        stop=300,
        step=5,
        value=100,
        label="relative permittivity, epsilon_r",
        show_value=True,
    )
    charge_control = mo.ui.slider(
        start=1,
        stop=3,
        step=1,
        value=1,
        label="defect charge magnitude, z",
        show_value=True,
    )
    surface_potential_control = mo.ui.slider(
        start=0.01,
        stop=0.25,
        step=0.01,
        value=0.16,
        label="core potential, phi_0 (V)",
        show_value=True,
    )
    show_linear_control = mo.ui.checkbox(
        value=True,
        label="show the low-potential Gouy-Chapman approximation",
    )
    profile_controls = mo.vstack(
        [
            mo.hstack(
                [temperature_control, log_concentration_control, epsilon_r_control],
                justify="start",
                align="center",
                wrap=True,
                gap=1.2,
            ),
            mo.hstack(
                [charge_control, surface_potential_control, show_linear_control],
                justify="start",
                align="center",
                wrap=True,
                gap=1.2,
            ),
        ]
    )
    return (
        charge_control,
        epsilon_r_control,
        log_concentration_control,
        profile_controls,
        show_linear_control,
        surface_potential_control,
        temperature_control,
    )


@app.cell
def _(mo, profile_controls):
    mo.vstack(
        [
            mo.md(r"""
            ## 1. From electrochemical equilibrium to a space-charge profile

            At equilibrium, every mobile defect has a flat electrochemical
            potential:

            \[
            \widetilde\mu_i(x)=\mu_i^0+k_BT\ln c_i(x)+z_i e\phi(x)
            =\text{constant}.
            \]

            Taking the bulk as the reference gives the Boltzmann distribution

            \[
            \frac{c_i(x)}{c_{i,\infty}}
            =\exp\!\left[-\frac{z_i e\phi(x)}{k_BT}\right].
            \]

            The chemical term \(k_BT\ln(c_i/c_{i,\infty})\) and electrical term
            \(z_i e\phi\) change in opposite directions. Their sum remains
            constant. For a mobile \(+ze/-ze\) pair,

            \[
            \rho(x)=ze[c_+(x)-c_-(x)],
            \qquad
            \frac{d^2\phi}{dx^2}=-\frac{\rho(x)}{\epsilon_0\epsilon_r}.
            \]

            This is the self-consistency loop: potential redistributes charged
            defects, and those defects create the charge density that bends the
            potential. With a positive core, \(\phi>0\): positive co-ions are
            depleted, negative counter-ions accumulate, and the surrounding
            space charge is negative.
            """),
            profile_controls,
        ]
    )
    return


@app.cell
def _(
    E_CHARGE_C,
    KB_J_PER_K,
    charge_control,
    debye_length_m,
    epsilon_r_control,
    gouy_chapman_profile,
    log_concentration_control,
    mott_schottky_profile,
    mott_schottky_width_m,
    ms_surface_charge_c_per_m2,
    np,
    surface_potential_control,
    temperature_control,
):
    temperature_k = float(temperature_control.value)
    bulk_concentration_cm3 = 10.0 ** float(log_concentration_control.value)
    epsilon_r = float(epsilon_r_control.value)
    charge_magnitude = float(charge_control.value)
    surface_potential_v = float(surface_potential_control.value)

    selected_debye_length_m = debye_length_m(
        temperature_k,
        bulk_concentration_cm3,
        epsilon_r,
        charge_magnitude,
    )
    selected_ms_width_m = mott_schottky_width_m(
        surface_potential_v,
        bulk_concentration_cm3,
        epsilon_r,
        charge_magnitude,
    )
    selected_ms_core_charge_c_per_m2 = ms_surface_charge_c_per_m2(
        surface_potential_v,
        bulk_concentration_cm3,
        epsilon_r,
        charge_magnitude,
    )
    surface_reduced_potential = (
        charge_magnitude * E_CHARGE_C * surface_potential_v
        / (KB_J_PER_K * temperature_k)
    )
    ms_surface_mobile_ratio = np.exp(-surface_reduced_potential)

    gc_distance_m = np.linspace(0.0, 6.0 * selected_debye_length_m, 700)
    ms_distance_m = np.linspace(0.0, 1.5 * selected_ms_width_m, 700)
    gc_profile = gouy_chapman_profile(
        gc_distance_m,
        surface_potential_v,
        temperature_k,
        bulk_concentration_cm3,
        epsilon_r,
        charge_magnitude,
    )
    ms_profile = mott_schottky_profile(
        ms_distance_m,
        surface_potential_v,
        temperature_k,
        bulk_concentration_cm3,
        epsilon_r,
        charge_magnitude,
    )
    return (
        bulk_concentration_cm3,
        charge_magnitude,
        epsilon_r,
        gc_distance_m,
        gc_profile,
        ms_distance_m,
        ms_profile,
        ms_surface_mobile_ratio,
        selected_debye_length_m,
        selected_ms_core_charge_c_per_m2,
        selected_ms_width_m,
        surface_potential_v,
        surface_reduced_potential,
        temperature_k,
    )


@app.cell
def _(
    gc_profile,
    mo,
    gc_distance_m,
    plt,
    selected_debye_length_m,
    show_linear_control,
    surface_reduced_potential,
):
    gc_figure, gc_axes = plt.subplots(1, 3, figsize=(14.4, 4.7), dpi=120)
    gc_distance_nm = gc_distance_m * 1.0e9
    gc_axes[0].plot(
        gc_distance_nm,
        gc_profile["potential_v"],
        color="#007C91",
        lw=3.0,
        label="exact Gouy-Chapman",
    )
    if bool(show_linear_control.value):
        gc_axes[0].plot(
            gc_distance_nm,
            gc_profile["linear_potential_v"],
            color="#EE9B00",
            lw=2.2,
            ls="--",
            label="small-potential limit",
        )
    gc_axes[0].axvline(
        selected_debye_length_m * 1.0e9,
        color="#666666",
        lw=1.2,
        ls=":",
        label=r"$\lambda_D$",
    )
    gc_axes[0].set(
        xlabel="distance from core, x (nm)",
        ylabel=r"electrostatic potential, $\phi$ (V)",
        title="Potential is screened into the bulk",
    )
    gc_axes[0].grid(alpha=0.22)
    gc_axes[0].legend(frameon=False, fontsize=8.5)

    gc_axes[1].semilogy(
        gc_distance_nm,
        gc_profile["positive_ratio"],
        color="#CC3311",
        lw=3.0,
        label=r"$c_+/c_{i,\infty}$",
    )
    gc_axes[1].semilogy(
        gc_distance_nm,
        gc_profile["negative_ratio"],
        color="#007C91",
        lw=3.0,
        label=r"$c_-/c_{i,\infty}$",
    )
    gc_axes[1].axhline(1.0, color="#777777", lw=1.0, ls=":")
    gc_axes[1].set(
        xlabel="distance from core, x (nm)",
        ylabel="concentration / bulk concentration",
        title="Co-ions deplete; counter-ions accumulate",
    )
    gc_axes[1].grid(which="both", alpha=0.22)
    gc_axes[1].legend(frameon=False)

    gc_axes[2].plot(
        gc_distance_nm,
        1.0e3 * gc_profile["chemical_positive_ev"],
        color="#CC3311",
        lw=3.0,
        label=r"chemical: $\Delta\mu_+$",
    )
    gc_axes[2].plot(
        gc_distance_nm,
        1.0e3 * gc_profile["electrical_positive_ev"],
        color="#007C91",
        lw=3.0,
        label=r"electrical: $+ze\phi$",
    )
    gc_axes[2].plot(
        gc_distance_nm,
        1.0e3 * gc_profile["electrochemical_positive_ev"],
        color="#222222",
        lw=1.8,
        ls="--",
        label=r"sum: $\Delta\widetilde\mu_+$",
    )
    gc_axes[2].set(
        xlabel="distance from core, x (nm)",
        ylabel="energy change (meV per defect)",
        title="Electrochemical potential stays flat",
    )
    gc_axes[2].grid(alpha=0.22)
    gc_axes[2].legend(frameon=False, fontsize=8.5)
    gc_figure.tight_layout()
    plt.close(gc_figure)

    gc_explanation = mo.md(
        rf"""
        ### Gouy–Chapman: all charged defects can move

        For the symmetric \(+ze/-ze\) pair, the exact planar solution is

        \[
        \frac{{\tanh[ze\phi(x)/(4k_BT)]}}
        {{\tanh[ze\phi_0/(4k_BT)]}}=e^{{-x/\lambda_D}},
        \qquad
        \lambda_D=\sqrt{{\frac{{\epsilon_0\epsilon_r k_BT}}
        {{2z^2e^2c_{{i,\infty}}}}}}.
        \]

        Bulk electroneutrality requires
        \(c_+(\infty)=c_-(\infty)=c_{{i,\infty}}\). Here
        \(\lambda_D=\mathbf{{{selected_debye_length_m * 1.0e9:.2f}\ nm}}\).
        A larger \(c_{{i,\infty}}\) gives more screening charge and a shorter
        layer; a larger \(\epsilon_r\) makes the layer wider.

        The important small-potential test is
        \(|ze\phi_0/(k_BT)|\ll1\). Its selected value is
        **{abs(surface_reduced_potential):.2f}**, so the dashed exponential is a
        limiting guide, not the curve used in the calculation. The solid curve
        always uses the exact Gouy–Chapman solution.

        **Continuum check.** The selected \(\lambda_D\) is
        {selected_debye_length_m * 1.0e9:.2f} nm. When a screening length becomes
        comparable to a lattice spacing (roughly below 1 nm), a continuum
        profile should be treated as qualitative and atomistic structure matters.
        """
    )
    mo.vstack([gc_explanation, gc_figure])
    return (gc_figure,)


@app.cell
def _(
    mo,
    ms_profile,
    np,
    plt,
    ms_distance_m,
    ms_surface_mobile_ratio,
    selected_ms_core_charge_c_per_m2,
    selected_ms_width_m,
):
    ms_figure, ms_axes = plt.subplots(1, 3, figsize=(14.4, 4.7), dpi=120)
    ms_distance_nm = ms_distance_m * 1.0e9
    ms_axes[0].plot(
        ms_distance_nm,
        ms_profile["potential_v"],
        color="#007C91",
        lw=3.0,
    )
    ms_axes[0].axvline(
        selected_ms_width_m * 1.0e9,
        color="#666666",
        lw=1.3,
        ls=":",
        label=r"depletion width $\lambda$",
    )
    ms_axes[0].set(
        xlabel="distance from core, x (nm)",
        ylabel=r"electrostatic potential, $\phi$ (V)",
        title="Frozen charge gives a parabola",
    )
    ms_axes[0].grid(alpha=0.22)
    ms_axes[0].legend(frameon=False)

    ms_axes[1].semilogy(
        ms_distance_nm,
        ms_profile["positive_ratio"],
        color="#CC3311",
        lw=3.0,
        label=r"mobile $c_+/c_{i,\infty}$",
    )
    ms_axes[1].semilogy(
        ms_distance_nm,
        ms_profile["negative_ratio"],
        color="#007C91",
        lw=2.5,
        ls="--",
        label=r"frozen $c_-/c_{i,\infty}$",
    )
    ms_axes[1].set(
        xlabel="distance from core, x (nm)",
        ylabel="concentration / bulk concentration",
        title="Mobile positive defects are depleted",
    )
    ms_axes[1].grid(which="both", alpha=0.22)
    ms_axes[1].legend(frameon=False)

    ms_axes[2].plot(
        ms_distance_nm,
        1.0e3 * ms_profile["chemical_positive_ev"],
        color="#CC3311",
        lw=3.0,
        label=r"chemical: $\Delta\mu_+$",
    )
    ms_axes[2].plot(
        ms_distance_nm,
        1.0e3 * ms_profile["electrical_positive_ev"],
        color="#007C91",
        lw=3.0,
        label=r"electrical: $+ze\phi$",
    )
    ms_axes[2].plot(
        ms_distance_nm,
        1.0e3 * ms_profile["electrochemical_positive_ev"],
        color="#222222",
        lw=1.8,
        ls="--",
        label=r"sum: $\Delta\widetilde\mu_+$",
    )
    ms_axes[2].set(
        xlabel="distance from core, x (nm)",
        ylabel="energy change (meV per defect)",
        title="The same equilibrium cancellation remains",
    )
    ms_axes[2].grid(alpha=0.22)
    ms_axes[2].legend(frameon=False, fontsize=8.5)
    ms_figure.tight_layout()
    plt.close(ms_figure)

    ms_explanation = mo.md(
        rf"""
        ## 2. Mott–Schottky: the majority dopants are frozen

        The negative defects remain at \(c_{{i,\infty}}\). When the positive
        mobile defects are strongly depleted, the space charge is approximately
        constant: \(\rho\approx-zec_{{i,\infty}}\). Poisson's equation then gives

        \[
        \phi(x)=\phi_0\left(\frac{{x}}{{\lambda}}-1\right)^2,
        \quad 0\le x\le\lambda,
        \qquad
        \lambda=\sqrt{{\frac{{2\epsilon_0\epsilon_r\phi_0}}
        {{ze c_{{i,\infty}}}}}}.
        \]

        Global charge compensation gives the second lecture relation

        \[
        Q_{{\rm core}}=-Q_{{\rm sc}}\approx ze c_{{i,\infty}}\lambda.
        \]

        Here \(\lambda=\mathbf{{{selected_ms_width_m * 1.0e9:.2f}\ nm}}\) and
        \(Q_{{\rm core}}={selected_ms_core_charge_c_per_m2:.3e}\) C/m². The
        plotted mobile concentration remains the positive Boltzmann value rather
        than being set to a literal zero. At the core,
        \(c_+(0)/c_{{i,\infty}}=\mathbf{{{ms_surface_mobile_ratio:.3f}}}\).
        The **Mott–Schottky approximation** is reliable when this ratio is much
        smaller than one, so the frozen negative dopants dominate the charge
        inside the depletion layer.
        """
    )
    mo.vstack([ms_explanation, ms_figure])
    return (ms_figure,)


@app.cell
def _(
    bulk_concentration_cm3,
    charge_magnitude,
    debye_length_m,
    epsilon_r,
    mo,
    mott_schottky_width_m,
    np,
    plt,
    surface_potential_v,
    temperature_k,
):
    concentration_sweep_cm3 = np.logspace(16.0, 21.0, 350)
    debye_sweep_nm = np.array(
        [
            debye_length_m(
                temperature_k,
                concentration,
                epsilon_r,
                charge_magnitude,
            )
            for concentration in concentration_sweep_cm3
        ]
    ) * 1.0e9
    ms_sweep_nm = np.array(
        [
            mott_schottky_width_m(
                surface_potential_v,
                concentration,
                epsilon_r,
                charge_magnitude,
            )
            for concentration in concentration_sweep_cm3
        ]
    ) * 1.0e9
    epsilon_sweep = np.logspace(np.log10(5.0), np.log10(300.0), 300)
    debye_epsilon_nm = np.array(
        [
            debye_length_m(
                temperature_k,
                bulk_concentration_cm3,
                relative_permittivity,
                charge_magnitude,
            )
            for relative_permittivity in epsilon_sweep
        ]
    ) * 1.0e9
    ms_epsilon_nm = np.array(
        [
            mott_schottky_width_m(
                surface_potential_v,
                bulk_concentration_cm3,
                relative_permittivity,
                charge_magnitude,
            )
            for relative_permittivity in epsilon_sweep
        ]
    ) * 1.0e9

    width_figure, width_axes = plt.subplots(1, 2, figsize=(12.8, 4.6), dpi=120)
    width_axes[0].loglog(
        concentration_sweep_cm3,
        debye_sweep_nm,
        color="#007C91",
        lw=3.0,
        label=r"Gouy-Chapman $\lambda_D$",
    )
    width_axes[0].loglog(
        concentration_sweep_cm3,
        ms_sweep_nm,
        color="#CC3311",
        lw=3.0,
        label=r"Mott-Schottky $\lambda$",
    )
    width_axes[0].axvline(
        bulk_concentration_cm3,
        color="#EE9B00",
        lw=1.5,
        ls="--",
        label="selected bulk concentration",
    )
    width_axes[0].set(
        xlabel=r"bulk concentration, $c_{i,\infty}$ (cm$^{-3}$)",
        ylabel="space-charge length (nm)",
        title="More defects screen over a shorter distance",
    )
    width_axes[0].grid(which="both", alpha=0.22)
    width_axes[0].legend(frameon=False)

    width_axes[1].loglog(
        epsilon_sweep,
        debye_epsilon_nm,
        color="#007C91",
        lw=3.0,
        label=r"Gouy-Chapman $\lambda_D$",
    )
    width_axes[1].loglog(
        epsilon_sweep,
        ms_epsilon_nm,
        color="#CC3311",
        lw=3.0,
        label=r"Mott-Schottky $\lambda$",
    )
    width_axes[1].axvline(
        epsilon_r,
        color="#EE9B00",
        lw=1.5,
        ls="--",
        label="selected permittivity",
    )
    width_axes[1].set(
        xlabel=r"relative permittivity, $\epsilon_r$",
        ylabel="space-charge length (nm)",
        title="Larger permittivity broadens the layer",
    )
    width_axes[1].grid(which="both", alpha=0.22)
    width_axes[1].legend(frameon=False)
    width_figure.tight_layout()
    plt.close(width_figure)

    mo.vstack(
        [
            mo.md(r"""
            ### Compare the two screening lengths

            Both lengths scale as \(c_{i,\infty}^{-1/2}\) and
            \(\epsilon_r^{1/2}\), but the interface enters differently:

            | case | screening charge | dependence on the core |
            |---|---|---|
            | Gouy–Chapman | two continuously redistributed mobile species | \(\lambda_D\) is independent of \(Q_{\rm core}\) and \(\phi_0\); changing \(\phi_0\) changes the profile amplitude |
            | Mott–Schottky | exposed frozen dopants in a depletion region | \(\lambda\propto\sqrt{\phi_0}\) and \(\lambda=Q_{\rm core}/(ze c_{i,\infty})\) |

            A higher bulk concentration shortens either layer because more
            charge is available per unit distance to screen the core.
            """),
            width_figure,
        ]
    )
    return (width_figure,)


@app.cell
def _(mo):
    stern_capacitance_control = mo.ui.slider(
        start=5,
        stop=80,
        step=1,
        value=20,
        label="Stern capacitance, C_s (microF/cm^2)",
        show_value=True,
    )

    transfer_coefficient_control = mo.ui.slider(
        start=0.1,
        stop=0.9,
        step=0.05,
        value=0.5,
        label="transfer coefficient, alpha",
        show_value=True,
    )
    reactant_charge_control = mo.ui.dropdown(
        options={
            "-2 (doubly negative R)": -2,
            "-1 (singly negative R)": -1,
            "0 (neutral R)": 0,
            "+1 (singly positive R)": 1,
            "+2 (doubly positive R)": 2,
        },
        value="+1 (singly positive R)",
        label="signed reactant charge number, z_R",
    )
    interface_controls = mo.hstack(
        [
            stern_capacitance_control,
            transfer_coefficient_control,
            reactant_charge_control,
        ],
        justify="start",
        align="center",
        wrap=True,
        gap=1.2,
    )
    return (
        interface_controls,
        reactant_charge_control,
        stern_capacitance_control,
        transfer_coefficient_control,
    )


@app.cell
def _(interface_controls, mo):
    mo.vstack(
        [
            mo.md(r"""
            ## 3. From charge to capacitance: Gouy–Chapman–Stern

            Gauss's law converts the Gouy–Chapman potential gradient into the
            core charge per area:

            \[
            Q_{\rm core}=\epsilon_0\epsilon_r
            \frac{2k_BT}{ze\lambda_D}
            \sinh\!\left(\frac{ze\phi_1}{2k_BT}\right),
            \qquad
            C_d=\frac{dQ_{\rm core}}{d\phi_1}
            =\frac{\epsilon_0\epsilon_r}{\lambda_D}
            \cosh\!\left(\frac{ze\phi_1}{2k_BT}\right).
            \]

            In this ideal symmetric model, \(\phi_0=0\) is the **point of zero
            charge (pZC)**. Gouy–Chapman alone predicts a rapidly increasing
            capacitance at large \(|\phi|\) because point defects can approach
            the core without a size limit. The **Stern layer** represents a
            finite closest approach. Its linear potential drop and the diffuse
            Gouy–Chapman layer carry the same charge:

            \[
            Q_{\rm core}=C_s(\phi_0-\phi_1)=Q_{\rm GC}(\phi_1),
            \qquad
            \frac{1}{C_{\rm tot}}=\frac{1}{C_s}+\frac{1}{C_d}.
            \]

            The profile drawing places \(x_1\) at 0.5 nm only to make the two
            regions visible. The physical compact-layer input is \(C_s\); an
            independent thickness would require a separate Stern permittivity.
            """),
            interface_controls,
        ]
    )
    return


@app.cell
def _(
    bulk_concentration_cm3,
    charge_magnitude,
    epsilon_r,
    gc_differential_capacitance_f_per_m2,
    gcs_state,
    gouy_chapman_profile,
    np,
    stern_capacitance_control,
    surface_potential_v,
    temperature_k,
):
    stern_capacitance_f_per_m2 = float(stern_capacitance_control.value) * 0.01
    stern_thickness_m = 0.5e-9
    selected_gcs_state = gcs_state(
        surface_potential_v,
        temperature_k,
        bulk_concentration_cm3,
        epsilon_r,
        charge_magnitude,
        stern_capacitance_f_per_m2,
    )
    gcs_potential_sweep_v = np.linspace(-0.25, 0.25, 301)
    gcs_states = [
        gcs_state(
            phi0,
            temperature_k,
            bulk_concentration_cm3,
            epsilon_r,
            charge_magnitude,
            stern_capacitance_f_per_m2,
        )
        for phi0 in gcs_potential_sweep_v
    ]
    gcs_phi1_sweep_v = np.array([state["phi1_v"] for state in gcs_states])
    gcs_charge_sweep_c_per_m2 = np.array(
        [state["charge_c_per_m2"] for state in gcs_states]
    )
    gcs_diffuse_capacitance_sweep = np.array(
        [state["diffuse_capacitance_f_per_m2"] for state in gcs_states]
    )
    gcs_total_capacitance_sweep = np.array(
        [state["total_capacitance_f_per_m2"] for state in gcs_states]
    )

    gcs_diffuse_distance_m = np.linspace(0.0, 6.0 * (
        gouy_chapman_profile(
            np.array([0.0]),
            selected_gcs_state["phi1_v"],
            temperature_k,
            bulk_concentration_cm3,
            epsilon_r,
            charge_magnitude,
        )["lambda_d_m"]
    ), 500)
    gcs_diffuse_profile = gouy_chapman_profile(
        gcs_diffuse_distance_m,
        selected_gcs_state["phi1_v"],
        temperature_k,
        bulk_concentration_cm3,
        epsilon_r,
        charge_magnitude,
    )
    stern_distance_m = np.linspace(0.0, stern_thickness_m, 80)
    stern_potential_v = np.linspace(
        surface_potential_v,
        selected_gcs_state["phi1_v"],
        stern_distance_m.size,
    )

    bare_gc_capacitance_selected = gc_differential_capacitance_f_per_m2(
        surface_potential_v,
        temperature_k,
        bulk_concentration_cm3,
        epsilon_r,
        charge_magnitude,
    )

    return (
        bare_gc_capacitance_selected,
        gcs_charge_sweep_c_per_m2,
        gcs_diffuse_capacitance_sweep,
        gcs_diffuse_distance_m,
        gcs_diffuse_profile,
        gcs_phi1_sweep_v,
        gcs_potential_sweep_v,
        gcs_states,
        gcs_total_capacitance_sweep,
        selected_gcs_state,
        stern_capacitance_f_per_m2,
        stern_distance_m,
        stern_potential_v,
        stern_thickness_m,
    )


@app.cell
def _(
    gcs_diffuse_capacitance_sweep,
    gcs_diffuse_distance_m,
    gcs_diffuse_profile,
    gcs_potential_sweep_v,
    gcs_total_capacitance_sweep,
    mo,
    plt,
    selected_gcs_state,
    stern_capacitance_f_per_m2,
    stern_distance_m,
    stern_potential_v,
    stern_thickness_m,
    surface_potential_v,
):
    gcs_figure, gcs_axes = plt.subplots(1, 2, figsize=(13.2, 4.8), dpi=120)
    gcs_axes[0].plot(
        stern_distance_m * 1.0e9,
        stern_potential_v,
        color="#EE9B00",
        lw=3.2,
        label="Stern layer",
    )
    gcs_axes[0].plot(
        (stern_thickness_m + gcs_diffuse_distance_m) * 1.0e9,
        gcs_diffuse_profile["potential_v"],
        color="#007C91",
        lw=3.2,
        label="diffuse layer",
    )
    gcs_axes[0].axvline(
        stern_thickness_m * 1.0e9,
        color="#666666",
        lw=1.2,
        ls=":",
        label=r"reaction plane $x_1$",
    )
    gcs_axes[0].scatter(
        [0.0, stern_thickness_m * 1.0e9],
        [surface_potential_v, selected_gcs_state["phi1_v"]],
        s=70,
        color=["#CC3311", "#007C91"],
        edgecolor="#222222",
        zorder=5,
    )
    gcs_axes[0].set(
        xlabel="distance from core, x (nm)",
        ylabel=r"electrostatic potential, $\phi$ (V)",
        title=r"The total drop splits into $\phi_0-\phi_1$ and $\phi_1$",
    )
    gcs_axes[0].grid(alpha=0.22)
    gcs_axes[0].legend(frameon=False)

    gcs_axes[1].plot(
        gcs_potential_sweep_v,
        gcs_diffuse_capacitance_sweep / 0.01,
        color="#007C91",
        lw=3.0,
        label=r"diffuse $C_d$",
    )
    gcs_axes[1].axhline(
        stern_capacitance_f_per_m2 / 0.01,
        color="#EE9B00",
        lw=2.5,
        ls="--",
        label=r"Stern $C_s$",
    )
    gcs_axes[1].plot(
        gcs_potential_sweep_v,
        gcs_total_capacitance_sweep / 0.01,
        color="#CC3311",
        lw=3.2,
        label=r"series total $C_{\rm tot}$",
    )
    gcs_axes[1].axvline(surface_potential_v, color="#555555", lw=1.2, ls=":")
    gcs_axes[1].set(
        xlabel=r"core potential relative to pZC, $\phi_0$ (V)",
        ylabel=r"differential capacitance ($\mu$F cm$^{-2}$)",
        title="The Stern layer limits the high-field capacitance",
    )
    gcs_axes[1].grid(alpha=0.22)
    gcs_axes[1].legend(frameon=False)
    gcs_figure.tight_layout()
    plt.close(gcs_figure)

    gcs_summary = mo.md(
        rf"""
        At the selected \(\phi_0={surface_potential_v:.3f}\) V, the calculation
        gives \(\phi_1=\mathbf{{{selected_gcs_state['phi1_v']:.4f}\ V}}\) and
        \(\phi_0-\phi_1=\mathbf{{{selected_gcs_state['stern_drop_v']:.4f}\ V}}\).
        The common charge is
        \(Q_{{\rm core}}={selected_gcs_state['charge_c_per_m2']:.3e}\) C/m².
        The diffuse and total capacitances are
        \(C_d={selected_gcs_state['diffuse_capacitance_f_per_m2'] / 0.01:.2f}\)
        and \(C_{{\rm tot}}={selected_gcs_state['total_capacitance_f_per_m2'] / 0.01:.2f}\)
        μF/cm². Because the two layers are in series, the total is always below
        either individual capacitance.
        """
    )
    mo.vstack([gcs_figure, gcs_summary])
    return (gcs_figure,)


@app.cell
def _(
    bulk_concentration_cm3,
    charge_magnitude,
    epsilon_r,
    frumkin_log10_anodic_current,
    np,
    reactant_charge_control,
    stern_capacitance_f_per_m2,
    surface_potential_v,
    temperature_k,
    transfer_coefficient_control,
):
    transfer_coefficient = float(transfer_coefficient_control.value)
    reactant_charge_number = float(reactant_charge_control.value)
    frumkin_potential_sweep_v = np.linspace(-0.25, 0.25, 301)
    frumkin_sweep = frumkin_log10_anodic_current(
        frumkin_potential_sweep_v,
        temperature_k,
        bulk_concentration_cm3,
        epsilon_r,
        charge_magnitude,
        stern_capacitance_f_per_m2,
        transfer_coefficient,
        reactant_charge_number,
    )
    selected_frumkin = frumkin_log10_anodic_current(
        np.array([surface_potential_v]),
        temperature_k,
        bulk_concentration_cm3,
        epsilon_r,
        charge_magnitude,
        stern_capacitance_f_per_m2,
        transfer_coefficient,
        reactant_charge_number,
    )
    return (
        frumkin_potential_sweep_v,
        frumkin_sweep,
        reactant_charge_number,
        selected_frumkin,
        transfer_coefficient,
    )


@app.cell
def _(
    frumkin_potential_sweep_v,
    frumkin_sweep,
    mo,
    plt,
    reactant_charge_number,
    selected_frumkin,
    surface_potential_v,
    transfer_coefficient,
):
    frumkin_figure, frumkin_axes = plt.subplots(1, 3, figsize=(14.4, 4.7), dpi=120)
    frumkin_axes[0].plot(
        frumkin_potential_sweep_v,
        frumkin_sweep["phi1_v"],
        color="#007C91",
        lw=3.0,
        label=r"diffuse drop $\phi_1$",
    )
    frumkin_axes[0].plot(
        frumkin_potential_sweep_v,
        frumkin_sweep["stern_drop_v"],
        color="#EE9B00",
        lw=3.0,
        label=r"Stern drop $\phi_0-\phi_1$",
    )
    frumkin_axes[0].plot(
        frumkin_potential_sweep_v,
        frumkin_potential_sweep_v,
        color="#777777",
        lw=1.2,
        ls=":",
        label=r"total $\phi_0$",
    )
    frumkin_axes[0].set(
        xlabel=r"driving potential, $\phi_0=E-E^{0\prime}$ (V)",
        ylabel="potential contribution (V)",
        title="GCS divides the applied potential",
    )
    frumkin_axes[0].grid(alpha=0.22)
    frumkin_axes[0].legend(frameon=False, fontsize=8.5)

    frumkin_axes[1].plot(
        frumkin_potential_sweep_v,
        frumkin_sweep["log10_local_ratio"],
        color="#007C91",
        lw=3.0,
    )
    frumkin_axes[1].axhline(0.0, color="#777777", lw=1.0, ls=":")
    frumkin_axes[1].set(
        xlabel=r"driving potential, $\phi_0=E-E^{0\prime}$ (V)",
        ylabel=r"$\log_{10}([R]_{x_1}/[R]_{\infty})$",
        title="The reaction-plane concentration changes",
    )
    frumkin_axes[1].grid(alpha=0.22)

    frumkin_axes[2].plot(
        frumkin_potential_sweep_v,
        frumkin_sweep["log10_naive_current"],
        color="#777777",
        lw=2.2,
        ls="--",
        label="use bulk concentration and total potential",
    )
    frumkin_axes[2].plot(
        frumkin_potential_sweep_v,
        frumkin_sweep["log10_corrected_current"],
        color="#CC3311",
        lw=3.2,
        label="Frumkin-corrected",
    )
    frumkin_axes[2].set(
        xlabel=r"driving potential, $\phi_0=E-E^{0\prime}$ (V)",
        ylabel=r"$\log_{10}$ normalized anodic current",
        title="The Frumkin effect can suppress or enhance",
    )
    frumkin_axes[2].grid(alpha=0.22)
    frumkin_axes[2].legend(frameon=False, fontsize=8.2)
    frumkin_figure.tight_layout()
    plt.close(frumkin_figure)

    selected_local_ratio = 10.0 ** float(selected_frumkin["log10_local_ratio"][0])
    selected_current_ratio = 10.0 ** float(
        selected_frumkin["log10_corrected_current"][0]
        - selected_frumkin["log10_naive_current"][0]
    )
    frumkin_text = mo.md(
        rf"""
        ## 4. The Frumkin effect: the reaction sees \(x=x_1\), not the bulk

        Following the lecture's anodic branch,

        \[
        I_a=FAk^0\exp\!\left[\frac{{(1-\alpha)F
        (E-E^{{0\prime}}-\phi_1)}}{{RT}}\right][R]_{{x=x_1}},
        \qquad
        \frac{{[R]_{{x=x_1}}}}{{[R]_{{x=\infty}}}}
        =\exp\!\left(-\frac{{z_R F\phi_1}}{{RT}}\right).
        \]

        The exchange current introduced in the lecture also uses reaction-plane
        concentrations,

        \[
        I_0=FAk^0[O]_{{x_1}}^{{1-\alpha}}[R]_{{x_1}}^\alpha,
        \]

        so a space-charge layer changes equilibrium kinetics even before a large
        overpotential is applied. The plotted anodic branch isolates the two
        Frumkin contributions:

        1. a **concentration effect** through \([R]_{{x=x_1}}\), and
        2. a **potential effect** because the reaction driving force contains
           \(E-E^{{0\prime}}-\phi_1\), not the full interfacial drop.

        To keep this teaching plot one-dimensional, we choose the formal
        potential to coincide with the pZC, so
        \(E-E^{{0\prime}}=\phi_0\). This is a reference choice for the plot, not
        a universal identity for every interface.

        Here \(z_R\) is a **signed** charge number. The lecture's notation
        \(R^{{-z}}\), with \(z>0\), corresponds to \(z_R=-z\). A positively
        charged proton-like reactant instead has \(z_R=+1\), the default used to
        reproduce the depletion picture in the Pt/BZY example. Relative to a
        curve that uses bulk concentration and the full potential,

        \[
        \log\!\left(\frac{{I_{{a,\rm Frumkin}}}}{{I_{{a,\rm naive}}}}\right)
        =-[(1-\alpha)+z_R]\frac{{F\phi_1}}{{RT}}.
        \]

        At \(\phi_0={surface_potential_v:.3f}\) V,
        \(\alpha={transfer_coefficient:.2f}\), and \(z_R={reactant_charge_number:+.0f}\),
        the local concentration ratio is **{selected_local_ratio:.3g}** and the
        corrected anodic current is **{selected_current_ratio:.3g} times** the
        naive curve. Either enhancement or suppression is possible because the
        concentration and potential contributions can compete.
        """
    )
    mo.vstack([frumkin_text, frumkin_figure])
    return (frumkin_figure,)


@app.cell
def _(
    E_CHARGE_C,
    EPSILON_0_F_PER_M,
    KB_J_PER_K,
    FARADAY_C_PER_MOL,
    GAS_CONSTANT_J_PER_MOL_K,
    bare_gc_capacitance_selected,
    bulk_concentration_cm3,
    charge_magnitude,
    epsilon_r,
    gc_distance_m,
    gc_profile,
    gc_surface_charge_c_per_m2,
    gcs_state,
    ms_distance_m,
    ms_profile,
    ms_surface_charge_c_per_m2,
    np,
    reactant_charge_number,
    selected_debye_length_m,
    selected_frumkin,
    selected_gcs_state,
    selected_ms_width_m,
    stern_capacitance_f_per_m2,
    surface_potential_v,
    temperature_k,
    transfer_coefficient,
):
    thermal_energy_j = KB_J_PER_K * temperature_k
    gc_boltzmann_residual = np.max(
        np.abs(
            np.log(gc_profile["positive_ratio"])
            + charge_magnitude * E_CHARGE_C * gc_profile["potential_v"] / thermal_energy_j
        )
    )
    gc_electrochemical_residual_ev = np.max(
        np.abs(gc_profile["electrochemical_positive_ev"])
    )
    gamma_selected = np.tanh(
        charge_magnitude * E_CHARGE_C * surface_potential_v
        / (4.0 * thermal_energy_j)
    )
    gc_tanh_identity = np.tanh(
        charge_magnitude * E_CHARGE_C * gc_profile["potential_v"]
        / (4.0 * thermal_energy_j)
    )
    gc_tanh_reference = gamma_selected * np.exp(
        -gc_distance_m / selected_debye_length_m
    )
    gc_solution_residual = np.max(np.abs(gc_tanh_identity - gc_tanh_reference))

    permittivity = EPSILON_0_F_PER_M * epsilon_r
    gc_surface_field_v_per_m = (
        2.0
        * KB_J_PER_K
        * temperature_k
        / (charge_magnitude * E_CHARGE_C * selected_debye_length_m)
        * np.sinh(
            charge_magnitude * E_CHARGE_C * surface_potential_v
            / (2.0 * thermal_energy_j)
        )
    )
    gc_gauss_residual = abs(
        permittivity * gc_surface_field_v_per_m
        - gc_surface_charge_c_per_m2(
            surface_potential_v,
            temperature_k,
            bulk_concentration_cm3,
            epsilon_r,
            charge_magnitude,
        )
    ) / max(abs(permittivity * gc_surface_field_v_per_m), 1.0e-300)

    ms_surface_residual = abs(ms_profile["potential_v"][0] - surface_potential_v)
    ms_far_residual = np.max(
        np.abs(ms_profile["potential_v"][ms_distance_m >= selected_ms_width_m])
    )
    ms_electrochemical_residual_ev = np.max(
        np.abs(ms_profile["electrochemical_positive_ev"])
    )
    expected_ms_curvature = (
        charge_magnitude
        * E_CHARGE_C
        * bulk_concentration_cm3
        * 1.0e6
        / permittivity
    )
    analytic_ms_curvature = 2.0 * surface_potential_v / selected_ms_width_m**2
    ms_poisson_residual = abs(
        analytic_ms_curvature - expected_ms_curvature
    ) / expected_ms_curvature
    ms_charge_from_width = ms_surface_charge_c_per_m2(
        surface_potential_v,
        bulk_concentration_cm3,
        epsilon_r,
        charge_magnitude,
    )
    ms_charge_from_gauss = (
        2.0 * permittivity * surface_potential_v / selected_ms_width_m
    )
    ms_charge_residual = abs(
        ms_charge_from_width - ms_charge_from_gauss
    ) / ms_charge_from_width

    capacitance_step_v = 1.0e-6
    gc_charge_plus = gc_surface_charge_c_per_m2(
        surface_potential_v + capacitance_step_v,
        temperature_k,
        bulk_concentration_cm3,
        epsilon_r,
        charge_magnitude,
    )
    gc_charge_minus = gc_surface_charge_c_per_m2(
        surface_potential_v - capacitance_step_v,
        temperature_k,
        bulk_concentration_cm3,
        epsilon_r,
        charge_magnitude,
    )
    numerical_gc_capacitance = (
        gc_charge_plus - gc_charge_minus
    ) / (2.0 * capacitance_step_v)
    gc_capacitance_residual = abs(
        numerical_gc_capacitance - bare_gc_capacitance_selected
    ) / bare_gc_capacitance_selected

    gcs_charge_residual = abs(
        stern_capacitance_f_per_m2 * selected_gcs_state["stern_drop_v"]
        - selected_gcs_state["charge_c_per_m2"]
    ) / max(abs(selected_gcs_state["charge_c_per_m2"]), 1.0e-300)
    gcs_voltage_residual = abs(
        surface_potential_v
        - selected_gcs_state["phi1_v"]
        - selected_gcs_state["stern_drop_v"]
    )
    gcs_series_capacitance_pass = bool(
        selected_gcs_state["total_capacitance_f_per_m2"]
        <= min(
            stern_capacitance_f_per_m2,
            selected_gcs_state["diffuse_capacitance_f_per_m2"],
        )
        * (1.0 + 1.0e-14)
    )
    frumkin_phi1_residual = abs(
        float(selected_frumkin["phi1_v"][0]) - selected_gcs_state["phi1_v"]
    )
    actual_frumkin_log10_factor = float(
        selected_frumkin["log10_corrected_current"][0]
        - selected_frumkin["log10_naive_current"][0]
    )
    expected_frumkin_log10_factor = (
        -((1.0 - transfer_coefficient) + reactant_charge_number)
        * FARADAY_C_PER_MOL
        * selected_gcs_state["phi1_v"]
        / (GAS_CONSTANT_J_PER_MOL_K * temperature_k * np.log(10.0))
    )
    frumkin_factor_residual = abs(
        actual_frumkin_log10_factor - expected_frumkin_log10_factor
    )
    positive_finite_values = np.array(
        [
            selected_debye_length_m,
            selected_ms_width_m,
            bare_gc_capacitance_selected,
            selected_gcs_state["diffuse_capacitance_f_per_m2"],
            selected_gcs_state["total_capacitance_f_per_m2"],
        ]
    )
    module04_validation = {
        "gc_boltzmann_residual": gc_boltzmann_residual,
        "gc_boltzmann_pass": gc_boltzmann_residual < 2.0e-13,
        "gc_electrochemical_residual_ev": gc_electrochemical_residual_ev,
        "gc_electrochemical_pass": gc_electrochemical_residual_ev < 1.0e-14,
        "gc_solution_residual": gc_solution_residual,
        "gc_solution_pass": gc_solution_residual < 2.0e-13,
        "gc_gauss_residual": gc_gauss_residual,
        "gc_gauss_pass": gc_gauss_residual < 2.0e-13,
        "gc_capacitance_residual": gc_capacitance_residual,
        "gc_capacitance_pass": gc_capacitance_residual < 2.0e-9,
        "ms_surface_residual": ms_surface_residual,
        "ms_far_residual": ms_far_residual,
        "ms_boundary_pass": max(ms_surface_residual, ms_far_residual) < 1.0e-14,
        "ms_electrochemical_residual_ev": ms_electrochemical_residual_ev,
        "ms_electrochemical_pass": ms_electrochemical_residual_ev < 1.0e-14,
        "ms_poisson_residual": ms_poisson_residual,
        "ms_poisson_pass": ms_poisson_residual < 2.0e-13,
        "ms_charge_residual": ms_charge_residual,
        "ms_charge_pass": ms_charge_residual < 2.0e-13,
        "gcs_charge_residual": gcs_charge_residual,
        "gcs_charge_pass": gcs_charge_residual < 2.0e-13,
        "gcs_voltage_residual": gcs_voltage_residual,
        "gcs_voltage_pass": gcs_voltage_residual < 1.0e-14,
        "gcs_series_capacitance_pass": gcs_series_capacitance_pass,
        "frumkin_phi1_residual": frumkin_phi1_residual,
        "frumkin_phi1_pass": frumkin_phi1_residual < 1.0e-14,
        "frumkin_factor_residual": frumkin_factor_residual,
        "frumkin_factor_pass": frumkin_factor_residual < 1.0e-13,
        "positive_finite_pass": bool(
            np.all(np.isfinite(positive_finite_values))
            and np.all(positive_finite_values > 0.0)
        ),
    }
    return (module04_validation,)


@app.cell
def _(mo, module04_validation):
    def _status(passed):
        return "PASS" if passed else "CHECK"

    mo.md(
        rf"""
        ## Numerical sanity checks

        | physical question | status | numerical result |
        |---|---:|---:|
        | do the Gouy–Chapman concentrations obey the Boltzmann law? | {_status(module04_validation['gc_boltzmann_pass'])} | max dimensionless residual {module04_validation['gc_boltzmann_residual']:.2e} |
        | do chemical and electrical energies cancel at equilibrium? | {_status(module04_validation['gc_electrochemical_pass'])} | max residual {module04_validation['gc_electrochemical_residual_ev']:.2e} eV/defect |
        | does the exact potential satisfy the lecture's tanh profile? | {_status(module04_validation['gc_solution_pass'])} | max residual {module04_validation['gc_solution_residual']:.2e} |
        | does the surface field give the same core charge by Gauss's law? | {_status(module04_validation['gc_gauss_pass'])} | relative mismatch {module04_validation['gc_gauss_residual']:.2e} |
        | is \(C_d\) really the slope \(dQ_{{\rm core}}/d\phi\)? | {_status(module04_validation['gc_capacitance_pass'])} | relative mismatch {module04_validation['gc_capacitance_residual']:.2e} |
        | does the Mott–Schottky parabola meet both boundary potentials? | {_status(module04_validation['ms_boundary_pass'])} | max residual {max(module04_validation['ms_surface_residual'], module04_validation['ms_far_residual']):.2e} V |
        | does its curvature match the frozen charge in Poisson's equation? | {_status(module04_validation['ms_poisson_pass'])} | relative mismatch {module04_validation['ms_poisson_residual']:.2e} |
        | does \(Q_{{\rm core}}=ze c_{{i,\infty}}\lambda\) agree with Gauss's law? | {_status(module04_validation['ms_charge_pass'])} | relative mismatch {module04_validation['ms_charge_residual']:.2e} |
        | is the mobile-defect electrochemical potential still flat? | {_status(module04_validation['ms_electrochemical_pass'])} | max residual {module04_validation['ms_electrochemical_residual_ev']:.2e} eV/defect |
        | do the Stern and diffuse layers carry the same charge? | {_status(module04_validation['gcs_charge_pass'])} | relative mismatch {module04_validation['gcs_charge_residual']:.2e} |
        | do the two potential drops add to \(\phi_0\)? | {_status(module04_validation['gcs_voltage_pass'])} | residual {module04_validation['gcs_voltage_residual']:.2e} V |
        | is a series capacitance no larger than either layer? | {_status(module04_validation['gcs_series_capacitance_pass'])} | series-capacitance ordering |
        | does the kinetics use the same reaction-plane \(\phi_1\) as GCS? | {_status(module04_validation['frumkin_phi1_pass'])} | residual {module04_validation['frumkin_phi1_residual']:.2e} V |
        | does the plotted Frumkin factor match its analytical expression? | {_status(module04_validation['frumkin_factor_pass'])} | log10 residual {module04_validation['frumkin_factor_residual']:.2e} |
        | are all characteristic lengths and capacitances positive and finite? | {_status(module04_validation['positive_finite_pass'])} | physical numerical values |

        **Why check these?** Each row protects one link in the reader's physical
        chain. The first five connect Boltzmann redistribution, the exact
        Gouy–Chapman profile, Gauss's law, and differential capacitance. The next
        four test the frozen-dopant Mott–Schottky approximation. The final five
        verify that the Stern and diffuse layers share one charge, split one
        voltage, and pass the same reaction-plane potential into the Frumkin
        kinetics.
        """
    )
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## Take-home map

    1. A different defect formation energy can charge an interface core.
    2. A flat \(\widetilde\mu_i\) gives the Boltzmann concentration profile, and
       Poisson's equation makes concentration, charge, and potential
       self-consistent.
    3. In Gouy–Chapman, all charged defects move. The exact nonlinear profile
       reduces to an exponential only when \(|ze\phi_0|\ll k_BT\).
    4. In Mott–Schottky, frozen dopants provide nearly constant depletion charge,
       so the potential is parabolic over a finite width \(\lambda\).
    5. \(\lambda_D\) is independent of core charge, whereas the
       Mott–Schottky \(\lambda\) grows with \(\phi_0\) or \(Q_{\rm core}\).
    6. Gauss's law gives \(Q_{\rm core}(\phi)\), and its derivative is the
       differential capacitance.
    7. The Stern and diffuse layers carry the same charge and act as two
       capacitances in series.
    8. A reaction at \(x=x_1\) experiences both the local concentration and the
       local potential. Their combined change is the **Frumkin effect**.

    **Model boundary.** The geometry is planar and one-dimensional. The
    Gouy–Chapman section uses an ideal symmetric \(+ze/-ze\) mobile pair. The
    Mott–Schottky section uses the depletion approximation with frozen negative
    dopants. The Stern layer is represented by a constant capacitance. The
    Frumkin plot aligns formal potential and pZC only to provide one transparent
    voltage coordinate. Continuum theory becomes qualitative when a screening
    length approaches atomic dimensions. Specific adsorption, finite-site
    crowding beyond Stern's closest-approach picture, and Marcus charge-transfer
    theory are deliberately left for later modules.
    """)
    return


if __name__ == "__main__":
    app.run()
