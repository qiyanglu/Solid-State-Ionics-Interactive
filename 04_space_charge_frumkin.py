# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "marimo>=0.23.14",
#     "matplotlib>=3.8",
#     "numpy>=1.26",
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

    return mo, np, plt


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
    model_details = mo.md(r"""
    # Space-Charge Layers and the Frumkin Effect

    **Guiding question.** How does an interface redistribute charged defects,
    and why does that redistribution change both capacitance and reaction rate?

    **Learning goals**

    1. Build a space-charge profile from electrochemical equilibrium and
       Poisson's equation.
    2. Distinguish the Gouy–Chapman and Mott–Schottky limiting pictures.
    3. Follow one interfacial charge through GCS capacitance to the Frumkin
       correction at the reaction plane.

    > **Predict before exploring.** If the bulk defect concentration increases,
    > should the space-charge layer extend farther into the crystal or become
    > more compact? What must happen to its capacitance?

    **Reader path and scope.** Sections 1–2 are the core space-charge lesson.
    Sections 3–4 are the interface extension: capacitance and reaction-plane
    kinetics. The geometry is planar and one-dimensional; $\phi_0$ is a boundary
    condition supplied by the charged core. See the shared
    [notation bridge](https://github.com/qiyanglu/Solid-State-Ionics-Interactive/blob/main/NOTATION.md).

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
    \rightarrow C_{\rm sc}
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
    | \(C_{\rm sc}\), \(C_s\) | diffuse-layer and Stern-layer differential capacitances per area |

    The bulk is the reference, \(\phi_\infty=0\), and the main profile sections
    use a positively charged core, as in the slides. Inputs are shown in
    **K, V, cm⁻³, nm, and μF/cm²**; calculations use SI internally. We use
    \(k_B,e\) for energies per defect and \(R,F\) for the equivalent molar
    kinetic equations.
    """)
    mo.vstack([
        mo.md(r"""
        # Space-charge layers: screening a charged interface

        **How does a charged interface rearrange defects in the nearby solid?**

        A charged core creates an electrostatic potential (phi(x)). Mobile
        charged defects respond through their electrochemical potential,

        $$\widetilde\mu_i=\mu_i^0+k_BT\ln c_i+z_i e\phi.$$

        At equilibrium this sum is spatially constant. Concentration therefore
        changes in exactly the way needed to oppose the electrical-energy
        change. We compare two one-dimensional limits:

        - **Gouy–Chapman:** positive and negative defects are both mobile.
        - **Mott–Schottky:** majority dopants are frozen and mobile defects are depleted.

        The core reader ends with the profile and screening length. GCS
        capacitance and the Frumkin reaction-plane correction remain available
        as advanced continuations.
        """),
        mo.accordion({"Model details — geometry, notation, units, and reader map": model_details}),
    ])
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
        log10_potential_factor = (
            -(1.0 - alpha) * molar_scale * phi1_values / np.log(10.0)
        )
        log10_concentration_factor = (
            -z_reactant * molar_scale * phi1_values / np.log(10.0)
        )
        log10_naive_current = (
            (1.0 - alpha) * molar_scale * potentials / np.log(10.0)
        )
        log10_total_frumkin_factor = (
            log10_potential_factor + log10_concentration_factor
        )
        log10_corrected_current = (
            log10_naive_current + log10_total_frumkin_factor
        )
        return {
            "phi1_v": phi1_values,
            "stern_drop_v": potentials - phi1_values,
            "log10_local_ratio": log10_concentration_factor,
            "log10_potential_factor": log10_potential_factor,
            "log10_concentration_factor": log10_concentration_factor,
            "log10_total_frumkin_factor": log10_total_frumkin_factor,
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
    profile_model_control = mo.ui.dropdown(
        options=["Gouy-Chapman", "Mott-Schottky", "Compare both"],
        value="Compare both",
        label="Space-charge model",
    )
    temperature_control = mo.ui.slider(
        start=300, stop=1200, step=25, value=800,
        label="Temperature (K)", show_value=True,
    )
    log_concentration_control = mo.ui.slider(
        start=16.0, stop=21.0, step=0.25, value=18.0,
        label="Bulk concentration exponent", show_value=True,
    )
    epsilon_r_control = mo.ui.slider(
        start=5, stop=300, step=5, value=100,
        label="Relative permittivity", show_value=True,
    )
    charge_control = mo.ui.slider(
        start=1, stop=3, step=1, value=1,
        label="Defect charge magnitude", show_value=True,
    )
    gc_surface_potential_control = mo.ui.slider(
        start=-0.50, stop=0.50, step=0.01, value=0.16,
        label="Gouy-Chapman core potential (V)", show_value=True,
    )
    ms_surface_potential_control = mo.ui.slider(
        start=0.01, stop=0.50, step=0.01, value=0.16,
        label="Mott-Schottky depletion potential (V)", show_value=True,
    )
    show_linear_control = mo.ui.checkbox(
        value=False,
        label="Show the low-potential Gouy-Chapman approximation",
    )
    return (
        charge_control,
        epsilon_r_control,
        gc_surface_potential_control,
        log_concentration_control,
        ms_surface_potential_control,
        profile_model_control,
        show_linear_control,
        temperature_control,
    )


@app.cell
def _(
    charge_control,
    epsilon_r_control,
    gc_surface_potential_control,
    log_concentration_control,
    mo,
    ms_surface_potential_control,
    profile_model_control,
    show_linear_control,
    temperature_control,
):
    derivation = mo.md(r"""
    At equilibrium,

    $$
    \frac{c_i(x)}{c_{i,\infty}}
    =\exp\!\left[-\frac{z_i e\phi(x)}{k_BT}\right],
    $$

    and the redistributed charge bends the potential through Poisson's equation,

    $$
    \rho(x)=\sum_i z_i e c_i(x),\qquad
    \frac{d^2\phi}{dx^2}=-\frac{\rho(x)}{\epsilon_0\epsilon_r}.
    $$

    Potential changes concentration, while the resulting charge density changes
    potential. The exact Gouy-Chapman solution and the Mott-Schottky depletion
    approximation close this loop in different ways.
    """)
    active_potential_control = (
        gc_surface_potential_control
        if profile_model_control.value == "Gouy-Chapman"
        else ms_surface_potential_control
    )
    core_profile_controls = mo.hstack(
        [profile_model_control, active_potential_control, log_concentration_control],
        justify="start", align="center", wrap=True, gap=1.2,
    )
    advanced_profile_controls = mo.hstack(
        [temperature_control, epsilon_r_control, charge_control, show_linear_control],
        justify="start", align="center", wrap=True, gap=1.2,
    )
    mo.vstack([
        mo.md(r"""
        ## 1. Compare two ways to screen the same positive core

        Choose a model, core potential, and bulk concentration. Ask first which
        species can move. Then look for the consequences in both the potential
        profile and the defect redistribution.
        """),
        core_profile_controls,
        mo.md(r"""
        The concentration control sets
        $c_{i,\infty}=10^x\,\mathrm{cm^{-3}}$. Gouy–Chapman permits either sign
        of $\phi_0$; the Mott–Schottky depletion approximation uses only the
        positive branch, $0.01\leq\phi_0\leq0.50\ \mathrm{V}$.
        """),
        mo.accordion({
            "Explore further — temperature, permittivity, charge, and linear guide":
            advanced_profile_controls,
            "Model details — Boltzmann distribution and Poisson equation": derivation,
        }),
    ])
    return

@app.cell
def _(
    E_CHARGE_C,
    KB_J_PER_K,
    charge_control,
    debye_length_m,
    epsilon_r_control,
    gc_surface_potential_control,
    gouy_chapman_profile,
    log_concentration_control,
    mott_schottky_profile,
    mott_schottky_width_m,
    ms_surface_charge_c_per_m2,
    ms_surface_potential_control,
    np,
    profile_model_control,
    temperature_control,
):
    temperature_k = float(temperature_control.value)
    bulk_concentration_cm3 = 10.0 ** float(log_concentration_control.value)
    epsilon_r = float(epsilon_r_control.value)
    charge_magnitude = float(charge_control.value)
    ms_surface_potential_v = float(ms_surface_potential_control.value)
    surface_potential_v = (
        float(gc_surface_potential_control.value)
        if profile_model_control.value == "Gouy-Chapman"
        else ms_surface_potential_v
    )

    selected_debye_length_m = debye_length_m(
        temperature_k,
        bulk_concentration_cm3,
        epsilon_r,
        charge_magnitude,
    )
    selected_ms_width_m = mott_schottky_width_m(
        ms_surface_potential_v,
        bulk_concentration_cm3,
        epsilon_r,
        charge_magnitude,
    )
    selected_ms_core_charge_c_per_m2 = ms_surface_charge_c_per_m2(
        ms_surface_potential_v,
        bulk_concentration_cm3,
        epsilon_r,
        charge_magnitude,
    )
    surface_reduced_potential = (
        charge_magnitude * E_CHARGE_C * surface_potential_v
        / (KB_J_PER_K * temperature_k)
    )
    ms_surface_mobile_ratio = np.exp(
        -charge_magnitude * E_CHARGE_C * ms_surface_potential_v
        / (KB_J_PER_K * temperature_k)
    )

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
        ms_surface_potential_v,
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
        ms_surface_potential_v,
        selected_debye_length_m,
        selected_ms_core_charge_c_per_m2,
        selected_ms_width_m,
        surface_potential_v,
        surface_reduced_potential,
        temperature_k,
    )

@app.cell
def _(
    gc_distance_m,
    gc_profile,
    mo,
    ms_distance_m,
    ms_profile,
    plt,
    profile_model_control,
    selected_debye_length_m,
    selected_ms_width_m,
    surface_reduced_potential,
):
    model = profile_model_control.value
    show_gc = model in ("Gouy-Chapman", "Compare both")
    show_ms = model in ("Mott-Schottky", "Compare both")
    gc_x_nm = gc_distance_m * 1.0e9
    ms_x_nm = ms_distance_m * 1.0e9

    profile_figure, (potential_axis, concentration_axis) = plt.subplots(
        1, 2, figsize=(12.8, 4.8), dpi=120, constrained_layout=True
    )
    if show_gc:
        potential_axis.plot(
            gc_x_nm, gc_profile["potential_v"], color="#4C7C86", lw=1.9,
            label="Gouy-Chapman",
        )
        potential_axis.axvline(
            selected_debye_length_m * 1.0e9, color="#4C7C86", lw=1.1, ls=":"
        )
        concentration_axis.semilogy(
            gc_x_nm, gc_profile["positive_ratio"], color="#B8734A", lw=1.8,
            label=r"GC positive co-ion, $c_+/c_{i,\infty}$",
        )
        concentration_axis.semilogy(
            gc_x_nm, gc_profile["negative_ratio"], color="#4C7C86", lw=1.8, ls="--",
            label=r"GC negative counter-ion, $c_-/c_{i,\infty}$",
        )
    if show_ms:
        potential_axis.plot(
            ms_x_nm, ms_profile["potential_v"], color="#7C6A91", lw=1.9,
            ls="--" if show_gc else "-", label="Mott-Schottky",
        )
        potential_axis.axvline(
            selected_ms_width_m * 1.0e9, color="#7C6A91", lw=1.1, ls=":"
        )
        concentration_axis.semilogy(
            ms_x_nm, ms_profile["positive_ratio"], color="#B8734A", lw=1.8,
            ls="-." if show_gc else "-", label=r"MS mobile positive defect",
        )
        concentration_axis.semilogy(
            ms_x_nm, ms_profile["negative_ratio"], color="#7C6A91", lw=1.7, ls=":",
            label=r"MS frozen negative dopant",
        )

    potential_axis.set(
        xlabel=r"Distance from core, $x$ (nm)",
        ylabel=r"Electrostatic potential, $\phi$ (V)",
        title="How far does the core potential extend?",
    )
    concentration_axis.axhline(1.0, color="#858B90", lw=1.0, ls=":")
    concentration_axis.set(
        xlabel=r"Distance from core, $x$ (nm)",
        ylabel=r"Concentration / bulk concentration",
        title="Which charged species rearrange?",
    )
    for axis in (potential_axis, concentration_axis):
        axis.grid(True, which="both", alpha=0.22)
        axis.legend(frameon=False, loc="best")
    plt.close(profile_figure)

    explanation = {
        "Gouy-Chapman": (
            "Both signs are mobile: positive co-ions leave the core region while "
            "negative counter-ions accumulate. The dotted marker is the Debye length."
        ),
        "Mott-Schottky": (
            "The negative dopants stay fixed. Depleting the mobile positive defects "
            "exposes nearly uniform negative charge, producing a parabolic potential."
        ),
        "Compare both": (
            "The same positive core is screened by two different charge inventories. "
            "GC continuously redistributes both signs; MS exposes frozen dopants over a depletion width."
        ),
    }[model]
    if show_gc and abs(surface_reduced_potential) > 5.0:
        validity_note = mo.callout(
            mo.md(r"""
            **Model caution.** Here $|ze\phi_0/(k_BT)|>5$. Ideal dilute
            Poisson–Boltzmann theory can then predict unrealistically large
            enrichment because it neglects finite site density, non-ideal
            activities, defect association, and field-dependent permittivity.
            """),
            kind="warn",
        )
    else:
        validity_note = mo.md("")
    mo.vstack([
        profile_figure,
        mo.md(
            explanation
            + f" Selected lengths: $\\lambda_D={selected_debye_length_m * 1.0e9:.2f}$ nm "
            + f"and $\\lambda={selected_ms_width_m * 1.0e9:.2f}$ nm."
        ),
        validity_note,
    ])
    return (profile_figure,)

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
        color="#4C7C86",
        lw=1.9,
        label="Exact Gouy-Chapman",
    )
    if bool(show_linear_control.value):
        gc_axes[0].plot(
            gc_distance_nm,
            gc_profile["linear_potential_v"],
            color="#C49345",
            lw=1.7,
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
        xlabel=r"Distance from core, $x$ (nm)",
        ylabel=r"Electrostatic potential, $\phi$ (V)",
        title="Potential is screened into the bulk",
    )
    gc_axes[0].grid(alpha=0.22)
    gc_axes[0].legend(frameon=False, fontsize=10)

    gc_axes[1].semilogy(
        gc_distance_nm,
        gc_profile["positive_ratio"],
        color="#B65C4A",
        lw=1.9,
        label=r"$c_+/c_{i,\infty}$",
    )
    gc_axes[1].semilogy(
        gc_distance_nm,
        gc_profile["negative_ratio"],
        color="#4C7C86",
        lw=1.9,
        label=r"$c_-/c_{i,\infty}$",
    )
    gc_axes[1].axhline(1.0, color="#858B90", lw=1.0, ls=":")
    gc_axes[1].set(
        xlabel=r"Distance from core, $x$ (nm)",
        ylabel="Concentration / bulk concentration",
        title="Co-ions deplete; counter-ions accumulate",
    )
    gc_axes[1].grid(which="both", alpha=0.22)
    gc_axes[1].legend(frameon=False)

    gc_axes[2].plot(
        gc_distance_nm,
        1.0e3 * gc_profile["chemical_positive_ev"],
        color="#B65C4A",
        lw=1.9,
        label=r"chemical: $\Delta\mu_+$",
    )
    gc_axes[2].plot(
        gc_distance_nm,
        1.0e3 * gc_profile["electrical_positive_ev"],
        color="#4C7C86",
        lw=1.9,
        label=r"electrical: $+ze\phi$",
    )
    gc_axes[2].plot(
        gc_distance_nm,
        1.0e3 * gc_profile["electrochemical_positive_ev"],
        color="#40464D",
        lw=1.8,
        ls="--",
        label=r"sum: $\Delta\widetilde\mu_+$",
    )
    gc_axes[2].set(
        xlabel=r"Distance from core, $x$ (nm)",
        ylabel="Energy change (meV per defect)",
        title="Electrochemical potential stays flat",
    )
    gc_axes[2].grid(alpha=0.22)
    gc_axes[2].legend(frameon=False, fontsize=10)
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
    mo.accordion({"Explore further — exact Gouy-Chapman profiles and electrochemical potential": mo.vstack([
        gc_explanation,
        gc_figure,
        mo.md(r"""
        Mobile co-ions and counter-ions redistribute
        until their chemical and electrical energy changes cancel everywhere.
        The Debye–Hückel exponential is only the low-potential guide; the plotted
        nonlinear profile is the exact planar Gouy–Chapman result.
        """),
    ])})
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
        color="#4C7C86",
        lw=1.9,
    )
    ms_axes[0].axvline(
        selected_ms_width_m * 1.0e9,
        color="#666666",
        lw=1.3,
        ls=":",
        label=r"depletion width $\lambda$",
    )
    ms_axes[0].set(
        xlabel=r"Distance from core, $x$ (nm)",
        ylabel=r"Electrostatic potential, $\phi$ (V)",
        title="Frozen charge gives a parabola",
    )
    ms_axes[0].grid(alpha=0.22)
    ms_axes[0].legend(frameon=False)

    ms_axes[1].semilogy(
        ms_distance_nm,
        ms_profile["positive_ratio"],
        color="#B65C4A",
        lw=1.9,
        label=r"mobile $c_+/c_{i,\infty}$",
    )
    ms_axes[1].semilogy(
        ms_distance_nm,
        ms_profile["negative_ratio"],
        color="#4C7C86",
        lw=1.8,
        ls="--",
        label=r"frozen $c_-/c_{i,\infty}$",
    )
    ms_axes[1].set(
        xlabel=r"Distance from core, $x$ (nm)",
        ylabel="Concentration / bulk concentration",
        title="Mobile positive defects are depleted",
    )
    ms_axes[1].grid(which="both", alpha=0.22)
    ms_axes[1].legend(frameon=False)

    ms_axes[2].plot(
        ms_distance_nm,
        1.0e3 * ms_profile["chemical_positive_ev"],
        color="#B65C4A",
        lw=1.9,
        label=r"chemical: $\Delta\mu_+$",
    )
    ms_axes[2].plot(
        ms_distance_nm,
        1.0e3 * ms_profile["electrical_positive_ev"],
        color="#4C7C86",
        lw=1.9,
        label=r"electrical: $+ze\phi$",
    )
    ms_axes[2].plot(
        ms_distance_nm,
        1.0e3 * ms_profile["electrochemical_positive_ev"],
        color="#40464D",
        lw=1.8,
        ls="--",
        label=r"sum: $\Delta\widetilde\mu_+$",
    )
    ms_axes[2].set(
        xlabel=r"Distance from core, $x$ (nm)",
        ylabel="Energy change (meV per defect)",
        title="The same equilibrium cancellation remains",
    )
    ms_axes[2].grid(alpha=0.22)
    ms_axes[2].legend(frameon=False, fontsize=10)
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
    mo.accordion({"Explore further — Mott-Schottky depletion and electrochemical potential": mo.vstack([
        ms_explanation,
        ms_figure,
        mo.md(r"""
        Frozen charge produces a parabolic potential, but
        the mobile positive-defect electrochemical potential still remains flat
        at equilibrium.
        """),
    ])})
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
    ms_surface_potential_v,
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
                ms_surface_potential_v,
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
                ms_surface_potential_v,
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
        color="#4C7C86",
        lw=1.9,
        label=r"Gouy-Chapman $\lambda_D$",
    )
    width_axes[0].loglog(
        concentration_sweep_cm3,
        ms_sweep_nm,
        color="#B65C4A",
        lw=1.9,
        label=r"Mott-Schottky $\lambda$",
    )
    width_axes[0].axvline(
        bulk_concentration_cm3,
        color="#C49345",
        lw=1.5,
        ls="--",
        label="selected bulk concentration",
    )
    width_axes[0].set(
        xlabel=r"Bulk concentration, $c_{i,\infty}$ (cm$^{-3}$)",
        ylabel="Space-charge length (nm)",
        title="More defects screen over a shorter distance",
    )
    width_axes[0].grid(which="both", alpha=0.22)
    width_axes[0].legend(frameon=False)

    width_axes[1].loglog(
        epsilon_sweep,
        debye_epsilon_nm,
        color="#4C7C86",
        lw=1.9,
        label=r"Gouy-Chapman $\lambda_D$",
    )
    width_axes[1].loglog(
        epsilon_sweep,
        ms_epsilon_nm,
        color="#B65C4A",
        lw=1.9,
        label=r"Mott-Schottky $\lambda$",
    )
    width_axes[1].axvline(
        epsilon_r,
        color="#C49345",
        lw=1.5,
        ls="--",
        label="selected permittivity",
    )
    width_axes[1].set(
        xlabel=r"Relative permittivity, $\epsilon_r$",
        ylabel="Space-charge length (nm)",
        title="Larger permittivity broadens the layer",
    )
    width_axes[1].grid(which="both", alpha=0.22)
    width_axes[1].legend(frameon=False)
    width_figure.tight_layout()
    plt.close(width_figure)

    mo.accordion({"Explore further — how concentration and permittivity set screening length": mo.vstack(
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

            **Prediction.** Concentration and permittivity set the spatial
            screening scale; only the Mott–Schottky depletion width also grows
            explicitly with the imposed core potential.
            """),
            width_figure,
            mo.md(r"""
            Raising the bulk concentration shortens both
            screening lengths, while raising permittivity broadens them. Only
            the Mott–Schottky width also grows with the imposed core potential.
            """),
        ]
    )})
    return (width_figure,)


@app.cell
def _(mo):
    stern_capacitance_control = mo.ui.slider(
        start=5,
        stop=80,
        step=1,
        value=20,
        label="Stern capacitance (microfarads per square centimeter)",
        show_value=True,
    )
    capacitance_view_control = mo.ui.dropdown(
        options=["Space-charge capacitance only", "Full GCS series capacitance"],
        value="Space-charge capacitance only",
        label="Capacitance view",
    )
    interface_controls = mo.hstack(
        [stern_capacitance_control, capacitance_view_control],
        justify="start",
        align="center",
        wrap=True,
        gap=1.2,
    )

    transfer_coefficient_control = mo.ui.slider(
        start=0.1,
        stop=0.9,
        step=0.05,
        value=0.5,
        label="Butler-Volmer transfer coefficient, alpha",
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
        label="Signed reactant charge number",
    )
    kinetic_controls = mo.hstack(
        [transfer_coefficient_control, reactant_charge_control],
        justify="start",
        align="center",
        wrap=True,
        gap=1.2,
    )
    return (
        capacitance_view_control,
        interface_controls,
        kinetic_controls,
        reactant_charge_control,
        stern_capacitance_control,
        transfer_coefficient_control,
    )


@app.cell
def _(interface_controls, mo):
    _gcs_intro = mo.vstack(
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
            C_{\rm sc}=\frac{dQ_{\rm core}}{d\phi_1}
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
            \frac{1}{C_{\rm tot}}=\frac{1}{C_s}+\frac{1}{C_{\rm sc}(\phi_1)}.
            \]

            At the pZC, \(\phi_0=\phi_1=0\), so

            \[
            C_{\rm sc,pZC}=\frac{\epsilon_0\epsilon_r}{\lambda_D},
            \qquad
            C_{\rm tot,pZC}
            =\left(\frac{1}{C_s}+\frac{1}{C_{\rm sc,pZC}}\right)^{-1}.
            \]

            Thus the pZC capacitance is close to \(C_s\) only when
            \(C_{\rm sc,pZC}\gg C_s\). In the ideal Gouy-Chapman model,
            \(C_{\rm sc}\) grows with \(|\phi_1|\), so the total capacitance approaches
            \(C_s\) at large field. The profile drawing places \(x_1\) at 0.5 nm
            only to make the two regions visible. The physical compact-layer
            input is \(C_s\); an independent thickness would require a separate
            Stern permittivity.
            """),
            interface_controls,
        ]
    )
    mo.accordion({"Advanced — Gouy–Chapman–Stern capacitance": _gcs_intro})
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
    gcs_potential_sweep_v = np.linspace(-0.50, 0.50, 401)
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
    bare_gc_capacitance_sweep = gc_differential_capacitance_f_per_m2(
        gcs_potential_sweep_v,
        temperature_k,
        bulk_concentration_cm3,
        epsilon_r,
        charge_magnitude,
    )
    pzc_diffuse_capacitance_f_per_m2 = gc_differential_capacitance_f_per_m2(
        0.0,
        temperature_k,
        bulk_concentration_cm3,
        epsilon_r,
        charge_magnitude,
    )
    pzc_total_capacitance_f_per_m2 = 1.0 / (
        1.0 / stern_capacitance_f_per_m2
        + 1.0 / pzc_diffuse_capacitance_f_per_m2
    )

    return (
        bare_gc_capacitance_selected,
        bare_gc_capacitance_sweep,
        gcs_charge_sweep_c_per_m2,
        gcs_diffuse_capacitance_sweep,
        gcs_diffuse_distance_m,
        gcs_diffuse_profile,
        gcs_phi1_sweep_v,
        gcs_potential_sweep_v,
        gcs_states,
        gcs_total_capacitance_sweep,
        pzc_diffuse_capacitance_f_per_m2,
        pzc_total_capacitance_f_per_m2,
        selected_gcs_state,
        stern_capacitance_f_per_m2,
        stern_distance_m,
        stern_potential_v,
        stern_thickness_m,
    )


@app.cell
def _(
    bare_gc_capacitance_sweep,
    capacitance_view_control,
    gcs_diffuse_capacitance_sweep,
    gcs_diffuse_distance_m,
    gcs_diffuse_profile,
    gcs_potential_sweep_v,
    gcs_total_capacitance_sweep,
    mo,
    plt,
    pzc_diffuse_capacitance_f_per_m2,
    pzc_total_capacitance_f_per_m2,
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
        color="#C49345",
        lw=2.0,
        label="Stern layer",
    )
    gcs_axes[0].plot(
        (stern_thickness_m + gcs_diffuse_distance_m) * 1.0e9,
        gcs_diffuse_profile["potential_v"],
        color="#4C7C86",
        lw=2.0,
        label="Space-charge layer",
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
        color=["#B65C4A", "#4C7C86"],
        edgecolor="#40464D",
        zorder=5,
    )
    gcs_axes[0].set(
        xlabel=r"Distance from core, $x$ (nm)",
        ylabel=r"Electrostatic potential, $\phi$ (V)",
        title=r"The total drop splits into $\phi_0-\phi_1$ and $\phi_1$",
    )
    gcs_axes[0].grid(alpha=0.22)
    gcs_axes[0].legend(frameon=False)

    if capacitance_view_control.value == "Space-charge capacitance only":
        gcs_axes[1].plot(
            gcs_potential_sweep_v,
            bare_gc_capacitance_sweep / 0.01,
            color="#4C7C86",
            lw=2.0,
            label=r"$C_{\rm sc}(\phi_1)$",
        )
        gcs_axes[1].scatter(
            [0.0], [pzc_diffuse_capacitance_f_per_m2 / 0.01],
            s=75, color="#4C7C86", edgecolor="white", zorder=5,
            label="pZC",
        )
        gcs_axes[1].set(
            xlabel=r"Space-charge potential, $\phi_1$ (V)",
            ylabel=r"$C_{\rm sc}$ ($\mu$F cm$^{-2}$)",
            title=r"At pZC, $C_{\rm sc}=\epsilon/\lambda_D$",
        )
        view_takeaway = (
            "This view isolates the Gouy–Chapman space-charge response. "
            "It does not yet include the compact Stern layer."
        )
    else:
        gcs_axes[1].plot(
            gcs_potential_sweep_v,
            gcs_diffuse_capacitance_sweep / 0.01,
            color="#4C7C86",
            lw=1.8,
            label=r"$C_{\rm sc}[\phi_1(\phi_0)]$",
        )
        gcs_axes[1].axhline(
            stern_capacitance_f_per_m2 / 0.01,
            color="#C49345",
            lw=1.7,
            ls="--",
            label=r"$C_s$",
        )
        gcs_axes[1].plot(
            gcs_potential_sweep_v,
            gcs_total_capacitance_sweep / 0.01,
            color="#B65C4A",
            lw=2.0,
            label=r"$C_{\rm tot}$",
        )
        gcs_axes[1].scatter(
            [0.0], [pzc_total_capacitance_f_per_m2 / 0.01],
            s=80, color="#B65C4A", edgecolor="white", zorder=5,
            label="total at pZC",
        )
        gcs_axes[1].set(
            xlabel=r"Applied core potential, $\phi_0$ (V)",
            ylabel=r"Differential capacitance ($\mu$F cm$^{-2}$)",
            title="Compact and space-charge layers add in series",
        )
        view_takeaway = (
            "The horizontal axis is the applied core potential $\\phi_0$; "
            "the space-charge curve is evaluated at the self-consistent "
            "reaction-plane potential $\\phi_1(\\phi_0)$."
        )
    gcs_axes[1].grid(alpha=0.22)
    gcs_axes[1].legend(frameon=False)
    gcs_figure.tight_layout()
    plt.close(gcs_figure)

    gcs_summary = mo.md(
        rf"""
        At the selected $\phi_0={surface_potential_v:.3f}$ V, the common charge
        is $Q_{{\rm core}}={selected_gcs_state['charge_c_per_m2']:.3e}$ C m$^{{-2}}$,
        with $\phi_1={selected_gcs_state['phi1_v']:.4f}$ V and
        $\phi_0-\phi_1={selected_gcs_state['stern_drop_v']:.4f}$ V.

        The selected space-charge and total capacitances are
        $C_{{\rm sc}}={selected_gcs_state['diffuse_capacitance_f_per_m2'] / 0.01:.2f}$
        and $C_{{\rm tot}}={selected_gcs_state['total_capacitance_f_per_m2'] / 0.01:.2f}$
        $\mu$F cm$^{{-2}}$. At the pZC,
        $C_{{\rm sc,pZC}}={pzc_diffuse_capacitance_f_per_m2 / 0.01:.2f}$,
        $C_s={stern_capacitance_f_per_m2 / 0.01:.2f}$, and
        $C_{{\rm tot,pZC}}={pzc_total_capacitance_f_per_m2 / 0.01:.2f}$
        $\mu$F cm$^{{-2}}$.

        The series total is near $C_{{\rm sc,pZC}}$ only when
        $C_s\gg C_{{\rm sc,pZC}}$; it is near $C_s$ only in the opposite limit.
        {view_takeaway}
        """
    )
    mo.accordion({
        "Advanced — GCS figure and interpretation": mo.vstack([gcs_figure, gcs_summary])
    })
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
    frumkin_potential_sweep_v = np.linspace(-0.50, 0.50, 401)
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
    kinetic_controls,
    mo,
    plt,
    reactant_charge_number,
    selected_frumkin,
    surface_potential_v,
    transfer_coefficient,
):
    frumkin_figure, frumkin_axes = plt.subplots(
        2,
        2,
        figsize=(13.8, 9.0),
        dpi=120,
    )
    frumkin_axes[0, 0].plot(
        frumkin_potential_sweep_v,
        frumkin_sweep["phi1_v"],
        color="#4C7C86",
        lw=1.9,
        label=r"diffuse drop $\phi_1$",
    )
    frumkin_axes[0, 0].plot(
        frumkin_potential_sweep_v,
        frumkin_sweep["stern_drop_v"],
        color="#C49345",
        lw=1.9,
        label=r"Stern drop $\phi_0-\phi_1$",
    )
    frumkin_axes[0, 0].plot(
        frumkin_potential_sweep_v,
        frumkin_potential_sweep_v,
        color="#858B90",
        lw=1.2,
        ls=":",
        label=r"total $\phi_0$",
    )
    frumkin_axes[0, 0].set(
        xlabel=r"Driving potential, $\phi_0=E-E^{0\prime}$ (V)",
        ylabel="Potential contribution (V)",
        title="GCS divides the applied potential",
    )
    frumkin_axes[0, 0].grid(alpha=0.22)
    frumkin_axes[0, 0].legend(frameon=False, fontsize=10)

    frumkin_axes[0, 1].plot(
        frumkin_potential_sweep_v,
        frumkin_sweep["log10_local_ratio"],
        color="#4C7C86",
        lw=1.9,
    )
    frumkin_axes[0, 1].axhline(0.0, color="#858B90", lw=1.0, ls=":")
    frumkin_axes[0, 1].set(
        xlabel=r"Driving potential, $\phi_0=E-E^{0\prime}$ (V)",
        ylabel=r"$\log_{10}([R]_{x_1}/[R]_{\infty})$",
        title="The reaction-plane concentration changes",
    )
    frumkin_axes[0, 1].grid(alpha=0.22)

    frumkin_axes[1, 0].plot(
        frumkin_potential_sweep_v,
        frumkin_sweep["log10_naive_current"],
        color="#858B90",
        lw=1.7,
        ls="--",
        label="bulk concentration + total potential",
    )
    frumkin_axes[1, 0].plot(
        frumkin_potential_sweep_v,
        frumkin_sweep["log10_corrected_current"],
        color="#B65C4A",
        lw=2.0,
        label="Frumkin-corrected",
    )
    frumkin_axes[1, 0].set(
        xlabel=r"Driving potential, $\phi_0=E-E^{0\prime}$ (V)",
        ylabel=r"$\log_{10}$ normalized anodic current",
        title="Current at the reaction plane",
    )
    frumkin_axes[1, 0].grid(alpha=0.22)
    frumkin_axes[1, 0].legend(frameon=False, fontsize=10)

    frumkin_axes[1, 1].plot(
        frumkin_potential_sweep_v,
        frumkin_sweep["log10_potential_factor"],
        color="#C49345",
        lw=1.9,
        label=rf"potential part, $1-\alpha={1.0 - transfer_coefficient:.2f}$",
    )
    frumkin_axes[1, 1].plot(
        frumkin_potential_sweep_v,
        frumkin_sweep["log10_concentration_factor"],
        color="#4C7C86",
        lw=1.9,
        label=rf"concentration part, $z_R={reactant_charge_number:+.0f}$",
    )
    frumkin_axes[1, 1].plot(
        frumkin_potential_sweep_v,
        frumkin_sweep["log10_total_frumkin_factor"],
        color="#B65C4A",
        lw=2.0,
        label="sum: Frumkin factor",
    )
    frumkin_axes[1, 1].axhline(0.0, color="#858B90", lw=1.0, ls=":")
    frumkin_axes[1, 1].set(
        xlabel=r"Driving potential, $\phi_0=E-E^{0\prime}$ (V)",
        ylabel=r"Contribution to $\log_{10}(I/I_{\rm naive})$",
        title=r"What $\alpha$ changes—and what it does not",
    )
    frumkin_axes[1, 1].grid(alpha=0.22)
    frumkin_axes[1, 1].legend(frameon=False, fontsize=10)
    frumkin_figure.tight_layout()
    plt.close(frumkin_figure)

    selected_local_ratio = 10.0 ** float(selected_frumkin["log10_local_ratio"][0])
    selected_current_ratio = 10.0 ** float(
        selected_frumkin["log10_total_frumkin_factor"][0]
    )
    frumkin_text = mo.md(
        rf"""
        ## 4. The Frumkin effect: the reaction sees \(x=x_1\), not the bulk

        **Butler–Volmer prerequisite and sign convention.** The plotted branch
        is anodic: increasing $E-E^{{0\prime}}$ drives oxidation in the positive
        current direction. The reactant charge $z_R$ is signed. At positive
        $\phi_1$, a positive reactant is depleted and a negative reactant is
        accumulated at the reaction plane.

        The **Butler-Volmer transfer coefficient** \(\alpha\) is a dimensionless
        kinetic parameter between 0 and 1. It describes how an interfacial
        overpotential changes the activation barriers of the two reaction
        directions. In the lecture's anodic branch, the potential sensitivity
        is proportional to \(1-\alpha\):

        \[
        I_a=FAk^0\exp\!\left[\frac{{(1-\alpha)F
        (E-E^{{0\prime}}-\phi_1)}}{{RT}}\right][R]_{{x=x_1}},
        \qquad
        \frac{{[R]_{{x=x_1}}}}{{[R]_{{x=\infty}}}}
        =\exp\!\left(-\frac{{z_R F\phi_1}}{{RT}}\right).
        \]

        The key separation is:

        - \(C_s\), the bulk concentration, and \(\epsilon_r\) determine the
          **equilibrium electrostatics** and GCS capacitance;
        - \(\alpha\) enters only the **charge-transfer kinetics**. It should not
          change the potential profile or capacitance.

        The lower-right panel isolates the \(\alpha\)-dependent potential
        contribution from the concentration contribution. Changing \(\alpha\)
        changes the kinetic potential factor and total Frumkin factor, while the
        GCS potential and capacitance remain unchanged.

        The exchange current introduced in the lecture also uses reaction-plane
        concentrations,

        \[
        I_0=FAk^0[O]_{{x_1}}^{{1-\alpha}}[R]_{{x_1}}^\alpha.
        \]

        To keep this teaching plot one-dimensional, the formal potential is set
        equal to the pZC reference, so \(E-E^{{0\prime}}=\phi_0\). This is a
        transparent reference choice, not a universal identity.

        Here \(z_R\) is a **signed** charge number. If a reaction is written
        \(R^{{-z}}\), with \(z>0\), then \(z_R=-z\). A positive \(\phi_1\)
        enriches a negative reactant and depletes a positive one. The default
        \(z_R=+1\) represents the proton-like depletion picture in the Pt/BZY
        example. Relative to a curve using bulk concentration and the full
        potential,

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

        **Figure takeaway.** The reaction rate samples both local reactant
        concentration and local potential at $x_1$; $\alpha$ changes kinetics
        but does not change the equilibrium GCS capacitance.
        """
    )
    mo.accordion({
        "Advanced — Butler–Volmer prerequisite and Frumkin correction":
        mo.vstack([frumkin_text, kinetic_controls, frumkin_figure])
    })
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
    gc_differential_capacitance_f_per_m2,
    gc_surface_charge_c_per_m2,
    gcs_diffuse_capacitance_sweep,
    gcs_phi1_sweep_v,
    gcs_state,
    gcs_total_capacitance_sweep,
    ms_distance_m,
    ms_profile,
    ms_surface_charge_c_per_m2,
    ms_surface_potential_v,
    np,
    reactant_charge_number,
    selected_debye_length_m,
    selected_frumkin,
    selected_gcs_state,
    frumkin_sweep,
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

    ms_surface_residual = abs(ms_profile["potential_v"][0] - ms_surface_potential_v)
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
    analytic_ms_curvature = 2.0 * ms_surface_potential_v / selected_ms_width_m**2
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
        2.0 * permittivity * ms_surface_potential_v / selected_ms_width_m
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
    gcs_charge_plus = gcs_state(
        surface_potential_v + capacitance_step_v,
        temperature_k,
        bulk_concentration_cm3,
        epsilon_r,
        charge_magnitude,
        stern_capacitance_f_per_m2,
    )["charge_c_per_m2"]
    gcs_charge_minus = gcs_state(
        surface_potential_v - capacitance_step_v,
        temperature_k,
        bulk_concentration_cm3,
        epsilon_r,
        charge_magnitude,
        stern_capacitance_f_per_m2,
    )["charge_c_per_m2"]
    numerical_gcs_capacitance = (
        gcs_charge_plus - gcs_charge_minus
    ) / (2.0 * capacitance_step_v)
    gcs_capacitance_residual = abs(
        numerical_gcs_capacitance
        - selected_gcs_state["total_capacitance_f_per_m2"]
    ) / selected_gcs_state["total_capacitance_f_per_m2"]
    pzc_diffuse_capacitance = permittivity / selected_debye_length_m
    calculated_pzc_space_charge_capacitance = gc_differential_capacitance_f_per_m2(
        0.0,
        temperature_k,
        bulk_concentration_cm3,
        epsilon_r,
        charge_magnitude,
    )
    gcs_pzc_space_charge_residual = abs(
        calculated_pzc_space_charge_capacitance - pzc_diffuse_capacitance
    ) / pzc_diffuse_capacitance
    expected_pzc_total_capacitance = 1.0 / (
        1.0 / stern_capacitance_f_per_m2
        + 1.0 / pzc_diffuse_capacitance
    )
    calculated_pzc_total_capacitance = gcs_state(
        0.0,
        temperature_k,
        bulk_concentration_cm3,
        epsilon_r,
        charge_magnitude,
        stern_capacitance_f_per_m2,
    )["total_capacitance_f_per_m2"]
    gcs_pzc_capacitance_residual = abs(
        calculated_pzc_total_capacitance - expected_pzc_total_capacitance
    ) / expected_pzc_total_capacitance
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
    positive_finite_values = np.concatenate([
        np.array([
            selected_debye_length_m,
            selected_ms_width_m,
            bare_gc_capacitance_selected,
            selected_gcs_state["diffuse_capacitance_f_per_m2"],
            selected_gcs_state["total_capacitance_f_per_m2"],
        ]),
        np.asarray(gcs_diffuse_capacitance_sweep),
        np.asarray(gcs_total_capacitance_sweep),
    ])
    expanded_finite_pass = bool(
        np.all(np.isfinite(gcs_phi1_sweep_v))
        and all(np.all(np.isfinite(values)) for values in frumkin_sweep.values())
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
        "gcs_capacitance_residual": gcs_capacitance_residual,
        "gcs_capacitance_pass": gcs_capacitance_residual < 2.0e-9,
        "gcs_pzc_space_charge_residual": gcs_pzc_space_charge_residual,
        "gcs_pzc_space_charge_pass": gcs_pzc_space_charge_residual < 2.0e-13,
        "gcs_pzc_capacitance_residual": gcs_pzc_capacitance_residual,
        "gcs_pzc_capacitance_pass": gcs_pzc_capacitance_residual < 2.0e-13,
        "frumkin_phi1_residual": frumkin_phi1_residual,
        "frumkin_phi1_pass": frumkin_phi1_residual < 1.0e-14,
        "frumkin_factor_residual": frumkin_factor_residual,
        "frumkin_factor_pass": frumkin_factor_residual < 1.0e-13,
        "positive_finite_pass": bool(
            expanded_finite_pass
            and np.all(np.isfinite(positive_finite_values))
            and np.all(positive_finite_values > 0.0)
        ),
    }
    return (module04_validation,)


@app.cell
def _(mo, module04_validation):
    def _status(passed):
        return "PASS" if passed else "CHECK"

    _checks = mo.md(
        rf"""
        ## Physical consistency checks

        | status | physical statement | why it matters |
        |---:|---|---|
        | {_status(module04_validation['gc_boltzmann_pass'] and module04_validation['gc_electrochemical_pass'] and module04_validation['gc_solution_pass'])} | the Gouy–Chapman concentration and potential profiles obey Boltzmann equilibrium | chemical and electrical energies balance throughout the diffuse layer |
        | {_status(module04_validation['gc_gauss_pass'] and module04_validation['gc_capacitance_pass'])} | the diffuse-layer charge agrees with Gauss's law and its slope gives \(C_{{\rm sc}}\) | charge, field, and capacitance are three views of the same interface |
        | {_status(module04_validation['ms_boundary_pass'] and module04_validation['ms_poisson_pass'] and module04_validation['ms_charge_pass'] and module04_validation['ms_electrochemical_pass'])} | the Mott–Schottky profile obeys its boundary potentials, frozen charge, and flat electrochemical potential | the depletion approximation remains self-consistent |
        | {_status(module04_validation['gcs_charge_pass'] and module04_validation['gcs_voltage_pass'] and module04_validation['gcs_series_capacitance_pass'] and module04_validation['gcs_capacitance_pass'] and module04_validation['gcs_pzc_space_charge_pass'] and module04_validation['gcs_pzc_capacitance_pass'])} | Stern and diffuse layers carry the same charge and add as series capacitances | the GCS voltage split and pZC limit agree |
        | {_status(module04_validation['frumkin_phi1_pass'] and module04_validation['frumkin_factor_pass'])} | the Frumkin response uses the same reaction-plane concentration and potential | equilibrium space charge and interfacial kinetics are connected consistently |
        | {_status(module04_validation['positive_finite_pass'])} | all lengths and capacitances are positive and finite | every plotted interfacial scale is physical |

        These checks follow the physical chain from redistribution to charge,
        capacitance, and the reaction-plane Frumkin effect.
        """
    )
    mo.accordion({"Physical consistency checks": _checks})
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## Take-home map

    1. **Equilibrium and Poisson must agree.** A flat $\widetilde\mu_i$ gives
       Boltzmann redistribution, while Poisson's equation makes concentration,
       charge, and potential self-consistent. Gouy–Chapman moves both signs;
       Mott–Schottky freezes the majority dopants.
    2. **Capacitance measures the charge response.** Gauss's law gives
       $Q_{\rm core}(\phi)$ and differentiation gives $C_{\rm sc}$. Stern and diffuse
       layers carry the same charge and add in series; near the pZC the total is
       close to $C_s$ only when $C_{\rm sc,pZC}\gg C_s$.
    3. **A reaction experiences the reaction plane.** The Frumkin effect combines
       the local concentration and local potential at $x=x_1$. The transfer
       coefficient $\alpha$ belongs to this kinetic step, not to equilibrium GCS
       capacitance.

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


@app.cell
def _(mo):
    mo.md(r"""
    **Continue:** [Module 05 — Stoichiometry Polarization](https://qiyanglu.github.io/Solid-State-Ionics-Interactive/05-stoichiometry-polarization/)
    """)
    return



if __name__ == "__main__":
    app.run()
