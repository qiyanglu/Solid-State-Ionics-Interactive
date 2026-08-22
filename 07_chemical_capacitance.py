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
    from scipy import optimize

    plt.rcParams.update(
        {
            "font.size": 15,
            "axes.titlesize": 17,
            "axes.labelsize": 15,
            "xtick.labelsize": 12,
            "ytick.labelsize": 12,
            "legend.fontsize": 11.5,
            "lines.linewidth": 2.0,
            "axes.facecolor": "#FCFCFA",
            "figure.facecolor": "white",
            "grid.color": "#C7CCD1",
            "grid.alpha": 0.25,
            "axes.titlepad": 10,
            "axes.labelpad": 7,
            "legend.frameon": False,
            "figure.dpi": 115,
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
    return mo, np, optimize, plt


@app.cell
def _(mo):
    mo.md(r"""
    # Module 07 · Chemical Capacitance: Thermodynamics, Defects, and Measurement

    **How can a solid store electrochemical charge by changing its neutral
    composition?**

    Chemical capacitance answers that question. It is not an extra electrostatic
    capacitor hidden inside a crystal. It is the differential storage associated
    with changing the concentration of a neutral chemical component.

    This reader follows one chain:

    \[
    G(c)\longrightarrow \mu(c)\longrightarrow
    \frac{\partial c}{\partial\mu}\longrightarrow C_{\rm chem}
    \longrightarrow D^\delta\text{ and a measured response}.
    \]

    **Learning goals**

    1. Distinguish dielectric charge separation from neutral chemical storage.
    2. Derive chemical capacitance from a chemical-potential curve.
    3. See how conductivity and storage combine to set chemical diffusivity.
    4. Recognize what thickness scaling, titration, and impedance can—and cannot—identify.
    """)
    return


@app.cell
def _(plt):
    _figure, _axes = plt.subplots(1, 2, figsize=(13.0, 4.2), constrained_layout=True)

    _ax = _axes[0]
    _ax.plot([0.20, 0.20], [0.15, 0.85], color="#526173", lw=8)
    _ax.plot([0.80, 0.80], [0.15, 0.85], color="#526173", lw=8)
    for _y in (0.30, 0.50, 0.70):
        _ax.text(0.14, _y, "+", color="#B8734A", fontsize=20, ha="center", va="center")
        _ax.text(0.86, _y, "−", color="#4C7C86", fontsize=20, ha="center", va="center")
    _ax.annotate("", xy=(0.68, 0.50), xytext=(0.32, 0.50),
                 arrowprops={"arrowstyle": "<->", "lw": 1.8, "color": "#7C6A91"})
    _ax.text(0.50, 0.92, "Dielectric capacitance", ha="center", weight="bold")
    _ax.text(0.50, 0.06, "opposite charge separates in space", ha="center")
    _ax.set(xlim=(0, 1), ylim=(0, 1))
    _ax.axis("off")

    _ax = _axes[1]
    _site_x = [0.18, 0.34, 0.50, 0.66, 0.82]
    for _x in _site_x:
        _ax.add_patch(plt.Circle((_x, 0.52), 0.065, facecolor="#DDE7E8", edgecolor="#4C7C86", lw=1.4))
    for _x in _site_x[:2]:
        _ax.add_patch(plt.Circle((_x, 0.52), 0.036, facecolor="#B8734A", edgecolor="white", lw=1.0))
    _ax.annotate("neutral species", xy=(0.34, 0.62), xytext=(0.56, 0.78),
                 arrowprops={"arrowstyle": "->", "color": "#B8734A"}, ha="center")
    _ax.annotate("", xy=(0.83, 0.27), xytext=(0.63, 0.27),
                 arrowprops={"arrowstyle": "->", "lw": 1.8, "color": "#5F8A6B"})
    _ax.text(0.50, 0.92, "Chemical capacitance", ha="center", weight="bold")
    _ax.text(0.50, 0.06, "neutral composition changes throughout a volume", ha="center")
    _ax.set(xlim=(0, 1), ylim=(0, 1))
    _ax.axis("off")
    plt.close(_figure)
    qualitative_capacitance_figure_07 = _figure
    return (qualitative_capacitance_figure_07,)


@app.cell
def _(mo, qualitative_capacitance_figure_07):
    mo.vstack([
        qualitative_capacitance_figure_07,
        mo.md(r"""
        A dielectric capacitor stores separated charge near interfaces. Chemical
        capacitance stores a change in **neutral composition** in an active
        volume. Real experiments may contain both, which is why a measured
        capacitance needs a physical model.
        """),
    ])
    return


@app.cell
def _(np):
    GAS_CONSTANT_J_PER_MOL_K = 8.314462618
    FARADAY_C_PER_MOL = 96485.33212

    def _positive_07(name, value):
        number = float(value)
        if not np.isfinite(number) or number <= 0.0:
            raise ValueError(f"{name} must be positive and finite")
        return number

    def dilute_storage_state_07(fraction, temperature_k, c_max_mol_per_m3, charge_number=1.0):
        """Ideal dilute storage with c = x c_max and mu/RT = ln(c/c_ref)."""
        x_value = _positive_07("fraction", fraction)
        if x_value >= 1.0:
            raise ValueError("fraction must be smaller than one")
        temperature = _positive_07("temperature_k", temperature_k)
        c_max = _positive_07("c_max_mol_per_m3", c_max_mol_per_m3)
        z_value = _positive_07("charge_number", abs(charge_number))
        concentration = x_value * c_max
        f_th = 1.0
        c_chem_v = z_value**2 * FARADAY_C_PER_MOL**2 * concentration / (
            GAS_CONSTANT_J_PER_MOL_K * temperature * f_th
        )
        return {
            "concentration_mol_per_m3": concentration,
            "mu_over_rt": np.log(x_value / 0.10),
            "f_th": f_th,
            "c_chem_v_f_per_m3": c_chem_v,
        }

    def lattice_gas_storage_state_07(fraction, temperature_k, c_max_mol_per_m3, charge_number=1.0):
        """Ideal lattice gas with finite site exclusion."""
        x_value = _positive_07("fraction", fraction)
        if x_value >= 1.0:
            raise ValueError("fraction must be smaller than one")
        temperature = _positive_07("temperature_k", temperature_k)
        c_max = _positive_07("c_max_mol_per_m3", c_max_mol_per_m3)
        z_value = _positive_07("charge_number", abs(charge_number))
        concentration = x_value * c_max
        f_th = 1.0 / (1.0 - x_value)
        c_chem_v = (
            z_value**2
            * FARADAY_C_PER_MOL**2
            * c_max
            * x_value
            * (1.0 - x_value)
            / (GAS_CONSTANT_J_PER_MOL_K * temperature)
        )
        return {
            "concentration_mol_per_m3": concentration,
            "mu_over_rt": np.log(x_value / (1.0 - x_value)),
            "f_th": f_th,
            "c_chem_v_f_per_m3": c_chem_v,
        }

    def mixed_transport_07(sigma_i_s_per_m, sigma_e_s_per_m, c_chem_v_f_per_m3):
        sigma_i = _positive_07("sigma_i_s_per_m", sigma_i_s_per_m)
        sigma_e = _positive_07("sigma_e_s_per_m", sigma_e_s_per_m)
        capacitance = _positive_07("c_chem_v_f_per_m3", c_chem_v_f_per_m3)
        sigma_amb = sigma_i * sigma_e / (sigma_i + sigma_e)
        return sigma_amb, sigma_amb / capacitance

    def slab_rc_07(length_m, area_m2, sigma_amb_s_per_m, c_chem_v_f_per_m3):
        length = _positive_07("length_m", length_m)
        area = _positive_07("area_m2", area_m2)
        sigma_amb = _positive_07("sigma_amb_s_per_m", sigma_amb_s_per_m)
        capacitance_v = _positive_07("c_chem_v_f_per_m3", c_chem_v_f_per_m3)
        resistance = length / (sigma_amb * area)
        capacitance = capacitance_v * area * length
        diffusivity = sigma_amb / capacitance_v
        return resistance, capacitance, resistance * capacitance, length**2 / diffusivity

    return (
        FARADAY_C_PER_MOL,
        GAS_CONSTANT_J_PER_MOL_K,
        dilute_storage_state_07,
        lattice_gas_storage_state_07,
        mixed_transport_07,
        slab_rc_07,
    )


@app.cell
def _(mo):
    storage_model_07 = mo.ui.dropdown(
        options=["Ideal dilute", "Ideal lattice gas"],
        value="Ideal lattice gas",
        label="Storage model",
    )
    stored_fraction_07 = mo.ui.slider(
        0.01, 0.99, value=0.20, step=0.01,
        label="Stored-site fraction, x", show_value=True,
    )
    temperature_07 = mo.ui.slider(
        300, 1200, value=800, step=25,
        label="Temperature, T (K)", show_value=True,
    )
    core_controls_07 = mo.hstack(
        [storage_model_07, stored_fraction_07, temperature_07],
        justify="start", align="center", wrap=True, gap=1.4,
    )
    return core_controls_07, storage_model_07, stored_fraction_07, temperature_07


@app.cell
def _(core_controls_07, mo):
    mo.vstack([
        mo.md(r"""
        ## 1. From chemical potential to differential storage

        If adding a neutral component changes charge by $dQ=zF V_{\rm act}\,dc$
        and voltage by $dU=d\mu/(zF)$, then

        \[
        \boxed{C_{\rm chem}=\frac{dQ}{dU}
        =z^2F^2V_{\rm act}\left(\frac{\partial c}{\partial\mu}\right)},
        \qquad
        c_{\rm chem}^V=\frac{C_{\rm chem}}{V_{\rm act}}.
        \]

        The derivative is the central idea: a shallow $\mu(c)$ curve means that
        a small voltage change stores a large composition change.
        """),
        core_controls_07,
    ])
    return


@app.cell
def _(
    FARADAY_C_PER_MOL,
    GAS_CONSTANT_J_PER_MOL_K,
    dilute_storage_state_07,
    lattice_gas_storage_state_07,
    np,
    storage_model_07,
    stored_fraction_07,
    temperature_07,
):
    c_max_mol_per_m3_07 = 1000.0
    charge_number_07 = 1.0
    selected_fraction_07 = float(stored_fraction_07.value)
    selected_temperature_k_07 = float(temperature_07.value)
    if storage_model_07.value == "Ideal dilute":
        selected_storage_state_07 = dilute_storage_state_07(
            selected_fraction_07, selected_temperature_k_07, c_max_mol_per_m3_07
        )
    else:
        selected_storage_state_07 = lattice_gas_storage_state_07(
            selected_fraction_07, selected_temperature_k_07, c_max_mol_per_m3_07
        )

    storage_fraction_grid_07 = np.linspace(0.002, 0.998, 600)
    dilute_mu_grid_07 = np.log(storage_fraction_grid_07 / 0.10)
    lattice_mu_grid_07 = np.log(
        storage_fraction_grid_07 / (1.0 - storage_fraction_grid_07)
    )
    dilute_c_normalized_07 = storage_fraction_grid_07
    lattice_c_normalized_07 = storage_fraction_grid_07 * (1.0 - storage_fraction_grid_07)
    capacitance_scale_07 = (
        FARADAY_C_PER_MOL**2
        * c_max_mol_per_m3_07
        / (GAS_CONSTANT_J_PER_MOL_K * selected_temperature_k_07)
    )
    return (
        c_max_mol_per_m3_07,
        capacitance_scale_07,
        charge_number_07,
        dilute_c_normalized_07,
        dilute_mu_grid_07,
        lattice_c_normalized_07,
        lattice_mu_grid_07,
        selected_fraction_07,
        selected_storage_state_07,
        selected_temperature_k_07,
        storage_fraction_grid_07,
    )


@app.cell
def _(
    capacitance_scale_07,
    dilute_c_normalized_07,
    dilute_mu_grid_07,
    lattice_c_normalized_07,
    lattice_mu_grid_07,
    mo,
    plt,
    selected_fraction_07,
    selected_storage_state_07,
    storage_fraction_grid_07,
    storage_model_07,
):
    _is_dilute = storage_model_07.value == "Ideal dilute"
    _mu_curve = dilute_mu_grid_07 if _is_dilute else lattice_mu_grid_07
    _c_curve = dilute_c_normalized_07 if _is_dilute else lattice_c_normalized_07
    _model_label = "ideal dilute" if _is_dilute else "ideal lattice gas"

    _figure, _axes = plt.subplots(1, 2, figsize=(13.0, 4.8), constrained_layout=True)
    _axes[0].plot(storage_fraction_grid_07, _mu_curve, color="#4C7C86")
    _axes[0].scatter(
        [selected_fraction_07], [selected_storage_state_07["mu_over_rt"]],
        s=85, color="#C49345", edgecolor="#40464D", zorder=4,
    )
    _axes[0].axhline(0.0, color="#858B90", lw=1.0, ls=":")
    _axes[0].set(
        xlabel=r"Stored-site fraction, $x$",
        ylabel=r"Reduced chemical potential, $(\mu-\mu^0)/(RT)$",
        title="How strongly does chemical potential rise?",
    )
    _axes[0].grid(alpha=0.24)

    _axes[1].plot(storage_fraction_grid_07, _c_curve, color="#B8734A")
    _axes[1].scatter(
        [selected_fraction_07],
        [selected_storage_state_07["c_chem_v_f_per_m3"] / capacitance_scale_07],
        s=85, color="#C49345", edgecolor="#40464D", zorder=4,
    )
    _axes[1].set(
        xlabel=r"Stored-site fraction, $x$",
        ylabel=r"$c_{\rm chem}^V RT/(F^2c_{\max})$",
        title="Differential storage is the inverse slope",
        ylim=(-0.02, 1.02 if _is_dilute else 0.28),
    )
    _axes[1].grid(alpha=0.24)
    plt.close(_figure)

    _summary = mo.md(
        f"""
        The selected **{_model_label}** state has
        $c={selected_storage_state_07['concentration_mol_per_m3']:.1f}$ mol m$^{{-3}}$,
        $f_{{\rm th}}={selected_storage_state_07['f_th']:.3g}$, and
        $c_{{\rm chem}}^V={selected_storage_state_07['c_chem_v_f_per_m3']:.3e}$ F m$^{{-3}}$.
        """
    )
    _items = [_figure, _summary]
    if _is_dilute and selected_fraction_07 > 0.10:
        _items.append(mo.callout(
            "The dilute model is being extended to a concentrated state. Use the lattice-gas model when finite site occupancy matters.",
            kind="warn",
        ))
    mo.vstack(_items)
    return


@app.cell
def _(mo):
    mo.md(r"""
    The same result can be written with the thermodynamic factor

    \[
    f_{\rm th}=\frac{c}{RT}\frac{\partial\mu}{\partial c},
    \qquad
    \boxed{c_{\rm chem}^V=\frac{z^2F^2c}{RT f_{\rm th}}}.
    \]

    For an ideal dilute solution, $f_{\rm th}=1$. For an ideal lattice gas,
    $f_{\rm th}=1/(1-x)$ and

    \[
    c_{\rm chem}^V=\frac{z^2F^2c_{\max}x(1-x)}{RT}.
    \]

    Site exclusion therefore makes the differential storage vanish as the host
    becomes empty **or** full. This module uses $f_{\rm th}$ for this derivative;
    no new transport factor is implied.
    """)
    return


@app.cell
def _(
    dilute_c_normalized_07,
    dilute_mu_grid_07,
    lattice_c_normalized_07,
    lattice_mu_grid_07,
    np,
    plt,
    selected_fraction_07,
    storage_fraction_grid_07,
    storage_model_07,
):
    _is_dilute = storage_model_07.value == "Ideal dilute"
    _x = storage_fraction_grid_07
    if _is_dilute:
        _g = _x * (np.log(_x / 0.10) - 1.0)
        _mu = dilute_mu_grid_07
        _capacity = dilute_c_normalized_07
    else:
        _g = _x * np.log(_x) + (1.0 - _x) * np.log(1.0 - _x)
        _mu = lattice_mu_grid_07
        _capacity = lattice_c_normalized_07
    _selected_index = int(np.argmin(np.abs(_x - selected_fraction_07)))

    _figure, _axes = plt.subplots(1, 3, figsize=(15.2, 4.5), constrained_layout=True)
    _axes[0].plot(_x, _g, color="#7C6A91")
    _axes[0].scatter([_x[_selected_index]], [_g[_selected_index]], s=75, color="#C49345", zorder=4)
    _axes[0].set(xlabel=r"$x$", ylabel=r"$G/(c_{\max}VRT)$", title=r"1 · Free energy, $G(x)$")
    _axes[1].plot(_x, _mu, color="#4C7C86")
    _axes[1].scatter([_x[_selected_index]], [_mu[_selected_index]], s=75, color="#C49345", zorder=4)
    _axes[1].axhline(0.0, color="#858B90", lw=1.0, ls=":")
    _axes[1].set(xlabel=r"$x$", ylabel=r"$(\mu-\mu^0)/(RT)$", title=r"2 · Slope, $\mu=\partial G/\partial n$")
    _axes[2].plot(_x, _capacity, color="#B8734A")
    _axes[2].scatter([_x[_selected_index]], [_capacity[_selected_index]], s=75, color="#C49345", zorder=4)
    _axes[2].set(xlabel=r"$x$", ylabel="normalized differential storage", title=r"3 · Inverse curvature, $C_{\rm chem}$")
    for _axis in _axes:
        _axis.grid(alpha=0.24)
    plt.close(_figure)
    linked_thermodynamic_figure_07 = _figure
    return (linked_thermodynamic_figure_07,)


@app.cell
def _(linked_thermodynamic_figure_07, mo):
    mo.vstack([
        mo.md("## 2. One free-energy curve produces all three views"),
        linked_thermodynamic_figure_07,
        mo.md(r"""
        Read left to right. The slope of $G$ is the chemical potential, and the
        inverse curvature controls how much composition changes for a small
        change in chemical potential. Chemical capacitance is therefore a
        thermodynamic susceptibility written in electrical units.
        """),
    ])
    return


@app.cell
def _(mo):
    log_sigma_i_07 = mo.ui.slider(
        -6.0, 2.0, value=-2.0, step=0.25,
        label="log10 ionic conductivity (S/m)", show_value=True,
    )
    log_sigma_ratio_07 = mo.ui.slider(
        -3.0, 4.0, value=2.0, step=0.25,
        label="log10 electronic/ionic conductivity", show_value=True,
    )
    transport_controls_07 = mo.hstack(
        [log_sigma_i_07, log_sigma_ratio_07],
        justify="start", align="center", wrap=True, gap=1.4,
    )
    return log_sigma_i_07, log_sigma_ratio_07, transport_controls_07


@app.cell
def _(
    log_sigma_i_07,
    log_sigma_ratio_07,
    mixed_transport_07,
    mo,
    np,
    plt,
    selected_storage_state_07,
    transport_controls_07,
):
    _sigma_i = 10.0 ** float(log_sigma_i_07.value)
    _ratio = 10.0 ** float(log_sigma_ratio_07.value)
    _sigma_e = _sigma_i * _ratio
    selected_sigma_amb_07, selected_diffusivity_m2_per_s_07 = mixed_transport_07(
        _sigma_i, _sigma_e, selected_storage_state_07["c_chem_v_f_per_m3"]
    )
    _ratio_grid = np.logspace(-3.0, 4.0, 450)
    _sigma_amb_grid = _sigma_i * _ratio_grid / (1.0 + _ratio_grid)
    _diffusivity_grid = _sigma_amb_grid / selected_storage_state_07["c_chem_v_f_per_m3"]

    _figure, _axes = plt.subplots(1, 2, figsize=(13.0, 4.6), constrained_layout=True)
    _axes[0].loglog(_ratio_grid, _sigma_amb_grid, color="#4C7C86", label=r"$\sigma_{\rm amb}$")
    _axes[0].axhline(_sigma_i, color="#858B90", ls="--", label=r"high-$\sigma_e$ limit: $\sigma_i$")
    _axes[0].scatter([_ratio], [selected_sigma_amb_07], s=75, color="#C49345", zorder=4)
    _axes[0].set(xlabel=r"$\sigma_e/\sigma_i$", ylabel=r"$\sigma_{\rm amb}$ (S m$^{-1}$)", title="Both carriers are required")
    _axes[1].loglog(_ratio_grid, _diffusivity_grid * 1.0e4, color="#B8734A")
    _axes[1].scatter([_ratio], [selected_diffusivity_m2_per_s_07 * 1.0e4], s=75, color="#C49345", zorder=4)
    _axes[1].set(xlabel=r"$\sigma_e/\sigma_i$", ylabel=r"$D^\delta$ (cm$^2$ s$^{-1}$)", title="Storage converts conductivity into diffusion")
    for _axis in _axes:
        _axis.grid(alpha=0.24, which="both")
    plt.close(_figure)

    mo.vstack([
        mo.md(r"""
        ## 3. Conductivity plus storage sets chemical diffusivity

        For a locally neutral mixed conductor,

        \[
        \sigma_{\rm amb}=\frac{\sigma_i\sigma_e}{\sigma_i+\sigma_e},
        \qquad
        \boxed{D^\delta=\frac{\sigma_{\rm amb}}{c_{\rm chem}^V}}.
        \]

        Conductivity says how easily charge moves; chemical capacitance says how
        much composition must move. Their ratio has units of diffusivity.
        """),
        transport_controls_07,
        _figure,
        mo.md(
            f"Selected: $\\sigma_{{\\rm amb}}={selected_sigma_amb_07:.3e}$ S m$^{{-1}}$ and "
            f"$D^\\delta={selected_diffusivity_m2_per_s_07 * 1.0e4:.3e}$ cm$^2$ s$^{{-1}}$. "
            "When electrons become much faster than ions, the ambipolar conductivity approaches the ionic conductivity rather than growing without limit."
        ),
    ])
    return selected_diffusivity_m2_per_s_07, selected_sigma_amb_07


@app.cell
def _(mo):
    log_thickness_um_07 = mo.ui.slider(
        -2.0, 3.0, value=1.0, step=0.10,
        label="log10 thickness (micrometers)", show_value=True,
    )
    interface_capacitance_areal_07 = mo.ui.slider(
        0.0, 200.0, value=40.0, step=5.0,
        label="Interface capacitance (microF/cm2)", show_value=True,
    )
    measurement_controls_07 = mo.hstack(
        [log_thickness_um_07, interface_capacitance_areal_07],
        justify="start", align="center", wrap=True, gap=1.4,
    )
    return interface_capacitance_areal_07, log_thickness_um_07, measurement_controls_07


@app.cell
def _(
    interface_capacitance_areal_07,
    log_thickness_um_07,
    measurement_controls_07,
    mo,
    np,
    plt,
    selected_diffusivity_m2_per_s_07,
    selected_sigma_amb_07,
    selected_storage_state_07,
    slab_rc_07,
):
    _area_m2 = 1.0e-4
    _length_um = 10.0 ** float(log_thickness_um_07.value)
    _length_m = _length_um * 1.0e-6
    _c_int_areal = float(interface_capacitance_areal_07.value)
    _length_grid_um = np.logspace(-2.0, 3.0, 350)
    _length_grid_m = _length_grid_um * 1.0e-6
    _chemical_areal = selected_storage_state_07["c_chem_v_f_per_m3"] * _length_grid_m * 100.0
    _measured_areal = _chemical_areal + _c_int_areal
    _selected_chemical_areal = selected_storage_state_07["c_chem_v_f_per_m3"] * _length_m * 100.0
    _selected_measured_areal = _selected_chemical_areal + _c_int_areal
    selected_measured_areal_07 = _selected_measured_areal
    selected_length_m_07 = _length_m
    selected_resistance_ohm_07, selected_capacitance_f_07, selected_rc_s_07, selected_diffusion_time_s_07 = slab_rc_07(
        _length_m,
        _area_m2,
        selected_sigma_amb_07,
        selected_storage_state_07["c_chem_v_f_per_m3"],
    )

    _figure, _axis = plt.subplots(figsize=(9.5, 5.0), constrained_layout=True)
    _axis.loglog(_length_grid_um, _chemical_areal, color="#4C7C86", label=r"bulk chemical: $c_{\rm chem}^V L$")
    _axis.loglog(_length_grid_um, _measured_areal, color="#B8734A", label=r"measured: $c_{\rm chem}^V L+C_{\rm int}/S$")
    if _c_int_areal > 0.0:
        _axis.axhline(_c_int_areal, color="#7C6A91", ls="--", label="interface intercept")
    _axis.scatter([_length_um], [_selected_measured_areal], s=80, color="#C49345", edgecolor="#40464D", zorder=4)
    _axis.set(
        xlabel=r"Active thickness, $L$ ($\mu$m)",
        ylabel=r"Capacitance per area ($\mu$F cm$^{-2}$)",
        title="Bulk chemical storage scales with active thickness",
    )
    _axis.grid(alpha=0.24, which="both")
    _axis.legend(loc="best")
    plt.close(_figure)

    mo.vstack([
        mo.md(r"""
        ## 4. How chemical capacitance appears in experiments

        A titration measures $dQ/dU$ directly. Impedance can also reveal a
        capacitance, but only through a model that separates bulk storage from
        interface, dielectric, and stray contributions. Thickness scaling is a
        useful test:

        \[
        C_{\rm meas}(L)=c_{\rm chem}^VSL+C_{\rm int}.
        \]
        """),
        measurement_controls_07,
        _figure,
        mo.md(
            f"At $L={_length_um:.3g}$ µm and $S=1$ cm², the model gives "
            f"$C_{{\\rm chem}}={selected_capacitance_f_07:.3e}$ F, "
            f"$R_{{\\rm amb}}={selected_resistance_ohm_07:.3e}$ Ω, and "
            f"$R_{{\\rm amb}}C_{{\\rm chem}}={selected_rc_s_07:.3e}$ s."
        ),
        mo.callout(
            mo.md(
                r"The exact slab identity $R_{\rm amb}C_{\rm chem}=L^2/D^\delta$ "
                "is the bridge to the transmission line in Module 08."
            ),
            kind="info",
        ),
    ])
    return (
        selected_capacitance_f_07,
        selected_diffusion_time_s_07,
        selected_length_m_07,
        selected_measured_areal_07,
        selected_rc_s_07,
        selected_resistance_ohm_07,
    )


@app.cell
def _(mo, np, plt):
    _eta = np.linspace(-8.0, 8.0, 500)
    _vacancy = np.exp(-0.5 * _eta)
    _hole = np.exp(0.5 * _eta)
    _storage = 1.0 / (1.0 / (4.0 * _vacancy) + 1.0 / _hole)

    _figure, _axes = plt.subplots(1, 2, figsize=(12.8, 4.5), constrained_layout=True)
    _axes[0].semilogy(_eta, 4.0 * _vacancy, color="#4C7C86", label=r"$4c_V$")
    _axes[0].semilogy(_eta, _hole, color="#B8734A", label=r"$c_h$")
    _axes[0].set(xlabel="Illustrative oxygen chemical-potential coordinate", ylabel="Relative concentration", title="Two defects must respond together")
    _axes[0].grid(alpha=0.24, which="both")
    _axes[0].legend()
    _axes[1].semilogy(_eta, _storage, color="#7C6A91", label=r"$[1/(4c_V)+1/c_h]^{-1}$")
    _axes[1].semilogy(_eta, np.minimum(4.0 * _vacancy, _hole), color="#858B90", ls="--", label="less abundant partner")
    _axes[1].set(xlabel="Illustrative oxygen chemical-potential coordinate", ylabel="Relative chemical storage", title="The minority partner limits storage")
    _axes[1].grid(alpha=0.24, which="both")
    _axes[1].legend()
    plt.close(_figure)

    _oxygen_reader = mo.vstack([
        mo.md(r"""
        ### Example A · A minimal oxygen mixed conductor

        For ideal oxygen vacancies and holes, a useful differential-storage form is

        \[
        c_{\rm chem}^V=\frac{F^2}{RT}
        \left(\frac{1}{4c_V}+\frac{1}{c_h}\right)^{-1}.
        \]

        The harmonic combination makes the less abundant charge-compensating
        defect the bottleneck. The curves below are an illustrative defect-chemistry
        coordinate, not a fitted Brouwer diagram.
        """),
        _figure,
    ])
    mo.accordion({"Explore further — defect chemistry controls storage": _oxygen_reader})
    return


@app.cell
def _(mo):
    mo.accordion({
        "Explore further — lithium insertion and the factor of two": mo.md(r"""
        ### Example B · An ideal lithium/electron pair

        For $\mathrm{Li}\rightleftharpoons\mathrm{Li^+}+e^-$ with
        $c_{\rm Li^+}=c_e=c$,

        \[
        \mu_{\rm Li}=\mu_{\rm Li^+}+\mu_e
        =\mu_{\rm Li}^0+2RT\ln(c/c_0),
        \qquad
        C_{\rm chem}=\frac{F^2V_{\rm act}c}{2RT}.
        \]

        The factor of two comes from this particular ideal, locally neutral
        one-ion/one-electron pair. It is not a universal prefactor for every
        insertion compound. Site exclusion, interactions, and multiple defects
        change $\partial\mu/\partial c$ and must be derived from their own free energy.
        """)
    })
    return


@app.cell
def _(mo):
    regular_interaction_07 = mo.ui.slider(
        0.0, 4.0, value=2.8, step=0.1,
        label=r"Regular-solution interaction, $\chi$", show_value=True,
    )
    return (regular_interaction_07,)


@app.cell
def _(mo, np, optimize, plt, regular_interaction_07):
    _chi = float(regular_interaction_07.value)
    _x = np.linspace(0.002, 0.998, 800)
    _g = _x * np.log(_x) + (1.0 - _x) * np.log(1.0 - _x) + _chi * _x * (1.0 - _x)
    _mu = np.log(_x / (1.0 - _x)) + _chi * (1.0 - 2.0 * _x)
    _curvature = 1.0 / _x + 1.0 / (1.0 - _x) - 2.0 * _chi
    _stable_capacity = np.where(_curvature > 0.0, np.minimum(1.0 / _curvature, 2.0), np.nan)
    _has_two_phase = _chi > 2.0
    _x_alpha = None
    _x_beta = None
    _spinodal_low = None
    _spinodal_high = None
    if _has_two_phase:
        _root = lambda value: np.log(value / (1.0 - value)) + _chi * (1.0 - 2.0 * value)
        _x_alpha = float(optimize.brentq(_root, 1.0e-8, 0.499999))
        _x_beta = 1.0 - _x_alpha
        _spinodal_low = 0.5 * (1.0 - np.sqrt(1.0 - 2.0 / _chi))
        _spinodal_high = 1.0 - _spinodal_low

    _figure, _axes = plt.subplots(1, 3, figsize=(15.2, 4.6), constrained_layout=True)
    _axes[0].plot(_x, _g, color="#7C6A91", label="homogeneous free energy")
    _axes[1].plot(_x, _mu, color="#4C7C86", label="homogeneous chemical potential")
    _axes[2].plot(_x, _stable_capacity, color="#B8734A", label="stable single-phase response (capped)")
    if _has_two_phase:
        _g_coex = _x_alpha * np.log(_x_alpha) + (1.0 - _x_alpha) * np.log(1.0 - _x_alpha) + _chi * _x_alpha * (1.0 - _x_alpha)
        _axes[0].plot([_x_alpha, _x_beta], [_g_coex, _g_coex], color="#5F8A6B", lw=2.6, label="common tangent")
        _axes[0].scatter([_x_alpha, _x_beta], [_g_coex, _g_coex], color="#5F8A6B", s=60, zorder=4)
        _axes[1].plot([_x_alpha, _x_beta], [0.0, 0.0], color="#5F8A6B", lw=2.6, label="two-phase equilibrium plateau")
        for _axis in _axes:
            _axis.axvspan(_spinodal_low, _spinodal_high, color="#B77A82", alpha=0.12)
        _broadened_peak = 1.45 * np.exp(-((_x - 0.5) / 0.12) ** 2)
        _axes[2].plot(_x, _broadened_peak, color="#5F8A6B", ls="--", label="broadened two-phase peak (schematic)")
    _axes[0].set(xlabel=r"$x$", ylabel=r"$G/(c_{\max}VRT)$", title="Interactions can create two minima")
    _axes[1].set(xlabel=r"$x$", ylabel=r"$(\mu-\mu^0)/(RT)$", title="Common tangent gives a voltage plateau")
    _axes[2].set(xlabel=r"$x$", ylabel="normalized differential response", title="A plateau appears as a capacity peak")
    for _axis in _axes:
        _axis.grid(alpha=0.24)
        _axis.legend(loc="best", fontsize=9.5)
    plt.close(_figure)

    _phase_statement = (
        f"For χ={_chi:.1f}, the common-tangent compositions are x={_x_alpha:.3f} and {_x_beta:.3f}."
        if _has_two_phase
        else f"For χ={_chi:.1f}, the homogeneous free energy remains convex and there is no two-phase interval."
    )
    _reader = mo.vstack([
        mo.md(r"""
        ### Example C · Regular-solution interactions and phase separation

        \[
        \frac{G}{c_{\max}VRT}=x\ln x+(1-x)\ln(1-x)+\chi x(1-x).
        \]

        When interactions make the homogeneous free energy non-convex, stable
        equilibrium follows its convex hull (the common tangent), not the
        negative-curvature branch. A two-phase voltage plateau therefore appears
        experimentally as a large, usually broadened differential-capacity peak.
        """),
        regular_interaction_07,
        _figure,
        mo.callout(
            _phase_statement
            + " Negative homogeneous curvature is a spinodal warning, not a stable negative chemical capacitance.",
            kind="warn" if _has_two_phase else "info",
        ),
    ])
    mo.accordion({"Explore further — phase separation without negative stable capacitance": _reader})
    regular_phase_data_07 = {
        "has_two_phase": _has_two_phase,
        "curvature": _curvature,
        "stable_capacity": _stable_capacity,
    }
    return (regular_phase_data_07,)


@app.cell
def _(
    FARADAY_C_PER_MOL,
    GAS_CONSTANT_J_PER_MOL_K,
    c_max_mol_per_m3_07,
    dilute_storage_state_07,
    interface_capacitance_areal_07,
    lattice_gas_storage_state_07,
    mixed_transport_07,
    np,
    regular_phase_data_07,
    selected_capacitance_f_07,
    selected_diffusion_time_s_07,
    selected_length_m_07,
    selected_measured_areal_07,
    selected_rc_s_07,
    selected_sigma_amb_07,
    selected_storage_state_07,
    selected_temperature_k_07,
    slab_rc_07,
    stored_fraction_07,
):
    _x = float(stored_fraction_07.value)
    _temperature = selected_temperature_k_07
    _volume = 2.3e-9
    _step = 1.0e-6

    def _numerical_dq_du(state_function):
        _minus = state_function(_x - _step, _temperature, c_max_mol_per_m3_07)
        _plus = state_function(_x + _step, _temperature, c_max_mol_per_m3_07)
        _q_minus = FARADAY_C_PER_MOL * _volume * _minus["concentration_mol_per_m3"]
        _q_plus = FARADAY_C_PER_MOL * _volume * _plus["concentration_mol_per_m3"]
        _u_minus = GAS_CONSTANT_J_PER_MOL_K * _temperature * _minus["mu_over_rt"] / FARADAY_C_PER_MOL
        _u_plus = GAS_CONSTANT_J_PER_MOL_K * _temperature * _plus["mu_over_rt"] / FARADAY_C_PER_MOL
        return (_q_plus - _q_minus) / (_u_plus - _u_minus)

    _dilute_state = dilute_storage_state_07(_x, _temperature, c_max_mol_per_m3_07)
    _lattice_state = lattice_gas_storage_state_07(_x, _temperature, c_max_mol_per_m3_07)
    _dilute_derivative_error = abs(
        _numerical_dq_du(dilute_storage_state_07)
        / (_dilute_state["c_chem_v_f_per_m3"] * _volume) - 1.0
    )
    _lattice_derivative_error = abs(
        _numerical_dq_du(lattice_gas_storage_state_07)
        / (_lattice_state["c_chem_v_f_per_m3"] * _volume) - 1.0
    )
    _dilute_fth_error = abs(_dilute_state["f_th"] - 1.0)
    _lattice_fth_error = abs(_lattice_state["f_th"] - 1.0 / (1.0 - _x))

    _interface_total_f = float(interface_capacitance_areal_07.value) * 1.0e-6
    _measured_total_f = selected_measured_areal_07 * 1.0e-6
    _interface_residual = abs((_measured_total_f - selected_capacitance_f_07) - _interface_total_f)
    _rc_identity_error = abs(selected_rc_s_07 / selected_diffusion_time_s_07 - 1.0)
    _diffusivity_identity_error = abs(
        selected_sigma_amb_07 / selected_storage_state_07["c_chem_v_f_per_m3"]
        * selected_diffusion_time_s_07 / selected_length_m_07**2 - 1.0
    )
    _sigma_i_check = 0.013
    _sigma_e_check = 1.0e8 * _sigma_i_check
    _sigma_amb_check, _ = mixed_transport_07(
        _sigma_i_check, _sigma_e_check, selected_storage_state_07["c_chem_v_f_per_m3"]
    )
    _high_electron_limit_error = abs(_sigma_amb_check / _sigma_i_check - 1.0)
    _positive_pass = bool(
        selected_storage_state_07["c_chem_v_f_per_m3"] > 0.0
        and selected_capacitance_f_07 > 0.0
        and selected_sigma_amb_07 > 0.0
    )
    _phase_mask_pass = bool(
        np.all(np.isnan(regular_phase_data_07["stable_capacity"][regular_phase_data_07["curvature"] <= 0.0]))
        and np.all(regular_phase_data_07["stable_capacity"][regular_phase_data_07["curvature"] > 0.0] >= 0.0)
    )
    module07_validation = {
        "dilute_derivative_pass": _dilute_derivative_error < 2.0e-9,
        "lattice_derivative_pass": _lattice_derivative_error < 2.0e-9,
        "dilute_fth_pass": _dilute_fth_error < 1.0e-14,
        "lattice_fth_pass": _lattice_fth_error < 1.0e-14,
        "positive_pass": _positive_pass,
        "interface_pass": _interface_residual < 1.0e-12,
        "diffusivity_pass": _diffusivity_identity_error < 1.0e-12,
        "rc_pass": _rc_identity_error < 1.0e-12,
        "high_electron_limit_pass": _high_electron_limit_error < 2.0e-8,
        "phase_mask_pass": _phase_mask_pass,
    }
    return (module07_validation,)


@app.cell
def _(mo, module07_validation):
    def _mark(passed):
        return "PASS" if passed else "CHECK"

    _checks = mo.md(
        rf"""
        ## Physical consistency checks

        | status | physical statement | why it matters |
        |---:|---|---|
        | {_mark(module07_validation['dilute_derivative_pass'])} | numerical $dQ/dU$ matches the ideal-dilute formula | capacitance really is the differential response |
        | {_mark(module07_validation['lattice_derivative_pass'])} | numerical $dQ/dU$ matches the lattice-gas formula | site exclusion is included consistently |
        | {_mark(module07_validation['dilute_fth_pass'])} | $f_{{\rm th}}=1$ for the ideal dilute model | the limiting thermodynamics is recovered |
        | {_mark(module07_validation['lattice_fth_pass'])} | $f_{{\rm th}}=1/(1-x)$ for the lattice gas | the chemical-potential slope and capacitance agree |
        | {_mark(module07_validation['positive_pass'])} | stable displayed capacitances, conductivities, and storage states are positive | every core state is physical |
        | {_mark(module07_validation['interface_pass'])} | the interface contribution is independent of thickness | the intercept and bulk slope remain distinct |
        | {_mark(module07_validation['diffusivity_pass'])} | $D^\delta=\sigma_{{\rm amb}}/c_{{\rm chem}}^V$ | transport and thermodynamics use the same state |
        | {_mark(module07_validation['rc_pass'])} | $R_{{\rm amb}}C_{{\rm chem}}=L^2/D^\delta$ | the slab and transmission-line time scales agree |
        | {_mark(module07_validation['high_electron_limit_pass'])} | $\sigma_e\gg\sigma_i$ gives $\sigma_{{\rm amb}}\to\sigma_i$ | the slower carrier remains the bottleneck |
        | {_mark(module07_validation['phase_mask_pass'])} | negative homogeneous curvature is excluded from stable-capacitance curves | the spinodal is not misread as negative stable storage |

        These checks protect the conceptual chain; students do not need them to
        operate the notebook.
        """
    )
    mo.accordion({"Physical consistency checks": _checks})
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## What to carry forward

    1. **Chemical capacitance is differential neutral storage.** It follows from
       $G(c)$ through $\mu(c)$ and $\partial c/\partial\mu$.
    2. **A capacitance value is not self-identifying.** Thickness scaling and a
       physical equivalent model help separate bulk chemical storage from
       interface and dielectric contributions.
    3. **Transport and storage form one time scale.** In the present local,
       single-phase model, $D^\delta=\sigma_{\rm amb}/c_{\rm chem}^V$ and
       $R_{\rm amb}C_{\rm chem}=L^2/D^\delta$.
    4. **Phase separation needs equilibrium construction.** A common tangent
       replaces the unstable homogeneous branch; negative curvature is a
       warning, not a stable negative capacitance.

    ### Sources and further reading

    - J. Maier, “Chemical resistance and chemical capacitance,”
      *Zeitschrift für Naturforschung B* **75** (2020),
      [doi:10.1515/znb-2019-0163](https://doi.org/10.1515/znb-2019-0163).
    - A. Schmid, G. M. Rupp, and J. Fleig, “How to determine the chemical
      capacitance of mixed conducting materials,” *Physical Chemistry Chemical
      Physics* **20** (2018),
      [doi:10.1039/C7CP07845E](https://doi.org/10.1039/C7CP07845E).
    - A. E. Bumberger, C. Steinbach, J. Ring, and J. Fleig, chemical-capacitance
      analysis of mixed conductors, *Chemistry of Materials* **34** (2022),
      [doi:10.1021/acs.chemmater.2c02614](https://doi.org/10.1021/acs.chemmater.2c02614).
    - L.-Q. Chen and J. C. Mauro, “Thermodynamics and capacitance,”
      *MRS Bulletin* **49** (2024),
      [doi:10.1557/s43577-024-00727-4](https://doi.org/10.1557/s43577-024-00727-4).

    **Model boundary.** The core assumes local equilibrium, one neutral storage
    coordinate, ideal dilute or ideal lattice-gas thermodynamics, constant
    temperature, and small-signal derivatives about a stable state. It omits
    elastic coupling, distributed phase-boundary motion, nonlinear reaction
    kinetics, and nonlocal electrostatics unless an Explore section states otherwise.

    **Previous:** [Module 06 — PITT and GITT](https://qiyanglu.github.io/Solid-State-Ionics-Interactive/06-pitt-gitt/)

    **Continue:** [Module 08 — Impedance, Warburg Diffusion, and Transmission Lines](https://qiyanglu.github.io/Solid-State-Ionics-Interactive/08-impedance-tlm/)
    """)
    return


if __name__ == "__main__":
    app.run()
