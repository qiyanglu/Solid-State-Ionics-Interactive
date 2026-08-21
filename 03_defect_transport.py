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
    mo.md(r"""
    # Defect transport: from atomic hopping to chemical diffusion

    **How can random atomic hops produce a predictable diffusivity?**

    A single defect wanders unpredictably, but many identical defects obey a
    simple statistical law. In one dimension we use the convention

    \[
    \Gamma=\nu e^{-\Delta H_{\rm mig}/(k_BT)},\qquad
    D=\frac{a^2\Gamma}{2}.
    \]

    We will build that result visually: one activated hop, many random hops,
    net exchange across a concentration gradient, and finally the coupled
    motion of \(\mathrm{Li^+}\) and electrons during chemical diffusion.
    Every spatial picture is one-dimensional.
    """)
    return


@app.cell
def _(np):
    KB_EV_PER_K = 8.617333262e-5
    KB_J_PER_K = 1.380649e-23
    E_CHARGE_C = 1.602176634e-19
    AVOGADRO_PER_MOL = 6.02214076e23
    FARADAY_C_PER_MOL = E_CHARGE_C * AVOGADRO_PER_MOL
    GAS_CONSTANT_J_PER_MOL_K = KB_J_PER_K * AVOGADRO_PER_MOL

    def require_positive(name, value):
        """Return a positive finite scalar."""
        number = float(value)
        if not np.isfinite(number) or number <= 0.0:
            raise ValueError(f"{name} must be positive and finite")
        return number

    def hop_frequency(temperature_k, migration_enthalpy_ev, attempt_frequency_hz):
        """Return Gamma = nu exp[-Delta H_mig/(k_B T)] in s^-1."""
        temperature = require_positive("temperature_k", temperature_k)
        barrier = require_positive("migration_enthalpy_ev", migration_enthalpy_ev)
        attempt = require_positive("attempt_frequency_hz", attempt_frequency_hz)
        return attempt * np.exp(-barrier / (KB_EV_PER_K * temperature))

    def hopping_diffusivity_1d(jump_distance_m, hop_frequency_hz):
        """Return D = a^2 Gamma / 2 in m^2/s for a 1D random walk."""
        distance = require_positive("jump_distance_m", jump_distance_m)
        frequency = require_positive("hop_frequency_hz", hop_frequency_hz)
        return 0.5 * distance**2 * frequency

    def biased_directional_rates(
        temperature_k,
        migration_enthalpy_ev,
        attempt_frequency_hz,
        jump_distance_m,
        charge_number,
        electric_field_v_per_m,
    ):
        """Return Gamma_+, Gamma_-, and zero-field Gamma in s^-1.

        Gamma is the total zero-field hop frequency used throughout this module.
        At zero field, half of the hops go toward +x and half toward -x.
        """
        temperature = require_positive("temperature_k", temperature_k)
        distance = require_positive("jump_distance_m", jump_distance_m)
        charge = float(charge_number)
        field = float(electric_field_v_per_m)
        if not np.isfinite(charge) or charge == 0.0:
            raise ValueError("charge_number must be finite and nonzero")
        if not np.isfinite(field):
            raise ValueError("electric_field_v_per_m must be finite")
        frequency = hop_frequency(
            temperature,
            migration_enthalpy_ev,
            attempt_frequency_hz,
        )
        half_bias = charge * field * distance / (
            2.0 * KB_EV_PER_K * temperature
        )
        return (
            0.5 * frequency * np.exp(half_bias),
            0.5 * frequency * np.exp(-half_bias),
            frequency,
        )

    def exact_hopping_drift_velocity(
        jump_distance_m,
        forward_rate_hz,
        backward_rate_hz,
    ):
        """Return a(Gamma_+ - Gamma_-) in m/s."""
        distance = require_positive("jump_distance_m", jump_distance_m)
        return distance * (float(forward_rate_hz) - float(backward_rate_hz))

    def nernst_einstein_drift_velocity(
        diffusivity_m2_per_s,
        temperature_k,
        charge_number,
        electric_field_v_per_m,
    ):
        """Return low-field particle drift z e D E/(k_B T) in m/s."""
        diffusivity = require_positive("diffusivity_m2_per_s", diffusivity_m2_per_s)
        temperature = require_positive("temperature_k", temperature_k)
        return (
            float(charge_number)
            * E_CHARGE_C
            * diffusivity
            * float(electric_field_v_per_m)
            / (KB_J_PER_K * temperature)
        )

    def simulate_zero_field_walks(walker_count=6000, step_count=320, seed=2026):
        """Return one fixed ensemble of unbiased one-dimensional lattice walks.

        The output is expressed only in hop number and lattice-site units.
        Physical distance and time are applied later through a and Gamma.
        """
        walkers = int(walker_count)
        steps = int(step_count)
        if walkers < 2 or steps < 2:
            raise ValueError("walker_count and step_count must be at least 2")
        rng = np.random.default_rng(int(seed))
        increments = rng.choice((-1.0, 1.0), size=(walkers, steps))
        positions_in_lattice_steps = np.concatenate(
            [np.zeros((walkers, 1)), np.cumsum(increments, axis=1)],
            axis=1,
        )
        hop_numbers = np.arange(steps + 1, dtype=float)
        mean_square_steps = np.mean(positions_in_lattice_steps**2, axis=0)
        return hop_numbers, positions_in_lattice_steps, mean_square_steps

    def evolve_periodic_master_equation(initial_occupancy, hop_frequency_hz, time_s):
        """Evolve dc_j/dt=(Gamma/2)(c_{j-1}-2c_j+c_{j+1}) in 1D."""
        profile = np.asarray(initial_occupancy, dtype=float)
        if profile.ndim != 1 or profile.size < 4:
            raise ValueError("initial_occupancy must be a 1D array with at least 4 sites")
        if np.any(~np.isfinite(profile)) or np.any(profile < 0.0):
            raise ValueError("initial_occupancy must be finite and nonnegative")
        frequency = require_positive("hop_frequency_hz", hop_frequency_hz)
        elapsed = float(time_s)
        if not np.isfinite(elapsed) or elapsed < 0.0:
            raise ValueError("time_s must be finite and nonnegative")
        modes = np.arange(profile.size)
        decay = np.exp(
            -2.0 * frequency * elapsed * np.sin(np.pi * modes / profile.size) ** 2
        )
        return np.real(np.fft.ifft(np.fft.fft(profile) * decay))

    def discrete_bond_fluxes(occupancy, hop_frequency_hz):
        """Return net 1D particle crossings per second on bonds j -> j+1."""
        profile = np.asarray(occupancy, dtype=float)
        frequency = require_positive("hop_frequency_hz", hop_frequency_hz)
        return 0.5 * frequency * (profile - np.roll(profile, -1))

    def fick_bond_fluxes(occupancy, jump_distance_m, diffusivity_m2_per_s):
        """Return -D dc/dx for c=occupancy/a on each 1D bond, in s^-1."""
        profile = np.asarray(occupancy, dtype=float)
        distance = require_positive("jump_distance_m", jump_distance_m)
        diffusivity = require_positive("diffusivity_m2_per_s", diffusivity_m2_per_s)
        concentration_per_m = profile / distance
        gradient_per_m2 = (
            np.roll(concentration_per_m, -1) - concentration_per_m
        ) / distance
        return -diffusivity * gradient_per_m2

    def molar_nernst_planck_flux(
        diffusivity_m2_per_s,
        concentration_mol_per_m3,
        concentration_gradient_mol_per_m4,
        potential_gradient_v_per_m,
        charge_number,
        temperature_k,
    ):
        """Return diffusion, electrical, and total 1D molar flux."""
        diffusivity = require_positive("diffusivity_m2_per_s", diffusivity_m2_per_s)
        concentration = require_positive(
            "concentration_mol_per_m3",
            concentration_mol_per_m3,
        )
        temperature = require_positive("temperature_k", temperature_k)
        diffusion_flux = -diffusivity * float(concentration_gradient_mol_per_m4)
        electrical_flux = (
            -float(charge_number)
            * FARADAY_C_PER_MOL
            * diffusivity
            * concentration
            * float(potential_gradient_v_per_m)
            / (GAS_CONSTANT_J_PER_MOL_K * temperature)
        )
        return diffusion_flux, electrical_flux, diffusion_flux + electrical_flux

    def ambipolar_internal_potential_gradient(
        ionic_diffusivity_m2_per_s,
        electronic_diffusivity_m2_per_s,
        concentration_mol_per_m3,
        concentration_gradient_mol_per_m4,
        temperature_k,
    ):
        """Return dphi/dx enforcing J_Li+ = J_e- in one dimension."""
        ionic = require_positive("ionic_diffusivity_m2_per_s", ionic_diffusivity_m2_per_s)
        electronic = require_positive(
            "electronic_diffusivity_m2_per_s",
            electronic_diffusivity_m2_per_s,
        )
        concentration = require_positive(
            "concentration_mol_per_m3",
            concentration_mol_per_m3,
        )
        temperature = require_positive("temperature_k", temperature_k)
        return (
            GAS_CONSTANT_J_PER_MOL_K
            * temperature
            / (FARADAY_C_PER_MOL * concentration)
            * (electronic - ionic)
            / (electronic + ionic)
            * float(concentration_gradient_mol_per_m4)
        )

    def ambipolar_fluxes(
        ionic_diffusivity_m2_per_s,
        electronic_diffusivity_m2_per_s,
        concentration_mol_per_m3,
        concentration_gradient_mol_per_m4,
        potential_gradient_v_per_m,
        temperature_k,
    ):
        """Return J_Li+ and J_e- in mol/(m^2 s)."""
        ionic = require_positive("ionic_diffusivity_m2_per_s", ionic_diffusivity_m2_per_s)
        electronic = require_positive(
            "electronic_diffusivity_m2_per_s",
            electronic_diffusivity_m2_per_s,
        )
        concentration = require_positive(
            "concentration_mol_per_m3",
            concentration_mol_per_m3,
        )
        temperature = require_positive("temperature_k", temperature_k)
        electrical_scale = (
            FARADAY_C_PER_MOL
            * concentration
            * float(potential_gradient_v_per_m)
            / (GAS_CONSTANT_J_PER_MOL_K * temperature)
        )
        gradient = float(concentration_gradient_mol_per_m4)
        ionic_flux = -ionic * (gradient + electrical_scale)
        electronic_flux = -electronic * (gradient - electrical_scale)
        return ionic_flux, electronic_flux

    def lithium_chemical_diffusivity(
        ionic_diffusivity_m2_per_s,
        electronic_diffusivity_m2_per_s,
    ):
        """Return D_Li^delta = 2 D_Li+ D_e-/(D_Li+ + D_e-) in m^2/s."""
        ionic = require_positive("ionic_diffusivity_m2_per_s", ionic_diffusivity_m2_per_s)
        electronic = require_positive(
            "electronic_diffusivity_m2_per_s",
            electronic_diffusivity_m2_per_s,
        )
        return 2.0 * ionic * electronic / (ionic + electronic)

    def conductivity_from_diffusivity(
        diffusivity_m2_per_s,
        concentration_mol_per_m3,
        temperature_k,
    ):
        """Return sigma=F^2 c D/(RT) in S/m for a monovalent carrier."""
        diffusivity = require_positive("diffusivity_m2_per_s", diffusivity_m2_per_s)
        concentration = require_positive(
            "concentration_mol_per_m3",
            concentration_mol_per_m3,
        )
        temperature = require_positive("temperature_k", temperature_k)
        return (
            FARADAY_C_PER_MOL**2
            * concentration
            * diffusivity
            / (GAS_CONSTANT_J_PER_MOL_K * temperature)
        )

    def lithium_chemical_diffusivity_from_conductivities(
        ionic_conductivity_s_per_m,
        electronic_conductivity_s_per_m,
        ionic_concentration_mol_per_m3,
        electronic_concentration_mol_per_m3,
        temperature_k,
    ):
        """Return D_Li^delta for the ideal dilute Li+/electron pair."""
        ionic_sigma = require_positive(
            "ionic_conductivity_s_per_m",
            ionic_conductivity_s_per_m,
        )
        electronic_sigma = require_positive(
            "electronic_conductivity_s_per_m",
            electronic_conductivity_s_per_m,
        )
        ionic_c = require_positive(
            "ionic_concentration_mol_per_m3",
            ionic_concentration_mol_per_m3,
        )
        electronic_c = require_positive(
            "electronic_concentration_mol_per_m3",
            electronic_concentration_mol_per_m3,
        )
        temperature = require_positive("temperature_k", temperature_k)
        harmonic_conductivity = (
            ionic_sigma * electronic_sigma / (ionic_sigma + electronic_sigma)
        )
        return (
            GAS_CONSTANT_J_PER_MOL_K
            * temperature
            / FARADAY_C_PER_MOL**2
            * harmonic_conductivity
            * (1.0 / ionic_c + 1.0 / electronic_c)
        )

    def characteristic_relaxation_time(length_m, diffusivity_m2_per_s):
        """Return the one-dimensional scaling time L^2/D in seconds."""
        length = require_positive("length_m", length_m)
        diffusivity = require_positive("diffusivity_m2_per_s", diffusivity_m2_per_s)
        return length**2 / diffusivity

    def slab_remaining_profile(normalized_position, time_over_l2_by_d, term_count=80):
        """Return normalized excess in a 1D slab held at zero at both faces."""
        coordinate = np.asarray(normalized_position, dtype=float)
        reduced_time = require_positive("time_over_l2_by_d", time_over_l2_by_d)
        odd = 2 * np.arange(int(term_count)) + 1
        terms = (
            4.0
            / np.pi
            * np.sin(np.pi * coordinate[:, None] * odd[None, :])
            / odd[None, :]
            * np.exp(-(np.pi * odd[None, :]) ** 2 * reduced_time)
        )
        return np.sum(terms, axis=1)

    def scaled_time_axis(times_s):
        """Return display-scaled times and a readable axis unit."""
        maximum = float(np.max(times_s))
        if maximum < 1.0e-6:
            return np.asarray(times_s) * 1.0e9, "ns"
        if maximum < 1.0e-3:
            return np.asarray(times_s) * 1.0e6, "µs"
        if maximum < 1.0:
            return np.asarray(times_s) * 1.0e3, "ms"
        return np.asarray(times_s), "s"

    return (
        AVOGADRO_PER_MOL,
        E_CHARGE_C,
        FARADAY_C_PER_MOL,
        GAS_CONSTANT_J_PER_MOL_K,
        KB_EV_PER_K,
        KB_J_PER_K,
        ambipolar_fluxes,
        ambipolar_internal_potential_gradient,
        biased_directional_rates,
        characteristic_relaxation_time,
        conductivity_from_diffusivity,
        discrete_bond_fluxes,
        evolve_periodic_master_equation,
        exact_hopping_drift_velocity,
        fick_bond_fluxes,
        hop_frequency,
        hopping_diffusivity_1d,
        lithium_chemical_diffusivity,
        lithium_chemical_diffusivity_from_conductivities,
        molar_nernst_planck_flux,
        nernst_einstein_drift_velocity,
        require_positive,
        scaled_time_axis,
        simulate_zero_field_walks,
        slab_remaining_profile,
    )


@app.cell
def _(mo):
    temperature = mo.ui.slider(
        start=400,
        stop=1600,
        step=25,
        value=900,
        label="Temperature, T (K)",
        show_value=True,
    )
    migration_barrier = mo.ui.slider(
        start=0.20,
        stop=1.50,
        step=0.05,
        value=0.70,
        label="Migration enthalpy (eV)",
        show_value=True,
    )
    log_attempt_frequency = mo.ui.slider(
        start=10.0,
        stop=14.0,
        step=0.25,
        value=13.0,
        label="Attempt-frequency exponent (base 10, per second)",
        show_value=True,
    )
    jump_distance = mo.ui.slider(
        start=0.20,
        stop=0.80,
        step=0.02,
        value=0.40,
        label="Jump distance, a (nm)",
        show_value=True,
    )
    hopping_controls = mo.hstack(
        [temperature, migration_barrier],
        justify="start",
        align="center",
        wrap=True,
        gap=1.5,
    )
    advanced_hopping_controls = mo.accordion(
        {
            "Explore further — jump geometry and attempt rate": mo.vstack(
                [log_attempt_frequency, jump_distance], gap=0.7
            )
        }
    )
    return (
        advanced_hopping_controls,
        hopping_controls,
        jump_distance,
        log_attempt_frequency,
        migration_barrier,
        temperature,
    )

@app.cell
def _(advanced_hopping_controls, hopping_controls, mo):
    mo.vstack(
        [
            mo.md(r"""
            ## 1. One activated hop

            A defect remains near a stable lattice site until thermal motion
            carries it over the migration barrier. Raise the temperature or
            lower the barrier and ask which part of the picture changes: the
            geometry stays fixed, while the hop frequency rises sharply.
            """),
            hopping_controls,
            advanced_hopping_controls,
        ]
    )
    return

@app.cell
def _(
    hopping_diffusivity_1d,
    hop_frequency,
    jump_distance,
    log_attempt_frequency,
    migration_barrier,
    np,
    simulate_zero_field_walks,
    temperature,
):
    temperature_k = float(temperature.value)
    migration_enthalpy_ev = float(migration_barrier.value)
    attempt_frequency_hz = 10.0 ** float(log_attempt_frequency.value)
    jump_distance_m = float(jump_distance.value) * 1.0e-9
    hop_frequency_hz = hop_frequency(
        temperature_k,
        migration_enthalpy_ev,
        attempt_frequency_hz,
    )
    defect_diffusivity_m2_per_s = hopping_diffusivity_1d(
        jump_distance_m,
        hop_frequency_hz,
    )
    walk_hop_numbers, walk_positions_in_steps, walk_msd_in_steps2 = (
        simulate_zero_field_walks(
            walker_count=6000,
            step_count=320,
            seed=2026,
        )
    )
    walk_times_s = walk_hop_numbers / hop_frequency_hz
    walk_positions_m = jump_distance_m * walk_positions_in_steps
    walk_msd_m2 = jump_distance_m**2 * walk_msd_in_steps2
    fit_start = walk_times_s.size // 5
    _fit_times = walk_times_s[fit_start:]
    _fit_msd = walk_msd_m2[fit_start:]
    msd_slope = float(np.dot(_fit_times, _fit_msd) / np.dot(_fit_times, _fit_times))
    extracted_diffusivity_m2_per_s = msd_slope / 2.0
    fitted_msd_m2 = msd_slope * walk_times_s
    residual_sum_squares = float(
        np.sum((_fit_msd - fitted_msd_m2[fit_start:]) ** 2)
    )
    total_sum_squares = float(
        np.sum((_fit_msd - np.mean(_fit_msd)) ** 2)
    )
    msd_fit_r_squared = 1.0 - residual_sum_squares / total_sum_squares
    return (
        attempt_frequency_hz,
        defect_diffusivity_m2_per_s,
        extracted_diffusivity_m2_per_s,
        fitted_msd_m2,
        jump_distance_m,
        migration_enthalpy_ev,
        msd_fit_r_squared,
        temperature_k,
        walk_hop_numbers,
        walk_msd_m2,
        walk_positions_m,
        walk_positions_in_steps,
        walk_times_s,
        hop_frequency_hz,
    )


@app.cell
def _(mo):
    mo.md(r"""
    **Predict before moving the controls.** Raising \(T\), lowering
    \(\Delta H_{\rm mig}\), or raising \(\nu\) should make the same number of
    hops occur in less physical time. Increasing \(a\) should enlarge every
    spatial step and increase \(D\) as \(a^2\).
    """)
    return


@app.cell
def _(
    defect_diffusivity_m2_per_s,
    extracted_diffusivity_m2_per_s,
    hop_frequency_hz,
    jump_distance_m,
    migration_enthalpy_ev,
    msd_fit_r_squared,
    np,
    plt,
    scaled_time_axis,
    walk_hop_numbers,
    walk_msd_m2,
    walk_positions_in_steps,
    walk_times_s,
):
    reaction_coordinate = np.linspace(0.0, 2.0, 500)
    barrier_energy_ev = 0.5 * migration_enthalpy_ev * (
        1.0 - np.cos(2.0 * np.pi * reaction_coordinate)
    )

    hop_figure, barrier_axis = plt.subplots(figsize=(11.5, 3.9), dpi=120)
    barrier_axis.plot(
        reaction_coordinate, barrier_energy_ev, color="#4C7C86", lw=1.7
    )
    barrier_axis.scatter(
        [0.0, 1.0, 2.0],
        [0.0, 0.0, 0.0],
        s=70,
        color="#C49345",
        edgecolor="#40464D",
        zorder=4,
        label="Equivalent sites",
    )
    barrier_axis.annotate(
        "",
        xy=(0.5, migration_enthalpy_ev),
        xytext=(0.5, 0.0),
        arrowprops={"arrowstyle": "<->", "color": "#B65C4A", "lw": 1.4},
    )
    barrier_axis.text(
        0.55,
        0.52 * migration_enthalpy_ev,
        r"$\Delta H_{\rm mig}$",
        color="#B65C4A",
    )
    barrier_axis.set(
        xlabel=r"Reaction coordinate, $x/a$",
        ylabel="Energy (eV)",
        title="A hop must cross the migration barrier",
        ylim=(-0.04 * migration_enthalpy_ev, 1.15 * migration_enthalpy_ev),
    )
    barrier_axis.legend(frameon=False, loc="upper right")
    barrier_axis.grid(alpha=0.2)
    hop_figure.tight_layout()
    plt.close(hop_figure)

    trajectory_figure, trajectory_axis = plt.subplots(
        figsize=(8.6, 4.5), dpi=120
    )
    for trajectory_index, (_color, _style) in enumerate(
        zip(("#4C7C86", "#B8734A", "#7C6A91"), ("-", "--", "-."))
    ):
        trajectory_axis.step(
            walk_hop_numbers,
            walk_positions_in_steps[trajectory_index],
            where="post",
            lw=1.5,
            color=_color,
            ls=_style,
            label=f"Walker {trajectory_index + 1}",
        )
    trajectory_axis.axhline(0.0, color="#73808C", lw=0.9, ls=":")
    trajectory_axis.set(
        xlabel="Hop number",
        ylabel=r"Position, $x/a$ (lattice steps)",
        title="A random walk is geometry before it is a clock",
    )
    trajectory_axis.grid(alpha=0.22)
    trajectory_axis.legend(frameon=False, loc="best")
    trajectory_figure.tight_layout()
    plt.close(trajectory_figure)

    selected_display_time, time_unit = scaled_time_axis(walk_times_s)
    msd_figure, msd_axis = plt.subplots(figsize=(8.6, 4.5), dpi=120)
    msd_axis.plot(
        selected_display_time,
        walk_msd_m2 * 1.0e18,
        color="#4C7C86",
        lw=1.7,
        marker="o",
        markevery=32,
        ms=3.5,
        label="Fixed random-walk ensemble",
    )
    msd_axis.plot(
        selected_display_time,
        2.0 * defect_diffusivity_m2_per_s * walk_times_s * 1.0e18,
        color="#B8734A",
        lw=1.5,
        ls="--",
        label=r"Analytical $2Dt$",
    )
    msd_axis.set(
        xlabel=f"Physical time ({time_unit})",
        ylabel=r"Mean-square displacement, $\langle x^2\rangle$ (nm$^2$)",
        title=r"The ensemble slope is $2D$",
    )
    msd_axis.grid(alpha=0.22)
    msd_axis.legend(frameon=False)
    msd_figure.tight_layout()
    plt.close(msd_figure)

    microscopic_summary = (
        rf"$D_{{\rm analytical}}={defect_diffusivity_m2_per_s * 1.0e4:.3e}$ "
        rf"cm² s⁻¹; $D_{{\rm fit}}={extracted_diffusivity_m2_per_s * 1.0e4:.3e}$ "
        rf"cm² s⁻¹; $D_{{\rm fit}}/D_{{\rm analytical}}="
        rf"{extracted_diffusivity_m2_per_s / defect_diffusivity_m2_per_s:.4f}$ "
        rf"($R^2={msd_fit_r_squared:.5f}$)."
    )
    return hop_figure, msd_figure, microscopic_summary, trajectory_figure

@app.cell
def _(hop_figure, msd_figure, microscopic_summary, mo, trajectory_figure):
    mo.vstack(
        [
            hop_figure,
            mo.md(r"""
            The barrier controls **how often** a hop occurs; it does not give a
            preferred direction between equivalent sites.

            ### Many hops give diffusivity

            The paths below use one fixed random realization. Moving a physical
            control therefore changes the clock or length scale, not the random
            sequence itself.
            """),
            trajectory_figure,
            mo.md(r"""
            Temperature and migration barrier do not change these paths in hop
            space. They change only how quickly the same hop sequence is
            traversed. The next figure applies the selected physical clock,
            $t=N_{\rm hop}/\Gamma$, and length scale, $x=a(x/a)$.
            """),
            msd_figure,
            mo.md(
                f"**At the selected state:** {microscopic_summary}  "
                r"Because $D\propto a^2\Gamma$, doubling $a$ gives $4D$, "
                r"while multiplying $\Gamma$ by ten gives $10D$."
            ),
        ]
    )
    return

@app.cell
def _(mo):
    master_time = mo.ui.slider(
        start=0.002,
        stop=0.080,
        step=0.002,
        value=0.012,
        label="Diffusion time, Dt / L²",
        show_value=True,
    )
    step_contrast = mo.ui.slider(
        start=0.10,
        stop=0.90,
        step=0.05,
        value=0.70,
        label="Initial left–right occupancy contrast",
        show_value=True,
    )
    master_controls = mo.hstack(
        [master_time, step_contrast],
        justify="start",
        align="center",
        wrap=True,
        gap=1.5,
    )
    return master_controls, master_time, step_contrast

@app.cell
def _(master_controls, mo):
    mo.vstack(
        [
            mo.md(r"""
            ## 2. Random exchange becomes Fick's law

            Put more defects on the left than on the right. Hops remain random,
            but more defects are available to cross from left to right. Across
            one bond,

            \[
            J_{j+1/2}=\frac{\Gamma}{2}(c_j-c_{j+1})
            =-\frac{D}{a^2}(c_{j+1}-c_j).
            \]

            The difference between two opposing random exchanges is therefore
            the finite-lattice form of Fick's law. The ends below are reflecting,
            so defects remain inside the sample.
            """),
            master_controls,
        ]
    )
    return

@app.cell
def _(
    defect_diffusivity_m2_per_s,
    jump_distance_m,
    master_time,
    np,
    step_contrast,
    hop_frequency_hz,
):
    master_site_count = 256
    master_position = (np.arange(master_site_count) + 0.5) / master_site_count
    contrast_value = float(step_contrast.value)
    left_occupancy = 0.5 * (1.0 + contrast_value)
    right_occupancy = 0.5 * (1.0 - contrast_value)
    initial_master_occupancy = np.where(
        master_position < 0.5,
        left_occupancy,
        right_occupancy,
    )

    fourier_time = float(master_time.value)
    mode_numbers = np.arange(1, 160, dtype=float)
    cosine_basis = np.cos(np.pi * mode_numbers[:, None] * master_position[None, :])
    mean_occupancy = float(np.mean(initial_master_occupancy))
    mode_coefficients = (2.0 / master_site_count) * (
        cosine_basis @ (initial_master_occupancy - mean_occupancy)
    )
    mode_decay = np.exp(-(np.pi * mode_numbers) ** 2 * fourier_time)
    evolved_master_occupancy = mean_occupancy + (
        mode_coefficients * mode_decay
    ) @ cosine_basis

    sample_length_m = master_site_count * jump_distance_m
    master_elapsed_s = (
        fourier_time * sample_length_m**2 / defect_diffusivity_m2_per_s
    )
    microscopic_bond_flux_per_s = 0.5 * hop_frequency_hz * (
        evolved_master_occupancy[:-1] - evolved_master_occupancy[1:]
    )
    fick_bond_flux_per_s = -(
        defect_diffusivity_m2_per_s / jump_distance_m**2
    ) * np.diff(evolved_master_occupancy)
    master_mass_relative_error = abs(
        np.sum(evolved_master_occupancy) - np.sum(initial_master_occupancy)
    ) / np.sum(initial_master_occupancy)
    master_flux_relative_error = float(
        np.max(np.abs(microscopic_bond_flux_per_s - fick_bond_flux_per_s))
        / max(np.max(np.abs(microscopic_bond_flux_per_s)), 1.0e-300)
    )
    return (
        evolved_master_occupancy,
        fick_bond_flux_per_s,
        initial_master_occupancy,
        master_elapsed_s,
        master_flux_relative_error,
        master_mass_relative_error,
        master_site_count,
        microscopic_bond_flux_per_s,
    )

@app.cell
def _(
    evolved_master_occupancy,
    initial_master_occupancy,
    master_elapsed_s,
    master_site_count,
    microscopic_bond_flux_per_s,
    mo,
    np,
    plt,
    scaled_time_axis,
):
    _master_position = (np.arange(master_site_count) + 0.5) / master_site_count
    bond_position = np.arange(1, master_site_count) / master_site_count
    master_display_time, master_time_unit = scaled_time_axis(
        np.array([master_elapsed_s])
    )
    master_figure, (profile_axis, _flux_axis) = plt.subplots(
        1, 2, figsize=(13.2, 4.8), dpi=120
    )
    profile_axis.plot(
        _master_position,
        initial_master_occupancy,
        color="#8C9196",
        lw=1.3,
        ls="--",
        label="Initial step",
    )
    profile_axis.plot(
        _master_position,
        evolved_master_occupancy,
        color="#4C7C86",
        lw=1.7,
        label="After diffusion",
    )
    profile_axis.set(
        xlabel=r"Position, $x/L$",
        ylabel="Site occupancy",
        title="Random hopping smooths the concentration step",
        ylim=(0.0, 1.0),
    )
    profile_axis.grid(alpha=0.22)
    profile_axis.legend(frameon=False)

    _flux_axis.plot(
        bond_position,
        microscopic_bond_flux_per_s,
        color="#B8734A",
        lw=1.6,
    )
    _flux_axis.axhline(0.0, color="#666D73", lw=1.0)
    _flux_axis.set(
        xlabel=r"Bond position, $x/L$",
        ylabel=r"Net exchange (site$^{-1}$ s$^{-1}$)",
        title="More defects cross from high to low concentration",
    )
    _flux_axis.grid(alpha=0.22)
    master_figure.tight_layout()
    plt.close(master_figure)

    master_summary = mo.md(
        f"""
        At this state, the profile has evolved for
        **{master_display_time[0]:.3g} {master_time_unit}**. The flux is zero at
        the reflecting ends and largest near the original interface. Both
        panels describe the same net exchange; no particle is assigned a
        deterministic force toward lower concentration.
        """
    )
    mo.vstack([master_figure, master_summary])
    return (master_figure,)

@app.cell
def _(mo):
    charge_selector = mo.ui.dropdown(
        options={"Positive defect, z = +1": 1, "Negative defect, z = -1": -1},
        value="Positive defect, z = +1",
        label="Defect charge",
    )
    field_sign = mo.ui.dropdown(
        options={"Toward +x": 1.0, "Toward -x": -1.0},
        value="Toward +x",
        label="Field direction",
    )
    log_electric_field = mo.ui.slider(
        start=2.0, stop=7.0, step=0.25, value=5.0,
        label="Field-magnitude exponent (base 10, volts per centimeter)", show_value=True,
    )
    relative_concentration_gradient = mo.ui.slider(
        start=-0.80, stop=0.80, step=0.05, value=0.30,
        label="Relative concentration gradient (per micrometer)", show_value=True,
    )
    field_controls = mo.hstack(
        [charge_selector, field_sign, log_electric_field],
        justify="start",
        align="center",
        wrap=True,
        gap=1.5,
    )
    return (
        charge_selector, field_controls, field_sign,
        log_electric_field, relative_concentration_gradient,
    )

@app.cell
def _(field_controls, mo):
    mo.vstack(
        [
            mo.md(r"""
            ## 3. An electric field biases the hops

            Define positive \(x\) to the right and
            \(\mathcal E=-d\phi/dx\). A positive field lowers one barrier and
            raises the other. The two directional rates are

            \[
            \Gamma_\pm=\frac{\Gamma}{2}
            \exp\!\left(\pm\frac{ze\mathcal Ea}{2k_BT}\right).
            \]

            Change the field direction and magnitude. The stable lattice sites
            do not move, but the energy difference between neighboring sites
            reverses with the field.
            """),
            field_controls,
        ]
    )
    return

@app.cell
def _(
    FARADAY_C_PER_MOL,
    GAS_CONSTANT_J_PER_MOL_K,
    biased_directional_rates,
    charge_selector,
    defect_diffusivity_m2_per_s,
    exact_hopping_drift_velocity,
    field_sign,
    jump_distance_m,
    log_electric_field,
    migration_enthalpy_ev,
    molar_nernst_planck_flux,
    nernst_einstein_drift_velocity,
    relative_concentration_gradient,
    temperature_k,
    attempt_frequency_hz,
):
    charge_number_value = int(charge_selector.value)
    electric_field_v_per_m = float(field_sign.value) * (
        100.0 * 10.0 ** float(log_electric_field.value)
    )
    forward_rate_hz, backward_rate_hz, _unbiased_hop_frequency_hz = (
        biased_directional_rates(
            temperature_k, migration_enthalpy_ev, attempt_frequency_hz,
            jump_distance_m, charge_number_value, electric_field_v_per_m,
        )
    )
    exact_drift_m_per_s = exact_hopping_drift_velocity(
        jump_distance_m, forward_rate_hz, backward_rate_hz
    )
    low_field_drift_m_per_s = nernst_einstein_drift_velocity(
        defect_diffusivity_m2_per_s, temperature_k, charge_number_value,
        electric_field_v_per_m,
    )
    hopping_bias_ratio = forward_rate_hz / backward_rate_hz
    field_work_per_hop_ev = (
        charge_number_value * electric_field_v_per_m * jump_distance_m
    )
    low_field_parameter = abs(field_work_per_hop_ev) / (
        8.617333262e-5 * temperature_k
    )

    demonstration_concentration_mol_per_m3 = 1000.0
    relative_gradient_per_m = float(relative_concentration_gradient.value) * 1.0e6
    demonstration_gradient_mol_per_m4 = (
        demonstration_concentration_mol_per_m3 * relative_gradient_per_m
    )
    demonstration_potential_gradient_v_per_m = (
        -GAS_CONSTANT_J_PER_MOL_K * temperature_k * relative_gradient_per_m
        / (charge_number_value * FARADAY_C_PER_MOL)
    )
    chemical_potential_gradient_j_per_mol_m = (
        GAS_CONSTANT_J_PER_MOL_K * temperature_k * relative_gradient_per_m
    )
    electrical_potential_gradient_j_per_mol_m = (
        charge_number_value * FARADAY_C_PER_MOL
        * demonstration_potential_gradient_v_per_m
    )
    electrochemical_gradient_j_per_mol_m = (
        chemical_potential_gradient_j_per_mol_m
        + electrical_potential_gradient_j_per_mol_m
    )
    (
        diffusion_flux_mol_per_m2_s,
        electrical_flux_mol_per_m2_s,
        total_np_flux_mol_per_m2_s,
    ) = molar_nernst_planck_flux(
        defect_diffusivity_m2_per_s, demonstration_concentration_mol_per_m3,
        demonstration_gradient_mol_per_m4,
        demonstration_potential_gradient_v_per_m, charge_number_value,
        temperature_k,
    )
    return (
        backward_rate_hz, charge_number_value,
        chemical_potential_gradient_j_per_mol_m,
        demonstration_concentration_mol_per_m3,
        demonstration_gradient_mol_per_m4,
        demonstration_potential_gradient_v_per_m,
        diffusion_flux_mol_per_m2_s, electric_field_v_per_m,
        electrical_flux_mol_per_m2_s,
        electrical_potential_gradient_j_per_mol_m,
        electrochemical_gradient_j_per_mol_m, exact_drift_m_per_s,
        field_work_per_hop_ev, forward_rate_hz, hopping_bias_ratio,
        low_field_drift_m_per_s, low_field_parameter, relative_gradient_per_m,
        total_np_flux_mol_per_m2_s,
    )

@app.cell
def _(
    charge_number_value,
    chemical_potential_gradient_j_per_mol_m,
    demonstration_potential_gradient_v_per_m,
    diffusion_flux_mol_per_m2_s,
    electric_field_v_per_m,
    electrical_flux_mol_per_m2_s,
    electrical_potential_gradient_j_per_mol_m,
    electrochemical_gradient_j_per_mol_m,
    exact_drift_m_per_s,
    field_work_per_hop_ev,
    forward_rate_hz,
    backward_rate_hz,
    hopping_bias_ratio,
    jump_distance_m,
    low_field_drift_m_per_s,
    low_field_parameter,
    migration_enthalpy_ev,
    mo,
    np,
    plt,
    relative_concentration_gradient,
    total_np_flux_mol_per_m2_s,
):
    field_coordinate = np.linspace(0.0, 2.0, 500)
    field_position_m = field_coordinate * jump_distance_m
    tilted_energy_ev = (
        0.5 * migration_enthalpy_ev
        * (1.0 - np.cos(2.0 * np.pi * field_coordinate))
        - charge_number_value * electric_field_v_per_m * field_position_m
    )
    field_figure, tilt_axis = plt.subplots(figsize=(11.5, 4.0), dpi=120)
    tilt_axis.plot(field_coordinate, tilted_energy_ev, color="#4C7C86", lw=1.7)
    site_energies = -charge_number_value * electric_field_v_per_m * (
        np.array([0.0, 1.0, 2.0]) * jump_distance_m
    )
    tilt_axis.scatter(
        [0.0, 1.0, 2.0], site_energies, s=65, color="#C49345",
        edgecolor="#40464D", zorder=4, label="Stable sites",
    )
    tilt_axis.set(
        xlabel=r"Position, $x/a$", ylabel="Energy (eV)",
        title="The electric field tilts the hopping landscape",
    )
    tilt_axis.grid(alpha=0.22)
    tilt_axis.legend(frameon=False)
    field_figure.tight_layout()
    plt.close(field_figure)

    gradient_values = 1.0e-6 * np.array([
        chemical_potential_gradient_j_per_mol_m,
        electrical_potential_gradient_j_per_mol_m,
        electrochemical_gradient_j_per_mol_m,
    ])
    flux_values = np.array([
        diffusion_flux_mol_per_m2_s, electrical_flux_mol_per_m2_s,
        total_np_flux_mol_per_m2_s,
    ])
    balance_figure, (_gradient_axis, _flux_axis) = plt.subplots(
        1, 2, figsize=(12.8, 4.5), dpi=120
    )
    _gradient_axis.bar(
        ["Chemical", "Electrical", "Total"], gradient_values,
        color=["#4C7C86", "#B8734A", "#C49345"],
    )
    _gradient_axis.axhline(0.0, color="#333333", lw=1.0)
    _gradient_axis.set(
        ylabel=r"Potential gradient (J mol$^{-1}$ $\mu$m$^{-1}$)",
        title="Potential gradients cancel",
    )
    _gradient_axis.grid(axis="y", alpha=0.22)
    _flux_axis.bar(
        ["Diffusion", "Electrical", "Total"], flux_values,
        color=["#4C7C86", "#B8734A", "#C49345"],
    )
    _flux_axis.axhline(0.0, color="#333333", lw=1.0)
    _flux_axis.set(
        ylabel=r"Molar flux (mol m$^{-2}$ s$^{-1}$)",
        title="The total flux is zero at equilibrium",
    )
    _flux_axis.grid(axis="y", alpha=0.22)
    balance_figure.tight_layout()
    plt.close(balance_figure)

    field_summary = mo.md(
        rf"""
        Neighboring sites differ by
        \(ze\mathcal{{E}}a={field_work_per_hop_ev:.3e}\) eV, so
        \(\Gamma_+/\Gamma_-={hopping_bias_ratio:.5g}\).
        The exact drift is **{exact_drift_m_per_s:.3e} m s^-1**; the low-field
        prediction is **{low_field_drift_m_per_s:.3e} m s^-1**.
        """
    )
    balance_details = mo.vstack([
        mo.md(r"""
        ### Chemical and electrical driving forces can cancel

        For a molar concentration \(c\),

        \[
        J=-D\frac{dc}{dx}-\frac{zFD}{RT}c\frac{d\phi}{dx}
          =-\frac{Dc}{RT}\frac{d\widetilde{\mu}}{dx},
        \qquad \widetilde{\mu}=\mu+zF\phi.
        \]

        A concentration gradient can therefore coexist with equilibrium:
        the chemical and electrical parts are nonzero, but their sum is zero.
        """),
        relative_concentration_gradient,
        balance_figure,
        mo.md(
            rf"The balancing field is "
            rf"\(\mathcal{{E}}={-demonstration_potential_gradient_v_per_m / 100.0:.3e}\) "
            r"V cm^-1 for the selected concentration gradient."
        ),
    ])
    mo.vstack([
        field_figure,
        field_summary,
        mo.accordion(
            {"Explore further - why equilibrium can have gradients": balance_details}
        ),
    ])
    return (field_figure,)

@app.cell
def _(mo):
    mo.accordion({
        "Model details - particle and molar flux notation": mo.md(r"""
        For number concentration \(c_N\),

        \[
        J_N=-D\frac{dc_N}{dx}+\frac{zeD}{k_BT}c_N\mathcal E .
        \]

        For molar concentration \(c\), use \(F=N_Ae\) and \(R=N_Ak_B\):

        \[
        J=-D\frac{dc}{dx}-\frac{zFD}{RT}c\frac{d\phi}{dx}.
        \]

        Particle-scale equations use \(k_B,e\); molar equations use \(R,F\).
        """)
    })
    return

@app.cell
def _(mo):
    mo.accordion({
        "Explore further - tracer, charge, and chemical diffusivity": mo.md(r"""
        The symbol \(D\) can describe different experiments:

        | notation | experiment | what moves? |
        |---|---|---|
        | \(D^*\) | isotope tracer profile | labeled atoms |
        | \(D^q\) | steady conductivity | charge-carrying ions |
        | \(D_{\rm Li}^{\delta}\) | chemical relaxation | \(\mathrm{Li^+}\) and \(e^-\) together |

        Correlation can make tracer and conductivity-derived values differ:
        \(D^*=H D^q\), where \(H\) is the Haven ratio. The core lesson below
        concerns chemical diffusion, where the composition itself changes.
        """)
    })
    return

@app.cell
def _(mo):
    log_ionic_diffusivity_cm2 = mo.ui.slider(
        start=-14.0,
        stop=-6.0,
        step=0.25,
        value=-10.0,
        label="Lithium-ion diffusivity exponent",
        show_value=True,
    )
    log_electronic_to_ionic_ratio = mo.ui.slider(
        start=-4.0,
        stop=4.0,
        step=0.25,
        value=3.0,
        label="Electron/ion diffusivity-ratio exponent",
        show_value=True,
    )
    ambipolar_controls = mo.hstack(
        [log_ionic_diffusivity_cm2, log_electronic_to_ionic_ratio],
        justify="start",
        align="center",
        wrap=True,
        gap=1.5,
    )
    return (
        ambipolar_controls,
        log_electronic_to_ionic_ratio,
        log_ionic_diffusivity_cm2,
    )


@app.cell
def _(ambipolar_controls, mo):
    mo.vstack(
        [
            mo.md(r"""
            ## 4. Chemical diffusion: two carriers, one bottleneck

            **If electrons become much faster than ions, can the chemical
            diffusivity grow without limit?**

            Consider the simplest neutral-composition reaction:

            \[
            \mathrm{Li}\rightleftharpoons \mathrm{Li^+}+e^- .
            \]

            In the bulk, local electroneutrality gives
            \(c_{\mathrm{Li^+}}=c_{e^-}=c_{\mathrm{Li}}\). Local equilibrium
            connects the neutral and charged chemical potentials:

            \[
            \mu_{\mathrm{Li}}
            =\widetilde{\mu}_{\mathrm{Li^+}}
            +\widetilde{\mu}_{e^-}.
            \]

            Open-circuit composition transport carries no net current. For this
            monovalent pair, the ion and electron therefore share one molar flux:

            \[
            J_{\mathrm{Li}}
            =J_{\mathrm{Li^+}}
            =J_{e^-}.
            \]

            If electrons would diffuse faster on their own, a small internal
            electric field slows them and speeds up \(\mathrm{Li^+}\). Only a
            minute charge displacement is needed to establish this field; the
            bulk remains locally electroneutral in the model.
            """),
            ambipolar_controls,
        ]
    )
    return


@app.cell
def _(
    FARADAY_C_PER_MOL,
    ambipolar_fluxes,
    ambipolar_internal_potential_gradient,
    conductivity_from_diffusivity,
    lithium_chemical_diffusivity,
    lithium_chemical_diffusivity_from_conductivities,
    log_electronic_to_ionic_ratio,
    log_ionic_diffusivity_cm2,
    temperature_k,
):
    ionic_diffusivity_cm2_per_s = 10.0 ** float(
        log_ionic_diffusivity_cm2.value
    )
    ionic_diffusivity_m2_per_s = ionic_diffusivity_cm2_per_s * 1.0e-4
    electronic_to_ionic_ratio = 10.0 ** float(
        log_electronic_to_ionic_ratio.value
    )
    electronic_diffusivity_m2_per_s = (
        ionic_diffusivity_m2_per_s * electronic_to_ionic_ratio
    )
    electronic_diffusivity_cm2_per_s = (
        electronic_diffusivity_m2_per_s * 1.0e4
    )

    ambipolar_concentration_mol_per_m3 = 1000.0
    ambipolar_relative_gradient_per_m = 2.0e5
    ambipolar_gradient_mol_per_m4 = (
        ambipolar_concentration_mol_per_m3
        * ambipolar_relative_gradient_per_m
    )
    internal_potential_gradient_v_per_m = (
        ambipolar_internal_potential_gradient(
            ionic_diffusivity_m2_per_s,
            electronic_diffusivity_m2_per_s,
            ambipolar_concentration_mol_per_m3,
            ambipolar_gradient_mol_per_m4,
            temperature_k,
        )
    )
    ionic_flux_mol_per_m2_s, electronic_flux_mol_per_m2_s = (
        ambipolar_fluxes(
            ionic_diffusivity_m2_per_s,
            electronic_diffusivity_m2_per_s,
            ambipolar_concentration_mol_per_m3,
            ambipolar_gradient_mol_per_m4,
            internal_potential_gradient_v_per_m,
            temperature_k,
        )
    )
    open_circuit_current_a_per_m2 = FARADAY_C_PER_MOL * (
        ionic_flux_mol_per_m2_s - electronic_flux_mol_per_m2_s
    )

    lithium_chemical_diffusivity_m2_per_s = lithium_chemical_diffusivity(
        ionic_diffusivity_m2_per_s,
        electronic_diffusivity_m2_per_s,
    )
    lithium_chemical_diffusivity_cm2_per_s = (
        lithium_chemical_diffusivity_m2_per_s * 1.0e4
    )
    analytic_common_flux_mol_per_m2_s = (
        -lithium_chemical_diffusivity_m2_per_s
        * ambipolar_gradient_mol_per_m4
    )

    ionic_conductivity_s_per_m = conductivity_from_diffusivity(
        ionic_diffusivity_m2_per_s,
        ambipolar_concentration_mol_per_m3,
        temperature_k,
    )
    electronic_conductivity_s_per_m = conductivity_from_diffusivity(
        electronic_diffusivity_m2_per_s,
        ambipolar_concentration_mol_per_m3,
        temperature_k,
    )
    conductivity_form_diffusivity_m2_per_s = (
        lithium_chemical_diffusivity_from_conductivities(
            ionic_conductivity_s_per_m,
            electronic_conductivity_s_per_m,
            ambipolar_concentration_mol_per_m3,
            ambipolar_concentration_mol_per_m3,
            temperature_k,
        )
    )
    return (
        ambipolar_concentration_mol_per_m3,
        ambipolar_gradient_mol_per_m4,
        ambipolar_relative_gradient_per_m,
        analytic_common_flux_mol_per_m2_s,
        conductivity_form_diffusivity_m2_per_s,
        electronic_conductivity_s_per_m,
        electronic_diffusivity_cm2_per_s,
        electronic_diffusivity_m2_per_s,
        electronic_flux_mol_per_m2_s,
        electronic_to_ionic_ratio,
        internal_potential_gradient_v_per_m,
        ionic_conductivity_s_per_m,
        ionic_diffusivity_cm2_per_s,
        ionic_diffusivity_m2_per_s,
        ionic_flux_mol_per_m2_s,
        lithium_chemical_diffusivity_cm2_per_s,
        lithium_chemical_diffusivity_m2_per_s,
        open_circuit_current_a_per_m2,
    )


@app.cell
def _(
    ambipolar_gradient_mol_per_m4,
    analytic_common_flux_mol_per_m2_s,
    electronic_diffusivity_m2_per_s,
    electronic_flux_mol_per_m2_s,
    electronic_to_ionic_ratio,
    internal_potential_gradient_v_per_m,
    ionic_diffusivity_m2_per_s,
    ionic_flux_mol_per_m2_s,
    lithium_chemical_diffusivity_cm2_per_s,
    mo,
    np,
    open_circuit_current_a_per_m2,
    plt,
):
    ratio_curve = np.logspace(-4.0, 4.0, 500)
    chemical_to_ionic_ratio_curve = 2.0 * ratio_curve / (1.0 + ratio_curve)

    uncoupled_ionic_flux = (
        -ionic_diffusivity_m2_per_s * ambipolar_gradient_mol_per_m4
    )
    uncoupled_electronic_flux = (
        -electronic_diffusivity_m2_per_s * ambipolar_gradient_mol_per_m4
    )
    flux_magnitudes = np.abs(
        [
            uncoupled_ionic_flux,
            uncoupled_electronic_flux,
            analytic_common_flux_mol_per_m2_s,
        ]
    )
    flux_floor = max(np.max(flux_magnitudes) * 1.0e-6, 1.0e-300)
    flux_magnitudes = np.maximum(flux_magnitudes, flux_floor)

    ambipolar_figure, (ratio_axis, coupling_axis) = plt.subplots(
        1,
        2,
        figsize=(13.8, 4.9),
        dpi=120,
    )
    ratio_axis.loglog(
        ratio_curve,
        chemical_to_ionic_ratio_curve,
        color="#4C7C86",
        lw=1.9,
        label=r"$2r/(1+r)$",
    )
    ratio_axis.scatter(
        [electronic_to_ionic_ratio],
        [
            2.0
            * electronic_to_ionic_ratio
            / (1.0 + electronic_to_ionic_ratio)
        ],
        s=95,
        color="#C49345",
        edgecolor="#40464D",
        zorder=5,
        label="selected values",
    )
    ratio_axis.axhline(1.0, color="#858B90", lw=1.0, ls=":")
    ratio_axis.set(
        xlabel=r"Mobility contrast, $r=D_{e^-}/D_{\rm Li^+}$",
        ylabel=r"$D_{\rm Li}^{\delta}/D_{\rm Li^+}$",
        title="The slower carrier limits chemical diffusion",
        xlim=(1.0e-4, 1.0e4),
        ylim=(1.0e-4, 10.0),
    )
    ratio_axis.grid(which="both", alpha=0.22)
    ratio_axis.legend(frameon=False)

    coupling_axis.bar(
        [
            r"$\mathrm{Li^+}$ alone",
            r"$e^-$ alone",
            "common flux\n" + r"$J_{\rm Li^+}=J_{e^-}$",
        ],
        flux_magnitudes,
        color=["#4C7C86", "#B65C4A", "#C49345"],
    )
    coupling_axis.set_yscale("log")
    coupling_axis.set(
        ylabel=r"|flux| (mol m$^{-2}$ s$^{-1}$)",
        title=r"The internal field makes $J_{\rm Li^+}=J_{e^-}$",
    )
    coupling_axis.tick_params(axis="x", rotation=12)
    coupling_axis.grid(axis="y", which="both", alpha=0.22)
    ambipolar_figure.tight_layout()
    plt.close(ambipolar_figure)

    ambipolar_summary = mo.md(
        rf"""
        The internal field forces the ion and electron to share one
        neutral-composition flux. The internal potential gradient is
        **{internal_potential_gradient_v_per_m / 100.0:.3e} V/cm**. It gives
        \(J_{{\rm Li^+}}={ionic_flux_mol_per_m2_s:.3e}\) and
        \(J_{{e^-}}={electronic_flux_mol_per_m2_s:.3e}\) mol/(m² s), while
        \(F(J_{{\rm Li^+}}-J_{{e^-}})
        ={open_circuit_current_a_per_m2:.3e}\) A/m².

        The common flux corresponds to
        \(D_{{\rm Li}}^{{\delta}}
        ={lithium_chemical_diffusivity_cm2_per_s:.3e}\) cm² s⁻¹.
        """
    )
    mo.vstack([ambipolar_figure, ambipolar_summary])
    return (ambipolar_figure,)


@app.cell
def _(
    conductivity_form_diffusivity_m2_per_s,
    electronic_conductivity_s_per_m,
    electronic_diffusivity_cm2_per_s,
    ionic_conductivity_s_per_m,
    ionic_diffusivity_cm2_per_s,
    lithium_chemical_diffusivity_cm2_per_s,
    mo,
):
    derivation = mo.md(
        rf"""
        ### Deriving the common chemical diffusivity

        For either carrier, the one-dimensional flux can be written in terms of
        its electrochemical-potential gradient:

        \[
        J_s=-\frac{{\sigma_s}}{{z_s^2F^2}}
        \frac{{d\widetilde{{\mu}}_s}}{{dx}}.
        \]

        For the monovalent Li-ion/electron pair, both carrier fluxes equal the
        common neutral-composition flux $J$. Differentiating the local-equilibrium
        relation then gives

        \[
        \frac{{d\mu_{{\rm Li}}}}{{dx}}
        =-F^2J\left(\frac{{1}}{{\sigma_{{\rm Li^+}}}}
        +\frac{{1}}{{\sigma_{{e^-}}}}\right),
        \]

        or

        \[
        J=-\frac{{\sigma_{{\rm Li^+}}\sigma_{{e^-}}}}
        {{F^2(\sigma_{{\rm Li^+}}+\sigma_{{e^-}})}}
        \frac{{d\mu_{{\rm Li}}}}{{dx}}.
        \]

        The model now makes its ideal-dilute assumption explicit:
        $c_{{\rm Li^+}}=c_{{e^-}}=c$ and
        $\mu_{{\rm Li}}=\mu_{{\rm Li}}^0+2RT\ln(c/c_0)$. Identifying
        $J=-D_{{\rm Li}}^\delta dc/dx$ gives

        \[
        D_{{\rm Li}}^\delta
        =\frac{{RT}}{{F^2}}
        \frac{{\sigma_{{e^-}}\sigma_{{\rm Li^+}}}}
             {{\sigma_{{e^-}}+\sigma_{{\rm Li^+}}}}
        \left(\frac{{1}}{{c_{{e^-}}}}
        +\frac{{1}}{{c_{{\rm Li^+}}}}\right).
        \]

        Using the Nernst–Einstein relation for each carrier reduces this to

        \[
        \boxed{{D_{{\rm Li}}^\delta
        =\frac{{2D_{{\rm Li^+}}D_{{e^-}}}}
        {{D_{{\rm Li^+}}+D_{{e^-}}}}}}.
        \]

        At the selected values,
        $D_{{\rm Li^+}}={ionic_diffusivity_cm2_per_s:.3e}$ cm² s⁻¹ and
        $D_{{e^-}}={electronic_diffusivity_cm2_per_s:.3e}$ cm² s⁻¹, with
        $\sigma_{{\rm Li^+}}={ionic_conductivity_s_per_m:.3e}$ S/m and
        $\sigma_{{e^-}}={electronic_conductivity_s_per_m:.3e}$ S/m.
        The conductivity form gives
        **{conductivity_form_diffusivity_m2_per_s * 1.0e4:.3e} cm² s⁻¹**, matching
        $D_{{\rm Li}}^\delta
        ={lithium_chemical_diffusivity_cm2_per_s:.3e}$ cm² s⁻¹.

        The factor of two belongs to this ideal, dilute, locally neutral pair.
        Outside that limit, the composition dependence of
        $\mu_{{\rm Li}}$ must be supplied separately.
        """
    )
    mo.accordion({"Model details - deriving the common diffusivity": derivation})
    return


@app.cell
def _(mo):
    log_sample_length = mo.ui.slider(
        start=-8.0,
        stop=-2.0,
        step=0.25,
        value=-5.0,
        label="Sample-length exponent (base 10, meters)",
        show_value=True,
    )
    reduced_profile_time = mo.ui.slider(
        start=0.01,
        stop=0.50,
        step=0.01,
        value=0.08,
        label="Reduced diffusion time",
        show_value=True,
    )
    relaxation_controls = mo.hstack(
        [log_sample_length, reduced_profile_time],
        justify="start",
        align="center",
        wrap=True,
        gap=1.5,
    )
    return log_sample_length, reduced_profile_time, relaxation_controls


@app.cell
def _(
    characteristic_relaxation_time,
    lithium_chemical_diffusivity_m2_per_s,
    log_sample_length,
    np,
    reduced_profile_time,
    slab_remaining_profile,
):
    selected_length_m = 10.0 ** float(log_sample_length.value)
    selected_relaxation_time_s = characteristic_relaxation_time(
        selected_length_m,
        lithium_chemical_diffusivity_m2_per_s,
    )
    length_curve_m = np.logspace(-8.0, -2.0, 400)
    relaxation_curve_s = (
        length_curve_m**2 / lithium_chemical_diffusivity_m2_per_s
    )
    slab_position = np.linspace(0.0, 1.0, 400)
    selected_slab_profile = slab_remaining_profile(
        slab_position,
        float(reduced_profile_time.value),
    )
    return (
        length_curve_m,
        relaxation_curve_s,
        selected_length_m,
        selected_relaxation_time_s,
        selected_slab_profile,
        slab_position,
    )


@app.cell
def _(
    length_curve_m,
    lithium_chemical_diffusivity_cm2_per_s,
    mo,
    np,
    plt,
    reduced_profile_time,
    relaxation_controls,
    relaxation_curve_s,
    selected_length_m,
    selected_relaxation_time_s,
    selected_slab_profile,
    slab_position,
):
    relaxation_figure, (time_axis, relaxation_profile_axis) = plt.subplots(
        1,
        2,
        figsize=(13.6, 4.8),
        dpi=120,
    )
    time_axis.loglog(
        length_curve_m,
        relaxation_curve_s,
        color="#4C7C86",
        lw=1.9,
        label=r"scaling time $t_D=L^2/D_{\rm Li}^{\delta}$",
    )
    time_axis.loglog(
        length_curve_m,
        relaxation_curve_s / np.pi**2,
        color="#B8734A",
        lw=1.8,
        ls="--",
        label=r"first mode $\tau^\delta=t_D/\pi^2$",
    )
    time_axis.scatter(
        [selected_length_m],
        [selected_relaxation_time_s],
        s=95,
        color="#C49345",
        edgecolor="#40464D",
        zorder=5,
        label="selected $t_D$",
    )
    time_axis.set(
        xlabel=r"Sample length, $L$ (m)",
        ylabel="Characteristic time (s)",
        title="Diffusion time grows as length squared",
    )
    time_axis.grid(which="both", alpha=0.22)
    time_axis.legend(frameon=False)

    relaxation_profile_axis.plot(
        slab_position,
        selected_slab_profile,
        color="#B65C4A",
        lw=1.9,
    )
    relaxation_profile_axis.set(
        xlabel=r"Position, $x/L$",
        ylabel="Remaining normalized composition change",
        title="One-dimensional slab relaxation",
        ylim=(-0.03, 1.03),
    )
    relaxation_profile_axis.grid(alpha=0.22)
    relaxation_figure.tight_layout()
    plt.close(relaxation_figure)

    if selected_relaxation_time_s < 1.0:
        time_text = f"{selected_relaxation_time_s:.3e} s"
    elif selected_relaxation_time_s < 3600.0:
        time_text = f"{selected_relaxation_time_s / 60.0:.2f} min"
    elif selected_relaxation_time_s < 86400.0:
        time_text = f"{selected_relaxation_time_s / 3600.0:.2f} h"
    else:
        time_text = f"{selected_relaxation_time_s / 86400.0:.2f} days"

    selected_first_mode_time_s = selected_relaxation_time_s / np.pi**2
    if selected_first_mode_time_s < 1.0:
        first_mode_text = f"{selected_first_mode_time_s:.3e} s"
    elif selected_first_mode_time_s < 3600.0:
        first_mode_text = f"{selected_first_mode_time_s / 60.0:.2f} min"
    elif selected_first_mode_time_s < 86400.0:
        first_mode_text = f"{selected_first_mode_time_s / 3600.0:.2f} h"
    else:
        first_mode_text = f"{selected_first_mode_time_s / 86400.0:.2f} days"

    relaxation_summary = mo.md(
        rf"""
        With \(D_{{\rm Li}}^\delta
        ={lithium_chemical_diffusivity_cm2_per_s:.3e}\) cm² s⁻¹ and
        \(L={selected_length_m:.3e}\) m, $t_D$ is **{time_text}** and
        $\tau^\delta$ is **{first_mode_text}**.
        The right panel shows the remaining composition change at Fourier
        number \(\theta={float(reduced_profile_time.value):.2f}\).
        """
    )
    relaxation_reader = mo.vstack([
        mo.md(r"""
        Chemical diffusivity also sets a sample-scale clock:

        \[
        t_D=\frac{L^2}{D_{\rm Li}^{\delta}},\qquad
        \tau^\delta=\frac{L^2}{\pi^2D_{\rm Li}^{\delta}}.
        \]

        This is why the same material responds rapidly as a thin film and
        slowly as a bulk sample.
        """),
        relaxation_controls,
        relaxation_figure,
        relaxation_summary,
    ])
    mo.accordion(
        {"Explore further - sample length and relaxation time": relaxation_reader}
    )
    return (relaxation_figure,)


@app.cell
def _(
    FARADAY_C_PER_MOL,
    KB_EV_PER_K,
    ambipolar_gradient_mol_per_m4,
    analytic_common_flux_mol_per_m2_s,
    attempt_frequency_hz,
    backward_rate_hz,
    biased_directional_rates,
    charge_number_value,
    conductivity_form_diffusivity_m2_per_s,
    defect_diffusivity_m2_per_s,
    diffusion_flux_mol_per_m2_s,
    electronic_diffusivity_m2_per_s,
    electronic_flux_mol_per_m2_s,
    electrical_flux_mol_per_m2_s,
    exact_drift_m_per_s,
    exact_hopping_drift_velocity,
    extracted_diffusivity_m2_per_s,
    field_work_per_hop_ev,
    forward_rate_hz,
    hop_frequency_hz,
    ionic_diffusivity_m2_per_s,
    ionic_flux_mol_per_m2_s,
    jump_distance_m,
    lithium_chemical_diffusivity_m2_per_s,
    low_field_drift_m_per_s,
    master_flux_relative_error,
    master_mass_relative_error,
    migration_enthalpy_ev,
    msd_fit_r_squared,
    nernst_einstein_drift_velocity,
    np,
    open_circuit_current_a_per_m2,
    temperature_k,
    total_np_flux_mol_per_m2_s,
):
    def _relative_error(value, reference):
        return abs(float(value) - float(reference)) / max(
            abs(float(reference)),
            1.0e-300,
        )

    zero_field_identity_error = _relative_error(
        defect_diffusivity_m2_per_s,
        0.5 * jump_distance_m**2 * hop_frequency_hz,
    )
    msd_diffusivity_error = _relative_error(
        extracted_diffusivity_m2_per_s,
        defect_diffusivity_m2_per_s,
    )
    detailed_balance_error = _relative_error(
        forward_rate_hz / backward_rate_hz,
        np.exp(field_work_per_hop_ev / (KB_EV_PER_K * temperature_k)),
    )
    _test_field_v_per_m = 1.0
    _test_forward_rate_hz, _test_backward_rate_hz, _ = biased_directional_rates(
        temperature_k,
        migration_enthalpy_ev,
        attempt_frequency_hz,
        jump_distance_m,
        charge_number_value,
        _test_field_v_per_m,
    )
    _test_exact_drift_m_per_s = exact_hopping_drift_velocity(
        jump_distance_m,
        _test_forward_rate_hz,
        _test_backward_rate_hz,
    )
    _test_low_field_drift_m_per_s = nernst_einstein_drift_velocity(
        defect_diffusivity_m2_per_s,
        temperature_k,
        charge_number_value,
        _test_field_v_per_m,
    )
    low_field_drift_error = _relative_error(
        _test_exact_drift_m_per_s,
        _test_low_field_drift_m_per_s,
    )

    flux_balance_scale = max(
        abs(diffusion_flux_mol_per_m2_s)
        + abs(electrical_flux_mol_per_m2_s),
        1.0e-300,
    )
    electrochemical_flux_relative = (
        abs(total_np_flux_mol_per_m2_s) / flux_balance_scale
    )
    ambipolar_flux_match_error = _relative_error(
        ionic_flux_mol_per_m2_s,
        electronic_flux_mol_per_m2_s,
    )
    ambipolar_current_relative = abs(open_circuit_current_a_per_m2) / max(
        FARADAY_C_PER_MOL
        * (abs(ionic_flux_mol_per_m2_s)
        + abs(electronic_flux_mol_per_m2_s)),
        1.0e-300,
    )
    analytic_chemical_flux_error = _relative_error(
        ionic_flux_mol_per_m2_s,
        analytic_common_flux_mol_per_m2_s,
    )
    conductivity_form_error = _relative_error(
        conductivity_form_diffusivity_m2_per_s,
        lithium_chemical_diffusivity_m2_per_s,
    )

    positive_values = np.array(
        [
            hop_frequency_hz,
            defect_diffusivity_m2_per_s,
            ionic_diffusivity_m2_per_s,
            electronic_diffusivity_m2_per_s,
            lithium_chemical_diffusivity_m2_per_s,
        ]
    )
    finite_values = np.array(
        [
            zero_field_identity_error,
            msd_diffusivity_error,
            detailed_balance_error,
            low_field_drift_error,
            master_mass_relative_error,
            master_flux_relative_error,
            electrochemical_flux_relative,
            ambipolar_flux_match_error,
            ambipolar_current_relative,
            analytic_chemical_flux_error,
            conductivity_form_error,
            msd_fit_r_squared,
            ambipolar_gradient_mol_per_m4,
        ]
    )
    transport_validation = {
        "zero_field_identity_error": zero_field_identity_error,
        "zero_field_identity_pass": zero_field_identity_error < 1.0e-14,
        "msd_diffusivity_error": msd_diffusivity_error,
        "msd_r_squared": msd_fit_r_squared,
        "msd_pass": (
            msd_diffusivity_error < 0.05
            and msd_fit_r_squared > 0.995
        ),
        "detailed_balance_error": detailed_balance_error,
        "detailed_balance_pass": detailed_balance_error < 1.0e-12,
        "low_field_drift_error": low_field_drift_error,
        "low_field_drift_pass": low_field_drift_error < 1.0e-7,
        "mass_error": master_mass_relative_error,
        "mass_pass": master_mass_relative_error < 1.0e-12,
        "fick_error": master_flux_relative_error,
        "fick_pass": master_flux_relative_error < 1.0e-12,
        "electrochemical_flux_relative": electrochemical_flux_relative,
        "electrochemical_flux_pass": electrochemical_flux_relative < 1.0e-12,
        "ambipolar_current_relative": ambipolar_current_relative,
        "ambipolar_current_pass": ambipolar_current_relative < 1.0e-12,
        "ambipolar_flux_match_error": ambipolar_flux_match_error,
        "ambipolar_flux_match_pass": ambipolar_flux_match_error < 1.0e-12,
        "analytic_chemical_flux_error": analytic_chemical_flux_error,
        "analytic_chemical_flux_pass": analytic_chemical_flux_error < 1.0e-12,
        "conductivity_form_error": conductivity_form_error,
        "conductivity_form_pass": conductivity_form_error < 1.0e-12,
        "positive_pass": bool(np.all(positive_values > 0.0)),
        "finite_pass": bool(np.all(np.isfinite(finite_values))),
    }
    return (transport_validation,)


@app.cell
def _(mo, transport_validation):
    def _check_mark(passed):
        return "PASS" if passed else "CHECK"

    _checks = mo.md(
        rf"""
        ## Physical consistency checks

        | status | physical statement | why it matters |
        |---:|---|---|
        | {_check_mark(transport_validation['zero_field_identity_pass'] and transport_validation['msd_pass'])} | atomic hops and the random walk give \(D=a^2\Gamma/2\) | microscopic motion connects to measurable diffusion |
        | {_check_mark(transport_validation['detailed_balance_pass'] and transport_validation['low_field_drift_pass'])} | biased hops obey detailed balance and approach Nernst–Einstein at low field | the field response has the correct equilibrium limit |
        | {_check_mark(transport_validation['mass_pass'] and transport_validation['fick_pass'])} | one-dimensional transport conserves defects and follows Fick's law | the concentration profile and flux describe the same motion |
        | {_check_mark(transport_validation['electrochemical_flux_pass'])} | chemical and electrical driving forces cancel at equilibrium | a flat electrochemical potential means zero total flux |
        | {_check_mark(transport_validation['ambipolar_flux_match_pass'] and transport_validation['ambipolar_current_pass'])} | Li⁺ and electrons carry a common chemical-diffusion flux with zero current | local neutrality couples the two carriers |
        | {_check_mark(transport_validation['analytic_chemical_flux_pass'] and transport_validation['conductivity_form_pass'])} | the two forms of \(D_{{\rm Li}}^\delta\) agree | transport and conductivity descriptions are consistent |
        | {_check_mark(transport_validation['positive_pass'] and transport_validation['finite_pass'])} | all rates and diffusivities are positive and finite | every displayed quantity is physical |

        These checks follow the main physical chain from atomic hopping to
        coupled chemical diffusion; they do not add new assumptions.
        """
    )
    mo.accordion({"Physical consistency checks": _checks})
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## What to carry forward

    Random hops produce \(D=a^2\Gamma/2\) in one dimension. A concentration
    difference creates a net flux because more defects exchange from the
    crowded side, and an electric field biases the two directional hop rates.

    For the ideal locally neutral pair
    \(\mathrm{Li}\rightleftharpoons\mathrm{Li^+}+e^-\), ions and electrons must
    share one flux. Their chemical diffusivity is therefore

    \[
    D_{\rm Li}^{\delta}
    =\frac{2D_{\rm Li^+}D_{e^-}}{D_{\rm Li^+}+D_{e^-}},
    \]

    so the slower carrier sets the bottleneck.
    """)
    return

@app.cell
def _(mo):
    mo.md(r"""
    **Continue:** [Module 04 — Space-Charge Layers and the Frumkin Effect](https://qiyanglu.github.io/Solid-State-Ionics-Interactive/04-space-charge-frumkin/)
    """)
    return



if __name__ == "__main__":
    app.run()
