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
            "font.size": 14,
            "axes.titlesize": 16,
            "axes.labelsize": 14,
            "xtick.labelsize": 12,
            "ytick.labelsize": 12,
            "legend.fontsize": 11,
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
    # Defect Transport: From Atomic Hopping to Chemical Diffusion

    **Guiding question.** How does a thermally activated atomic jump become a
    measurable transport coefficient—and why can a composition profile relax
    only when ions and electrons move together?

    **Learning goals**

    1. Connect an activated hop rate to a one-dimensional diffusivity and
       Fickian flux.
    2. Read transport as motion down an electrochemical-potential gradient.
    3. Explain why neutral chemical diffusion couples ionic and electronic
       carriers and is limited by the slower one.

    > **Predict before exploring.** If the electron diffusivity becomes much
    > larger than the Li-ion diffusivity, can $D_{\rm Li}^{\delta}$ grow without
    > limit, or does it approach a bottleneck value?

    **Model scope.** All spatial pictures are one-dimensional. Particle-scale
    equations use $k_B,e$; molar equations use $R,F$. The electric field is
    $\mathcal E=-\partial\phi/\partial x$. See the shared
    [notation bridge](https://github.com/qiyanglu/Solid-State-Ionics-Interactive/blob/main/NOTATION.md).

    **Flux and current bridge.** Positive flux points toward $+x$:

    \[
    J_N\;[\mathrm{particles\,m^{-2}\,s^{-1}}],\qquad
    J=J_N/N_A\;[\mathrm{mol\,m^{-2}\,s^{-1}}],
    \]
    \[
    j=zFJ\;[\mathrm{A\,m^{-2}}],\qquad I=jS\;[\mathrm A].
    \]

    This module follows one connected argument:

    \[
    \text{activated hops}
    \rightarrow \langle x^2\rangle
    \rightarrow D
    \rightarrow J
    \rightarrow \widetilde{\mu}
    \rightarrow D_{\rm Li}^{\delta}
    \rightarrow \tau^\delta .
    \]

    We begin with an ideal one-dimensional lattice, then add concentration and
    electrical driving forces, and finally couple \(\mathrm{Li^+}\) and
    \(e^-\) under local electroneutrality. Every spatial equation and simulation
    stays in one dimension, following the lecture derivations.

    Lecture-facing controls and results use **K, eV, nm, s, V/cm, and cm²/s**.
    Molar fluxes are reported in mol/(m² s), with \(R\) and \(F\) used for
    molar equations.
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

        Gamma is the total zero-field hop frequency used in the lecture slides.
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

    def simulate_zero_field_walks(
        jump_distance_m,
        hop_frequency_hz,
        walker_count=2000,
        step_count=320,
        seed=2026,
    ):
        """Simulate unbiased 1D walks with mean interval 1/Gamma.

        Each observation interval contains one jump, chosen toward +x or -x
        with equal probability.
        """
        distance = require_positive("jump_distance_m", jump_distance_m)
        frequency = require_positive("hop_frequency_hz", hop_frequency_hz)
        walkers = int(walker_count)
        steps = int(step_count)
        if walkers < 2 or steps < 2:
            raise ValueError("walker_count and step_count must be at least 2")
        rng = np.random.default_rng(int(seed))
        increments = rng.choice((-1.0, 1.0), size=(walkers, steps))
        lattice_positions = np.concatenate(
            [np.zeros((walkers, 1)), np.cumsum(increments, axis=1)],
            axis=1,
        )
        positions_m = distance * lattice_positions
        times_s = np.arange(steps + 1, dtype=float) / frequency
        msd_m2 = np.mean(positions_m**2, axis=0)
        return times_s, positions_m, msd_m2

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
        """Return the dilute-limit D_Li^delta expression used in Lecture 5."""
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
        label="Migration enthalpy, ΔH_mig (eV)",
        show_value=True,
    )
    log_attempt_frequency = mo.ui.slider(
        start=10.0,
        stop=14.0,
        step=0.25,
        value=13.0,
        label="log10 attempt frequency, ν (s^-1)",
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
        [temperature, migration_barrier, log_attempt_frequency, jump_distance],
        justify="start",
        align="center",
        wrap=True,
        gap=1.5,
    )
    return (
        hopping_controls,
        jump_distance,
        log_attempt_frequency,
        migration_barrier,
        temperature,
    )


@app.cell
def _(hopping_controls, mo):
    mo.vstack(
        [
            mo.md(r"""
            ## 1. One activated hop

            Use the same notation as the lecture slides:

            \[
            \Gamma=\nu\exp\!\left(-\frac{\Delta H_{\rm mig}}{k_BT}\right),
            \qquad
            D=\frac{1}{2}a^2\Gamma \quad \text{(one dimension)}.
            \]

            Here \(\Gamma\) is the total hop frequency and \(1/\Gamma\) is the
            mean time between hops. Each hop goes either left or right with
            equal probability when no electric field is present. Therefore

            **Simulation convention.** The embedded random walk places one hop
            in each fixed observation interval $1/\Gamma$. It is a transparent
            demonstration of the spatial statistics, not a continuous-time
            Poisson clock for the waiting times.

            \[
            \langle x^2\rangle=2Dt.
            \]

            This module stays in one dimension throughout, matching the
            derivation used in class.
            """),
            hopping_controls,
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
    walk_times_s, walk_positions_m, walk_msd_m2 = simulate_zero_field_walks(
        jump_distance_m,
        hop_frequency_hz,
    )
    fit_start = walk_times_s.size // 5
    msd_slope, msd_intercept = np.polyfit(
        walk_times_s[fit_start:],
        walk_msd_m2[fit_start:],
        1,
    )
    extracted_diffusivity_m2_per_s = msd_slope / 2.0
    fitted_msd_m2 = msd_slope * walk_times_s + msd_intercept
    residual_sum_squares = float(
        np.sum((walk_msd_m2[fit_start:] - fitted_msd_m2[fit_start:]) ** 2)
    )
    total_sum_squares = float(
        np.sum(
            (
                walk_msd_m2[fit_start:]
                - np.mean(walk_msd_m2[fit_start:])
            )
            ** 2
        )
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
        walk_msd_m2,
        walk_positions_m,
        walk_times_s,
        hop_frequency_hz,
    )


@app.cell
def _(
    defect_diffusivity_m2_per_s,
    extracted_diffusivity_m2_per_s,
    jump_distance_m,
    migration_enthalpy_ev,
    msd_fit_r_squared,
    np,
    plt,
    scaled_time_axis,
    walk_msd_m2,
    walk_positions_m,
    walk_times_s,
    hop_frequency_hz,
):
    plt.rcParams.update(
        {
            "font.size": 14,
            "axes.titlesize": 16,
            "axes.labelsize": 14,
            "legend.fontsize": 11,
        }
    )
    display_time, time_unit = scaled_time_axis(walk_times_s)
    reaction_coordinate = np.linspace(0.0, 2.0, 500)
    barrier_energy_ev = 0.5 * migration_enthalpy_ev * (
        1.0 - np.cos(2.0 * np.pi * reaction_coordinate)
    )

    microscopic_figure, microscopic_axes = plt.subplots(
        1,
        3,
        figsize=(15.0, 4.8),
        dpi=120,
    )
    barrier_axis, trajectory_axis, msd_axis = microscopic_axes
    barrier_axis.plot(
        reaction_coordinate,
        barrier_energy_ev,
        color="#4C7C86",
        lw=3.0,
    )
    barrier_axis.scatter(
        [0.0, 1.0, 2.0],
        [0.0, 0.0, 0.0],
        s=80,
        color="#C49345",
        edgecolor="#40464D",
        zorder=4,
        label="equivalent sites",
    )
    barrier_axis.annotate(
        "",
        xy=(0.5, migration_enthalpy_ev),
        xytext=(0.5, 0.0),
        arrowprops={"arrowstyle": "<->", "color": "#B65C4A", "lw": 1.8},
    )
    barrier_axis.text(
        0.55,
        0.52 * migration_enthalpy_ev,
        r"$H_{\rm mig}$",
        color="#B65C4A",
    )
    barrier_axis.set(
        xlabel="position / a",
        ylabel="energy (eV)",
        title="Thermally activated hops",
        ylim=(-0.04 * migration_enthalpy_ev, 1.15 * migration_enthalpy_ev),
    )
    barrier_axis.legend(frameon=False, loc="upper right")
    barrier_axis.grid(alpha=0.2)

    for trajectory_index in range(8):
        trajectory_axis.step(
            display_time,
            walk_positions_m[trajectory_index] / jump_distance_m,
            where="post",
            lw=1.2,
            alpha=0.75,
        )
    trajectory_axis.set(
        xlabel=f"time ({time_unit})",
        ylabel="position / a",
        title="Unbiased paths wander",
    )
    trajectory_axis.grid(alpha=0.22)

    msd_axis.plot(
        display_time,
        walk_msd_m2 * 1.0e18,
        color="#4C7C86",
        lw=2.8,
        label="deterministic-seed simulation",
    )
    msd_axis.plot(
        display_time,
        2.0 * defect_diffusivity_m2_per_s * walk_times_s * 1.0e18,
        color="#B8734A",
        lw=2.0,
        ls="--",
        label=r"$2Dt$",
    )
    msd_axis.set(
        xlabel=f"time ({time_unit})",
        ylabel=r"$\langle x^2\rangle$ (nm$^2$)",
        title="The ensemble reveals D",
    )
    msd_axis.grid(alpha=0.22)
    msd_axis.legend(frameon=False)
    microscopic_figure.tight_layout()
    plt.close(microscopic_figure)

    microscopic_summary = (
        f"Gamma = {hop_frequency_hz:.3e} s^-1; "
        f"D = a^2 Gamma / 2 = {defect_diffusivity_m2_per_s * 1.0e4:.3e} cm^2/s; "
        f"MSD fit gives {extracted_diffusivity_m2_per_s * 1.0e4:.3e} cm^2/s "
        f"(R^2 = {msd_fit_r_squared:.5f})."
    )
    return microscopic_figure, microscopic_summary


@app.cell
def _(microscopic_figure, microscopic_summary, mo):
    mo.vstack(
        [
            microscopic_figure,
            mo.md(f"**Current result.** {microscopic_summary}"),
            mo.md(r"""
            **Figure takeaway.** A single path is noisy and has no preferred direction. Diffusivity is
            an ensemble property: many random paths produce a linear mean-square
            displacement. Raising \(T\) or lowering \(H_{\rm mig}\) increases
            \(\Gamma\), so the same number of lattice steps occurs in less time.
            """),
        ]
    )
    return

@app.cell
def _(mo):
    master_time = mo.ui.slider(
        start=0.0,
        stop=20.0,
        step=0.5,
        value=4.0,
        label="Reduced time, Γt",
        show_value=True,
    )
    step_contrast = mo.ui.slider(
        start=0.10,
        stop=0.90,
        step=0.05,
        value=0.70,
        label="Initial concentration-step contrast",
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
            ## 2. From hopping to Fick's law

            Place more defects on the left than on the right. A concentration
            gradient does **not** exert a mechanical force on each defect.
            Individual hops remain random. There are simply more possible
            left-to-right departures than right-to-left departures, so their
            difference is a net flux.

            For site occupancy \(c_j\), the zero-field master equation is

            \[
            \frac{dc_j}{dt}
            =\frac{\Gamma}{2}(c_{j-1}-2c_j+c_{j+1}).
            \]

            With \(D=a^2\Gamma/2\), its long-wavelength limit is
            \(\partial c/\partial t=D\,\partial^2c/\partial x^2\), and the
            bond flux is exactly the discrete form of \(J=-D\,dc/dx\).
            """),
            master_controls,
        ]
    )
    return


@app.cell
def _(
    defect_diffusivity_m2_per_s,
    discrete_bond_fluxes,
    evolve_periodic_master_equation,
    fick_bond_fluxes,
    jump_distance_m,
    master_time,
    np,
    step_contrast,
    hop_frequency_hz,
):
    master_site_count = 256
    contrast_value = float(step_contrast.value)
    left_occupancy = 0.5 * (1.0 + contrast_value)
    right_occupancy = 0.5 * (1.0 - contrast_value)
    initial_master_occupancy = np.where(
        np.arange(master_site_count) < master_site_count // 2,
        left_occupancy,
        right_occupancy,
    )
    master_elapsed_s = float(master_time.value) / hop_frequency_hz
    evolved_master_occupancy = evolve_periodic_master_equation(
        initial_master_occupancy,
        hop_frequency_hz,
        master_elapsed_s,
    )
    microscopic_bond_flux_per_s = discrete_bond_fluxes(
        evolved_master_occupancy,
        hop_frequency_hz,
    )
    fick_bond_flux_per_s = fick_bond_fluxes(
        evolved_master_occupancy,
        jump_distance_m,
        defect_diffusivity_m2_per_s,
    )
    master_mass_relative_error = abs(
        np.sum(evolved_master_occupancy)
        - np.sum(initial_master_occupancy)
    ) / np.sum(initial_master_occupancy)
    master_flux_relative_error = float(
        np.max(
            np.abs(microscopic_bond_flux_per_s - fick_bond_flux_per_s)
        )
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
    master_flux_relative_error,
    master_mass_relative_error,
    master_site_count,
    microscopic_bond_flux_per_s,
    mo,
    np,
    plt,
    scaled_time_axis,
):
    master_position = np.arange(master_site_count) / master_site_count
    master_display_time, master_time_unit = scaled_time_axis(
        np.array([master_elapsed_s])
    )
    master_figure, (profile_axis, _flux_axis) = plt.subplots(
        1,
        2,
        figsize=(13.2, 4.8),
        dpi=120,
    )
    profile_axis.plot(
        master_position,
        initial_master_occupancy,
        color="#999999",
        lw=2.0,
        ls="--",
        label="initial step",
    )
    profile_axis.plot(
        master_position,
        evolved_master_occupancy,
        color="#4C7C86",
        lw=3.0,
        label="master-equation solution",
    )
    profile_axis.set(
        xlabel="position / periodic cell length",
        ylabel="site occupancy",
        title="Random hopping smooths a concentration step",
        ylim=(0.0, 1.0),
    )
    profile_axis.grid(alpha=0.22)
    profile_axis.legend(frameon=False)

    flux_scale = max(np.max(np.abs(microscopic_bond_flux_per_s)), 1.0)
    _flux_axis.plot(
        master_position,
        microscopic_bond_flux_per_s / flux_scale,
        color="#B8734A",
        lw=2.5,
    )
    _flux_axis.axhline(0.0, color="#666D73", lw=1.0)
    _flux_axis.set(
        xlabel="bond position / periodic cell length",
        ylabel="net bond flux / max |flux|",
        title="Opposing random exchanges leave a net flux",
    )
    _flux_axis.grid(alpha=0.22)
    master_figure.tight_layout()
    plt.close(master_figure)

    master_summary = mo.md(
        f"""
        **Figure takeaway.** The selected time is **{master_display_time[0]:.3g} {master_time_unit}**.
        Total defect number changes by only
        **{master_mass_relative_error:.2e}** (relative), while the microscopic
        bond flux and discrete Fick flux differ by
        **{master_flux_relative_error:.2e}** (relative).
        Because this teaching lattice is periodic, the concentration step also
        wraps across the cell edge. That creates a second, oppositely directed
        flux peak at the boundary.
        """
    )
    mo.vstack([master_figure, master_summary])
    return (master_figure,)


@app.cell
def _(mo):
    charge_selector = mo.ui.dropdown(
        options={
            "positive defect, z = +1": 1,
            "negative defect, z = -1": -1,
        },
        value="positive defect, z = +1",
        label="Defect charge",
    )
    field_sign = mo.ui.dropdown(
        options={
            "E points toward +x": 1.0,
            "E points toward -x": -1.0,
        },
        value="E points toward +x",
        label="Field direction",
    )
    log_electric_field = mo.ui.slider(
        start=-2.0,
        stop=6.0,
        step=0.25,
        value=1.0,
        label="log10 |E| when balance is off (V/cm)",
        show_value=True,
    )
    relative_concentration_gradient = mo.ui.slider(
        start=-0.80,
        stop=0.80,
        step=0.05,
        value=0.30,
        label="(1/c) dc/dx (per µm)",
        show_value=True,
    )
    enforce_electrochemical_balance = mo.ui.checkbox(
        value=True,
        label="Set electrical gradient to cancel chemical gradient",
    )
    field_controls = mo.hstack(
        [
            charge_selector,
            field_sign,
            log_electric_field,
            relative_concentration_gradient,
            enforce_electrochemical_balance,
        ],
        justify="start",
        align="center",
        wrap=True,
        gap=1.5,
    )
    return (
        charge_selector,
        enforce_electrochemical_balance,
        field_controls,
        field_sign,
        log_electric_field,
        relative_concentration_gradient,
    )


@app.cell
def _(field_controls, mo):
    mo.vstack(
        [
            mo.md(r"""
            ## 3. Add an electric field: electrochemical potential

            Define positive \(x\) to the right and \(E=-d\phi/dx\). A species
            with charge \(q=ze\) has symmetric-barrier rates

            \[
            \Gamma_+=\frac{\Gamma}{2}
            \exp\!\left(\frac{zeEa}{2k_BT}\right),\qquad
            \Gamma_-=\frac{\Gamma}{2}
            \exp\!\left(-\frac{zeEa}{2k_BT}\right).
            \]

            The factors \(1/2\) mean that, without a field, half of all hops go
            in each direction. The exact drift is
            \(v=a(\Gamma_+-\Gamma_-)\). At low field,
            \(|zeEa|\ll k_BT\), it becomes the Nernst–Einstein result

            \[
            v=\frac{zeD}{k_BT}E.
            \]

            The same field is used in the hopping landscape and in the
            macroscopic flux example below. With the cancellation toggle on,
            that field is calculated from the chosen concentration gradient so
            that \(d\widetilde{\mu}/dx=0\). With it off, the field controls are
            used directly.
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
    enforce_electrochemical_balance,
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
    demonstration_concentration_mol_per_m3 = 1000.0
    relative_gradient_per_m = (
        float(relative_concentration_gradient.value) * 1.0e6
    )
    demonstration_gradient_mol_per_m4 = (
        demonstration_concentration_mol_per_m3 * relative_gradient_per_m
    )
    if enforce_electrochemical_balance.value:
        demonstration_potential_gradient_v_per_m = (
            -GAS_CONSTANT_J_PER_MOL_K
            * temperature_k
            * relative_gradient_per_m
            / (charge_number_value * FARADAY_C_PER_MOL)
        )
    else:
        selected_electric_field_v_per_m = float(field_sign.value) * (
            100.0 * 10.0 ** float(log_electric_field.value)
        )
        demonstration_potential_gradient_v_per_m = -selected_electric_field_v_per_m

    electric_field_v_per_m = -demonstration_potential_gradient_v_per_m
    forward_rate_hz, backward_rate_hz, _unbiased_hop_frequency_hz = biased_directional_rates(
        temperature_k,
        migration_enthalpy_ev,
        attempt_frequency_hz,
        jump_distance_m,
        charge_number_value,
        electric_field_v_per_m,
    )
    exact_drift_m_per_s = exact_hopping_drift_velocity(
        jump_distance_m,
        forward_rate_hz,
        backward_rate_hz,
    )
    low_field_drift_m_per_s = nernst_einstein_drift_velocity(
        defect_diffusivity_m2_per_s,
        temperature_k,
        charge_number_value,
        electric_field_v_per_m,
    )
    hopping_bias_ratio = forward_rate_hz / backward_rate_hz
    field_work_per_hop_ev = (
        charge_number_value * electric_field_v_per_m * jump_distance_m
    )
    low_field_parameter = abs(field_work_per_hop_ev) / (
        8.617333262e-5 * temperature_k
    )

    chemical_potential_gradient_j_per_mol_m = (
        GAS_CONSTANT_J_PER_MOL_K
        * temperature_k
        * relative_gradient_per_m
    )
    electrical_potential_gradient_j_per_mol_m = (
        charge_number_value
        * FARADAY_C_PER_MOL
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
        defect_diffusivity_m2_per_s,
        demonstration_concentration_mol_per_m3,
        demonstration_gradient_mol_per_m4,
        demonstration_potential_gradient_v_per_m,
        charge_number_value,
        temperature_k,
    )
    return (
        backward_rate_hz,
        charge_number_value,
        chemical_potential_gradient_j_per_mol_m,
        demonstration_concentration_mol_per_m3,
        demonstration_gradient_mol_per_m4,
        demonstration_potential_gradient_v_per_m,
        diffusion_flux_mol_per_m2_s,
        electric_field_v_per_m,
        electrical_flux_mol_per_m2_s,
        electrical_potential_gradient_j_per_mol_m,
        electrochemical_gradient_j_per_mol_m,
        exact_drift_m_per_s,
        field_work_per_hop_ev,
        forward_rate_hz,
        hopping_bias_ratio,
        low_field_drift_m_per_s,
        low_field_parameter,
        relative_gradient_per_m,
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
    total_np_flux_mol_per_m2_s,
):
    field_coordinate = np.linspace(0.0, 2.0, 500)
    field_position_m = field_coordinate * jump_distance_m
    tilted_energy_ev = (
        0.5
        * migration_enthalpy_ev
        * (1.0 - np.cos(2.0 * np.pi * field_coordinate))
        - charge_number_value * electric_field_v_per_m * field_position_m
    )
    field_figure, (tilt_axis, gradient_axis, _flux_axis) = plt.subplots(
        1,
        3,
        figsize=(15.2, 4.9),
        dpi=120,
    )
    tilt_axis.plot(field_coordinate, tilted_energy_ev, color="#4C7C86", lw=3.0)
    tilt_axis.scatter(
        [0.0, 1.0, 2.0],
        [
            0.0,
            -charge_number_value * electric_field_v_per_m * jump_distance_m,
            -2.0 * charge_number_value * electric_field_v_per_m * jump_distance_m,
        ],
        s=65,
        color="#C49345",
        edgecolor="#40464D",
        zorder=4,
    )
    tilt_axis.set(
        xlabel="position / a",
        ylabel="energy (eV)",
        title="The field tilts the hopping landscape",
    )
    tilt_axis.grid(alpha=0.22)

    gradient_values = np.array(
        [
            chemical_potential_gradient_j_per_mol_m,
            electrical_potential_gradient_j_per_mol_m,
            electrochemical_gradient_j_per_mol_m,
        ]
    ) * 1.0e-6
    gradient_axis.bar(
        ["chemical", "electrical", "total"],
        gradient_values,
        color=["#4C7C86", "#B8734A", "#C49345"],
    )
    gradient_axis.axhline(0.0, color="#333333", lw=1.0)
    gradient_axis.set(
        ylabel="potential change (J/mol per µm)",
        title=r"$d\widetilde{\mu}/dx=d\mu/dx+zF\,d\phi/dx$",
    )
    gradient_axis.tick_params(axis="x", rotation=18)
    gradient_axis.grid(axis="y", alpha=0.22)

    flux_values = np.array(
        [
            diffusion_flux_mol_per_m2_s,
            electrical_flux_mol_per_m2_s,
            total_np_flux_mol_per_m2_s,
        ]
    )
    _flux_axis.bar(
        ["diffusion", "electrical", "total"],
        flux_values,
        color=["#4C7C86", "#B8734A", "#C49345"],
    )
    _flux_axis.axhline(0.0, color="#333333", lw=1.0)
    _flux_axis.set(
        ylabel=r"molar flux (mol m$^{-2}$ s$^{-1}$)",
        title="Nonzero parts can cancel exactly",
    )
    _flux_axis.tick_params(axis="x", rotation=18)
    _flux_axis.grid(axis="y", alpha=0.22)
    field_figure.tight_layout()
    plt.close(field_figure)

    balance_scale = max(
        abs(diffusion_flux_mol_per_m2_s)
        + abs(electrical_flux_mol_per_m2_s),
        1.0e-300,
    )
    field_summary = mo.md(
        rf"""
        **Figure takeaway.** The microscopic and macroscopic descriptions use the same electric driving force. \(ze\mathcal{{E}}a={field_work_per_hop_ev:.3e}\) eV,
        \(|zeEa|/(k_BT)={low_field_parameter:.3e}\), and
        \(\Gamma_+/\Gamma_-={hopping_bias_ratio:.5g}\). The exact drift is
        **{exact_drift_m_per_s:.3e} m/s**; the low-field prediction is
        **{low_field_drift_m_per_s:.3e} m/s**.
        At the default equilibrium field, the tilt is intentionally tiny beside
        the **{migration_enthalpy_ev:.2f} eV** migration barrier; that contrast
        is the physical point.

        **One field, two views.** The field is
        **{electric_field_v_per_m / 100.0:.3e} V/cm**, so
        \(d\phi/dx={demonstration_potential_gradient_v_per_m / 100.0:.3e}\) V/cm.
        At equilibrium the chemical and electrical contributions oppose one
        another, so the total flux vanishes even though each contribution is
        nonzero.
        """
    )
    mo.vstack([field_figure, field_summary])
    return (field_figure,)


@app.cell
def _(mo):
    mo.md(r"""
    For number concentration \(c_N\), the one-dimensional low-field equation is

    \[
    J_N=-D\frac{dc_N}{dx}+\frac{zeD}{k_BT}c_N\mathcal E .
    \]

    For molar concentration \(c\), use \(F=N_Ae\) and \(R=N_Ak_B\):

    \[
    J=-D\frac{dc}{dx}-\frac{zFD}{RT}c\frac{d\phi}{dx}
      =-\frac{Dc}{RT}\frac{d(\mu+zF\phi)}{dx}.
    \]

    Thus \(\widetilde{\mu}=\mu+zF\phi\) is the electrochemical potential.
    Equilibrium requires \(d\widetilde{\mu}/dx=0\), not separately
    \(d\mu/dx=0\) and \(d\phi/dx=0\). Keep the cancellation toggle on to see
    a nonzero chemical gradient and nonzero electrical gradient produce zero
    total flux.
    """)
    return

@app.cell
def _(mo):
    mo.md(r"""
    ## 4. Three diffusivities used in the lectures

    The symbol \(D\) does not always describe the same experiment.

    | notation | experiment | what moves? | composition changes? |
    |---|---|---|---|
    | \(D^*\) | isotope tracer profile | labeled atoms | no |
    | \(D^q\) | steady conductivity | charge-carrying ions | no |
    | \(D_{\rm Li}^{\delta}\) | chemical relaxation | \(\mathrm{Li^+}\) and \(e^-\) together | yes |

    The first sections of this notebook calculate the microscopic diffusivity
    of an ideal mobile defect from its hops. A tracer or conductivity
    measurement can contain additional information about which atoms carry
    those hops and how their motions are correlated. Chemical diffusion is
    different because the material's composition changes.

    Tracer and conductivity-derived diffusion are related by the **Haven
    ratio**,

    \[
    D^*=H D^q.
    \]

    Here $H$ summarizes correlations between successive ionic motions; $H=1$
    is the uncorrelated limit, not a rule for all solids.
    """)
    return


@app.cell
def _(mo):
    log_ionic_diffusivity_cm2 = mo.ui.slider(
        start=-14.0,
        stop=-6.0,
        step=0.25,
        value=-10.0,
        label="log10 D_Li⁺ (cm²/s)",
        show_value=True,
    )
    log_electronic_to_ionic_ratio = mo.ui.slider(
        start=-4.0,
        stop=4.0,
        step=0.25,
        value=3.0,
        label="log10(D_e⁻ / D_Li⁺)",
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
            ## 5. Chemical diffusion: Li⁺ and electrons move together

            Follow the Li example used in Lecture 5:

            \[
            \mathrm{Li}\rightleftharpoons \mathrm{Li^+}+e^- .
            \]

            In the bulk, local charge neutrality gives
            \(c_{\mathrm{Li^+}}=c_{e^-}=c_{\mathrm{Li}}\). Local equilibrium gives

            \[
            \mu_{\mathrm{Li}}
            =\widetilde{\mu}_{\mathrm{Li^+}}
            +\widetilde{\mu}_{e^-}.
            \]

            During one-dimensional chemical diffusion, the steady fluxes are
            equal:

            \[
            J_{\mathrm{Li}}
            =J_{\mathrm{Li^+}}
            =J_{e^-}.
            \]

            If electrons would diffuse faster on their own, a small internal
            electric field slows them and speeds up \(\mathrm{Li^+}\). The bulk
            remains charge neutral; the separated-charge sketch in the lecture
            is only an exaggerated intermediate picture.
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
        lw=3.0,
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
        label="current state",
    )
    ratio_axis.axhline(1.0, color="#858B90", lw=1.0, ls=":")
    ratio_axis.set(
        xlabel=r"mobility contrast, $r=D_{e^-}/D_{\rm Li^+}$",
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
        **Figure takeaway.** The internal field forces the ion and electron
        to share one neutral-composition flux. The internal potential gradient is
        **{internal_potential_gradient_v_per_m / 100.0:.3e} V/cm**. It gives
        \(J_{{\rm Li^+}}={ionic_flux_mol_per_m2_s:.3e}\) and
        \(J_{{e^-}}={electronic_flux_mol_per_m2_s:.3e}\) mol/(m² s), while
        \(F(J_{{\rm Li^+}}-J_{{e^-}})
        ={open_circuit_current_a_per_m2:.3e}\) A/m².

        The common flux corresponds to
        \(D_{{\rm Li}}^{{\delta}}
        ={lithium_chemical_diffusivity_cm2_per_s:.3e}\) cm²/s.
        """
    )
    mo.vstack([ambipolar_figure, ambipolar_summary])
    return (ambipolar_figure,)


@app.cell
def _(mo):
    _general_chemical_diffusion = mo.md(r"""
    For the monovalent pair used here, equal molar fluxes and local equilibrium
    give

    \[
    J=-\frac{\sigma_i\sigma_e}{F^2(\sigma_i+\sigma_e)}
    \frac{\partial\mu_{\rm Li}}{\partial x}
    =-D_{\rm Li}^{\delta}\frac{\partial c}{\partial x},
    \]

    so the general composition-dependent result is

    \[
    \boxed{D_{\rm Li}^{\delta}=
    \frac{\sigma_i\sigma_e}{F^2(\sigma_i+\sigma_e)}
    \frac{\partial\mu_{\rm Li}}{\partial c}}.
    \]

    The transport prefactor contains the two conductivity pathways; the
    derivative $\partial\mu_{\rm Li}/\partial c$ is the thermodynamic
    susceptibility for molar $\mu_{\rm Li}$ (J mol$^{-1}$) and molar $c$
    (mol m$^{-3}$). Together they give m$^2$ s$^{-1}$. For the ideal dilute
    pair, $\mu_{\rm Li}=\mu_{\rm Li}^0+2RT\ln(c/c_0)$, so the expression
    reduces to the Lecture 5 formula shown next. The compact harmonic-mean
    result is therefore a special ideal limit, not a universal MIEC identity.
    """)
    mo.accordion({"General chemical diffusion (non-ideal bridge)": _general_chemical_diffusion})
    return


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
    mo.md(
        rf"""
        ### The Lecture 5 result

        The flux law used in the slides is

        \[
        J_i=-\frac{{\sigma_i}}{{z_i^2F^2}}
        \frac{{d\widetilde{{\mu}}_i}}{{dx}}.
        \]

        Combining local equilibrium with
        \(J_{{\rm Li^+}}=J_{{e^-}}\) gives, in the dilute limit,

        \[
        D_{{\rm Li}}^\delta
        =\frac{{RT}}{{F^2}}
        \frac{{\sigma_{{e^-}}\sigma_{{\rm Li^+}}}}
             {{\sigma_{{e^-}}+\sigma_{{\rm Li^+}}}}
        \left(\frac{{1}}{{c_{{e^-}}}}
        +\frac{{1}}{{c_{{\rm Li^+}}}}\right).
        \]

        For \(c_{{e^-}}=c_{{\rm Li^+}}\), the Nernst-Einstein relation reduces
        this to

        \[
        D_{{\rm Li}}^\delta
        =\frac{{2D_{{\rm Li^+}}D_{{e^-}}}}
        {{D_{{\rm Li^+}}+D_{{e^-}}}}.
        \]

        Here \(D_{{\rm Li^+}}={ionic_diffusivity_cm2_per_s:.3e}\) cm²/s,
        \(D_{{e^-}}={electronic_diffusivity_cm2_per_s:.3e}\) cm²/s,
        \(\sigma_{{\rm Li^+}}={ionic_conductivity_s_per_m:.3e}\) S/m, and
        \(\sigma_{{e^-}}={electronic_conductivity_s_per_m:.3e}\) S/m.
        The conductivity expression gives
        **{conductivity_form_diffusivity_m2_per_s * 1.0e4:.3e} cm²/s**, matching
        \(D_{{\rm Li}}^\delta
        ={lithium_chemical_diffusivity_cm2_per_s:.3e}\) cm²/s.

        We stop at the dilute-limit Li expression developed in the slides.
        """
    )
    return


@app.cell
def _(mo):
    log_sample_length = mo.ui.slider(
        start=-8.0,
        stop=-2.0,
        step=0.25,
        value=-5.0,
        label="log10 sample length, L (m)",
        show_value=True,
    )
    reduced_profile_time = mo.ui.slider(
        start=0.01,
        stop=0.50,
        step=0.01,
        value=0.08,
        label=r"Fourier number, theta = D_Li^delta t / L²",
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
def _(mo, relaxation_controls):
    mo.vstack(
        [
            mo.md(r"""
            ## 6. Chemical diffusivity sets the relaxation time

            Two related clocks must not be given the same name:

            \[
            t_D=\frac{L^2}{D_{\rm Li}^{\delta}},\qquad
            \tau^\delta=\frac{L^2}{\pi^2D_{\rm Li}^{\delta}},\qquad
            \theta=\frac{D_{\rm Li}^{\delta}t}{L^2}.
            \]

            $t_D$ is the direct diffusion scaling time. $\tau^\delta$ is the
            slowest-mode time for this slab relaxation. The profile control is
            the Fourier number $\theta$, not an unspecified “reduced time.”

            Increasing the length by a factor of ten increases both diffusion
            times by a factor of one hundred. This is why the same material can
            respond quickly as a thin film and slowly as a bulk sample.
            """),
            relaxation_controls,
        ]
    )
    return


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
        lw=3.0,
        label=r"scaling time $t_D=L^2/D_{\rm Li}^{\delta}$",
    )
    time_axis.loglog(
        length_curve_m,
        relaxation_curve_s / np.pi**2,
        color="#B8734A",
        lw=2.5,
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
        xlabel="sample length, L (m)",
        ylabel="characteristic time (s)",
        title="Diffusion time grows as length squared",
    )
    time_axis.grid(which="both", alpha=0.22)
    time_axis.legend(frameon=False)

    relaxation_profile_axis.plot(
        slab_position,
        selected_slab_profile,
        color="#B65C4A",
        lw=3.0,
    )
    relaxation_profile_axis.set(
        xlabel="position, x / L",
        ylabel="remaining normalized composition change",
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
        **Figure takeaway.** With \(D_{{\rm Li}}^\delta
        ={lithium_chemical_diffusivity_cm2_per_s:.3e}\) cm²/s and
        \(L={selected_length_m:.3e}\) m, $t_D$ is **{time_text}** and
        $\tau^\delta$ is **{first_mode_text}**.
        The right panel shows the remaining composition change at Fourier
        number \(\theta={float(reduced_profile_time.value):.2f}\).
        """
    )
    mo.vstack([relaxation_figure, relaxation_summary])
    return (relaxation_figure,)


@app.cell
def _(
    FARADAY_C_PER_MOL,
    KB_EV_PER_K,
    ambipolar_gradient_mol_per_m4,
    analytic_common_flux_mol_per_m2_s,
    backward_rate_hz,
    conductivity_form_diffusivity_m2_per_s,
    defect_diffusivity_m2_per_s,
    diffusion_flux_mol_per_m2_s,
    electronic_diffusivity_m2_per_s,
    electronic_flux_mol_per_m2_s,
    electrical_flux_mol_per_m2_s,
    exact_drift_m_per_s,
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
    msd_fit_r_squared,
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
    low_field_drift_error = _relative_error(
        exact_drift_m_per_s,
        low_field_drift_m_per_s,
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
            msd_diffusivity_error < 0.08
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

        These checks follow the notebook's main chain from atomic hopping to
        coupled chemical diffusion; they do not add new assumptions.
        """
    )
    mo.accordion({"Physical consistency checks": _checks})
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## Take-home map

    1. **Hops become diffusion.** In one dimension,
       $\Gamma_{\rm hop}=\nu e^{-\Delta H_{\rm mig}/k_BT}$ gives
       $D=a^2\Gamma_{\rm hop}/2$, $\langle x^2\rangle=2Dt$, and the
       long-wavelength flux $J=-D\,dc/dx$.
    2. **Electrochemical potential unifies the driving forces.** A concentration
       gradient and $\mathcal E=-\partial\phi/\partial x$ enter through
       $\widetilde\mu=\mu+zF\phi$; equilibrium means their contributions cancel.
    3. **Chemical diffusion moves a neutral composition.** For the dilute ideal
       Li pair, equal Li-ion/electron fluxes give
       $D_{\rm Li}^{\delta}=2D_{\rm Li^+}D_{e^-}/(D_{\rm Li^+}+D_{e^-})$.
       This compact factor of two is specific to the stated ideal pair, not a
       universal formula for every defect reaction.

    **Model boundary.** Every spatial equation and simulation in this notebook
    is one-dimensional. We use ideal, dilute concentrations for the Li chemical
    diffusion derivation, exactly where the lecture obtains the compact
    expression above. The Haven ratio introduces correlation without
    calculating it here; interfaces are left for later treatment.
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
