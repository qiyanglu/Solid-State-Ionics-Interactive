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
    # Defect Transport: From Atomic Hopping to Chemical Diffusion

    **Guiding question.** How does a thermally activated atomic jump become a
    measurable transport coefficient—and why can a composition profile relax
    only when ions and electrons move together?

    This module follows one connected argument:

    \[
    \text{activated hops}
    \rightarrow \langle x^2\rangle
    \rightarrow D
    \rightarrow J
    \rightarrow \widetilde{\mu}
    \rightarrow D_{\rm chem}
    \rightarrow \tau .
    \]

    We begin with an ideal one-dimensional lattice, then add concentration and
    electrical driving forces, and finally couple \(\mathrm{Li^+}\) and
    \(e^-\) under local electroneutrality. The simple models expose the physics;
    correlation factors, space charge, electrode polarization, and specific
    PITT/GITT protocols are deliberately left for later modules.

    Canonical units are **K, eV per particle, m, s, V/m, mol/m³, mol/(m² s),
    J/mol, S/m, and m²/s**. Particle and molar equations are kept separate.
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

    def hopping_rate_per_neighbor(temperature_k, migration_enthalpy_ev, attempt_frequency_hz):
        """Return the zero-field rate to each neighbor, w0, in s^-1."""
        temperature = require_positive("temperature_k", temperature_k)
        barrier = require_positive("migration_enthalpy_ev", migration_enthalpy_ev)
        attempt = require_positive("attempt_frequency_hz", attempt_frequency_hz)
        return attempt * np.exp(-barrier / (KB_EV_PER_K * temperature))

    def hopping_diffusivity_1d(jump_distance_m, rate_per_neighbor_hz):
        """Return D=a^2 w0 in m^2/s for two 1D neighbors, each at rate w0."""
        distance = require_positive("jump_distance_m", jump_distance_m)
        rate = require_positive("rate_per_neighbor_hz", rate_per_neighbor_hz)
        return distance**2 * rate

    def biased_neighbor_rates(
        temperature_k,
        migration_enthalpy_ev,
        attempt_frequency_hz,
        jump_distance_m,
        charge_number,
        electric_field_v_per_m,
    ):
        """Return w+, w-, and w0 in s^-1 for E pointing along positive x.

        The electrostatic energy change over a +x hop is -z e E a. Because
        energies are supplied in eV, z E a is already the corresponding eV
        value. A symmetric transition state splits the bias equally.
        """
        temperature = require_positive("temperature_k", temperature_k)
        distance = require_positive("jump_distance_m", jump_distance_m)
        charge = float(charge_number)
        field = float(electric_field_v_per_m)
        if not np.isfinite(charge) or charge == 0.0:
            raise ValueError("charge_number must be finite and nonzero")
        if not np.isfinite(field):
            raise ValueError("electric_field_v_per_m must be finite")
        base_rate = hopping_rate_per_neighbor(
            temperature,
            migration_enthalpy_ev,
            attempt_frequency_hz,
        )
        half_bias = charge * field * distance / (
            2.0 * KB_EV_PER_K * temperature
        )
        return (
            base_rate * np.exp(half_bias),
            base_rate * np.exp(-half_bias),
            base_rate,
        )

    def exact_hopping_drift_velocity(jump_distance_m, forward_rate_hz, backward_rate_hz):
        """Return a(w+ - w-) in m/s."""
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
        rate_per_neighbor_hz,
        walker_count=2000,
        step_count=320,
        seed=2026,
    ):
        """Simulate unbiased 1D walks at observation times spaced by 1/(2w0).

        Each observation interval contains one jump. The total leaving rate is
        2w0 because w0 is the rate to each of the two neighbors.
        """
        distance = require_positive("jump_distance_m", jump_distance_m)
        rate = require_positive("rate_per_neighbor_hz", rate_per_neighbor_hz)
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
        times_s = np.arange(steps + 1, dtype=float) / (2.0 * rate)
        msd_m2 = np.mean(positions_m**2, axis=0)
        return times_s, positions_m, msd_m2

    def evolve_periodic_master_equation(initial_occupancy, rate_per_neighbor_hz, time_s):
        """Evolve dc_j/dt=w0(c_{j-1}-2c_j+c_{j+1}) on a periodic lattice."""
        profile = np.asarray(initial_occupancy, dtype=float)
        if profile.ndim != 1 or profile.size < 4:
            raise ValueError("initial_occupancy must be a 1D array with at least 4 sites")
        if np.any(~np.isfinite(profile)) or np.any(profile < 0.0):
            raise ValueError("initial_occupancy must be finite and nonnegative")
        rate = require_positive("rate_per_neighbor_hz", rate_per_neighbor_hz)
        elapsed = float(time_s)
        if not np.isfinite(elapsed) or elapsed < 0.0:
            raise ValueError("time_s must be finite and nonnegative")
        modes = np.arange(profile.size)
        decay = np.exp(
            -4.0 * rate * elapsed * np.sin(np.pi * modes / profile.size) ** 2
        )
        return np.real(np.fft.ifft(np.fft.fft(profile) * decay))

    def discrete_bond_fluxes(occupancy, rate_per_neighbor_hz):
        """Return net particle crossings per second on bonds j -> j+1."""
        profile = np.asarray(occupancy, dtype=float)
        return float(rate_per_neighbor_hz) * (profile - np.roll(profile, -1))

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
        """Return diffusion, electrical, and total molar flux in mol/(m^2 s)."""
        diffusivity = require_positive("diffusivity_m2_per_s", diffusivity_m2_per_s)
        concentration = require_positive(
            "concentration_mol_per_m3",
            concentration_mol_per_m3,
        )
        temperature = require_positive("temperature_k", temperature_k)
        chemical_flux = -diffusivity * float(concentration_gradient_mol_per_m4)
        electrical_flux = (
            -float(charge_number)
            * FARADAY_C_PER_MOL
            * diffusivity
            * concentration
            * float(potential_gradient_v_per_m)
            / (GAS_CONSTANT_J_PER_MOL_K * temperature)
        )
        return chemical_flux, electrical_flux, chemical_flux + electrical_flux

    def ambipolar_internal_potential_gradient(
        ionic_diffusivity_m2_per_s,
        electronic_diffusivity_m2_per_s,
        concentration_mol_per_m3,
        concentration_gradient_mol_per_m4,
        temperature_k,
    ):
        """Return dphi/dx in V/m enforcing J_i=J_e for Li <-> Li+ + e-."""
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
        """Return Li+ and electron molar fluxes in mol/(m^2 s)."""
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

    def ideal_ambipolar_diffusivity(
        ionic_diffusivity_m2_per_s,
        electronic_diffusivity_m2_per_s,
    ):
        """Return 2 Di De/(Di+De) in m^2/s for Li <-> Li+ + e-."""
        ionic = require_positive("ionic_diffusivity_m2_per_s", ionic_diffusivity_m2_per_s)
        electronic = require_positive(
            "electronic_diffusivity_m2_per_s",
            electronic_diffusivity_m2_per_s,
        )
        return 2.0 * ionic * electronic / (ionic + electronic)

    def chemical_diffusivity(
        ionic_diffusivity_m2_per_s,
        electronic_diffusivity_m2_per_s,
        thermodynamic_factor,
    ):
        """Return D_chem=2 Di De/(Di+De) Theta in m^2/s."""
        factor = require_positive("thermodynamic_factor", thermodynamic_factor)
        return ideal_ambipolar_diffusivity(
            ionic_diffusivity_m2_per_s,
            electronic_diffusivity_m2_per_s,
        ) * factor

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

    def characteristic_relaxation_time(length_m, diffusivity_m2_per_s):
        """Return the scaling time L^2/D in seconds."""
        length = require_positive("length_m", length_m)
        diffusivity = require_positive("diffusivity_m2_per_s", diffusivity_m2_per_s)
        return length**2 / diffusivity

    def slab_remaining_profile(normalized_position, time_over_l2_by_d, term_count=80):
        """Return normalized excess in a slab held at zero excess at both faces."""
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
        biased_neighbor_rates,
        characteristic_relaxation_time,
        chemical_diffusivity,
        conductivity_from_diffusivity,
        discrete_bond_fluxes,
        evolve_periodic_master_equation,
        exact_hopping_drift_velocity,
        fick_bond_fluxes,
        hopping_diffusivity_1d,
        hopping_rate_per_neighbor,
        ideal_ambipolar_diffusivity,
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
        label="Migration enthalpy, H_mig (eV)",
        show_value=True,
    )
    log_attempt_frequency = mo.ui.slider(
        start=10.0,
        stop=14.0,
        step=0.25,
        value=13.0,
        label="log10 attempt frequency, nu0 (s^-1)",
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

            Let \(w_0=\nu_0\exp[-H_{\rm mig}/(k_BT)]\) be the rate **to each
            neighbor**. In one dimension there are two neighbors, so the total
            leaving rate is \(2w_0\). This convention gives

            \[
            \langle x^2\rangle=2Dt,\qquad D=a^2w_0.
            \]

            The lecture slides use \(\Gamma\) for the **total** 1D hop
            frequency. The two notations are therefore related by

            \[
            \Gamma=2w_0,
            \qquad
            D=\frac{1}{2}a^2\Gamma=a^2w_0.
            \]

            This is one physical result written with two rate conventions, not
            two different diffusivities.
            """),
            hopping_controls,
        ]
    )
    return


@app.cell
def _(
    hopping_diffusivity_1d,
    hopping_rate_per_neighbor,
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
    zero_field_rate_hz = hopping_rate_per_neighbor(
        temperature_k,
        migration_enthalpy_ev,
        attempt_frequency_hz,
    )
    defect_diffusivity_m2_per_s = hopping_diffusivity_1d(
        jump_distance_m,
        zero_field_rate_hz,
    )
    walk_times_s, walk_positions_m, walk_msd_m2 = simulate_zero_field_walks(
        jump_distance_m,
        zero_field_rate_hz,
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
        zero_field_rate_hz,
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
    zero_field_rate_hz,
):
    plt.rcParams.update(
        {
            "font.size": 12,
            "axes.titlesize": 15,
            "axes.labelsize": 13,
            "legend.fontsize": 10,
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
        color="#007C91",
        lw=3.0,
    )
    barrier_axis.scatter(
        [0.0, 1.0, 2.0],
        [0.0, 0.0, 0.0],
        s=80,
        color="#EE9B00",
        edgecolor="#222222",
        zorder=4,
        label="equivalent sites",
    )
    barrier_axis.annotate(
        "",
        xy=(0.5, migration_enthalpy_ev),
        xytext=(0.5, 0.0),
        arrowprops={"arrowstyle": "<->", "color": "#CC3311", "lw": 1.8},
    )
    barrier_axis.text(
        0.55,
        0.52 * migration_enthalpy_ev,
        r"$H_{\rm mig}$",
        color="#CC3311",
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
        color="#007C91",
        lw=2.8,
        label="deterministic-seed simulation",
    )
    msd_axis.plot(
        display_time,
        2.0 * defect_diffusivity_m2_per_s * walk_times_s * 1.0e18,
        color="#D55E00",
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
        f"w0 = {zero_field_rate_hz:.3e} s^-1 per neighbor; "
        f"D = a^2 w0 = {defect_diffusivity_m2_per_s:.3e} m^2/s; "
        f"MSD fit gives {extracted_diffusivity_m2_per_s:.3e} m^2/s "
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
            A single path is noisy and has no preferred direction. Diffusivity is
            an ensemble property: many random paths produce a linear mean-square
            displacement. Raising \(T\) or lowering \(H_{\rm mig}\) increases
            \(w_0\), so the same number of lattice steps occurs in less time.
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
        label="Reduced time, w0 t",
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
            =w_0(c_{j-1}-2c_j+c_{j+1}).
            \]

            With \(D=a^2w_0\), its long-wavelength limit is
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
    zero_field_rate_hz,
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
    master_elapsed_s = float(master_time.value) / zero_field_rate_hz
    evolved_master_occupancy = evolve_periodic_master_equation(
        initial_master_occupancy,
        zero_field_rate_hz,
        master_elapsed_s,
    )
    microscopic_bond_flux_per_s = discrete_bond_fluxes(
        evolved_master_occupancy,
        zero_field_rate_hz,
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
        color="#007C91",
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
        color="#D55E00",
        lw=2.5,
    )
    _flux_axis.axhline(0.0, color="#555555", lw=1.0)
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
        The selected time is **{master_display_time[0]:.3g} {master_time_unit}**.
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
        start=3.0,
        stop=9.0,
        step=0.25,
        value=8.0,
        label="log10 |E| when balance is off (V/m)",
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
            w_+=w_0\exp\!\left(\frac{zeEa}{2k_BT}\right),\qquad
            w_-=w_0\exp\!\left(-\frac{zeEa}{2k_BT}\right).
            \]

            The exact hopping drift is \(v=a(w_+-w_-)\). At low field,
            \(|zeEa|\ll k_BT\), it becomes the Nernst–Einstein result

            \[
            v=\frac{zeD}{k_BT}E.
            \]

            The same field is used in the hopping landscape and in the
            macroscopic flux example below. With the cancellation toggle on,
            that field is calculated from the chosen concentration gradient so
            that \(\nabla\widetilde{\mu}=0\). With it off, the field controls are
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
    biased_neighbor_rates,
    charge_selector,
    defect_diffusivity_m2_per_s,
    enforce_electrochemical_balance,
    exact_hopping_drift_velocity,
    field_sign,
    jump_distance_m,
    log_electric_field,
    migration_enthalpy_ev,
    master_mass_relative_error,
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
            10.0 ** float(log_electric_field.value)
        )
        demonstration_potential_gradient_v_per_m = -selected_electric_field_v_per_m

    electric_field_v_per_m = -demonstration_potential_gradient_v_per_m
    forward_rate_hz, backward_rate_hz, biased_base_rate_hz = biased_neighbor_rates(
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
        biased_base_rate_hz,
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
    tilt_axis.plot(field_coordinate, tilted_energy_ev, color="#007C91", lw=3.0)
    tilt_axis.scatter(
        [0.0, 1.0, 2.0],
        [
            0.0,
            -charge_number_value * electric_field_v_per_m * jump_distance_m,
            -2.0 * charge_number_value * electric_field_v_per_m * jump_distance_m,
        ],
        s=65,
        color="#EE9B00",
        edgecolor="#222222",
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
        color=["#007C91", "#D55E00", "#EE9B00"],
    )
    gradient_axis.axhline(0.0, color="#333333", lw=1.0)
    gradient_axis.set(
        ylabel="potential change (J/mol per µm)",
        title=r"$\nabla\widetilde{\mu}=\nabla\mu+zF\nabla\phi$",
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
        color=["#007C91", "#D55E00", "#EE9B00"],
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
        **Microscopic field response.** \(zeEa={field_work_per_hop_ev:.3e}\) eV,
        \(|zeEa|/(k_BT)={low_field_parameter:.3e}\), and
        \(w_+/w_-={hopping_bias_ratio:.5g}\). The exact drift is
        **{exact_drift_m_per_s:.3e} m/s**; the low-field prediction is
        **{low_field_drift_m_per_s:.3e} m/s**.
        At the default equilibrium field, the tilt is intentionally tiny beside
        the **{migration_enthalpy_ev:.2f} eV** migration barrier; that contrast
        is the physical point.

        **One field, two views.** The field is
        **{electric_field_v_per_m:.3e} V/m**, so
        \(d\phi/dx={demonstration_potential_gradient_v_per_m:.3e}\) V/m.
        The residual
        total flux is **{total_np_flux_mol_per_m2_s:.3e} mol/(m² s)**, or
        **{abs(total_np_flux_mol_per_m2_s) / balance_scale:.2e}** of the two
        opposing contributions.
        """
    )
    mo.vstack([field_figure, field_summary])
    return (field_figure,)


@app.cell
def _(mo):
    mo.md(r"""
    For number concentration \(c_N\) (particles/m³), the low-field equation is

    \[
    J_N=-D\nabla c_N+\frac{zeD}{k_BT}c_NE .
    \]

    For molar concentration \(c\) (mol/m³), use \(F=N_Ae\) and \(R=N_Ak_B\):

    \[
    J=-D\nabla c-\frac{zFD}{RT}c\nabla\phi
      =-\frac{Dc}{RT}\nabla(\mu+zF\phi).
    \]

    Thus \(\widetilde{\mu}=\mu+zF\phi\) is the electrochemical potential.
    Equilibrium requires \(\nabla\widetilde{\mu}=0\), not separately
    \(\nabla\mu=0\) and \(\nabla\phi=0\). Keep the cancellation toggle on to see
    a nonzero chemical gradient and nonzero electrical gradient produce zero
    total flux.
    """)
    return

@app.cell
def _(mo):
    mo.md(r"""
    ## 4. “Which diffusivity?” is part of the physics

    The symbol \(D\) is not universal.

    | Quantity | What it follows | Why it can differ |
    |---|---|---|
    | defect diffusivity \(D\) | motion of an ideal defect population | set here by the jump rate |
    | tracer diffusivity \(D^*\) | labeled atoms | remembers atom–defect exchange geometry |
    | conductivity or self-diffusivity \(D^q\) | ion motion inferred from \(\sigma\), as named in the lectures | the concentration and correlations used in Nernst–Einstein matter |
    | chemical diffusivity \(D_{\rm chem}\) | relaxation of composition | couples carriers and thermodynamic response |

    In an ideal uncorrelated lattice some of these coincide. In real solids,
    tracer correlation factors and collective (Haven-type) correlations often
    separate them. In particular, a conductivity-derived diffusivity for the
    lattice ions is not automatically the same as the diffusivity of the mobile
    defects that enable those ions to move. The next derivation concerns
    **chemical** diffusion and keeps that distinction visible.
    """)
    return


@app.cell
def _(mo):
    log_ionic_diffusivity = mo.ui.slider(
        start=-16.0,
        stop=-8.0,
        step=0.25,
        value=-12.0,
        label="log10 ionic diffusivity, Di (m²/s)",
        show_value=True,
    )
    log_electronic_to_ionic_ratio = mo.ui.slider(
        start=-4.0,
        stop=4.0,
        step=0.25,
        value=3.0,
        label="log10(De/Di)",
        show_value=True,
    )
    thermodynamic_factor_control = mo.ui.slider(
        start=0.20,
        stop=4.00,
        step=0.10,
        value=1.00,
        label="Thermodynamic factor, Theta",
        show_value=True,
    )
    ambipolar_controls = mo.hstack(
        [
            log_ionic_diffusivity,
            log_electronic_to_ionic_ratio,
            thermodynamic_factor_control,
        ],
        justify="start",
        align="center",
        wrap=True,
        gap=1.5,
    )
    return (
        ambipolar_controls,
        log_electronic_to_ionic_ratio,
        log_ionic_diffusivity,
        thermodynamic_factor_control,
    )


@app.cell
def _(ambipolar_controls, mo):
    mo.vstack(
        [
            mo.md(r"""
            ## 5. Ambipolar diffusion: ions and electrons must cooperate

            Consider the neutral insertion reaction

            \[
            \mathrm{Li}\rightleftharpoons \mathrm{Li^+}+e^-,
            \qquad c_i=c_e=c
            \]

            under local electroneutrality. For \(+x\) gradients,

            \[
            J_i=-D_i\left[\frac{dc}{dx}+\frac{F}{RT}c\frac{d\phi}{dx}\right],
            \qquad
            J_e=-D_e\left[\frac{dc}{dx}-\frac{F}{RT}c\frac{d\phi}{dx}\right].
            \]

            Open circuit means \(j=F(J_i-J_e)=0\), hence \(J_i=J_e\). The
            internal field prevents charge separation: it slows the faster
            species and accelerates the slower one.
            """),
            ambipolar_controls,
        ]
    )
    return


@app.cell
def _(
    FARADAY_C_PER_MOL,
    GAS_CONSTANT_J_PER_MOL_K,
    ambipolar_fluxes,
    ambipolar_internal_potential_gradient,
    chemical_diffusivity,
    conductivity_from_diffusivity,
    ideal_ambipolar_diffusivity,
    log_electronic_to_ionic_ratio,
    log_ionic_diffusivity,
    temperature_k,
    thermodynamic_factor_control,
):
    ionic_diffusivity_m2_per_s = 10.0 ** float(log_ionic_diffusivity.value)
    electronic_to_ionic_ratio = 10.0 ** float(
        log_electronic_to_ionic_ratio.value
    )
    electronic_diffusivity_m2_per_s = (
        ionic_diffusivity_m2_per_s * electronic_to_ionic_ratio
    )
    thermodynamic_factor_value = float(thermodynamic_factor_control.value)
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
    ionic_flux_mol_per_m2_s, electronic_flux_mol_per_m2_s = ambipolar_fluxes(
        ionic_diffusivity_m2_per_s,
        electronic_diffusivity_m2_per_s,
        ambipolar_concentration_mol_per_m3,
        ambipolar_gradient_mol_per_m4,
        internal_potential_gradient_v_per_m,
        temperature_k,
    )
    open_circuit_current_a_per_m2 = FARADAY_C_PER_MOL * (
        ionic_flux_mol_per_m2_s - electronic_flux_mol_per_m2_s
    )
    ideal_chemical_diffusivity_m2_per_s = ideal_ambipolar_diffusivity(
        ionic_diffusivity_m2_per_s,
        electronic_diffusivity_m2_per_s,
    )
    selected_chemical_diffusivity_m2_per_s = chemical_diffusivity(
        ionic_diffusivity_m2_per_s,
        electronic_diffusivity_m2_per_s,
        thermodynamic_factor_value,
    )
    analytic_common_flux_mol_per_m2_s = (
        -ideal_chemical_diffusivity_m2_per_s
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
    neutral_chemical_potential_derivative_j_m3_per_mol2 = (
        2.0
        * GAS_CONSTANT_J_PER_MOL_K
        * temperature_k
        * thermodynamic_factor_value
        / ambipolar_concentration_mol_per_m3
    )
    conductivity_form_chemical_diffusivity_m2_per_s = (
        ionic_conductivity_s_per_m
        * electronic_conductivity_s_per_m
        / (ionic_conductivity_s_per_m + electronic_conductivity_s_per_m)
        / FARADAY_C_PER_MOL**2
        * neutral_chemical_potential_derivative_j_m3_per_mol2
    )
    return (
        ambipolar_concentration_mol_per_m3,
        ambipolar_gradient_mol_per_m4,
        ambipolar_relative_gradient_per_m,
        analytic_common_flux_mol_per_m2_s,
        conductivity_form_chemical_diffusivity_m2_per_s,
        electronic_conductivity_s_per_m,
        electronic_diffusivity_m2_per_s,
        electronic_flux_mol_per_m2_s,
        electronic_to_ionic_ratio,
        ideal_chemical_diffusivity_m2_per_s,
        internal_potential_gradient_v_per_m,
        ionic_conductivity_s_per_m,
        ionic_diffusivity_m2_per_s,
        ionic_flux_mol_per_m2_s,
        neutral_chemical_potential_derivative_j_m3_per_mol2,
        open_circuit_current_a_per_m2,
        selected_chemical_diffusivity_m2_per_s,
        thermodynamic_factor_value,
    )


@app.cell
def _(
    ambipolar_gradient_mol_per_m4,
    analytic_common_flux_mol_per_m2_s,
    electronic_diffusivity_m2_per_s,
    electronic_flux_mol_per_m2_s,
    electronic_to_ionic_ratio,
    ideal_chemical_diffusivity_m2_per_s,
    internal_potential_gradient_v_per_m,
    ionic_diffusivity_m2_per_s,
    ionic_flux_mol_per_m2_s,
    mo,
    np,
    open_circuit_current_a_per_m2,
    plt,
    selected_chemical_diffusivity_m2_per_s,
    thermodynamic_factor_value,
):
    ratio_curve = np.logspace(-4.0, 4.0, 500)
    ideal_ratio_curve = 2.0 * ratio_curve / (1.0 + ratio_curve)
    selected_ratio_curve = thermodynamic_factor_value * ideal_ratio_curve

    uncoupled_ionic_flux = -ionic_diffusivity_m2_per_s * ambipolar_gradient_mol_per_m4
    uncoupled_electronic_flux = (
        -electronic_diffusivity_m2_per_s * ambipolar_gradient_mol_per_m4
    )
    flux_magnitudes = np.maximum(
        np.abs(
            [
                uncoupled_ionic_flux,
                uncoupled_electronic_flux,
                ionic_flux_mol_per_m2_s,
            ]
        ),
        1.0e-300,
    )

    ambipolar_figure, (ratio_axis, coupling_axis) = plt.subplots(
        1,
        2,
        figsize=(13.2, 4.9),
        dpi=120,
    )
    ratio_axis.loglog(
        ratio_curve,
        ideal_ratio_curve,
        color="#007C91",
        lw=3.0,
        label=r"ideal: $2r/(1+r)$",
    )
    if abs(thermodynamic_factor_value - 1.0) > 1.0e-12:
        ratio_axis.loglog(
            ratio_curve,
            selected_ratio_curve,
            color="#D55E00",
            lw=2.2,
            ls="--",
            label=rf"with $\Theta={thermodynamic_factor_value:.2f}$",
        )
    ratio_axis.scatter(
        [electronic_to_ionic_ratio],
        [
            selected_chemical_diffusivity_m2_per_s
            / ionic_diffusivity_m2_per_s
        ],
        s=95,
        color="#EE9B00",
        edgecolor="#222222",
        zorder=5,
        label="current state",
    )
    ratio_axis.axhline(1.0, color="#777777", lw=1.0, ls=":")
    ratio_axis.set(
        xlabel=r"mobility contrast, $D_e/D_i$",
        ylabel=r"$D_{\rm chem}/D_i$",
        title="The slower carrier limits ideal ambipolar motion",
        xlim=(1.0e-4, 1.0e4),
        ylim=(1.0e-4, 10.0),
    )
    ratio_axis.grid(which="both", alpha=0.22)
    ratio_axis.legend(frameon=False)

    coupling_axis.bar(
        ["ion alone", "electron alone", "common coupled flux\nJi = Je (Θ = 1)"],
        flux_magnitudes,
        color=["#007C91", "#CC3311", "#EE9B00"],
    )
    coupling_axis.set_yscale("log")
    coupling_axis.set(
        ylabel=r"|flux| (mol m$^{-2}$ s$^{-1}$)",
        title=r"The internal field makes $J_i=J_e$",
    )
    coupling_axis.tick_params(axis="x", rotation=15)
    coupling_axis.grid(axis="y", which="both", alpha=0.22)
    ambipolar_figure.tight_layout()
    plt.close(ambipolar_figure)

    ambipolar_summary = mo.md(
        rf"""
        The internal potential gradient is
        **{internal_potential_gradient_v_per_m:.3e} V/m**. It gives
        \(J_i={ionic_flux_mol_per_m2_s:.3e}\) and
        \(J_e={electronic_flux_mol_per_m2_s:.3e}\) mol/(m² s), while
        \(j=F(J_i-J_e)={open_circuit_current_a_per_m2:.3e}\) A/m².
        The analytic ideal common flux,
        \(-D_{{\rm chem,ideal}}dc/dx\), is
        **{analytic_common_flux_mol_per_m2_s:.3e} mol/(m² s)**.

        Here \(D_{{\rm chem,ideal}}={ideal_chemical_diffusivity_m2_per_s:.3e}\)
        m²/s and the selected thermodynamic factor gives
        \(D_{{\rm chem}}={selected_chemical_diffusivity_m2_per_s:.3e}\) m²/s.
        The flux bars show the ideal \(\Theta=1\) coupling step; the next
        section then adds the thermodynamic factor explicitly.
        """
    )
    mo.vstack([ambipolar_figure, ambipolar_summary])
    return (ambipolar_figure,)


@app.cell
def _(mo):
    mo.md(r"""
    Solving \(J_i=J_e\) gives

    \[
    \frac{d\phi}{dx}
    =\frac{RT}{Fc}\frac{D_e-D_i}{D_i+D_e}\frac{dc}{dx},
    \qquad
    J_i=J_e=-\frac{2D_iD_e}{D_i+D_e}\frac{dc}{dx}.
    \]

    Therefore

    \[
    \boxed{D_{\rm chem,ideal}=\frac{2D_iD_e}{D_i+D_e}}.
    \]

    The factor of two is **not universal**. It is specific to this 1:1
    \(\mathrm{Li}\rightleftharpoons\mathrm{Li^+}+e^-\) model and to defining
    the neutral chemical potential as the sum of the ionic and electronic
    electrochemical potentials.
    """)
    return


@app.cell
def _(
    conductivity_form_chemical_diffusivity_m2_per_s,
    electronic_conductivity_s_per_m,
    ionic_conductivity_s_per_m,
    mo,
    neutral_chemical_potential_derivative_j_m3_per_mol2,
    selected_chemical_diffusivity_m2_per_s,
    thermodynamic_factor_value,
):
    mo.md(
        rf"""
        ## 6. Thermodynamic factor: transport also feels free-energy curvature

        The general conductivity form for this monovalent pair is

        \[
        D_{{\rm chem}}
        =\frac{{\sigma_i\sigma_e}}{{\sigma_i+\sigma_e}}
         \frac{{1}}{{F^2}}
         \frac{{\partial\mu_{{\rm neutral}}}}{{\partial c}}.
        \]

        Unit check:
        \(\sigma\) is S/m, \(F\) is C/mol, and
        \(\partial\mu/\partial c\) is J m³/mol²; the product is m²/s.
        With ideal Nernst–Einstein conductivities

        \[
        \sigma_i=\frac{{F^2cD_i}}{{RT}},\qquad
        \sigma_e=\frac{{F^2cD_e}}{{RT}},
        \]

        define

        \[
        \Theta=\frac{{c}}{{2RT}}
        \frac{{\partial\mu_{{\rm neutral}}}}{{\partial c}}.
        \]

        Then \(\Theta=1\) for this ideal neutral solution and

        \[
        D_{{\rm chem}}=
        \frac{{2D_iD_e}}{{D_i+D_e}}\,\Theta.
        \]

        At the selected state,
        \(\sigma_i={ionic_conductivity_s_per_m:.3e}\) S/m,
        \(\sigma_e={electronic_conductivity_s_per_m:.3e}\) S/m, and
        \(\partial\mu_{{\rm neutral}}/\partial c=
        {neutral_chemical_potential_derivative_j_m3_per_mol2:.3e}\)
        J m³/mol². The conductivity expression gives
        **{conductivity_form_chemical_diffusivity_m2_per_s:.3e} m²/s**, matching
        the displayed \(D_{{\rm chem}}=
        {selected_chemical_diffusivity_m2_per_s:.3e}\) m²/s at
        \(\Theta={thermodynamic_factor_value:.2f}\).

        This is the bridge back to Module 01: the slope or curvature of free
        energy determines how strongly a composition gradient changes chemical
        potential. Kinetics supplies mobility; thermodynamics supplies the
        restoring force.
        """
    )
    return


@app.cell
def _(mo):
    log_sample_length = mo.ui.slider(
        start=-9.0,
        stop=-3.0,
        step=0.25,
        value=-5.0,
        label="log10 sample thickness, L (m)",
        show_value=True,
    )
    log_time_over_tau = mo.ui.slider(
        start=-3.0,
        stop=1.0,
        step=0.25,
        value=-1.5,
        label="log10(t / [L²/Dchem])",
        show_value=True,
    )
    length_controls = mo.hstack(
        [log_sample_length, log_time_over_tau],
        justify="start",
        align="center",
        wrap=True,
        gap=1.5,
    )
    return length_controls, log_sample_length, log_time_over_tau


@app.cell
def _(length_controls, mo):
    mo.vstack(
        [
            mo.md(r"""
            ## 7. From diffusivity to a relaxation time

            A diffusion length \(L\) introduces the scaling

            \[
            \tau\sim\frac{L^2}{D_{\rm chem}}.
            \]

            Boundary conditions add numerical factors such as \(\pi^2\), but
            the \(L^2\) dependence is the essential reason a thin film can
            equilibrate in milliseconds while a millimeter body takes days.
            This length–time bridge underlies later interpretation of
            polarization and intermittent-titration experiments.
            """),
            length_controls,
        ]
    )
    return


@app.cell
def _(
    characteristic_relaxation_time,
    log_sample_length,
    log_time_over_tau,
    np,
    selected_chemical_diffusivity_m2_per_s,
    slab_remaining_profile,
):
    selected_length_m = 10.0 ** float(log_sample_length.value)
    selected_relaxation_time_s = characteristic_relaxation_time(
        selected_length_m,
        selected_chemical_diffusivity_m2_per_s,
    )
    length_curve_m = np.logspace(-9.0, -3.0, 400)
    relaxation_curve_s = (
        length_curve_m**2 / selected_chemical_diffusivity_m2_per_s
    )
    selected_time_over_tau = 10.0 ** float(log_time_over_tau.value)
    slab_coordinate = np.linspace(0.0, 1.0, 401)
    slab_excess_profile = slab_remaining_profile(
        slab_coordinate,
        selected_time_over_tau,
    )
    return (
        length_curve_m,
        relaxation_curve_s,
        selected_length_m,
        selected_relaxation_time_s,
        selected_time_over_tau,
        slab_coordinate,
        slab_excess_profile,
    )


@app.cell
def _(
    length_curve_m,
    mo,
    np,
    plt,
    relaxation_curve_s,
    selected_length_m,
    selected_relaxation_time_s,
    selected_time_over_tau,
    slab_coordinate,
    slab_excess_profile,
):
    def _readable_length(length_m):
        if length_m < 1.0e-6:
            return f"{length_m * 1.0e9:.3g} nm"
        if length_m < 1.0e-3:
            return f"{length_m * 1.0e6:.3g} µm"
        return f"{length_m * 1.0e3:.3g} mm"

    def _readable_time(time_s):
        if time_s < 1.0e-3:
            return f"{time_s * 1.0e6:.3g} µs"
        if time_s < 1.0:
            return f"{time_s * 1.0e3:.3g} ms"
        if time_s < 3600.0:
            return f"{time_s:.3g} s"
        if time_s < 86400.0:
            return f"{time_s / 3600.0:.3g} h"
        return f"{time_s / 86400.0:.3g} days"

    length_figure, (time_axis, slab_axis) = plt.subplots(
        1,
        2,
        figsize=(13.2, 4.9),
        dpi=120,
    )
    time_axis.loglog(
        length_curve_m,
        relaxation_curve_s,
        color="#007C91",
        lw=3.0,
    )
    time_axis.scatter(
        [selected_length_m],
        [selected_relaxation_time_s],
        s=100,
        color="#EE9B00",
        edgecolor="#222222",
        zorder=5,
    )
    time_axis.set(
        xlabel="sample thickness, L (m)",
        ylabel=r"$L^2/D_{\rm chem}$ (s)",
        title="Diffusion time grows as length squared",
    )
    time_axis.grid(which="both", alpha=0.22)

    slab_axis.plot(
        slab_coordinate,
        slab_excess_profile,
        color="#D55E00",
        lw=3.0,
    )
    slab_axis.fill_between(
        slab_coordinate,
        0.0,
        slab_excess_profile,
        color="#D55E00",
        alpha=0.12,
    )
    slab_axis.set(
        xlabel="position / L",
        ylabel="normalized concentration excess",
        title=rf"Slab relaxation at $t/\tau={selected_time_over_tau:.3g}$",
        xlim=(0.0, 1.0),
        ylim=(0.0, 1.08),
    )
    slab_axis.grid(alpha=0.22)
    length_figure.tight_layout()
    plt.close(length_figure)

    length_summary = mo.md(
        rf"""
        For **{_readable_length(selected_length_m)}**,
        \(L^2/D_{{\rm chem}}\) is
        **{_readable_time(selected_relaxation_time_s)}**. The profile uses
        fixed surface composition at both faces; a different geometry or
        boundary condition changes the prefactor, not the \(L^2/D\) scaling.
        """
    )
    mo.vstack([length_figure, length_summary])
    return (length_figure,)

@app.cell
def _(
    FARADAY_C_PER_MOL,
    GAS_CONSTANT_J_PER_MOL_K,
    ambipolar_fluxes,
    ambipolar_gradient_mol_per_m4,
    ambipolar_concentration_mol_per_m3,
    ambipolar_internal_potential_gradient,
    analytic_common_flux_mol_per_m2_s,
    attempt_frequency_hz,
    backward_rate_hz,
    biased_neighbor_rates,
    charge_number_value,
    chemical_diffusivity,
    conductivity_form_chemical_diffusivity_m2_per_s,
    defect_diffusivity_m2_per_s,
    discrete_bond_fluxes,
    electronic_diffusivity_m2_per_s,
    electronic_flux_mol_per_m2_s,
    evolved_master_occupancy,
    exact_hopping_drift_velocity,
    extracted_diffusivity_m2_per_s,
    field_work_per_hop_ev,
    fick_bond_fluxes,
    forward_rate_hz,
    hopping_diffusivity_1d,
    ideal_chemical_diffusivity_m2_per_s,
    ionic_diffusivity_m2_per_s,
    ionic_flux_mol_per_m2_s,
    jump_distance_m,
    migration_enthalpy_ev,
    master_mass_relative_error,
    molar_nernst_planck_flux,
    msd_fit_r_squared,
    nernst_einstein_drift_velocity,
    np,
    selected_chemical_diffusivity_m2_per_s,
    temperature_k,
    thermodynamic_factor_value,
    zero_field_rate_hz,
):
    def _relative_error(value, reference):
        return abs(float(value) - float(reference)) / max(
            abs(float(reference)),
            1.0e-300,
        )

    zero_field_identity_error = _relative_error(
        defect_diffusivity_m2_per_s,
        jump_distance_m**2 * zero_field_rate_hz,
    )
    msd_diffusivity_error = _relative_error(
        extracted_diffusivity_m2_per_s,
        defect_diffusivity_m2_per_s,
    )
    detailed_balance_reference = np.exp(
        field_work_per_hop_ev / (8.617333262e-5 * temperature_k)
    )
    detailed_balance_error = _relative_error(
        forward_rate_hz / backward_rate_hz,
        detailed_balance_reference,
    )

    low_field_test_e_v_per_m = (
        1.0e-4
        * 8.617333262e-5
        * temperature_k
        / jump_distance_m
    )
    (
        low_field_forward_rate_hz,
        low_field_backward_rate_hz,
        low_field_base_rate_hz,
    ) = biased_neighbor_rates(
        temperature_k,
        migration_enthalpy_ev,
        attempt_frequency_hz,
        jump_distance_m,
        1,
        low_field_test_e_v_per_m,
    )
    low_field_exact_velocity_m_per_s = exact_hopping_drift_velocity(
        jump_distance_m,
        low_field_forward_rate_hz,
        low_field_backward_rate_hz,
    )
    low_field_reference_diffusivity_m2_per_s = hopping_diffusivity_1d(
        jump_distance_m,
        low_field_base_rate_hz,
    )
    low_field_ne_velocity_m_per_s = nernst_einstein_drift_velocity(
        low_field_reference_diffusivity_m2_per_s,
        temperature_k,
        1,
        low_field_test_e_v_per_m,
    )
    low_field_ne_error = _relative_error(
        low_field_exact_velocity_m_per_s,
        low_field_ne_velocity_m_per_s,
    )

    validation_master_flux = discrete_bond_fluxes(
        evolved_master_occupancy,
        zero_field_rate_hz,
    )
    validation_fick_flux = fick_bond_fluxes(
        evolved_master_occupancy,
        jump_distance_m,
        defect_diffusivity_m2_per_s,
    )
    master_fick_error = float(
        np.max(np.abs(validation_master_flux - validation_fick_flux))
        / max(np.max(np.abs(validation_master_flux)), 1.0e-300)
    )

    validation_concentration = 1000.0
    validation_relative_gradient = 2.5e5
    validation_gradient = (
        validation_concentration * validation_relative_gradient
    )
    validation_potential_gradient = (
        -GAS_CONSTANT_J_PER_MOL_K
        * temperature_k
        * validation_relative_gradient
        / (charge_number_value * FARADAY_C_PER_MOL)
    )
    validation_diffusion_flux, validation_electrical_flux, validation_total_flux = (
        molar_nernst_planck_flux(
            defect_diffusivity_m2_per_s,
            validation_concentration,
            validation_gradient,
            validation_potential_gradient,
            charge_number_value,
            temperature_k,
        )
    )
    electrochemical_cancellation_error = abs(validation_total_flux) / max(
        abs(validation_diffusion_flux) + abs(validation_electrical_flux),
        1.0e-300,
    )

    validation_internal_gradient = ambipolar_internal_potential_gradient(
        ionic_diffusivity_m2_per_s,
        electronic_diffusivity_m2_per_s,
        ambipolar_concentration_mol_per_m3,
        ambipolar_gradient_mol_per_m4,
        temperature_k,
    )
    validation_ionic_flux, validation_electronic_flux = ambipolar_fluxes(
        ionic_diffusivity_m2_per_s,
        electronic_diffusivity_m2_per_s,
        ambipolar_concentration_mol_per_m3,
        ambipolar_gradient_mol_per_m4,
        validation_internal_gradient,
        temperature_k,
    )
    ambipolar_flux_match_error = _relative_error(
        validation_ionic_flux,
        validation_electronic_flux,
    )
    ambipolar_current_relative = abs(
        FARADAY_C_PER_MOL
        * (validation_ionic_flux - validation_electronic_flux)
    ) / max(
        FARADAY_C_PER_MOL
        * (abs(validation_ionic_flux) + abs(validation_electronic_flux)),
        1.0e-300,
    )
    analytic_ambipolar_flux_error = _relative_error(
        validation_ionic_flux,
        -ideal_chemical_diffusivity_m2_per_s
        * ambipolar_gradient_mol_per_m4,
    )
    displayed_flux_consistency_error = max(
        _relative_error(ionic_flux_mol_per_m2_s, analytic_common_flux_mol_per_m2_s),
        _relative_error(electronic_flux_mol_per_m2_s, analytic_common_flux_mol_per_m2_s),
    )
    theta_multiplication_error = _relative_error(
        selected_chemical_diffusivity_m2_per_s,
        ideal_chemical_diffusivity_m2_per_s * thermodynamic_factor_value,
    )
    conductivity_form_error = _relative_error(
        conductivity_form_chemical_diffusivity_m2_per_s,
        selected_chemical_diffusivity_m2_per_s,
    )

    positive_quantities = np.array(
        [
            zero_field_rate_hz,
            forward_rate_hz,
            backward_rate_hz,
            defect_diffusivity_m2_per_s,
            ionic_diffusivity_m2_per_s,
            electronic_diffusivity_m2_per_s,
            ideal_chemical_diffusivity_m2_per_s,
            selected_chemical_diffusivity_m2_per_s,
        ]
    )
    finite_quantities = np.array(
        [
            zero_field_identity_error,
            msd_diffusivity_error,
            detailed_balance_error,
            low_field_ne_error,
            master_fick_error,
            electrochemical_cancellation_error,
            ambipolar_flux_match_error,
            ambipolar_current_relative,
            analytic_ambipolar_flux_error,
            displayed_flux_consistency_error,
            theta_multiplication_error,
            conductivity_form_error,
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
        "low_field_ne_error": low_field_ne_error,
        "low_field_ne_pass": low_field_ne_error < 1.0e-7,
        "master_mass_error": master_mass_relative_error,
        "master_mass_pass": master_mass_relative_error < 1.0e-12,
        "master_fick_error": master_fick_error,
        "master_fick_pass": master_fick_error < 1.0e-12,
        "cancellation_error": electrochemical_cancellation_error,
        "cancellation_pass": electrochemical_cancellation_error < 1.0e-12,
        "ambipolar_current_relative": ambipolar_current_relative,
        "ambipolar_current_pass": ambipolar_current_relative < 1.0e-12,
        "ambipolar_flux_match_error": ambipolar_flux_match_error,
        "ambipolar_flux_match_pass": ambipolar_flux_match_error < 1.0e-12,
        "analytic_ambipolar_error": analytic_ambipolar_flux_error,
        "analytic_ambipolar_pass": analytic_ambipolar_flux_error < 1.0e-12,
        "displayed_flux_error": displayed_flux_consistency_error,
        "displayed_flux_pass": displayed_flux_consistency_error < 1.0e-12,
        "theta_error": theta_multiplication_error,
        "theta_pass": theta_multiplication_error < 1.0e-14,
        "conductivity_form_error": conductivity_form_error,
        "conductivity_form_pass": conductivity_form_error < 1.0e-12,
        "positive_pass": bool(np.all(positive_quantities > 0.0)),
        "finite_pass": bool(
            np.all(np.isfinite(positive_quantities))
            and np.all(np.isfinite(finite_quantities))
            and np.all(np.isfinite(evolved_master_occupancy))
        ),
    }
    return (transport_validation,)


@app.cell
def _(mo, transport_validation):
    def _check_mark(passed):
        return "PASS" if passed else "CHECK"

    mo.md(
        rf"""
        ## 8. Numerical sanity checks

        | Check | Result | Error or diagnostic |
        |---|---:|---:|
        | zero-field random-walk identity \(D=a^2w_0\) | {_check_mark(transport_validation['zero_field_identity_pass'])} | relative error {transport_validation['zero_field_identity_error']:.2e} |
        | simulated MSD is linear and recovers \(D\) | {_check_mark(transport_validation['msd_pass'])} | \(D\) error {transport_validation['msd_diffusivity_error']:.2e}, \(R^2={transport_validation['msd_r_squared']:.6f}\) |
        | biased rates obey detailed balance \(w_+/w_-=\exp(zeEa/k_BT)\) | {_check_mark(transport_validation['detailed_balance_pass'])} | relative error {transport_validation['detailed_balance_error']:.2e} |
        | exact hopping drift approaches Nernst–Einstein at low field | {_check_mark(transport_validation['low_field_ne_pass'])} | relative error {transport_validation['low_field_ne_error']:.2e} |
        | master equation conserves defect number | {_check_mark(transport_validation['master_mass_pass'])} | relative error {transport_validation['master_mass_error']:.2e} |
        | microscopic bond flux equals discrete Fick flux | {_check_mark(transport_validation['master_fick_pass'])} | relative error {transport_validation['master_fick_error']:.2e} |
        | opposing chemical and electrical gradients give zero total flux | {_check_mark(transport_validation['cancellation_pass'])} | scaled residual {transport_validation['cancellation_error']:.2e} |
        | ambipolar solution has zero open-circuit current | {_check_mark(transport_validation['ambipolar_current_pass'])} | scaled residual {transport_validation['ambipolar_current_relative']:.2e} |
        | ambipolar solution gives \(J_i=J_e\) | {_check_mark(transport_validation['ambipolar_flux_match_pass'])} | relative mismatch {transport_validation['ambipolar_flux_match_error']:.2e} |
        | common flux gives \(2D_iD_e/(D_i+D_e)\) | {_check_mark(transport_validation['analytic_ambipolar_pass'])} | relative error {transport_validation['analytic_ambipolar_error']:.2e} |
        | displayed ionic and electronic fluxes match the analytic flux | {_check_mark(transport_validation['displayed_flux_pass'])} | maximum relative error {transport_validation['displayed_flux_error']:.2e} |
        | \(\Theta\) multiplies the ideal chemical diffusivity | {_check_mark(transport_validation['theta_pass'])} | relative error {transport_validation['theta_error']:.2e} |
        | conductivity form matches \(D_{{\rm chem,ideal}}\Theta\) | {_check_mark(transport_validation['conductivity_form_pass'])} | relative error {transport_validation['conductivity_form_error']:.2e} |
        | rates and diffusivities remain positive | {_check_mark(transport_validation['positive_pass'])} | strict positivity |
        | all validation quantities remain finite | {_check_mark(transport_validation['finite_pass'])} | finite arrays and scalars |

        The stochastic MSD check uses a fixed random seed and a tolerance that
        reflects finite ensemble sampling. The conservation, detailed-balance,
        flux-cancellation, and ambipolar checks test identities independently of
        what is visually resolvable in a plot.

        **Why check these?** The first six rows connect the atomic jump model to
        diffusion and drift. The next row verifies the central equilibrium idea
        that chemical and electrical driving forces can cancel. The ambipolar
        rows verify that ions and electrons move together without net current.
        The final rows check the thermodynamic factor and guard against
        nonphysical numerical values.
        """
    )
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## Take-home map

    1. A migration barrier sets a per-neighbor hopping rate.
    2. Random hops generate \(\langle x^2\rangle=2Dt\); a concentration
       imbalance turns symmetric exchange into Fick flux.
    3. An electric field biases forward and backward barriers, producing the
       Nernst–Einstein drift in the low-field limit.
    4. Charged species respond to gradients of
       \(\widetilde{\mu}=\mu+zF\phi\). Chemical and electrical pieces can be
       nonzero while their sum—and therefore the flux—is zero.
    5. During neutral composition relaxation, local electroneutrality couples
       ionic and electronic motion. The internal field slows the faster carrier
       and accelerates the slower one.
    6. \(D_{\rm chem}\) combines kinetic coefficients with a thermodynamic
       factor, and the experimental time scale grows as \(L^2/D_{\rm chem}\).

    **Model boundaries.** The lattice is one-dimensional and ideal; transition
    states are symmetric; hops are uncorrelated; the field is spatially uniform
    in the hopping example; the Nernst–Planck equations use ideal dilute
    activities until \(\Theta\) is introduced; and the ambipolar example assumes
    monovalent carriers, local equilibrium, and local electroneutrality.
    """)
    return


if __name__ == "__main__":
    app.run()
