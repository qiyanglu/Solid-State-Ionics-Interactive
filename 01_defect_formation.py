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
    import math

    import marimo as mo
    import matplotlib.pyplot as plt
    import numpy as np
    from scipy import optimize, special

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
    return math, mo, np, optimize, plt, special


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
    # Defect formation thermodynamics

    **Why do defects exist at equilibrium when each one costs energy?**

    A perfect crystal has low formation energy, but a crystal containing defects
    can be arranged in many more ways. Equilibrium is the compromise between
    those two tendencies:

    \[
    G(n)=n\Delta g_f^0-T S_{\rm config}(n).
    \]

    Use the lattice picture first. Then we will count its configurations and see
    the same equilibrium appear as both a minimum of \(G\) and a zero of the
    defect chemical potential.
    """)
    return


@app.cell
def _(np, special):
    KB_EV_PER_K = 8.617333262e-5

    def validate_site_count(site_count):
        """Return an integer number of equivalent sites."""
        count = int(site_count)
        if count != site_count or count < 2:
            raise ValueError("site_count must be an integer of at least 2")
        return count

    def formation_free_energy_ev(
        temperature_k,
        formation_enthalpy_ev,
        formation_entropy_kb,
    ):
        """Return Delta g_f^0 in eV per defect.

        formation_entropy_kb is the non-configurational formation entropy
        (for example, vibrational entropy) in units of k_B per defect.
        """
        if not np.isfinite(temperature_k) or temperature_k <= 0.0:
            raise ValueError("temperature_k must be positive and finite")
        if not np.isfinite(formation_enthalpy_ev):
            raise ValueError("formation_enthalpy_ev must be finite")
        if not np.isfinite(formation_entropy_kb):
            raise ValueError("formation_entropy_kb must be finite")
        return (
            float(formation_enthalpy_ev)
            - KB_EV_PER_K * float(temperature_k) * float(formation_entropy_kb)
        )

    def log_multiplicity(site_count, defect_count):
        """Return ln[binomial(N, n)] for integer finite-lattice states."""
        count = validate_site_count(site_count)
        defects = np.asarray(defect_count, dtype=float)
        if np.any(~np.isfinite(defects)):
            raise ValueError("defect_count must be finite")
        if np.any(defects < 0.0) or np.any(defects > count):
            raise ValueError("defect_count must lie between 0 and N")
        if np.any(defects != np.rint(defects)):
            raise ValueError("defect_count must contain integers")
        return (
            special.gammaln(count + 1.0)
            - special.gammaln(defects + 1.0)
            - special.gammaln(count - defects + 1.0)
        )

    def stirling_entropy_kb_per_site(defect_fraction):
        """Return s_config/k_B per site, including its endpoint limits."""
        fraction = np.asarray(defect_fraction, dtype=float)
        if np.any(~np.isfinite(fraction)):
            raise ValueError("defect_fraction must be finite")
        if np.any(fraction < 0.0) or np.any(fraction > 1.0):
            raise ValueError("defect_fraction must lie in [0, 1]")
        entropy = np.zeros_like(fraction)
        interior = (fraction > 0.0) & (fraction < 1.0)
        values = fraction[interior]
        entropy[interior] = -(
            values * np.log(values)
            + (1.0 - values) * np.log1p(-values)
        )
        return entropy

    def stirling_free_energy_ev_per_site(
        defect_fraction,
        temperature_k,
        formation_free_energy_ev_value,
    ):
        """Return G/N in eV per site in the Stirling limit."""
        fraction = np.asarray(defect_fraction, dtype=float)
        entropy_kb = stirling_entropy_kb_per_site(fraction)
        return (
            fraction * float(formation_free_energy_ev_value)
            - KB_EV_PER_K * float(temperature_k) * entropy_kb
        )

    def exact_free_energy_ev_per_site(
        site_count,
        defect_count,
        temperature_k,
        formation_free_energy_ev_value,
    ):
        """Return exact finite-N G/N in eV per site."""
        count = validate_site_count(site_count)
        defects = np.asarray(defect_count, dtype=float)
        return (
            defects * float(formation_free_energy_ev_value)
            - KB_EV_PER_K
            * float(temperature_k)
            * log_multiplicity(count, defects)
        ) / count

    def chemical_potential_ev(
        defect_fraction,
        temperature_k,
        formation_free_energy_ev_value,
    ):
        """Return d(G/N)/dx in eV per defect in the Stirling limit."""
        fraction = np.asarray(defect_fraction, dtype=float)
        if np.any(fraction <= 0.0) or np.any(fraction >= 1.0):
            raise ValueError("chemical potential requires 0 < x < 1")
        return (
            float(formation_free_energy_ev_value)
            + KB_EV_PER_K
            * float(temperature_k)
            * (np.log(fraction) - np.log1p(-fraction))
        )

    def equilibrium_fraction(temperature_k, formation_free_energy_ev_value):
        """Return the exact thermodynamic-limit logistic equilibrium."""
        reduced_energy = float(formation_free_energy_ev_value) / (
            KB_EV_PER_K * float(temperature_k)
        )
        return float(special.expit(-reduced_energy))

    def dilute_fraction(temperature_k, formation_free_energy_ev_value):
        """Return exp[-Delta g_f^0/(k_B T)], without asserting validity."""
        reduced_energy = float(formation_free_energy_ev_value) / (
            KB_EV_PER_K * float(temperature_k)
        )
        return float(np.exp(np.clip(-reduced_energy, -745.0, 700.0)))

    def finite_equilibrium_index(
        site_count,
        temperature_k,
        formation_free_energy_ev_value,
    ):
        """Return the discrete n minimizing exact finite-N free energy."""
        count = validate_site_count(site_count)
        states = np.arange(count + 1)
        energies = exact_free_energy_ev_per_site(
            count,
            states,
            temperature_k,
            formation_free_energy_ev_value,
        )
        return int(np.argmin(energies))

    return (
        KB_EV_PER_K,
        chemical_potential_ev,
        dilute_fraction,
        equilibrium_fraction,
        exact_free_energy_ev_per_site,
        finite_equilibrium_index,
        formation_free_energy_ev,
        log_multiplicity,
        stirling_entropy_kb_per_site,
        stirling_free_energy_ev_per_site,
        validate_site_count,
    )


@app.cell
def _(mo):
    temperature = mo.ui.slider(
        start=300,
        stop=1800,
        step=25,
        value=1000,
        label="Temperature, T (K)",
        show_value=True,
    )
    formation_enthalpy = mo.ui.slider(
        start=0.10,
        stop=2.50,
        step=0.05,
        value=0.45,
        label="Formation enthalpy (eV per defect)",
        show_value=True,
    )
    formation_entropy = mo.ui.slider(
        start=-2.0,
        stop=10.0,
        step=0.25,
        value=3.0,
        label="Formation entropy (kB per defect)",
        show_value=True,
    )
    lattice_sites = mo.ui.slider(
        start=16,
        stop=400,
        step=8,
        value=200,
        label="Finite lattice size, N",
        show_value=True,
    )
    show_components = mo.ui.checkbox(
        value=True,
        label="Show energetic and entropic contributions",
    )
    show_exact_points = mo.ui.checkbox(
        value=True,
        label="Show exact finite-lattice states",
    )
    controls = mo.hstack(
        [temperature, formation_enthalpy, formation_entropy],
        justify="start",
        align="center",
        wrap=True,
        gap=1.5,
    )
    advanced_controls = mo.accordion(
        {
            "Explore finite-size and display options": mo.vstack(
                [lattice_sites, show_components, show_exact_points],
                gap=0.7,
            )
        }
    )
    return (
        advanced_controls,
        controls,
        formation_enthalpy,
        formation_entropy,
        lattice_sites,
        show_components,
        show_exact_points,
        temperature,
    )


@app.cell
def _(advanced_controls, controls, mo):
    mo.vstack(
        [
            mo.md("## Explore the competition"),
            controls,
            mo.md(r"""
            Change one quantity at a time. Higher temperature or a more
            positive formation entropy lowers the effective defect cost and
            should move equilibrium toward a larger defect fraction.

            The default is intentionally easy to see; it is a teaching state,
            not a claim that oxide defects are usually this concentrated.
            """),
            advanced_controls,
        ]
    )
    return


@app.cell
def _(
    KB_EV_PER_K,
    dilute_fraction,
    equilibrium_fraction,
    exact_free_energy_ev_per_site,
    finite_equilibrium_index,
    formation_enthalpy,
    formation_entropy,
    formation_free_energy_ev,
    lattice_sites,
    log_multiplicity,
    math,
    np,
    special,
    stirling_entropy_kb_per_site,
    temperature,
    validate_site_count,
):
    temperature_k = float(temperature.value)
    formation_enthalpy_ev = float(formation_enthalpy.value)
    formation_entropy_kb = float(formation_entropy.value)
    site_count = validate_site_count(int(lattice_sites.value))

    delta_g0_ev = formation_free_energy_ev(
        temperature_k,
        formation_enthalpy_ev,
        formation_entropy_kb,
    )
    kbt_ev = KB_EV_PER_K * temperature_k
    reduced_formation_energy = delta_g0_ev / kbt_ev
    equilibrium_x = equilibrium_fraction(temperature_k, delta_g0_ev)
    dilute_x = dilute_fraction(temperature_k, delta_g0_ev)

    discrete_n = np.arange(site_count + 1)
    discrete_x = discrete_n / site_count
    log_omega = log_multiplicity(site_count, discrete_n)
    entropy_exact_kb_per_site = log_omega / site_count
    entropy_stirling_on_discrete = stirling_entropy_kb_per_site(discrete_x)
    free_energy_exact_ev_per_site = exact_free_energy_ev_per_site(
        site_count,
        discrete_n,
        temperature_k,
        delta_g0_ev,
    )
    finite_n_eq = finite_equilibrium_index(
        site_count,
        temperature_k,
        delta_g0_ev,
    )
    finite_x_eq = finite_n_eq / site_count
    omega_at_finite_min = math.comb(site_count, finite_n_eq)

    _equilibrium_logit = -reduced_formation_energy
    _focus_logit = np.linspace(
        max(-700.0, _equilibrium_logit - 3.5),
        min(700.0, _equilibrium_logit + 3.5),
        601,
    )
    focus_x = special.expit(_focus_logit)
    formation_term_ev = focus_x * delta_g0_ev
    entropy_term_ev = -kbt_ev * stirling_entropy_kb_per_site(focus_x)
    total_free_energy_ev = formation_term_ev + entropy_term_ev
    chemical_potential_values_ev = delta_g0_ev + kbt_ev * _focus_logit
    return (
        chemical_potential_values_ev,
        delta_g0_ev,
        dilute_x,
        discrete_n,
        discrete_x,
        entropy_exact_kb_per_site,
        entropy_term_ev,
        equilibrium_x,
        finite_n_eq,
        finite_x_eq,
        focus_x,
        formation_enthalpy_ev,
        formation_entropy_kb,
        formation_term_ev,
        free_energy_exact_ev_per_site,
        log_omega,
        omega_at_finite_min,
        reduced_formation_energy,
        site_count,
        temperature_k,
        total_free_energy_ev,
    )


@app.cell
def _(
    delta_g0_ev,
    dilute_x,
    equilibrium_x,
    finite_n_eq,
    finite_x_eq,
    formation_enthalpy_ev,
    formation_entropy_kb,
    mo,
    omega_at_finite_min,
    reduced_formation_energy,
    temperature_k,
):
    _dilute_error = abs(dilute_x - equilibrium_x) / equilibrium_x
    if omega_at_finite_min < 1_000_000:
        _omega_text = f"{omega_at_finite_min:,}"
    else:
        _omega_text = (
            rf"$10^{{{len(str(omega_at_finite_min)) - 1:d}}}$"
            " (order of magnitude)"
        )
    _dilute_status = (
        "accurate in this state"
        if equilibrium_x < 0.01
        else "not a dilute approximation in this state"
    )

    _state_table = mo.md(
        rf"""
        ### Current thermodynamic state

        | Quantity | Value |
        |---|---:|
        | \(T\) | {temperature_k:.0f} K |
        | \(\Delta h_f\) | {formation_enthalpy_ev:.3f} eV/defect |
        | \(\Delta s_f^0\) | {formation_entropy_kb:.2f} \(k_B\)/defect |
        | \(\Delta g_f^0(T)\) | {delta_g0_ev:.4f} eV/defect |
        | \(\Delta g_f^0/(k_BT)\) | {reduced_formation_energy:.3f} |
        | thermodynamic \(x_{{\rm eq}}\) | {equilibrium_x:.4e} |
        | most probable finite-\(N\) macrostate | \(n={finite_n_eq}\), \(x={finite_x_eq:.4e}\) |
        | \(\Omega(N,n_{{\rm min}})\) | {_omega_text} |
        | dilute \(\exp[-\Delta g_f^0/(k_BT)]\) | {dilute_x:.4e} |

        The dilute relative error is {_dilute_error:.2e}; it is
        **{_dilute_status}**. The finite-\(N\) mode is discrete, so its most
        probable macrostate can be \(n=0\) in a small classroom lattice even
        when the ensemble-mean occupancy is nonzero.
        """
    )
    mo.accordion({"Explore further — current state": _state_table})
    return


@app.cell
def _(finite_n_eq, np, plt, site_count):
    _columns = int(np.ceil(np.sqrt(site_count)))
    _rows = int(np.ceil(site_count / _columns))
    _indices = np.arange(site_count)
    _x_positions = _indices % _columns
    _y_positions = _rows - 1 - _indices // _columns
    _defect_mask = np.zeros(site_count, dtype=bool)
    if finite_n_eq > 0:
        _rng = np.random.default_rng(2024)
        _defect_indices = _rng.choice(
            site_count,
            size=finite_n_eq,
            replace=False,
        )
        _defect_mask[_defect_indices] = True

    plt.rcParams.update(
        {
            "font.size": 13,
            "axes.titlesize": 15,
            "axes.labelsize": 13,
            "legend.fontsize": 10.5,
        }
    )
    _lattice_fig, _lattice_axis = plt.subplots(figsize=(11.5, 4.8), dpi=120)
    _lattice_axis.scatter(
        _x_positions[~_defect_mask],
        _y_positions[~_defect_mask],
        s=115,
        color="#4C7C86",
        edgecolor="white",
        linewidth=0.8,
        label="occupied site",
        zorder=2,
    )
    if np.any(_defect_mask):
        _lattice_axis.scatter(
            _x_positions[_defect_mask],
            _y_positions[_defect_mask],
            s=130,
            marker="s",
            facecolor="white",
            edgecolor="#B8734A",
            linewidth=1.7,
            label="defect",
            zorder=3,
        )
    else:
        _lattice_axis.scatter(
            [],
            [],
            s=130,
            marker="s",
            facecolor="white",
            edgecolor="#B8734A",
            linewidth=1.7,
            label="defect (none in this finite-N mode)",
        )
    _lattice_axis.set_aspect("equal")
    _lattice_axis.set_xlim(-1.0, _columns)
    _lattice_axis.set_ylim(-1.0, _rows)
    _lattice_axis.axis("off")
    _lattice_axis.set_title("A finite lattice at its most probable composition")
    _lattice_axis.legend(
        loc="center left",
        bbox_to_anchor=(1.01, 0.5),
        frameon=False,
    )
    _lattice_fig.tight_layout()
    plt.close(_lattice_fig)
    lattice_figure = _lattice_fig
    return (lattice_figure,)


@app.cell
def _(lattice_figure):
    lattice_figure
    return


@app.cell
def _(finite_n_eq, finite_x_eq, mo, site_count):
    _finite_note = (
        "This finite lattice cannot display less than one defect, so its most "
        "probable composition can be the perfect lattice even when the "
        "thermodynamic mean is nonzero."
        if finite_n_eq == 0
        else
        "The open squares are defects placed at reproducibly randomized sites. "
        "Moving them among equivalent sites changes the configuration without "
        "changing the formation-energy term."
    )
    mo.md(
        rf"""
        \(N={site_count}\), \(n={finite_n_eq}\), and
        \(x=n/N={finite_x_eq:.4g}\). {_finite_note}

        This lattice is a **configurational schematic**, not an atomistic
        simulation.

        ## From arrangements to entropy

        For \(n\) defects on \(N\) equivalent sites,

        \[
        \Omega(N,n)=\binom{{N}}{{n}},\qquad
        S_{{\rm config}}=k_B\ln\Omega.
        \]

        A perfect lattice has only one arrangement. Mixed occupied and defect
        sites have many, creating an entropic reason for defects to appear.
        """
    )
    return


@app.cell
def _(
    discrete_x,
    entropy_exact_kb_per_site,
    finite_x_eq,
    log_omega,
    np,
    plt,
    site_count,
    stirling_entropy_kb_per_site,
):
    _entropy_x = np.linspace(0.0, 1.0, 801)
    _entropy_stirling = stirling_entropy_kb_per_site(_entropy_x)

    _entropy_fig, (_multiplicity_axis, _entropy_axis) = plt.subplots(
        1,
        2,
        figsize=(13.0, 5.2),
        dpi=120,
    )
    _multiplicity_axis.plot(
        discrete_x,
        log_omega,
        "o",
        ms=3.2,
        color="#4C7C86",
        alpha=0.75,
        label=r"exact $\ln\binom{N}{n}$",
    )
    _multiplicity_axis.plot(
        _entropy_x,
        site_count * _entropy_stirling,
        color="#B8734A",
        lw=1.9,
        label="Stirling limit",
    )
    _multiplicity_axis.axvline(finite_x_eq, color="#40464D", ls=":", lw=1.2)
    _multiplicity_axis.set(
        xlabel="Defect fraction, x = n/N",
        ylabel=r"$\ln\Omega$",
        title="How many arrangements are possible?",
        xlim=(0.0, 1.0),
    )
    _multiplicity_axis.grid(alpha=0.25)
    _multiplicity_axis.legend(frameon=False)

    _entropy_axis.plot(
        discrete_x,
        entropy_exact_kb_per_site,
        "o",
        ms=3.2,
        color="#4C7C86",
        alpha=0.75,
        label=r"exact $S_{\rm config}/(Nk_B)$",
    )
    _entropy_axis.plot(
        _entropy_x,
        _entropy_stirling,
        color="#B8734A",
        lw=1.9,
        label="Stirling limit",
    )
    _entropy_axis.axvline(finite_x_eq, color="#40464D", ls=":", lw=1.2)
    _entropy_axis.set(
        xlabel="Defect fraction, x = n/N",
        ylabel=r"Configurational entropy per site, $s_{\rm config}/k_B$",
        title="Mixing creates configurational entropy",
        xlim=(0.0, 1.0),
    )
    _entropy_axis.grid(alpha=0.25)
    _entropy_axis.legend(frameon=False)
    _entropy_fig.tight_layout()
    plt.close(_entropy_fig)
    entropy_figure = _entropy_fig
    return (entropy_figure,)


@app.cell
def _(entropy_figure):
    entropy_figure
    return


@app.cell
def _(mo):
    mo.md(r"""
    A perfect lattice has only one configuration, so
    \(\ln\Omega=0\) and \(S_{\rm config}=0\). Adding defects initially
    creates many new choices. Entropy is maximized near equal populations of
    occupied and defect sites. The continuous Stirling line approaches the
    marker-displayed finite-\(N\) states, whose allowed compositions remain
    separated by $1/N$.
    """)
    return


@app.cell
def _(
    chemical_potential_values_ev,
    discrete_x,
    entropy_term_ev,
    equilibrium_x,
    finite_n_eq,
    finite_x_eq,
    focus_x,
    formation_term_ev,
    free_energy_exact_ev_per_site,
    np,
    plt,
    show_components,
    show_exact_points,
    total_free_energy_ev,
):
    _energy_fig, (_free_energy_axis, _mu_axis) = plt.subplots(
        1,
        2,
        figsize=(13.2, 5.8),
        dpi=120,
        gridspec_kw={"width_ratios": [1.45, 1.0]},
    )
    _energy_scale = 1000.0
    if show_components.value:
        _free_energy_axis.plot(
            focus_x,
            _energy_scale * formation_term_ev,
            color="#B65C4A",
            lw=1.7,
            label=r"formation: $x\Delta g_f^0$",
        )
        _free_energy_axis.plot(
            focus_x,
            _energy_scale * entropy_term_ev,
            color="#7C6A91",
            lw=1.7,
            label=r"entropy: $-Ts_{\rm config}$",
        )
    _free_energy_axis.plot(
        focus_x,
        _energy_scale * total_free_energy_ev,
        color="#4C7C86",
        lw=2.0,
        label=r"total: $G/N$",
        zorder=4,
    )
    if show_exact_points.value:
        _finite_mask = (
            (discrete_x > 0.0)
            & (discrete_x < 1.0)
            & (discrete_x >= focus_x[0])
            & (discrete_x <= focus_x[-1])
        )
        _free_energy_axis.plot(
            discrete_x[_finite_mask],
            _energy_scale * free_energy_exact_ev_per_site[_finite_mask],
            "o",
            ms=4.0,
            mfc="white",
            mec="#40464D",
            alpha=0.8,
            label="exact finite-N states",
            zorder=5,
        )
    _equilibrium_energy = float(
        np.interp(equilibrium_x, focus_x, total_free_energy_ev)
    )
    _free_energy_axis.scatter(
        [equilibrium_x],
        [_energy_scale * _equilibrium_energy],
        s=105,
        color="#C49345",
        edgecolor="#40464D",
        zorder=7,
        label=r"thermodynamic $x_{\rm eq}$",
    )
    if finite_n_eq > 0 and focus_x[0] <= finite_x_eq <= focus_x[-1]:
        _finite_energy = free_energy_exact_ev_per_site[finite_n_eq]
        _free_energy_axis.scatter(
            [finite_x_eq],
            [_energy_scale * _finite_energy],
            s=80,
            marker="D",
            color="#40464D",
            zorder=7,
            label="most probable finite-N macrostate",
        )

    _use_log_x = equilibrium_x < 0.05
    if _use_log_x:
        _free_energy_axis.set_xscale("log")
        _mu_axis.set_xscale("log")
    _free_energy_axis.set_xlim(focus_x[0], focus_x[-1])
    _free_energy_axis.set_xlabel("Defect fraction, x")
    _free_energy_axis.set_ylabel("Free energy per site (meV/site)")
    _free_energy_axis.set_title(
        "Free energy balances cost and entropy"
    )
    _free_energy_axis.grid(alpha=0.25)
    _free_energy_axis.legend(frameon=False, fontsize=11)

    _mu_axis.plot(
        focus_x,
        chemical_potential_values_ev,
        color="#4C7C86",
        lw=1.9,
    )
    _mu_axis.axhline(0.0, color="#666D73", ls="--", lw=1.4)
    _mu_axis.axvline(equilibrium_x, color="#C49345", ls=":", lw=2.0)
    _mu_axis.scatter(
        [equilibrium_x],
        [0.0],
        s=100,
        color="#C49345",
        edgecolor="#40464D",
        zorder=5,
    )
    _mu_axis.set_xlim(focus_x[0], focus_x[-1])
    _mu_axis.set_xlabel("Defect fraction, x")
    _mu_axis.set_ylabel(r"Chemical potential, $\mu_D$ (eV/defect)")
    _mu_axis.set_title(r"The same equilibrium has $\mu_D=0$")
    _mu_axis.grid(alpha=0.25)
    _energy_fig.tight_layout()
    plt.close(_energy_fig)
    free_energy_figure = _energy_fig
    return (free_energy_figure,)


@app.cell
def _(free_energy_figure):
    free_energy_figure
    return


@app.cell
def _(mo):
    mo.md(r"""
    Read the left panel first. The formation term
    penalizes defects, while $-Ts_{\rm config}$ lowers the free energy. Their
    solid-line sum has a marked minimum. Nothing in the calculation
    assumes that equilibrium concentration in advance.

    Differentiating that same Stirling-limit free energy gives

    \[
    \mu_D(x)=\frac{1}{N}\frac{\partial G}{\partial x}
    =\Delta g_f^0+k_BT\ln\frac{x}{1-x}.
    \]

    Here $\mu_D=0$ means equilibrium relative to the chosen perfect-lattice +
    defect reference reaction. It is the stationarity condition for exchanging
    one occupied site with one defect site, not an assertion that every
    absolute chemical potential is zero.

    The right panel shows that the minimum and the zero of \(\mu_D\) are the
    same state. Solving \(\mu_D=0\) gives

    \[
    \frac{x_{\rm eq}}{1-x_{\rm eq}}
    =\exp\left(-\frac{\Delta g_f^0}{k_BT}\right),
    \qquad
    \boxed{x_{\rm eq}=
    \frac{1}{1+\exp[\Delta g_f^0/(k_BT)]}}.
    \]
    """)
    return


@app.cell
def _(
    discrete_n,
    equilibrium_x,
    finite_x_eq,
    log_omega,
    np,
    plt,
    reduced_formation_energy,
    site_count,
    special,
):
    _driving_force = np.linspace(0.0, 16.0, 321)
    _thermodynamic_fraction = special.expit(-_driving_force)
    _dilute_fraction_curve = np.exp(-_driving_force)
    _finite_fraction_curve = np.empty_like(_driving_force)
    for _index, _energy_ratio in enumerate(_driving_force):
        _reduced_free_energy = (
            discrete_n * _energy_ratio - log_omega
        ) / site_count
        _finite_fraction_curve[_index] = (
            int(np.argmin(_reduced_free_energy)) / site_count
        )

    _comparison_fig, _comparison_axis = plt.subplots(
        figsize=(11.8, 5.8),
        dpi=120,
    )
    _comparison_axis.semilogy(
        _driving_force,
        _thermodynamic_fraction,
        color="#4C7C86",
        lw=1.9,
        label="Stirling / thermodynamic limit",
    )
    _comparison_axis.semilogy(
        _driving_force,
        _dilute_fraction_curve,
        color="#B8734A",
        lw=1.7,
        ls="--",
        label=r"dilute $\exp[-\Delta g_f^0/(k_BT)]$",
    )
    _positive_finite = _finite_fraction_curve > 0.0
    _comparison_axis.step(
        _driving_force[_positive_finite],
        _finite_fraction_curve[_positive_finite],
        where="mid",
        color="#40464D",
        lw=2.0,
        label=f"most probable finite-N macrostate (N = {site_count})",
    )
    _comparison_axis.axvspan(
        np.log(20.0),
        16.0,
        color="#5F8A6B",
        alpha=0.08,
        label="dilute relative error below 5%",
    )
    if 0.0 <= reduced_formation_energy <= 16.0:
        _comparison_axis.scatter(
            [reduced_formation_energy],
            [equilibrium_x],
            s=100,
            color="#C49345",
            edgecolor="#40464D",
            zorder=6,
            label="selected state",
        )
        if finite_x_eq > 0.0:
            _comparison_axis.scatter(
                [reduced_formation_energy],
                [finite_x_eq],
                s=70,
                marker="D",
                color="#40464D",
                zorder=6,
            )
    _comparison_axis.set(
        xlim=(0.0, 16.0),
        ylim=(1.0e-7, 1.2),
        xlabel=r"Formation driving force, $\Delta g_f^0/(k_BT)$",
        ylabel="Defect fraction",
        title="Exact, thermodynamic, and dilute limits",
    )
    _comparison_axis.grid(which="both", alpha=0.25)
    _comparison_axis.legend(
        loc="upper right",
        frameon=False,
        fontsize=11,
    )
    _comparison_fig.tight_layout()
    plt.close(_comparison_fig)
    approximation_figure = _comparison_fig
    return (approximation_figure,)


@app.cell
def _(approximation_figure, mo):
    mo.vstack([
        approximation_figure,
        mo.md(r"""
        The dilute exponential converges to the exact
        thermodynamic fraction only at low $x$. A finite lattice adds discrete
        composition steps, and its most-probable macrostate can reach $n=0$ even
        while the ensemble mean remains positive.
        """),
    ])
    return


@app.cell
def _(
    KB_EV_PER_K,
    chemical_potential_ev,
    delta_g0_ev,
    equilibrium_x,
    finite_x_eq,
    log_multiplicity,
    np,
    optimize,
    site_count,
    special,
    stirling_entropy_kb_per_site,
    stirling_free_energy_ev_per_site,
    temperature_k,
):
    _equilibrium_logit = float(
        np.log(equilibrium_x) - np.log1p(-equilibrium_x)
    )

    def _objective_in_logit(logit_value):
        _fraction = special.expit(logit_value)
        return float(
            stirling_free_energy_ev_per_site(
                _fraction,
                temperature_k,
                delta_g0_ev,
            )
        )

    _minimum_result = optimize.minimize_scalar(
        _objective_in_logit,
        bounds=(_equilibrium_logit - 8.0, _equilibrium_logit + 8.0),
        method="bounded",
        options={"xatol": 1.0e-11},
    )
    _minimum_x = float(special.expit(_minimum_result.x))

    _root_logit = optimize.brentq(
        lambda _logit: (
            delta_g0_ev + KB_EV_PER_K * temperature_k * _logit
        ),
        _equilibrium_logit - 8.0,
        _equilibrium_logit + 8.0,
        xtol=1.0e-13,
    )
    _root_x = float(special.expit(_root_logit))
    _mu_at_equilibrium = float(
        chemical_potential_ev(
            equilibrium_x,
            temperature_k,
            delta_g0_ev,
        )
    )

    _large_count = 1_000_000
    _large_fraction = 0.2
    _large_defects = int(_large_count * _large_fraction)
    _large_exact_entropy = float(
        log_multiplicity(_large_count, _large_defects) / _large_count
    )
    _large_stirling_entropy = float(
        stirling_entropy_kb_per_site(_large_fraction)
    )
    _large_entropy_relative_error = abs(
        _large_exact_entropy - _large_stirling_entropy
    ) / _large_stirling_entropy

    _dilute_reference_ratio = 12.0
    _dilute_reference = float(np.exp(-_dilute_reference_ratio))
    _exact_reference = float(special.expit(-_dilute_reference_ratio))
    _dilute_reference_relative_error = abs(
        _dilute_reference - _exact_reference
    ) / _exact_reference

    _selected_finite_rounding = abs(finite_x_eq - equilibrium_x)
    validation = {
        "minimum_x": _minimum_x,
        "minimum_error": abs(_minimum_x - equilibrium_x),
        "minimum_pass": bool(
            _minimum_result.success
            and abs(_minimum_result.x - _equilibrium_logit) < 2.0e-6
        ),
        "root_x": _root_x,
        "root_error": abs(_root_x - equilibrium_x),
        "mu_at_equilibrium": abs(_mu_at_equilibrium),
        "root_pass": bool(
            abs(_root_x - equilibrium_x) < 1.0e-11
            and abs(_mu_at_equilibrium) < 1.0e-11
        ),
        "finite_rounding": _selected_finite_rounding,
        "finite_rounding_pass": bool(
            _selected_finite_rounding <= 1.0 / site_count + 1.0e-14
        ),
        "large_entropy_relative_error": _large_entropy_relative_error,
        "large_entropy_pass": bool(_large_entropy_relative_error < 2.0e-5),
        "dilute_reference_relative_error": (
            _dilute_reference_relative_error
        ),
        "dilute_pass": bool(
            _dilute_reference_relative_error < 1.0e-5
        ),
        "all_finite": bool(
            np.all(
                np.isfinite(
                    [
                        _minimum_x,
                        _root_x,
                        _large_entropy_relative_error,
                        _dilute_reference_relative_error,
                    ]
                )
            )
        ),
    }
    return (validation,)


@app.cell
def _(mo, validation):
    def _mark(passed):
        return "PASS" if passed else "CHECK"

    _checks = mo.md(
        rf"""
        ## Physical consistency checks

        | status | physical statement | why it matters |
        |---:|---|---|
        | {_mark(validation['minimum_pass'])} | the minimum of \(G(x)\) occurs at the analytical \(x_{{\rm eq}}\) | equilibrium is the lowest-free-energy composition |
        | {_mark(validation['root_pass'])} | \(\mu_D=0\) at the same composition | the chemical-potential and free-energy views agree |
        | {_mark(validation['finite_rounding_pass'])} | the most probable finite lattice lies within one site of \(x_{{\rm eq}}\) | a finite lattice changes only in steps of \(1/N\) |
        | {_mark(validation['large_entropy_pass'])} | the exact entropy approaches Stirling's expression for large \(N\) | the smooth thermodynamic curve has the correct large-system limit |
        | {_mark(validation['dilute_pass'])} | the dilute exponential agrees at low defect fraction | the familiar approximation is used only in its proper limit |
        | {_mark(validation['all_finite'])} | all displayed quantities are well defined | every plotted point has physical meaning |

        These checks connect the different descriptions of the same
        equilibrium. They are not additional assumptions in the model.
        """
    )
    mo.accordion({"Physical consistency checks": _checks})
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## What to carry forward

    The exact finite lattice contains only
    \(x=0,1/N,\ldots,1\), while Stirling's approximation turns those states
    into a smooth thermodynamic curve. Its minimum gives

    \[
    x_{\rm eq}=\frac{1}{1+\exp[\Delta g_f^0/(k_BT)]}.
    \]

    Only when \(x_{\rm eq}\ll1\) may \(1-x_{\rm eq}\) be replaced by one, giving
    the familiar dilute exponential.

    The central lesson is simple: multiplicity creates configurational entropy,
    and that entropy competes with the energy cost of forming defects. The
    minimum of \(G\) and the zero of \(\mu_D\) are two views of the same
    equilibrium.

    **Model boundary.** Sites are equivalent and defects do not interact.
    Charged defects, coupled reactions, pressure, and electroneutrality begin in
    the next module.
    """)
    return


@app.cell
def _(mo):
    mo.md(r"""
    **Continue:** [Module 02 — Brouwer Diagram Explorer](https://qiyanglu.github.io/Solid-State-Ionics-Interactive/02-brouwer-sto/)
    """)
    return


if __name__ == "__main__":
    app.run()
