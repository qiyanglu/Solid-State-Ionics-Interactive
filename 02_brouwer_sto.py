import marimo

__generated_with = "0.23.14"
app = marimo.App(width="full")


@app.cell
def _():
    import marimo as mo
    import matplotlib.pyplot as plt
    import numpy as np
    from scipy import optimize

    return mo, np, optimize, plt


@app.cell
def _(mo):
    mo.md(r"""
    # Brouwer diagram explorer: weakly acceptor-doped SrTiO$_3$

    **Question.** Can Brouwer regimes and slopes emerge from equilibrium
    thermodynamics and electroneutrality, without putting a power law into the
    calculation?

    We retain only oxygen vacancies $V=[V_O^{\bullet\bullet}]$, electrons
    $n=[e']$, holes $p=[h^\bullet]$, and a fixed concentration
    $A=[A_{Ti}']$ of fully ionized acceptors. Concentrations are in cm$^{-3}$,
    temperature is in K, and $p_{O_2}$ is in bar.

    $$O_O^x \rightleftharpoons V_O^{\bullet\bullet}+2e'+\tfrac12O_2,$$

    $$K_{red}=Vn^2p_{O_2}^{1/2}, \qquad K_{eh}=np,$$

    $$\boxed{2V+p=A+n}. $$

    The last equation is **exact electroneutrality** for this four-species model:
    the two positive effective charges on each oxygen vacancy plus holes balance
    ionized acceptors plus electrons.

    The quantitative SrTiO$_3$ equilibrium constants are

    $$K_{red}(T)=6.616\times10^{68}
    \exp\!\left[-\frac{5.581}{k_BT}\right]
    \ \mathrm{cm^{-9}\,bar^{1/2}},$$

    $$K_{eh}(T)=7.49\times10^{45}
    \exp\!\left[-\frac{3.459}{k_BT}\right]
    \ \mathrm{cm^{-6}},$$

    with $k_B=8.617333262\times10^{-5}$ eV K$^{-1}$. Every plotted point comes
    from these two mass-action laws plus exact charge neutrality. Limiting
    Brouwer balances are used only afterward, as interpretation and checks.
    """)
    return


@app.cell
def _(np, optimize):
    KB_EV_PER_K = 8.617333262e-5

    def log_equilibrium_constants(temperature_k):
        """Return ln(K_red) and ln(K_eh) for temperature_k in K.

        K_red is referenced to concentrations in cm^-3 and oxygen pressure in
        bar, so its units are cm^-9 bar^(1/2). K_eh has units cm^-6.
        """
        if not np.isfinite(temperature_k) or temperature_k <= 0.0:
            raise ValueError("temperature_k must be positive and finite")
        log_k_red = (
            np.log(6.616)
            + 68.0 * np.log(10.0)
            - 5.581 / (KB_EV_PER_K * temperature_k)
        )
        log_k_eh = (
            np.log(7.49)
            + 45.0 * np.log(10.0)
            - 3.459 / (KB_EV_PER_K * temperature_k)
        )
        return log_k_red, log_k_eh

    def solve_equilibrium(oxygen_pressure_bar, temperature_k, acceptor_cm3):
        """Solve one positive equilibrium point in canonical units.

        Parameters are oxygen pressure in bar, temperature in K, and fixed
        acceptor concentration in cm^-3. Returned V, n, and p are in cm^-3.

        The mass-action laws eliminate V and p exactly. The remaining charge
        balance is solved for ln(n) with a bracketed root finder. Its positive
        side decreases monotonically while its negative side increases, so the
        positive solution is unique.
        """
        if not np.isfinite(oxygen_pressure_bar) or oxygen_pressure_bar <= 0.0:
            raise ValueError("oxygen_pressure_bar must be positive and finite")
        if not np.isfinite(acceptor_cm3) or acceptor_cm3 <= 0.0:
            raise ValueError("acceptor_cm3 must be positive and finite")

        log_k_red, log_k_eh = log_equilibrium_constants(temperature_k)
        log_pressure = np.log(oxygen_pressure_bar)
        log_acceptor = np.log(acceptor_cm3)

        def log_charge_balance(log_n):
            log_v = log_k_red - 2.0 * log_n - 0.5 * log_pressure
            log_p = log_k_eh - log_n
            log_positive = np.logaddexp(np.log(2.0) + log_v, log_p)
            log_negative = np.logaddexp(log_acceptor, log_n)
            return log_positive - log_negative

        log_n = optimize.brentq(
            log_charge_balance,
            np.log(1.0e-80),
            np.log(1.0e40),
            xtol=1.0e-12,
            rtol=4.0 * np.finfo(float).eps,
        )
        log_v = log_k_red - 2.0 * log_n - 0.5 * log_pressure
        log_p = log_k_eh - log_n
        return np.exp(log_v), np.exp(log_n), np.exp(log_p)

    def solve_curve(oxygen_pressures_bar, temperature_k, acceptor_cm3):
        """Solve a 1-D pressure grid; no Brouwer power law is used."""
        pressures = np.asarray(oxygen_pressures_bar, dtype=float)
        if pressures.ndim != 1 or pressures.size == 0:
            raise ValueError("oxygen_pressures_bar must be a non-empty 1-D array")
        if not np.all(np.isfinite(pressures)) or np.any(pressures <= 0.0):
            raise ValueError("all oxygen pressures must be positive and finite")

        vacancies = np.empty_like(pressures)
        electrons = np.empty_like(pressures)
        holes = np.empty_like(pressures)
        for index, pressure in enumerate(pressures):
            vacancies[index], electrons[index], holes[index] = solve_equilibrium(
                pressure, temperature_k, acceptor_cm3
            )
        return {
            "V": vacancies,
            "n": electrons,
            "p": holes,
            "A": np.full_like(pressures, acceptor_cm3),
        }

    return log_equilibrium_constants, solve_curve


@app.cell
def _(mo):
    temperature = mo.ui.slider(
        start=700,
        stop=1500,
        step=1,
        value=973,
        label="Temperature, T (K)",
        show_value=True,
    )
    log_acceptor = mo.ui.slider(
        start=13.0,
        stop=21.0,
        step=0.1,
        value=18.0,
        label=r"log$_{10}$(A / cm$^{-3}$)",
        show_value=True,
    )
    show_guides = mo.ui.checkbox(
        value=True,
        label="Show post-solution limiting-slope guides",
    )
    controls = mo.hstack(
        [temperature, log_acceptor, show_guides],
        justify="start",
        align="center",
        wrap=True,
        gap=2.0,
    )
    return controls, log_acceptor, show_guides, temperature


@app.cell
def _(controls, mo):
    mo.vstack(
        [
            mo.md("## Explore the exact equilibrium"),
            controls,
            mo.md(
                "The pressure window is fixed at "
                r"$10^{-25}\leq p_{O_2}\leq 1$ bar so temperatures and acceptor "
                "levels can be compared directly.\n\n"
                "**Suggested comparisons:** keep the 973 K, $10^{18}$ cm$^{-3}$ "
                "default for the acceptor plateau; raise $T$ toward 1500 K to "
                "bring intrinsic reduction into view; or try 850 K and "
                r"$\log_{10}A\approx13.3$ to reveal oxidizing compensation."
            ),
        ]
    )
    return


@app.cell
def _(
    log_acceptor,
    log_equilibrium_constants,
    np,
    solve_curve,
    temperature,
):
    oxygen_pressure = np.logspace(-25.0, 0.0, 301)
    acceptor = 10.0 ** log_acceptor.value
    concentrations = solve_curve(
        oxygen_pressure,
        float(temperature.value),
        acceptor,
    )
    log_pressure = np.log10(oxygen_pressure)
    local_slopes = {
        species: np.gradient(np.log10(values), log_pressure)
        for species, values in concentrations.items()
    }
    log_k_red, log_k_eh = log_equilibrium_constants(float(temperature.value))
    thermodynamic_state = {
        "log10_k_red": float(log_k_red / np.log(10.0)),
        "log10_k_eh": float(log_k_eh / np.log(10.0)),
        "log10_intrinsic_carrier": float(0.5 * log_k_eh / np.log(10.0)),
    }
    return (
        acceptor,
        concentrations,
        local_slopes,
        log_pressure,
        oxygen_pressure,
        thermodynamic_state,
    )


@app.cell
def _(acceptor, mo, temperature, thermodynamic_state):
    mo.md(
        rf"""
        ### Current thermodynamic state

        | Quantity | Value |
        |---|---:|
        | $T$ | {temperature.value:.0f} K |
        | $A$ | {acceptor:.3e} cm$^{{-3}}$ |
        | $\log_{{10}}(K_{{red}}/[\mathrm{{cm^{{-9}}\,bar^{{1/2}}}}])$ | {thermodynamic_state['log10_k_red']:.3f} |
        | $\log_{{10}}(K_{{eh}}/[\mathrm{{cm^{{-6}}}}])$ | {thermodynamic_state['log10_k_eh']:.3f} |
        | $\log_{{10}}(\sqrt{{K_{{eh}}}}/[\mathrm{{cm^{{-3}}}}])$ | {thermodynamic_state['log10_intrinsic_carrier']:.3f} |

        $\sqrt{{K_{{eh}}}}$ is the concentration at the electron–hole crossover
        $n=p$; it is a derived scale, not a constraint supplied to the solver.
        """
    )
    return


@app.cell
def _(acceptor, concentrations, log_pressure, np):
    vacancies = concentrations["V"]
    electrons = concentrations["n"]
    holes = concentrations["p"]
    dominance_tolerance = 0.10

    regime_masks = {
        "reducing": (
            (np.abs(2.0 * vacancies / electrons - 1.0) < dominance_tolerance)
            & (acceptor / electrons < dominance_tolerance)
            & (holes / electrons < dominance_tolerance)
        ),
        "acceptor": (
            (np.abs(2.0 * vacancies / acceptor - 1.0) < dominance_tolerance)
            & (electrons / acceptor < dominance_tolerance)
            & (holes / acceptor < dominance_tolerance)
        ),
        "oxidizing": (
            (np.abs(holes / acceptor - 1.0) < dominance_tolerance)
            & (electrons / acceptor < dominance_tolerance)
            & (2.0 * vacancies / acceptor < dominance_tolerance)
        ),
    }

    log_carrier_ratio = np.log10(electrons) - np.log10(holes)
    if log_carrier_ratio[0] * log_carrier_ratio[-1] <= 0.0:
        crossover_log_pressure = float(
            np.interp(0.0, log_carrier_ratio[::-1], log_pressure[::-1])
        )
    else:
        crossover_log_pressure = None
    return crossover_log_pressure, dominance_tolerance, regime_masks


@app.cell
def _(
    acceptor,
    concentrations,
    local_slopes,
    log_equilibrium_constants,
    log_pressure,
    np,
    oxygen_pressure,
    regime_masks,
    temperature,
):
    def median_slope(species, mask):
        if np.count_nonzero(mask) < 8:
            return None
        return float(np.median(local_slopes[species][mask]))

    def sampled_span(mask):
        indices = np.flatnonzero(mask)
        if indices.size < 8:
            return None
        return float(log_pressure[indices[0]]), float(log_pressure[indices[-1]])

    _log_k_red, _log_k_eh = log_equilibrium_constants(float(temperature.value))
    v_values = concentrations["V"]
    n_values = concentrations["n"]
    p_values = concentrations["p"]

    red_log_residual = (
        np.log(v_values)
        + 2.0 * np.log(n_values)
        + 0.5 * np.log(oxygen_pressure)
        - _log_k_red
    ) / np.log(10.0)
    eh_log_residual = (
        np.log(n_values) + np.log(p_values) - _log_k_eh
    ) / np.log(10.0)
    charge_residual = np.abs(
        2.0 * v_values + p_values - acceptor - n_values
    ) / (2.0 * v_values + p_values + acceptor + n_values)

    sanity = {
        "positive": bool(
            all(np.all(values > 0.0) for values in concentrations.values())
        ),
        "finite": bool(
            all(np.all(np.isfinite(values)) for values in concentrations.values())
        ),
        "red_residual": float(np.max(np.abs(red_log_residual))),
        "eh_residual": float(np.max(np.abs(eh_log_residual))),
        "charge_residual": float(np.max(charge_residual)),
    }
    slope_summary = {
        "reducing": {
            "span": sampled_span(regime_masks["reducing"]),
            "V": median_slope("V", regime_masks["reducing"]),
            "n": median_slope("n", regime_masks["reducing"]),
            "p": median_slope("p", regime_masks["reducing"]),
        },
        "acceptor": {
            "span": sampled_span(regime_masks["acceptor"]),
            "V": median_slope("V", regime_masks["acceptor"]),
            "n": median_slope("n", regime_masks["acceptor"]),
            "p": median_slope("p", regime_masks["acceptor"]),
        },
        "oxidizing": {
            "span": sampled_span(regime_masks["oxidizing"]),
            "V": median_slope("V", regime_masks["oxidizing"]),
            "n": median_slope("n", regime_masks["oxidizing"]),
            "p": median_slope("p", regime_masks["oxidizing"]),
        },
    }
    return sanity, slope_summary


@app.cell
def _(log_pressure, np):
    def add_slope_guide(axis, y_values, mask, slope, label):
        """Draw a short comparison segment without changing calculated data."""
        indices = np.flatnonzero(mask)
        if indices.size < 8:
            return
        center_index = int(indices[indices.size // 2])
        available_span = log_pressure[indices[-1]] - log_pressure[indices[0]]
        half_width = min(1.25, max(0.35, 0.22 * available_span))
        x_center = log_pressure[center_index]
        x_guide = np.array([x_center - half_width, x_center + half_width])
        y_center = np.log10(y_values[center_index])
        y_guide = y_center + slope * (x_guide - x_center)
        axis.plot(
            x_guide,
            y_guide,
            color="black",
            lw=1.7,
            ls=(0, (4, 2)),
            zorder=6,
        )
        axis.annotate(
            label,
            xy=(x_guide[-1], y_guide[-1]),
            xytext=(5, 3),
            textcoords="offset points",
            fontsize=10,
            color="black",
        )

    return (add_slope_guide,)


@app.cell
def _(
    acceptor,
    add_slope_guide,
    concentrations,
    crossover_log_pressure,
    log_pressure,
    np,
    plt,
    regime_masks,
    show_guides,
    temperature,
):
    plt.rcParams.update(
        {
            "font.size": 12,
            "axes.titlesize": 16,
            "axes.labelsize": 14,
            "legend.fontsize": 11,
            "lines.linewidth": 2.8,
        }
    )
    figure, axis = plt.subplots(figsize=(12.0, 7.0), dpi=120)

    if crossover_log_pressure is not None:
        axis.axvspan(
            log_pressure[0],
            crossover_log_pressure,
            color="#4477AA",
            alpha=0.035,
        )
        axis.axvspan(
            crossover_log_pressure,
            log_pressure[-1],
            color="#CC6677",
            alpha=0.035,
        )
    elif concentrations["n"][0] > concentrations["p"][0]:
        axis.axvspan(log_pressure[0], log_pressure[-1], color="#4477AA", alpha=0.035)
    else:
        axis.axvspan(log_pressure[0], log_pressure[-1], color="#CC6677", alpha=0.035)

    acceptor_indices = np.flatnonzero(regime_masks["acceptor"])
    if acceptor_indices.size:
        axis.axvspan(
            log_pressure[acceptor_indices[0]],
            log_pressure[acceptor_indices[-1]],
            color="#228833",
            alpha=0.075,
            label=r"strict $2V\approx A$ window",
        )

    axis.plot(
        log_pressure,
        np.log10(concentrations["V"]),
        color="#EE7733",
        label=r"$V_O^{\bullet\bullet}$ ($V$)",
    )
    axis.plot(
        log_pressure,
        np.log10(concentrations["n"]),
        color="#0077BB",
        label=r"$e'$ ($n$)",
    )
    axis.plot(
        log_pressure,
        np.log10(concentrations["p"]),
        color="#CC3311",
        label=r"$h^\bullet$ ($p$)",
    )
    axis.plot(
        log_pressure,
        np.log10(concentrations["A"]),
        color="#555555",
        ls="--",
        lw=2.2,
        label=r"fixed $A_{Ti}'$ ($A$)",
    )

    if crossover_log_pressure is not None:
        axis.axvline(crossover_log_pressure, color="#666666", ls=":", lw=1.3)
        axis.text(
            crossover_log_pressure,
            0.04,
            r" $n=p$",
            transform=axis.get_xaxis_transform(),
            rotation=90,
            va="bottom",
            ha="left",
            fontsize=10,
            color="#555555",
        )
        axis.text(
            0.02,
            0.96,
            "electron-rich (reducing) side",
            transform=axis.transAxes,
            color="#225588",
            va="top",
            fontsize=11,
        )
        axis.text(
            0.98,
            0.96,
            "hole-rich (oxidizing) side",
            transform=axis.transAxes,
            color="#992233",
            ha="right",
            va="top",
            fontsize=11,
        )

    if show_guides.value:
        add_slope_guide(
            axis,
            concentrations["n"],
            regime_masks["reducing"],
            -1.0 / 6.0,
            r"$n:\;-1/6$",
        )
        add_slope_guide(
            axis,
            concentrations["n"],
            regime_masks["acceptor"],
            -1.0 / 4.0,
            r"$n:\;-1/4$",
        )
        add_slope_guide(
            axis,
            concentrations["p"],
            regime_masks["acceptor"],
            1.0 / 4.0,
            r"$p:\;+1/4$",
        )
        add_slope_guide(
            axis,
            concentrations["V"],
            regime_masks["acceptor"],
            0.0,
            r"$V:\;0$",
        )
        add_slope_guide(
            axis,
            concentrations["V"],
            regime_masks["oxidizing"],
            -1.0 / 2.0,
            r"$V:\;-1/2$",
        )
        add_slope_guide(
            axis,
            concentrations["p"],
            regime_masks["oxidizing"],
            0.0,
            r"$p:\;0$",
        )

    plotted_logs = np.concatenate(
        [np.log10(values) for values in concentrations.values()]
    )
    axis.set_ylim(
        np.floor(plotted_logs.min()) - 0.5,
        np.ceil(plotted_logs.max()) + 0.7,
    )
    axis.set_xlim(-25.0, 0.0)
    axis.set_xlabel(r"$\log_{10}(p_{O_2}\,/\,\mathrm{bar})$")
    axis.set_ylabel(r"$\log_{10}(c\,/\,\mathrm{cm}^{-3})$")
    axis.set_title(
        rf"Exact defect equilibrium in SrTiO$_3$:  T = {temperature.value:.0f} K,  "
        rf"A = $10^{{{np.log10(acceptor):.1f}}}$ cm$^{{-3}}$"
    )
    axis.xaxis.set_major_locator(plt.MultipleLocator(5.0))
    axis.xaxis.set_minor_locator(plt.MultipleLocator(1.0))
    axis.grid(which="major", color="#B0B0B0", alpha=0.40)
    axis.grid(which="minor", color="#D0D0D0", alpha=0.18)
    axis.legend(loc="center left", bbox_to_anchor=(1.01, 0.5), frameon=False)
    figure.tight_layout()
    plt.close(figure)
    return (figure,)


@app.cell
def _(figure):
    figure
    return


@app.cell
def _(mo):
    mo.md(r"""
    **How to read the shading.** The very light blue/red tint identifies the
    majority carrier ($n>p$ or $p>n$); it is **not** a charge-compensation
    approximation. The green band marks points that pass the strict
    acceptor-compensation test $2V\approx A$. The vertical dotted line is the
    derived crossover $n=p$. Short black dashed segments are reference slopes
    added after solving; the colored curves are always exact numerical results.
    """)
    return


@app.cell
def _(
    crossover_log_pressure,
    dominance_tolerance,
    mo,
    sanity,
    slope_summary,
):
    def residual_mark(value, tolerance):
        return "✓" if value < tolerance else "⚠"

    def mass_action_result(value, tolerance=1e-10, display_floor=1e-12):
        """Report roundoff-level log residuals consistently for students."""
        mark = residual_mark(value, tolerance)
        if value < display_floor:
            return rf"{mark} numerical zero ($< {display_floor:.0e}$)"
        return f"{mark} {value:.2e}"

    def slope_cell(value, expected, tolerance=0.035):
        if value is None:
            return "not sampled"
        mark = "✓" if abs(value - expected) < tolerance else "△"
        return f"{mark} {value:+.3f}"

    def span_cell(span):
        if span is None:
            return "not sampled"
        return f"{span[0]:.2f} to {span[1]:.2f}"

    positivity_mark = "✓" if sanity["positive"] else "⚠"
    finite_mark = "✓" if sanity["finite"] else "⚠"
    if crossover_log_pressure is None:
        crossover_text = "outside the plotted pressure window"
    else:
        crossover_text = f"$p_{{O_2}}={10.0 ** crossover_log_pressure:.3e}$ bar"

    mo.md(
        rf"""
        ## Numerical sanity checks

        | Check | Result |
        |---|---:|
        | All concentrations positive | {positivity_mark} |
        | All concentrations finite | {finite_mark} |
        | max $\lvert\log_{{10}}(Vn^2p_{{O_2}}^{{1/2}}/K_{{red}})\rvert$ (target: 0) | {mass_action_result(sanity['red_residual'])} |
        | max $\lvert\log_{{10}}(np/K_{{eh}})\rvert$ (target: 0) | {mass_action_result(sanity['eh_residual'])} |
        | max scaled charge residual | {residual_mark(sanity['charge_residual'], 1e-12)} {sanity['charge_residual']:.2e} |

        The two mass-action targets are zero because each logarithm contains an
        equilibrium ratio whose unlogged target is one. Values below $10^{{-12}}$
        are reported uniformly as **numerical zero**: an exact floating-point
        cancellation and a roundoff-level value such as $10^{{-14}}$ carry the
        same physical meaning here.

        **Why check these?** Positivity and finiteness ensure that every plotted
        point is a physical concentration. The two mass-action rows verify the
        defect reactions, and the charge row verifies electroneutrality. The
        slope rows then ask whether a true limiting balance is actually visible
        in the selected pressure window; they do not supply slopes to the solver.

        The electron–hole crossover is {crossover_text}.

        ### Regime coverage and measured slopes

        A limiting balance is sampled only where every neglected charge term is
        below {100.0 * dominance_tolerance:.0f}% of the retained scale and the
        proposed balance closes within the same tolerance. At least eight grid
        points are required. Pressure spans are reported as
        $\log_{{10}}(p_{{O_2}}/\mathrm{{bar}})$.

        | Limiting balance | Sampled span |
        |---|---:|
        | $2V\simeq n$ (intrinsic reduction) | {span_cell(slope_summary['reducing']['span'])} |
        | $2V\simeq A$ (acceptor compensation) | {span_cell(slope_summary['acceptor']['span'])} |
        | $p\simeq A$ (oxidizing compensation) | {span_cell(slope_summary['oxidizing']['span'])} |

        | Balance | Species | Expected | Median numerical slope |
        |---|---:|---:|---:|
        | $2V\simeq n$ | $V$ | $-1/6$ | {slope_cell(slope_summary['reducing']['V'], -1/6)} |
        | $2V\simeq n$ | $n$ | $-1/6$ | {slope_cell(slope_summary['reducing']['n'], -1/6)} |
        | $2V\simeq n$ plus $np=K_{{eh}}$ | $p$ | $+1/6$ | {slope_cell(slope_summary['reducing']['p'], 1/6)} |
        | $2V\simeq A$ | $V$ | $0$ | {slope_cell(slope_summary['acceptor']['V'], 0)} |
        | $2V\simeq A$ | $n$ | $-1/4$ | {slope_cell(slope_summary['acceptor']['n'], -1/4)} |
        | $2V\simeq A$ plus $np=K_{{eh}}$ | $p$ | $+1/4$ | {slope_cell(slope_summary['acceptor']['p'], 1/4)} |
        | $p\simeq A$ | $p$ | $0$ | {slope_cell(slope_summary['oxidizing']['p'], 0)} |
        | $p\simeq A$ plus $np=K_{{eh}}$ | $n$ | $0$ | {slope_cell(slope_summary['oxidizing']['n'], 0)} |
        | $p\simeq A$ | $V$ | $-1/2$ | {slope_cell(slope_summary['oxidizing']['V'], -1/2)} |

        ✓ means the measured slope is within 0.035 of its asymptotic value; △
        means the dominance test is met but convergence to the ideal slope is
        not yet tight. “Not sampled” is a physical result, not a failed test:
        the selected $T$, $A$, and pressure window decide which limits exist on
        screen.
        """
    )
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## Where the power laws come from

    The numerical solver above does not know any of these approximations. They
    are derived afterward by deciding which terms dominate exact charge
    neutrality:

    - **Strong reduction:** if $2V\simeq n$, substitution into
      $K_{red}=Vn^2p_{O_2}^{1/2}$ gives
      $V,n\propto p_{O_2}^{-1/6}$. Since $np=K_{eh}$,
      $p\propto p_{O_2}^{+1/6}$.
    - **Acceptor-compensated plateau:** if $2V\simeq A$, then $V$ has slope zero,
      $n\propto p_{O_2}^{-1/4}$, and
      $p=K_{eh}/n\propto p_{O_2}^{+1/4}$.
    - **Strong oxidation:** if $p\simeq A$, then $p$ and
      $n=K_{eh}/p$ approach slope-zero plateaus, while
      $V\propto p_{O_2}^{-1/2}$.

    The bends between these limits are not patched together. They appear because
    all four terms in $2V+p=A+n$ are retained while the balance changes smoothly.
    Likewise, the point $n=p$ separates electron-rich and hole-rich behavior but
    is not itself a charge-compensation approximation.

    **Assumptions.** Defects and carriers are treated as dilute ideal species;
    activities are represented by concentrations in the stated units; oxygen
    activity is represented by pressure in bar; acceptors are fixed and fully
    ionized. Cation defects, defect complexes, band-density limits, non-ideal
    activities, and changes in site density are deliberately outside this first
    teaching model.

    **Numerical method.** At each pressure, the mass-action equations are used
    exactly to write $V=K_{red}/(n^2p_{O_2}^{1/2})$ and $p=K_{eh}/n$. A bracketed
    root solve in $\ln n$ then enforces electroneutrality. Working in log space
    preserves positivity and avoids overflow across many decades.
    """)
    return


if __name__ == "__main__":
    app.run()
