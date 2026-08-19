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
            "lines.solid_capstyle": "round",
            "figure.dpi": 115,
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
    # Impedance Spectroscopy: From a Semicircle to a Transmission Line

    **Guiding question.** What physical process is able to respond at each
    frequency, and how does that response appear in the complex impedance?

    Electrochemical impedance spectroscopy (EIS) applies a small sinusoidal
    voltage and measures the sinusoidal current. Sweeping the frequency separates
    processes by their characteristic times: a fast bulk response, a slower
    interface, or still slower chemical diffusion. The shape is therefore a
    compressed map of **dynamics**, not a collection of arcs to name by eye.

    This notebook follows three steps:

    1. build the complex-number language of EIS from voltage and current waves;
    2. see the Warburg response emerge from one-dimensional diffusion;
    3. join ionic and electronic transport with chemical storage in a
       two-rail transmission line model (TLM).

    The convention throughout is

    $$e^{i\omega t},\qquad Z=\frac{\widehat V}{\widehat I}=Z'+iZ'',$$

    so a capacitor has $Z''<0$. A Nyquist plot displays $Z'$ horizontally and
    **$-Z''$ vertically**, putting passive capacitive responses above the axis.
    """)
    return


@app.cell
def _(np):
    def _positive_07(name, value):
        number = float(value)
        if not np.isfinite(number) or number <= 0.0:
            raise ValueError(f"{name} must be positive and finite")
        return number

    def waveform_data_07(frequency_hz, current_lead_deg, sample_count=500):
        """Return two cycles of normalized voltage and current."""
        frequency = _positive_07("frequency_hz", frequency_hz)
        lead_rad = np.deg2rad(float(current_lead_deg))
        phase = np.linspace(0.0, 4.0 * np.pi, int(sample_count))
        time_s = phase / (2.0 * np.pi * frequency)
        voltage = np.cos(phase)
        current = 0.72 * np.cos(phase + lead_rad)
        impedance_phase_deg = -float(current_lead_deg)
        return time_s, voltage, current, impedance_phase_deg

    def rc_impedance_07(frequency_hz, series_resistance_ohm, branches):
        """Series resistance followed by any number of parallel R-C branches."""
        frequencies = np.asarray(frequency_hz, dtype=float)
        omega = 2.0 * np.pi * frequencies
        impedance = np.full(frequencies.shape, complex(series_resistance_ohm), dtype=complex)
        for resistance_ohm, capacitance_f in branches:
            resistance = _positive_07("resistance_ohm", resistance_ohm)
            capacitance = _positive_07("capacitance_f", capacitance_f)
            impedance += resistance / (1.0 + 1j * omega * resistance * capacitance)
        return impedance

    def _stable_tanh_positive_real_07(value):
        """Evaluate tanh without overflow when the real part is large and positive."""
        values = np.asarray(value, dtype=complex)
        flat_values = values.reshape(-1)
        flat_result = np.empty(flat_values.shape, dtype=complex)
        large = flat_values.real > 20.0
        flat_result[~large] = np.tanh(flat_values[~large])
        decay = np.exp(-2.0 * flat_values[large])
        flat_result[large] = (1.0 - decay) / (1.0 + decay)
        return flat_result.reshape(values.shape)


    def warburg_impedance_07(reduced_omega, boundary):
        """Dimensionless semi-infinite or finite-length Warburg impedance."""
        omega_tilde = np.asarray(reduced_omega, dtype=float)
        if np.any(omega_tilde <= 0.0):
            raise ValueError("reduced_omega must be positive")
        q_value = np.sqrt(1j * omega_tilde)
        tanh_q = _stable_tanh_positive_real_07(q_value)
        if boundary == "semi-infinite":
            return 1.0 / q_value
        if boundary == "open":
            return tanh_q / q_value
        if boundary == "blocked":
            return 1.0 / (q_value * tanh_q)
        raise ValueError(f"unknown Warburg boundary: {boundary}")

    def warburg_profile_07(reduced_omega, position, boundary):
        """Concentration phasor normalized to unit amplitude at x = 0."""
        omega_tilde = _positive_07("reduced_omega", reduced_omega)
        xi = np.asarray(position, dtype=float)
        q_value = np.sqrt(1j * omega_tilde)
        if boundary == "open":
            return np.sinh(q_value * (1.0 - xi)) / np.sinh(q_value)
        if boundary == "blocked":
            return np.cosh(q_value * (1.0 - xi)) / np.cosh(q_value)
        if boundary == "semi-infinite":
            return np.exp(-q_value * xi)
        raise ValueError(f"unknown Warburg boundary: {boundary}")

    def warburg_scales_07(
        length_um,
        diffusivity_cm2_per_s,
        concentration_mol_per_m3,
        area_cm2,
        temperature_k,
        charge_number,
    ):
        """Return the diffusion frequency and dilute finite-length resistance."""
        length_m = _positive_07("length_um", length_um) * 1.0e-6
        diffusivity_m2_per_s = _positive_07(
            "diffusivity_cm2_per_s", diffusivity_cm2_per_s
        ) * 1.0e-4
        concentration = _positive_07("concentration_mol_per_m3", concentration_mol_per_m3)
        area_m2 = _positive_07("area_cm2", area_cm2) * 1.0e-4
        temperature = _positive_07("temperature_k", temperature_k)
        electrons = _positive_07("charge_number", charge_number)
        gas_constant = 8.314462618
        faraday = 96485.33212
        tau_diffusion_s = length_m**2 / diffusivity_m2_per_s
        frequency_diffusion_hz = 1.0 / (2.0 * np.pi * tau_diffusion_s)
        resistance_diffusion_ohm = (
            gas_constant
            * temperature
            * length_m
            / (
                electrons**2
                * faraday**2
                * concentration
                * area_m2
                * diffusivity_m2_per_s
            )
        )
        return {
            "length_m": length_m,
            "diffusivity_m2_per_s": diffusivity_m2_per_s,
            "tau_diffusion_s": tau_diffusion_s,
            "frequency_diffusion_hz": frequency_diffusion_hz,
            "resistance_diffusion_ohm": resistance_diffusion_ohm,
        }

    def _tlm_dynamic_factors_07(reduced_omega):
        """Stable endpoint factors for the continuous dual-rail TLM."""
        k_value = np.sqrt(1j * float(reduced_omega))
        if abs(k_value) < 1.0e-5:
            k_squared = k_value**2
            mean_factor = k_squared / 2.0 - k_squared**2 / 24.0
            difference_factor = 1.0 + k_squared / 12.0 - k_squared**2 / 720.0
        else:
            mean_factor = k_value * np.tanh(k_value / 2.0)
            difference_factor = (k_value / 2.0) / np.tanh(k_value / 2.0)
        return k_value, mean_factor, difference_factor

    def tlm_parameters_07(parallel_resistance_ohm, chemical_capacitance_f, conductivity_ratio):
        """Map classroom controls to the two rail resistances and chemical time."""
        r_parallel = _positive_07("parallel_resistance_ohm", parallel_resistance_ohm)
        capacitance = _positive_07("chemical_capacitance_f", chemical_capacitance_f)
        ratio = _positive_07("conductivity_ratio", conductivity_ratio)
        t_e = ratio / (1.0 + ratio)
        t_i = 1.0 / (1.0 + ratio)
        r_e = r_parallel / t_e
        r_i = r_parallel / t_i
        tau_chemical = (r_e + r_i) * capacitance
        frequency_chemical = 1.0 / (2.0 * np.pi * tau_chemical)
        return {
            "R_parallel_ohm": r_parallel,
            "C_chemical_f": capacitance,
            "sigma_e_over_sigma_i": ratio,
            "t_e": t_e,
            "t_i": t_i,
            "R_e_ohm": r_e,
            "R_i_ohm": r_i,
            "R_sum_ohm": r_e + r_i,
            "tau_chemical_s": tau_chemical,
            "frequency_chemical_hz": frequency_chemical,
        }

    def _tlm_endpoint_maps_07(reduced_omega, parameters):
        """Map modal coefficients [U, Q, m, delta] to both endpoint states."""
        _, factor_d, factor_h = _tlm_dynamic_factors_07(reduced_omega)
        r_e = parameters["R_e_ohm"]
        r_i = parameters["R_i_ohm"]
        r_sum = parameters["R_sum_ohm"]
        r_parallel = parameters["R_parallel_ohm"]
        alpha_e = r_e / r_sum
        alpha_i = r_i / r_sum

        left = np.array(
            [
                [1.0, 0.0, alpha_e, 0.5 * alpha_e],
                [1.0, 0.0, -alpha_i, -0.5 * alpha_i],
                [0.0, alpha_i / r_parallel, factor_d / r_sum, factor_h / r_sum],
                [0.0, alpha_e / r_parallel, -factor_d / r_sum, -factor_h / r_sum],
            ],
            dtype=complex,
        )
        right = np.array(
            [
                [1.0, -1.0, alpha_e, -0.5 * alpha_e],
                [1.0, -1.0, -alpha_i, 0.5 * alpha_i],
                [0.0, alpha_i / r_parallel, -factor_d / r_sum, factor_h / r_sum],
                [0.0, alpha_e / r_parallel, factor_d / r_sum, -factor_h / r_sum],
            ],
            dtype=complex,
        )
        return left, right

    def tlm_solution_07(reduced_omega, parameters, contact_case):
        """Solve one continuous TLM frequency for selected ideal contacts."""
        omega_tilde = _positive_07("reduced_omega", reduced_omega)
        left, right = _tlm_endpoint_maps_07(omega_tilde, parameters)
        if contact_case == "electron contacts, ions blocked":
            boundary_matrix = np.vstack((left[0], left[3], right[0], right[3]))
        elif contact_case == "cross-selective contacts":
            boundary_matrix = np.vstack((left[0], left[3], right[2], right[1]))
        elif contact_case == "both carriers reversible":
            boundary_matrix = np.vstack((left[0], left[1], right[0], right[1]))
        else:
            raise ValueError(f"unknown TLM contact case: {contact_case}")
        if contact_case == "both carriers reversible":
            target = np.array([1.0, 1.0, 0.0, 0.0], dtype=complex)
        else:
            target = np.array([1.0, 0.0, 0.0, 0.0], dtype=complex)
        coefficients = np.linalg.solve(boundary_matrix, target)
        total_current = coefficients[1] / parameters["R_parallel_ohm"]
        impedance = 1.0 / total_current
        residual = np.max(np.abs(boundary_matrix @ coefficients - target))
        return {
            "coefficients": coefficients,
            "impedance_ohm": impedance,
            "boundary_residual": residual,
            "left_state": left @ coefficients,
            "right_state": right @ coefficients,
        }

    def tlm_spectrum_07(reduced_omega, parameters, contact_case):
        omega_array = np.asarray(reduced_omega, dtype=float)
        return np.array(
            [
                tlm_solution_07(value, parameters, contact_case)["impedance_ohm"]
                for value in omega_array
            ]
        )

    def _sinhc_complex_07(value):
        values = np.asarray(value, dtype=complex)
        result = np.empty(values.shape, dtype=complex)
        small = np.abs(values) < 1.0e-6
        result[small] = 1.0 + values[small] ** 2 / 6.0 + values[small] ** 4 / 120.0
        result[~small] = np.sinh(values[~small]) / values[~small]
        return result

    def tlm_profile_07(position, reduced_omega, parameters, contact_case):
        """Return continuous voltage-equivalent potentials and rail currents."""
        xi = np.asarray(position, dtype=float)
        solution = tlm_solution_07(reduced_omega, parameters, contact_case)
        u_value, q_value, mean_value, difference_value = solution["coefficients"]
        k_value, _, _ = _tlm_dynamic_factors_07(reduced_omega)
        zeta = xi - 0.5
        r_e = parameters["R_e_ohm"]
        r_i = parameters["R_i_ohm"]
        r_sum = parameters["R_sum_ohm"]
        r_parallel = parameters["R_parallel_ohm"]
        alpha_e = r_e / r_sum
        alpha_i = r_i / r_sum

        common_voltage = u_value - xi * q_value
        chemical_voltage = (
            mean_value * np.cosh(k_value * zeta) / np.cosh(k_value / 2.0)
            - difference_value
            * zeta
            * _sinhc_complex_07(k_value * zeta)
            / _sinhc_complex_07(np.array([k_value / 2.0]))[0]
        )
        weighted_current = (
            -mean_value
            * k_value**2
            * zeta
            * _sinhc_complex_07(k_value * zeta)
            / np.cosh(k_value / 2.0)
            + difference_value
            * np.cosh(k_value * zeta)
            / _sinhc_complex_07(np.array([k_value / 2.0]))[0]
        ) / r_sum
        total_current = np.full(xi.shape, q_value / r_parallel, dtype=complex)
        u_e = common_voltage + alpha_e * chemical_voltage
        u_i = common_voltage - alpha_i * chemical_voltage
        i_e = alpha_i * total_current + weighted_current
        i_i = alpha_e * total_current - weighted_current
        return {
            "position": xi,
            "u_e": u_e,
            "u_i": u_i,
            "u_chemical": chemical_voltage,
            "I_e": i_e,
            "I_i": i_i,
            "I_total": total_current,
            "solution": solution,
        }

    return (
        rc_impedance_07,
        tlm_parameters_07,
        tlm_profile_07,
        tlm_solution_07,
        tlm_spectrum_07,
        warburg_impedance_07,
        warburg_profile_07,
        warburg_scales_07,
        waveform_data_07,
    )


@app.cell
def _(mo):
    mo.md(r"""
    ## 1. A sinusoid becomes one complex number

    In the linear small-signal limit, a sinusoidal input produces a sinusoidal
    response at the same frequency:

    $$
    \Delta V(t)=\Re\!\left[\widehat V e^{i\omega t}\right],\qquad
    \Delta I(t)=\Re\!\left[\widehat I e^{i\omega t}\right].
    $$

    The complex ratio $Z=\widehat V/\widehat I$ stores two measurements at once:
    its magnitude is the voltage/current amplitude ratio and its phase is the
    voltage phase minus the current phase. The perturbation must be small enough
    that the material is approximately linear around its operating point.
    """)
    return


@app.cell
def _(mo):
    waveform_frequency_07 = mo.ui.slider(
        start=-1.0, stop=2.0, step=0.1, value=0.5, label=r"$\log_{10}(f/\mathrm{Hz})$"
    )
    waveform_lead_07 = mo.ui.slider(
        start=0, stop=90, step=5, value=45, label="current lead (degrees)"
    )
    mo.hstack([waveform_frequency_07, waveform_lead_07], justify="start", gap=2.0)
    return waveform_frequency_07, waveform_lead_07


@app.cell
def _(waveform_data_07, waveform_frequency_07, waveform_lead_07):
    waveform_result_07 = waveform_data_07(
        10.0 ** waveform_frequency_07.value, waveform_lead_07.value
    )
    return (waveform_result_07,)


@app.cell
def _(np, plt, waveform_result_07):
    _time_s, _voltage, _current, _z_phase = waveform_result_07
    _period_s = (_time_s[-1] - _time_s[0]) / 2.0
    _time_cycles = _time_s / _period_s
    _figure, (_axis_wave, _axis_phasor) = plt.subplots(
        1, 2, figsize=(12.5, 4.4), constrained_layout=True
    )
    _axis_wave.plot(_time_cycles, _voltage, lw=2.7, color="#4C7C86", label=r"$\Delta V/V_a$")
    _axis_wave.plot(_time_cycles, _current, lw=2.7, color="#B8734A", label=r"$\Delta I/I_a$")
    _axis_wave.axhline(0.0, color="#73808C", lw=0.9)
    _axis_wave.set(xlabel="time / period", ylabel="normalized signal", title="The current leads the voltage")
    _axis_wave.grid(True)
    _axis_wave.legend(loc="upper right")

    _lead_rad = np.deg2rad(-_z_phase)
    _axis_phasor.arrow(0.0, 0.0, 1.0, 0.0, width=0.012, color="#4C7C86", length_includes_head=True)
    _axis_phasor.arrow(
        0.0,
        0.0,
        0.78 * np.cos(_lead_rad),
        0.78 * np.sin(_lead_rad),
        width=0.012,
        color="#B8734A",
        length_includes_head=True,
    )
    _axis_phasor.text(1.03, 0.0, r"$\widehat V$", color="#4C7C86", va="center")
    _axis_phasor.text(
        0.82 * np.cos(_lead_rad),
        0.82 * np.sin(_lead_rad),
        r"$\widehat I$",
        color="#B8734A",
        ha="left",
        va="bottom",
    )
    _axis_phasor.text(
        0.03,
        0.04,
        rf"$\arg Z={_z_phase:.0f}^\circ$",
        transform=_axis_phasor.transAxes,
        color="#526173",
    )
    _axis_phasor.set(xlim=(-0.15, 1.25), ylim=(-0.15, 1.05), aspect="equal", title="Phasors rotate together")
    _axis_phasor.set_xlabel("real component")
    _axis_phasor.set_ylabel("imaginary component")
    _axis_phasor.grid(True)
    _figure
    return


@app.cell
def _(mo):
    mo.md(r"""
    ### Three ideal elements are enough to read the first Nyquist plot

    | element | impedance | Nyquist signature |
    |---|---:|---|
    | resistor | $Z_R=R$ | point on the real axis |
    | capacitor | $Z_C=1/(i\omega C)$ | vertical capacitive line |
    | parallel $R\parallel C$ | $Z=R/(1+i\omega RC)$ | semicircle |

    For a parallel $R\parallel C$ element, the relaxation time is
    $\tau=RC$ and the top of its ideal semicircle occurs at
    $\omega\tau=1$. Two processes are resolved as two arcs only when their time
    scales are sufficiently separated. An arc is therefore evidence of a
    relaxation time, not a unique microscopic label.

    The capacitor expression above uses capacitance $C$. This corrects an
    obvious $R/C$ symbol slip in the lecture graphic while preserving its
    intended convention.
    """)
    return


@app.cell
def _(mo):
    rc_series_07 = mo.ui.slider(start=0, stop=200, step=5, value=25, label=r"$R_s$ ($\Omega$)")
    rc_resistance_1_07 = mo.ui.slider(
        start=100, stop=1500, step=50, value=600, label=r"$R_1$ ($\Omega$)"
    )
    rc_log_tau_1_07 = mo.ui.slider(
        start=-4.0, stop=1.0, step=0.1, value=-1.7, label=r"$\log_{10}(\tau_1/\mathrm{s})$"
    )
    rc_show_second_07 = mo.ui.checkbox(value=True, label="show a second relaxation")
    rc_resistance_ratio_07 = mo.ui.slider(
        start=0.25, stop=2.0, step=0.05, value=1.4, label=r"$R_2/R_1$"
    )
    rc_log_separation_07 = mo.ui.slider(
        start=0.0,
        stop=4.0,
        step=0.1,
        value=2.0,
        label=r"$\log_{10}(\tau_2/\tau_1)$",
    )
    mo.vstack(
        [
            mo.hstack([rc_series_07, rc_resistance_1_07, rc_log_tau_1_07], justify="start", gap=1.5),
            mo.hstack(
                [rc_show_second_07, rc_resistance_ratio_07, rc_log_separation_07],
                justify="start",
                gap=1.5,
            ),
        ]
    )
    return (
        rc_log_separation_07,
        rc_log_tau_1_07,
        rc_resistance_1_07,
        rc_resistance_ratio_07,
        rc_series_07,
        rc_show_second_07,
    )


@app.cell
def _(
    np,
    rc_impedance_07,
    rc_log_separation_07,
    rc_log_tau_1_07,
    rc_resistance_1_07,
    rc_resistance_ratio_07,
    rc_series_07,
    rc_show_second_07,
):
    rc_frequency_07 = np.logspace(-4.5, 5.5, 600)
    rc_tau_1_07 = 10.0 ** rc_log_tau_1_07.value
    rc_capacitance_1_07 = rc_tau_1_07 / rc_resistance_1_07.value
    rc_tau_2_07 = rc_tau_1_07 * 10.0 ** rc_log_separation_07.value
    rc_resistance_2_07 = rc_resistance_1_07.value * rc_resistance_ratio_07.value
    rc_capacitance_2_07 = rc_tau_2_07 / rc_resistance_2_07
    rc_branches_07 = [(rc_resistance_1_07.value, rc_capacitance_1_07)]
    if rc_show_second_07.value:
        rc_branches_07.append((rc_resistance_2_07, rc_capacitance_2_07))
    rc_spectrum_07 = rc_impedance_07(rc_frequency_07, rc_series_07.value, rc_branches_07)
    return (
        rc_capacitance_1_07,
        rc_capacitance_2_07,
        rc_frequency_07,
        rc_spectrum_07,
        rc_tau_1_07,
        rc_tau_2_07,
    )


@app.cell
def _(
    np,
    plt,
    rc_capacitance_1_07,
    rc_capacitance_2_07,
    rc_frequency_07,
    rc_show_second_07,
    rc_spectrum_07,
    rc_tau_1_07,
    rc_tau_2_07,
):
    _figure, (_axis_nyquist, _axis_bode) = plt.subplots(
        1, 2, figsize=(13.2, 5.0), constrained_layout=True
    )
    _axis_nyquist.plot(rc_spectrum_07.real, -rc_spectrum_07.imag, lw=2.8, color="#4C7C86")
    _frequency_markers = np.array([1.0 / (2.0 * np.pi * rc_tau_1_07)])
    if rc_show_second_07.value:
        _frequency_markers = np.append(_frequency_markers, 1.0 / (2.0 * np.pi * rc_tau_2_07))
    for _marker_number, _marker_frequency in enumerate(_frequency_markers, start=1):
        _marker_index = int(np.argmin(np.abs(np.log(rc_frequency_07 / _marker_frequency))))
        _axis_nyquist.scatter(
            rc_spectrum_07.real[_marker_index],
            -rc_spectrum_07.imag[_marker_index],
            s=75,
            color="#B8734A" if _marker_number == 1 else "#7C6A91",
            zorder=4,
            label=rf"$\omega\tau_{_marker_number}=1$",
        )
    _axis_nyquist.set(
        xlabel=r"$Z'$ ($\Omega$)",
        ylabel=r"$-Z''$ ($\Omega$)",
        title="Nyquist: time-scale separation makes arcs visible",
    )
    _axis_nyquist.set_aspect("equal", adjustable="datalim")
    _axis_nyquist.grid(True)
    _axis_nyquist.legend(loc="best")

    _axis_bode.semilogx(rc_frequency_07, np.abs(rc_spectrum_07), lw=2.7, color="#4C7C86")
    _axis_bode_phase = _axis_bode.twinx()
    _axis_bode_phase.semilogx(
        rc_frequency_07, np.angle(rc_spectrum_07, deg=True), lw=2.3, color="#B8734A"
    )
    _axis_bode.set(xlabel="frequency (Hz)", ylabel=r"$|Z|$ ($\Omega$)", title="Bode: magnitude and phase retain frequency")
    _axis_bode_phase.set_ylabel(r"phase of $Z$ (degrees)", color="#B8734A")
    _axis_bode_phase.tick_params(axis="y", colors="#B8734A")
    _axis_bode.grid(True)
    _axis_bode.text(
        0.03,
        0.05,
        rf"$C_1={1e6 * rc_capacitance_1_07:.1f}\ \mu$F"
        + (rf"\n$C_2={1e6 * rc_capacitance_2_07:.1f}\ \mu$F" if rc_show_second_07.value else ""),
        transform=_axis_bode.transAxes,
        color="#526173",
    )
    _figure
    return


@app.cell
def _(mo):
    mo.md(r"""
    ### A ceramic interpretation — with a caution

    In the lecture's brick-layer picture, a high-frequency arc is often
    associated with grain interiors and a lower-frequency arc with grain
    boundaries. The assignment is plausible when the two regions have distinct
    $R C$ times and their fitted capacitances scale sensibly with geometry.
    When the times overlap, the Nyquist curve need not reveal two complete
    semicircles. Frequency range, capacitance, thickness scaling, and independent
    microstructural knowledge should support the assignment; arc order alone is
    not enough.
    """)
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## 2. Warburg impedance: diffusion written in the frequency domain

    Consider a neutral composition variable $c$ in a one-dimensional slab of
    length $L$ and cross-sectional area $S$:

    $$
    \frac{\partial c}{\partial t}=D^\delta\frac{\partial^2c}{\partial x^2}
    \quad\Longrightarrow\quad
    i\omega\widehat{\Delta c}=D^\delta\frac{d^2\widehat{\Delta c}}{dx^2}.
    $$

    Here $D^\delta$ is the same **chemical diffusivity** used in Modules 05 and
    06. A concentration wave penetrates a distance of order
    $\sqrt{D^\delta/\omega}$. At high frequency it samples only a thin region;
    at low frequency it can reach the far boundary.

    For a semi-infinite sample,

    $$
    Z_W=\frac{W}{\sqrt{i\omega}}
       =\frac{W(1-i)}{\sqrt{2\omega}},
    $$

    so $Z'=-Z''$: the famous $45^\circ$ line follows from diffusion, rather
    than being inserted as an equivalent-circuit element.

    With
    $\widetilde\omega=\omega L^2/D^\delta$ and
    $R_D=RTL/(n^2F^2c_0SD^\delta)$ for a dilute host, the two finite-length
    boundary conditions are

    $$
    \widetilde Z_{\rm open}
      =\frac{\tanh\sqrt{i\widetilde\omega}}{\sqrt{i\widetilde\omega}},
    \qquad
    \widetilde Z_{\rm blocked}
      =\frac{\coth\sqrt{i\widetilde\omega}}{\sqrt{i\widetilde\omega}}.
    $$

    In this section, **open** means the far-face composition perturbation is
    clamped to zero; **blocked** means its flux is zero. These diffusion labels
    should not be confused with electrical open circuit.
    """)
    return


@app.cell
def _(mo):
    warburg_boundary_07 = mo.ui.dropdown(
        options={
            "Open far boundary: fixed composition": "open",
            "Blocked far boundary: zero flux": "blocked",
        },
        value="Open far boundary: fixed composition",
        label="profile boundary",
    )
    warburg_log_omega_07 = mo.ui.slider(
        start=-3.0,
        stop=3.0,
        step=0.1,
        value=0.0,
        label=r"$\log_{10}\widetilde\omega$",
    )
    warburg_log_diffusivity_07 = mo.ui.slider(
        start=-12.0,
        stop=-5.0,
        step=0.25,
        value=-8.0,
        label=r"$\log_{10}(D^\delta/\mathrm{cm^2\,s^{-1}})$",
    )
    warburg_length_07 = mo.ui.slider(
        start=10, stop=500, step=10, value=100, label=r"$L$ ($\mu$m)"
    )
    mo.hstack(
        [
            warburg_boundary_07,
            warburg_log_omega_07,
            warburg_log_diffusivity_07,
            warburg_length_07,
        ],
        justify="start",
        gap=1.4,
    )
    return (
        warburg_boundary_07,
        warburg_length_07,
        warburg_log_diffusivity_07,
        warburg_log_omega_07,
    )


@app.cell
def _(
    np,
    warburg_boundary_07,
    warburg_impedance_07,
    warburg_length_07,
    warburg_log_diffusivity_07,
    warburg_log_omega_07,
    warburg_profile_07,
    warburg_scales_07,
):
    warburg_reduced_omega_07 = np.logspace(-4.0, 4.0, 500)
    warburg_spectra_07 = {
        _boundary: warburg_impedance_07(warburg_reduced_omega_07, _boundary)
        for _boundary in ("semi-infinite", "open", "blocked")
    }
    warburg_selected_omega_07 = 10.0 ** warburg_log_omega_07.value
    warburg_position_07 = np.linspace(0.0, 1.0, 260)
    warburg_selected_profile_07 = warburg_profile_07(
        warburg_selected_omega_07, warburg_position_07, warburg_boundary_07.value
    )
    warburg_phases_07 = np.array([0.0, 0.5, 1.0, 1.5]) * np.pi
    warburg_scale_data_07 = warburg_scales_07(
        warburg_length_07.value,
        10.0 ** warburg_log_diffusivity_07.value,
        concentration_mol_per_m3=160.0,
        area_cm2=0.5,
        temperature_k=800.0,
        charge_number=1.0,
    )
    return (
        warburg_phases_07,
        warburg_position_07,
        warburg_reduced_omega_07,
        warburg_scale_data_07,
        warburg_selected_omega_07,
        warburg_selected_profile_07,
        warburg_spectra_07,
    )


@app.cell
def _(
    np,
    plt,
    warburg_boundary_07,
    warburg_phases_07,
    warburg_position_07,
    warburg_reduced_omega_07,
    warburg_scale_data_07,
    warburg_selected_omega_07,
    warburg_selected_profile_07,
    warburg_spectra_07,
):
    _figure, (_axis_warburg, _axis_profile) = plt.subplots(
        1, 2, figsize=(13.4, 5.2), constrained_layout=True
    )
    _colors = {
        "semi-infinite": "#7A8793",
        "open": "#4C7C86",
        "blocked": "#B8734A",
    }
    _labels = {
        "semi-infinite": "semi-infinite",
        "open": "finite, fixed composition",
        "blocked": "finite, zero flux",
    }
    for _boundary, _impedance in warburg_spectra_07.items():
        _plot_mask = (-_impedance.imag <= 4.0) & (_impedance.real <= 2.2)
        _axis_warburg.plot(
            _impedance.real[_plot_mask],
            -_impedance.imag[_plot_mask],
            lw=2.6,
            color=_colors[_boundary],
            label=_labels[_boundary],
        )
    _selected_impedance = warburg_spectra_07[warburg_boundary_07.value]
    _selected_index = int(
        np.argmin(np.abs(np.log(warburg_reduced_omega_07 / warburg_selected_omega_07)))
    )
    _axis_warburg.scatter(
        _selected_impedance.real[_selected_index],
        -_selected_impedance.imag[_selected_index],
        s=85,
        color="#C49345",
        edgecolor="white",
        linewidth=0.8,
        zorder=5,
        label="selected frequency",
    )
    _axis_warburg.plot([0.0, 1.2], [0.0, 1.2], ls="--", lw=1.3, color="#9AA3AB")
    _axis_warburg.annotate("45° semi-infinite limit", xy=(0.65, 0.65), xytext=(0.92, 0.28), arrowprops={"arrowstyle": "->", "color": "#73808C"}, color="#526173")
    _axis_warburg.set(
        xlim=(-0.04, 1.55),
        ylim=(-0.04, 2.0),
        xlabel=r"$\widetilde Z'$",
        ylabel=r"$-\widetilde Z''$",
        title="The far boundary changes only the low-frequency end",
    )
    _axis_warburg.grid(True)
    _axis_warburg.legend(loc="upper left")

    _profile_colors = ["#4C7C86", "#B8734A", "#7C6A91", "#5F8A6B"]
    for _phase, _color in zip(warburg_phases_07, _profile_colors):
        _snapshot = np.real(warburg_selected_profile_07 * np.exp(1j * _phase))
        _axis_profile.plot(
            warburg_position_07,
            _snapshot,
            lw=2.4,
            color=_color,
            label=rf"$\omega t={_phase / np.pi:.1f}\pi$",
        )
    _axis_profile.axhline(0.0, color="#73808C", lw=0.9)
    _axis_profile.set(
        xlabel=r"position $x/L$",
        ylabel=r"$\Delta c(x,t)/|\widehat{\Delta c}(0)|$",
        title=f"Concentration wave: {warburg_boundary_07.value} far boundary",
    )
    _axis_profile.grid(True)
    _axis_profile.legend(loc="best", ncols=2)
    _axis_profile.text(
        0.03,
        0.05,
        rf"$f_D={warburg_scale_data_07['frequency_diffusion_hz']:.3g}$ Hz"
        + "\n"
        + rf"$R_D={warburg_scale_data_07['resistance_diffusion_ohm']:.3g}\ \Omega$",
        transform=_axis_profile.transAxes,
        color="#526173",
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.82},
    )
    _figure
    return


@app.cell
def _(mo):
    mo.md(r"""
    ### Read the two finite-length limits

    $$
    \begin{aligned}
    \widetilde Z_{\rm open}&\longrightarrow 1-
      \frac{i\widetilde\omega}{3},\\[2mm]
    \widetilde Z_{\rm blocked}&\longrightarrow
      \frac{1}{3}+\frac{1}{i\widetilde\omega},
    \end{aligned}
    \qquad \widetilde\omega\ll1.
    $$

    A composition reservoir allows a steady diffusive flux, so the open case
    ends at a finite resistance. A blocking wall cannot sustain a dc flux, so
    material accumulates and the response becomes capacitive. At
    $\widetilde\omega\gg1$, neither experiment feels the far face and both
    recover the same semi-infinite $45^\circ$ response.

    The dimensionless curves contain the shape. The controls show how $L$ and
    $D^\delta$ move that shape along the laboratory frequency axis through
    $f_D=D^\delta/(2\pi L^2)$. The displayed resistance scale is an illustrative
    dilute example with $T=800$ K, $c_0=160$ mol m$^{-3}$, $S=0.5$ cm$^2$,
    and $n=1$; those fixed values affect $R_D$ but not the dimensionless curve.
    """)
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## 3. Why a mixed conductor needs a transmission line

    A homogeneous MIEC has two conducting pathways. The electronic rail carries
    $I_e$ with resistance per length $r_e=1/(\sigma_e S)$; the ionic rail
    carries $I_i$ with $r_i=1/(\sigma_i S)$. Local stoichiometry can change when
    charge transfers between the rails, represented in the small-signal limit by
    the **chemical capacitance per length** $c_{\rm chem}$.

    The continuous TLM used here is

    $$
    \frac{du_e}{dx}=-r_e I_e,\qquad
    \frac{du_i}{dx}=-r_i I_i,
    $$

    $$
    \frac{dI_e}{dx}=-i\omega c_{\rm chem}(u_e-u_i),\qquad
    \frac{dI_i}{dx}=+i\omega c_{\rm chem}(u_e-u_i).
    $$

    $u_e$ and $u_i$ are **voltage-equivalent electrochemical potentials** in
    volts. Their difference drives local chemical storage. Adding the last two
    equations gives $d(I_e+I_i)/dx=0$: total current is conserved even though
    its division between carriers changes with position.

    A useful time scale is

    $$
    \tau_{\rm chem}=(R_e+R_i)C_{\rm chem},\qquad
    f_{\rm chem}=\frac{1}{2\pi\tau_{\rm chem}},
    $$

    where $R_e=r_eL$, $R_i=r_iL$, and
    $C_{\rm chem}=c_{\rm chem}L$. It is a scaling frequency, not a promise that
    every boundary condition has a peak exactly there.
    """)
    return


@app.cell
def _(np, plt):
    _figure, _axis = plt.subplots(figsize=(12.6, 3.2), constrained_layout=True)
    _x_nodes = np.linspace(0.12, 0.88, 6)
    _y_e, _y_i = 0.72, 0.25
    _axis.plot([0.05, 0.95], [_y_e, _y_e], color="#B8734A", lw=4.0)
    _axis.plot([0.05, 0.95], [_y_i, _y_i], color="#4C7C86", lw=4.0)
    for _node in _x_nodes:
        _axis.plot([_node, _node], [_y_i + 0.06, _y_e - 0.06], color="#7C6A91", lw=1.8)
        _axis.plot([_node - 0.018, _node + 0.018], [0.50, 0.50], color="#7C6A91", lw=2.0)
        _axis.plot([_node - 0.018, _node + 0.018], [0.46, 0.46], color="#7C6A91", lw=2.0)
    _axis.text(0.50, 0.82, r"electronic rail: $r_e\,dx$", ha="center", color="#8F5638", fontsize=15)
    _axis.text(0.50, 0.10, r"ionic rail: $r_i\,dx$", ha="center", color="#3D6972", fontsize=15)
    _axis.text(0.50, 0.53, r"distributed $c_{\rm chem}\,dx$", ha="center", color="#665777", fontsize=14)
    _axis.text(0.02, 0.48, "left contacts", ha="left", va="center", color="#526173")
    _axis.text(0.98, 0.48, "right contacts", ha="right", va="center", color="#526173")
    _axis.set(xlim=(0.0, 1.0), ylim=(0.0, 1.0), title="A continuous two-rail transmission line")
    _axis.axis("off")
    _figure
    return


@app.cell
def _(mo):
    tlm_contact_case_07 = mo.ui.dropdown(
        options={
            "Electron-reversible contacts; ions blocked at both faces": "electron contacts, ions blocked",
            "Current collector on left; ion electrolyte on right": "cross-selective contacts",
            "Both carriers reversible at both faces": "both carriers reversible",
        },
        value="Electron-reversible contacts; ions blocked at both faces",
        label="contact boundary conditions",
    )
    tlm_log_ratio_07 = mo.ui.slider(
        start=-3.0,
        stop=3.0,
        step=0.25,
        value=0.0,
        label=r"$\log_{10}(\sigma_e/\sigma_i)$",
    )
    tlm_log_resistance_07 = mo.ui.slider(
        start=0.0,
        stop=4.0,
        step=0.25,
        value=2.0,
        label=r"$\log_{10}(R_\parallel/\Omega)$",
    )
    tlm_log_capacitance_07 = mo.ui.slider(
        start=-7.0,
        stop=-1.0,
        step=0.25,
        value=-4.0,
        label=r"$\log_{10}(C_{\rm chem}/\mathrm{F})$",
    )
    tlm_log_selected_frequency_07 = mo.ui.slider(
        start=-3.0,
        stop=3.0,
        step=0.1,
        value=0.0,
        label=r"$\log_{10}(f/f_{\rm chem})$",
    )
    tlm_profile_phase_07 = mo.ui.slider(
        start=0, stop=330, step=30, value=60, label=r"snapshot phase $\omega t$ (degrees)"
    )
    mo.vstack(
        [
            tlm_contact_case_07,
            mo.hstack(
                [tlm_log_ratio_07, tlm_log_resistance_07, tlm_log_capacitance_07],
                justify="start",
                gap=1.4,
            ),
            mo.hstack(
                [tlm_log_selected_frequency_07, tlm_profile_phase_07],
                justify="start",
                gap=1.4,
            ),
        ]
    )
    return (
        tlm_contact_case_07,
        tlm_log_capacitance_07,
        tlm_log_ratio_07,
        tlm_log_resistance_07,
        tlm_log_selected_frequency_07,
        tlm_profile_phase_07,
    )


@app.cell
def _(
    np,
    tlm_contact_case_07,
    tlm_log_capacitance_07,
    tlm_log_ratio_07,
    tlm_log_resistance_07,
    tlm_log_selected_frequency_07,
    tlm_parameters_07,
    tlm_profile_07,
    tlm_spectrum_07,
):
    tlm_parameter_data_07 = tlm_parameters_07(
        10.0 ** tlm_log_resistance_07.value,
        10.0 ** tlm_log_capacitance_07.value,
        10.0 ** tlm_log_ratio_07.value,
    )
    tlm_frequency_ratio_07 = np.logspace(-3.0, 3.0, 420)
    tlm_selected_omega_07 = 10.0 ** tlm_log_selected_frequency_07.value
    tlm_spectrum_data_07 = tlm_spectrum_07(
        tlm_frequency_ratio_07, tlm_parameter_data_07, tlm_contact_case_07.value
    )
    tlm_position_07 = np.linspace(0.0, 1.0, 260)
    tlm_profile_data_07 = tlm_profile_07(
        tlm_position_07,
        tlm_selected_omega_07,
        tlm_parameter_data_07,
        tlm_contact_case_07.value,
    )
    return (
        tlm_frequency_ratio_07,
        tlm_parameter_data_07,
        tlm_position_07,
        tlm_profile_data_07,
        tlm_selected_omega_07,
        tlm_spectrum_data_07,
    )


@app.cell
def _(
    np,
    plt,
    tlm_contact_case_07,
    tlm_frequency_ratio_07,
    tlm_parameter_data_07,
    tlm_selected_omega_07,
    tlm_spectrum_data_07,
):
    _normalized_impedance = tlm_spectrum_data_07 / tlm_parameter_data_07["R_parallel_ohm"]
    _minus_imaginary = -_normalized_impedance.imag
    _real_extent = max(1.0, float(np.nanmax(_normalized_impedance.real)))
    _finite_capacitive = _minus_imaginary[
        np.isfinite(_minus_imaginary) & (_minus_imaginary >= 0.0)
    ]
    _imaginary_extent = float(np.max(_finite_capacitive))
    _nyquist_cap = max(0.02, 1.18 * _imaginary_extent)
    if tlm_contact_case_07.value == "cross-selective contacts":
        _nyquist_cap = min(8.0, max(1.4, 2.2 * _real_extent))
    _visible = (
        np.isfinite(_normalized_impedance.real)
        & np.isfinite(_minus_imaginary)
        & (_minus_imaginary >= -1.0e-9)
        & (_minus_imaginary <= _nyquist_cap)
    )
    _selected_index = int(
        np.argmin(np.abs(np.log(tlm_frequency_ratio_07 / tlm_selected_omega_07)))
    )

    _figure, (_axis_nyquist, _axis_bode) = plt.subplots(
        1, 2, figsize=(13.4, 5.2), constrained_layout=True
    )
    _axis_nyquist.plot(
        _normalized_impedance.real[_visible],
        _minus_imaginary[_visible],
        lw=2.8,
        color="#4C7C86",
    )
    if _visible[_selected_index]:
        _axis_nyquist.scatter(
            _normalized_impedance.real[_selected_index],
            _minus_imaginary[_selected_index],
            s=90,
            color="#C49345",
            edgecolor="white",
            linewidth=0.8,
            zorder=5,
            label="selected frequency",
        )
    if np.any(~_visible & (_minus_imaginary > _nyquist_cap)):
        _axis_nyquist.annotate(
            "low-frequency branch continues upward",
            xy=(float(_normalized_impedance.real[_visible][0]), _nyquist_cap * 0.96),
            xytext=(0.44, 0.72),
            textcoords="axes fraction",
            arrowprops={"arrowstyle": "->", "color": "#73808C"},
            color="#526173",
        )
    _axis_nyquist.set(
        ylim=(-0.03 * _nyquist_cap, _nyquist_cap),
        xlabel=r"$Z'/R_\parallel$",
        ylabel=r"$-Z''/R_\parallel$",
        title=tlm_contact_case_07.value.capitalize()
        + "\nNyquist response selected by the contacts",
    )
    _axis_nyquist.grid(True)
    if _visible[_selected_index]:
        _axis_nyquist.legend(loc="best")

    _axis_bode.loglog(
        tlm_frequency_ratio_07,
        np.abs(_normalized_impedance),
        lw=2.7,
        color="#4C7C86",
    )
    _axis_phase = _axis_bode.twinx()
    _axis_phase.semilogx(
        tlm_frequency_ratio_07,
        np.angle(_normalized_impedance, deg=True),
        lw=2.3,
        color="#B8734A",
    )
    _axis_bode.axvline(1.0, color="#9AA3AB", lw=1.2, ls="--")
    _axis_bode.set(
        xlabel=r"$f/f_{\rm chem}=\omega\tau_{\rm chem}$",
        ylabel=r"$|Z|/R_\parallel$",
        title="Frequency keeps the time-scale information",
    )
    _axis_phase.set_ylabel(r"phase of $Z$ (degrees)", color="#B8734A")
    _axis_phase.tick_params(axis="y", colors="#B8734A")
    _axis_bode.grid(True, which="both")
    _axis_bode.text(
        0.03,
        0.05,
        rf"$R_e={tlm_parameter_data_07['R_e_ohm']:.3g}\ \Omega$"
        + "\n"
        + rf"$R_i={tlm_parameter_data_07['R_i_ohm']:.3g}\ \Omega$"
        + "\n"
        + rf"$f_{{\rm chem}}={tlm_parameter_data_07['frequency_chemical_hz']:.3g}$ Hz",
        transform=_axis_bode.transAxes,
        color="#526173",
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.82},
    )
    _figure
    return


@app.cell
def _(
    np,
    plt,
    tlm_parameter_data_07,
    tlm_position_07,
    tlm_profile_data_07,
    tlm_profile_phase_07,
    tlm_selected_omega_07,
):
    _phase_factor = np.exp(1j * np.deg2rad(tlm_profile_phase_07.value))
    _u_e_snapshot = np.real(tlm_profile_data_07["u_e"] * _phase_factor)
    _u_i_snapshot = np.real(tlm_profile_data_07["u_i"] * _phase_factor)
    _chemical_snapshot = np.real(tlm_profile_data_07["u_chemical"] * _phase_factor)
    _current_scale = max(
        np.max(np.abs(tlm_profile_data_07["I_total"])),
        1.0 / tlm_parameter_data_07["R_parallel_ohm"],
    )
    _i_e_snapshot = np.real(tlm_profile_data_07["I_e"] * _phase_factor) / _current_scale
    _i_i_snapshot = np.real(tlm_profile_data_07["I_i"] * _phase_factor) / _current_scale
    _i_total_snapshot = np.real(tlm_profile_data_07["I_total"] * _phase_factor) / _current_scale

    _figure, (_axis_potential, _axis_current) = plt.subplots(
        1, 2, figsize=(13.4, 4.9), constrained_layout=True
    )
    _axis_potential.plot(tlm_position_07, _u_e_snapshot, lw=2.8, color="#B8734A", label=r"$u_e$")
    _axis_potential.plot(tlm_position_07, _u_i_snapshot, lw=2.8, color="#4C7C86", label=r"$u_i$")
    _axis_potential.plot(
        tlm_position_07,
        _chemical_snapshot,
        lw=2.0,
        ls="--",
        color="#7C6A91",
        label=r"$u_e-u_i$",
    )
    _axis_potential.axhline(0.0, color="#73808C", lw=0.9)
    _axis_potential.set(
        xlabel=r"position $x/L$",
        ylabel="instantaneous voltage-equivalent potential (V)",
        title="Potentials show where chemical storage is driven",
    )
    _axis_potential.grid(True)
    _axis_potential.legend(loc="best")

    _axis_current.plot(tlm_position_07, _i_e_snapshot, lw=2.8, color="#B8734A", label=r"$I_e$")
    _axis_current.plot(tlm_position_07, _i_i_snapshot, lw=2.8, color="#4C7C86", label=r"$I_i$")
    _axis_current.plot(
        tlm_position_07,
        _i_total_snapshot,
        lw=2.1,
        ls="--",
        color="#5F8A6B",
        label=r"$I_e+I_i$",
    )
    _axis_current.axhline(0.0, color="#73808C", lw=0.9)
    _axis_current.set(
        xlabel=r"position $x/L$",
        ylabel="instantaneous current / common scale",
        title="Carrier currents exchange, but their sum stays constant",
    )
    _axis_current.grid(True)
    _axis_current.legend(loc="best")
    _figure.suptitle(
        rf"Internal TLM state at $f/f_{{\rm chem}}={tlm_selected_omega_07:.3g}$ and $\omega t={tlm_profile_phase_07.value}^\circ$",
        fontsize=17,
    )
    _figure
    return


@app.cell
def _(mo, tlm_contact_case_07):
    _interpretations = {
        "electron contacts, ions blocked": (
            "Electrons can cross both interfaces, whereas ions accumulate internally. "
            "The dc limit remains electronically conducting, but chemical redistribution "
            "adds a frequency-dependent polarization."
        ),
        "cross-selective contacts": (
            "The left face passes electrons and the right face passes ions. No carrier has "
            "a complete dc path by itself, so the low-frequency response becomes blocking "
            "and capacitive. This is the selective-contact geometry behind Modules 06 and many "
            "chemical-diffusion measurements."
        ),
        "both carriers reversible": (
            "Both rails are pinned to the same reservoir voltage at each face. Then "
            "$u_e-u_i=0$, chemical storage is not excited, and the response is simply "
            "$R_e\\parallel R_i$."
        ),
    }
    mo.callout(
        mo.md("**What this contact choice means.** " + _interpretations[tlm_contact_case_07.value]),
        kind="info",
    )
    return


@app.cell
def _(
    np,
    rc_impedance_07,
    tlm_parameters_07,
    tlm_profile_07,
    tlm_solution_07,
    tlm_spectrum_07,
    warburg_impedance_07,
    waveform_data_07,
):
    _check_frequency = np.logspace(-5.0, 5.0, 800)
    _check_rc = rc_impedance_07(_check_frequency, 0.0, [(1.0, 1.0)])
    rc_circle_error_07 = np.max(
        np.abs((_check_rc.real - 0.5) ** 2 + _check_rc.imag**2 - 0.25)
    )
    _peak_index = int(np.argmax(-_check_rc.imag))
    rc_peak_error_07 = abs(2.0 * np.pi * _check_frequency[_peak_index] - 1.0)

    _waveform_check = waveform_data_07(3.0, 37.0)
    phasor_sign_error_07 = abs(_waveform_check[3] + 37.0)

    _omega_high = 1.0e4
    _warburg_semi_high = warburg_impedance_07(np.array([_omega_high]), "semi-infinite")[0]
    warburg_high_frequency_error_07 = max(
        abs(warburg_impedance_07(np.array([_omega_high]), "open")[0] / _warburg_semi_high - 1.0),
        abs(warburg_impedance_07(np.array([_omega_high]), "blocked")[0] / _warburg_semi_high - 1.0),
    )
    _omega_low = 1.0e-6
    _warburg_open_low = warburg_impedance_07(np.array([_omega_low]), "open")[0]
    _warburg_blocked_low = warburg_impedance_07(np.array([_omega_low]), "blocked")[0]
    warburg_low_frequency_error_07 = max(
        abs(_warburg_open_low.real - 1.0),
        abs(_warburg_blocked_low.real - 1.0 / 3.0),
        abs((-_warburg_blocked_low.imag) * _omega_low - 1.0),
    )
    _warburg_passivity_grid = np.logspace(-5.0, 5.0, 500)
    warburg_passivity_margin_07 = min(
        np.min(warburg_impedance_07(_warburg_passivity_grid, _boundary).real)
        for _boundary in ("semi-infinite", "open", "blocked")
    )

    _tlm_check_parameters = tlm_parameters_07(100.0, 1.0e-4, 100.0)
    _tlm_check_frequency = np.logspace(-3.0, 3.0, 220)
    _tlm_reversible = tlm_spectrum_07(
        _tlm_check_frequency, _tlm_check_parameters, "both carriers reversible"
    )
    tlm_reversible_error_07 = np.max(
        np.abs(_tlm_reversible / _tlm_check_parameters["R_parallel_ohm"] - 1.0)
    )
    _tlm_cases = (
        "electron contacts, ions blocked",
        "cross-selective contacts",
        "both carriers reversible",
    )
    tlm_passivity_margin_07 = min(
        np.min(tlm_spectrum_07(_tlm_check_frequency, _tlm_check_parameters, _case).real)
        for _case in _tlm_cases
    )
    tlm_boundary_error_07 = max(
        tlm_solution_07(10.0 ** _frequency, _tlm_check_parameters, _case)["boundary_residual"]
        for _case in _tlm_cases
        for _frequency in (-3.0, 0.0, 3.0)
    )
    _tlm_check_profile = tlm_profile_07(
        np.linspace(0.0, 1.0, 251),
        1.7,
        _tlm_check_parameters,
        "cross-selective contacts",
    )
    tlm_current_conservation_error_07 = np.ptp(_tlm_check_profile["I_total"]) / max(
        abs(_tlm_check_profile["I_total"][0]), 1.0e-30
    )
    tlm_finiteness_07 = all(
        np.all(np.isfinite(_tlm_check_profile[_key]))
        for _key in ("u_e", "u_i", "u_chemical", "I_e", "I_i", "I_total")
    )
    return (
        phasor_sign_error_07,
        rc_circle_error_07,
        rc_peak_error_07,
        tlm_boundary_error_07,
        tlm_current_conservation_error_07,
        tlm_finiteness_07,
        tlm_passivity_margin_07,
        tlm_reversible_error_07,
        warburg_high_frequency_error_07,
        warburg_low_frequency_error_07,
        warburg_passivity_margin_07,
    )


@app.cell
def _(
    mo,
    phasor_sign_error_07,
    rc_circle_error_07,
    rc_peak_error_07,
    tlm_boundary_error_07,
    tlm_current_conservation_error_07,
    tlm_finiteness_07,
    tlm_passivity_margin_07,
    tlm_reversible_error_07,
    warburg_high_frequency_error_07,
    warburg_low_frequency_error_07,
    warburg_passivity_margin_07,
):
    _checks = [
        (
            "Phasor sign convention",
            phasor_sign_error_07 < 1.0e-12,
            "A leading capacitive current must give a negative impedance phase for $e^{i\\omega t}$.",
        ),
        (
            "Ideal $R\\parallel C$ semicircle",
            max(rc_circle_error_07, rc_peak_error_07) < 2.0e-2,
            "The geometric arc and its $\\omega RC=1$ apex must come from the same circuit equation.",
        ),
        (
            "Finite Warburg limits",
            max(warburg_high_frequency_error_07, warburg_low_frequency_error_07) < 2.0e-6,
            "Both boundaries share the high-frequency 45 degree limit and approach their own low-frequency forms.",
        ),
        (
            "Passive Warburg response",
            warburg_passivity_margin_07 > 0.0,
            "Diffusion dissipates energy, so the real part of its impedance cannot be negative.",
        ),
        (
            "TLM boundary conditions",
            tlm_boundary_error_07 < 2.0e-10,
            "Every spectrum must satisfy the selected passing or blocking conditions at both faces.",
        ),
        (
            "TLM total-current conservation",
            abs(tlm_current_conservation_error_07) < 1.0e-12,
            "Current may transfer between rails, but $I_e+I_i$ must be independent of position.",
        ),
        (
            "Reversible-contact limit",
            tlm_reversible_error_07 < 2.0e-10,
            "Pinning both rails removes chemical polarization and leaves $R_e\\parallel R_i$.",
        ),
        (
            "Finite, passive TLM solution",
            tlm_finiteness_07 and tlm_passivity_margin_07 > -1.0e-10,
            "The internal state must stay finite and a passive line cannot generate power.",
        ),
    ]
    _rows = []
    for _name, _passed, _why in _checks:
        _rows.append(
            f"| {'PASS' if _passed else 'CHECK'} | {_name} | {_why} |"
        )
    _heading = [
        "## 4. Physical consistency checks",
        "",
        "These checks connect the sign convention, circuit response, diffusion boundaries, "
        "and transmission-line interpretation.",
        "",
        "| status | check | why it matters |",
        "|---:|---|---|",
    ]
    mo.md("\n".join(_heading + _rows))
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## 5. Take-home map

    $$
    \boxed{
    \text{small sinusoid}
    \rightarrow Z(\omega)
    \rightarrow \text{time scale}
    \rightarrow \text{transport + storage + boundaries}
    }
    $$

    - A semicircle follows from a relaxation; its diameter and characteristic
      frequency carry different information.
    - A Warburg response is the frequency-domain solution of chemical diffusion.
      Its low-frequency end identifies what the far boundary allows.
    - In a MIEC, the TLM keeps electronic conduction, ionic conduction, chemical
      storage, and contact selectivity in one constrained physical model.
    - The same interior material can show a different spectrum when the contacts
      change. Always state geometry, sign convention, and boundary conditions
      before assigning a feature.

    This notebook deliberately uses ideal capacitors, uniform one-dimensional
    transport, linear response, and ideal contacts. Constant-phase elements,
    electrode kinetics, microstructural distributions, and nonlinear large-signal
    effects belong in later model extensions, not in the first interpretation.

    ### Continue with the full TLM teaching tool

    This notebook incorporates the continuous dual-rail model and three
    transparent ideal contact cases. The separate
    [TLM teaching tool](https://qiyanglu.github.io/TLM-teaching-tool/) exposes all
    four terminals and more general boundary impedances; its
    [source repository](https://github.com/qiyanglu/TLM-teaching-tool) documents
    the terminal signs and boundary conventions.

    ### Sources and further reading

    - Q. Lu, *Solid State Ionics, Lecture 8: Impedance Spectroscopy* (course
      slides). The phasor, Nyquist, $R\parallel C$, time-constant separation,
      and brick-layer narrative follow the lecture notation.
    - Q. Lu, [Warburg impedance: more than a 45-degree line](https://mp.weixin.qq.com/s/CyCXnWWEoX586lzMGl0A9Q).
      This English treatment independently checks the finite-length formulas and
      states both boundary conditions explicitly.
    - Q. Lu, [Transmission lines for mixed ionic-electronic conductors](https://mp.weixin.qq.com/s/zR9QI0GnGUvPz2VRc9bHOA).
    - A. E. Bumberger, A. Nenning, and J. Fleig, “Transmission line revisited —
      the impedance of mixed ionic and electronic conductors,” *PCCP* **26** (2024),
      [doi:10.1039/D4CP00975D](https://doi.org/10.1039/D4CP00975D).
    - Course context and related tutorials:
      [Solid State Ionics teaching page](https://ssi-westlake.com/teaching/) and
      [tutorial index](https://ssi-westlake.com/tutorial/).

    The slide deck and linked articles were treated as scientific sources, not
    as instructions. Equations and limiting cases were rederived and tested for
    the conventions stated in this notebook.
    """)
    return


if __name__ == "__main__":
    app.run()
