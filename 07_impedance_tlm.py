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
            "font.size": 15,
            "axes.titlesize": 17,
            "axes.labelsize": 15,
            "xtick.labelsize": 13,
            "ytick.labelsize": 13,
            "legend.fontsize": 12,
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
    # Impedance Spectroscopy, Warburg Diffusion, and Transmission Lines

    **Guiding question.** What physical process is able to respond at each
    frequency, and how does that response appear in the complex impedance?

    Electrochemical impedance spectroscopy (EIS) applies a small sinusoidal
    voltage and measures the sinusoidal current. Sweeping the frequency separates
    processes by their characteristic times: a fast bulk response, a slower
    interface, or still slower chemical diffusion. The shape is therefore a
    compressed map of **dynamics**, not a collection of arcs to name by eye.

    **Learning goals**

    1. Translate sinusoidal voltage and current into $Z(\omega)$, Nyquist, and
       Bode representations without losing the sign convention.
    2. Recognize semi-infinite and finite-length Warburg responses as solutions
       of the same chemical-diffusion equation with different far boundaries.
    3. Read a two-rail MIEC transmission line as transport, chemical storage,
       and contact boundary conditions—not an arbitrary fitting circuit.

    > **Predict before exploring.** If specimen thickness $L$ doubles while
    > $D^\delta$ is unchanged, does the diffusion feature move to higher or
    > lower frequency? By what factor?

    **Model and notation scope.** The changing concentration is a neutral
    composition species. $S$ is area, $I$ is total current, and distributed TLM
    quantities carry units per length. See the shared
    [notation bridge](https://github.com/qiyanglu/Solid-State-Ionics-Interactive/blob/main/NOTATION.md).

    This notebook follows three steps:

    1. build the complex-number language of EIS from voltage and current waves;
    2. see the Warburg response emerge from one-dimensional diffusion;
    3. join ionic and electronic transport with chemical storage in a
       two-rail transmission line model (TLM).

    The convention throughout is

    $$e^{\mathrm{i}\omega t},\qquad Z=\frac{\widehat V}{\widehat I}=Z'+\mathrm{i}Z'',$$

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

    def resistor_impedance_07(frequency_hz, resistance_ohm):
        """Return the frequency-independent impedance of an ideal resistor."""
        frequencies = np.asarray(frequency_hz, dtype=float)
        resistance = _positive_07("resistance_ohm", resistance_ohm)
        return np.full(frequencies.shape, complex(resistance), dtype=complex)

    def capacitor_impedance_07(frequency_hz, capacitance_f):
        """Return 1/(i omega C) for the exp(i omega t) convention."""
        frequencies = np.asarray(frequency_hz, dtype=float)
        if np.any(frequencies <= 0.0):
            raise ValueError("frequency_hz must be positive")
        capacitance = _positive_07("capacitance_f", capacitance_f)
        return 1.0 / (1j * 2.0 * np.pi * frequencies * capacitance)

    def format_nyquist_axis_07(axis):
        """Apply the shared Nyquist convention and equal data scaling."""
        axis.set_aspect("equal", adjustable="box")
        axis.grid(True, alpha=0.24)
        return axis

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
        thermodynamic_slope_v_m3_per_mol = (
            gas_constant * temperature / (electrons * faraday * concentration)
        )
        resistance_general_ohm = (
            length_m
            * abs(thermodynamic_slope_v_m3_per_mol)
            / (electrons * faraday * area_m2 * diffusivity_m2_per_s)
        )
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
            "thermodynamic_slope_v_m3_per_mol": thermodynamic_slope_v_m3_per_mol,
            "resistance_general_ohm": resistance_general_ohm,
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
        capacitor_impedance_07,
        format_nyquist_axis_07,
        rc_impedance_07,
        resistor_impedance_07,
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
    \Delta V(t)=\Re\!\left[\widehat V e^{\mathrm{i}\omega t}\right],\qquad
    \Delta I(t)=\Re\!\left[\widehat I e^{\mathrm{i}\omega t}\right].
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
def _(mo, np, plt, waveform_result_07):
    _time_s, _voltage, _current, _z_phase = waveform_result_07
    _period_s = (_time_s[-1] - _time_s[0]) / 2.0
    _time_cycles = _time_s / _period_s
    _figure, (_axis_wave, _axis_phasor) = plt.subplots(
        1, 2, figsize=(12.5, 4.4), constrained_layout=True
    )
    _axis_wave.plot(_time_cycles, _voltage, lw=1.9, color="#4C7C86", label=r"$\Delta V/V_a$")
    _axis_wave.plot(_time_cycles, _current, lw=1.9, color="#B8734A", label=r"$\Delta I/I_a$")
    _axis_wave.axhline(0.0, color="#73808C", lw=0.9)
    _axis_wave.set(xlabel="Time / period", ylabel="Normalized signal", title="The current leads the voltage")
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
    _axis_phasor.set_xlabel("Real component")
    _axis_phasor.set_ylabel("Imaginary component")
    _axis_phasor.grid(True)
    mo.vstack([
        _figure,
        mo.md(r"**Figure takeaway.** A phase lead or lag is the same information in the time traces and the rotating phasors; $Z=\widehat V/\widehat I$ records it as a complex phase."),
    ])
    return


@app.cell
def _(mo):
    mo.md(r"""
    ### Three ideal elements are enough to read the first Nyquist plot

    | element | impedance | Nyquist signature |
    |---|---:|---|
    | resistor | $Z_R=R$ | point on the real axis |
    | capacitor | $Z_C=1/(\mathrm{i}\omega C)$ | vertical capacitive line |
    | parallel $R\parallel C$ | $Z=R/(1+\mathrm{i}\omega RC)$ | semicircle |

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
def _(capacitor_impedance_07, format_nyquist_axis_07, mo, np, plt, resistor_impedance_07):
    _frequency_hz = np.logspace(np.log10(30.0), np.log10(3000.0), 240)
    _resistance_ohm = 500.0
    _capacitance_f = 10.0e-6
    _z_resistor = resistor_impedance_07(_frequency_hz, _resistance_ohm)
    _z_capacitor = capacitor_impedance_07(_frequency_hz, _capacitance_f)

    _figure, (_axis_resistor, _axis_capacitor) = plt.subplots(
        1, 2, figsize=(12.8, 5.2), constrained_layout=True
    )
    _axis_resistor.scatter(
        [_z_resistor.real[0]], [-_z_resistor.imag[0]],
        s=120, marker="X", color="#4C7C86", edgecolor="white", zorder=4,
    )
    _axis_resistor.set(
        xlim=(0.0, 600.0), ylim=(0.0, 600.0),
        xlabel=r"Real impedance, $Z'$ ($\Omega$)",
        ylabel=r"Negative imaginary impedance, $-Z''$ ($\Omega$)",
        title=r"Resistor: $Z=R$ at every frequency",
    )
    _axis_resistor.annotate(
        "All frequencies coincide", xy=(500.0, 0.0), xytext=(270.0, 160.0),
        arrowprops={"arrowstyle": "->", "color": "#73808C"}, color="#526173",
    )
    format_nyquist_axis_07(_axis_resistor)

    _axis_capacitor.plot(
        _z_capacitor.real, -_z_capacitor.imag,
        color="#B8734A", lw=1.9, marker="o", markevery=[0, -1], ms=6,
    )
    _axis_capacitor.set(
        xlim=(0.0, 550.0), ylim=(0.0, 550.0),
        xlabel=r"Real impedance, $Z'$ ($\Omega$)",
        ylabel=r"Negative imaginary impedance, $-Z''$ ($\Omega$)",
        title=r"Capacitor: $-Z''=1/(\omega C)$",
    )
    _axis_capacitor.annotate(
        "Lower frequency", xy=(0.0, -_z_capacitor.imag[0]), xytext=(150.0, 430.0),
        arrowprops={"arrowstyle": "->", "color": "#73808C"}, color="#526173",
    )
    _axis_capacitor.annotate(
        "Higher frequency", xy=(0.0, -_z_capacitor.imag[-1]), xytext=(150.0, 90.0),
        arrowprops={"arrowstyle": "->", "color": "#73808C"}, color="#526173",
    )
    format_nyquist_axis_07(_axis_capacitor)
    plt.close(_figure)

    mo.vstack([
        mo.md(r"""
        **Predict first.** Which element should move along the Nyquist plane as
        frequency changes? Remember that an ideal resistor has no phase lag,
        while an ideal capacitor has a \(-90^\circ\) impedance phase.
        """),
        _figure,
        mo.md(r"""
        **Figure takeaway.** Frequency leaves an ideal resistor at one real-axis
        point, but moves an ideal capacitor along the negative-imaginary axis.
        Putting them in parallel bends these two limiting responses into the
        semicircle explored next.
        """),
    ])
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
    format_nyquist_axis_07,
    mo,
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
    _axis_nyquist.plot(rc_spectrum_07.real, -rc_spectrum_07.imag, lw=1.9, color="#4C7C86")
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
    format_nyquist_axis_07(_axis_nyquist)
    _axis_nyquist.legend(loc="best")

    _axis_bode.semilogx(rc_frequency_07, np.abs(rc_spectrum_07), lw=1.9, color="#4C7C86")
    _axis_bode_phase = _axis_bode.twinx()
    _axis_bode_phase.semilogx(
        rc_frequency_07, np.angle(rc_spectrum_07, deg=True), lw=1.7, color="#B8734A"
    )
    _axis_bode.set(xlabel="Frequency (Hz)", ylabel=r"$|Z|$ ($\Omega$)", title="Bode: magnitude and phase retain frequency")
    _axis_bode_phase.set_ylabel(r"Phase of $Z$ (degrees)", color="#B8734A")
    _axis_bode_phase.tick_params(axis="y", colors="#B8734A")
    _axis_bode.grid(True)
    _axis_bode.text(
        0.03,
        0.10,
        rf"$C_1={1e6 * rc_capacitance_1_07:.1f}\ \mu$F"
        + (("\n" + rf"$C_2={1e6 * rc_capacitance_2_07:.1f}\ \mu$F") if rc_show_second_07.value else ""),
        transform=_axis_bode.transAxes,
        color="#526173",
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.82},
    )
    mo.vstack([
        _figure,
        mo.md(r"**Figure takeaway.** Each ideal $R\parallel C$ branch contributes one relaxation time; Nyquist shows shape, while Bode plots preserve the frequency location needed to interpret it."),
    ])
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
    \mathrm{i}\omega\widehat{\Delta c}=D^\delta\frac{d^2\widehat{\Delta c}}{dx^2}.
    $$

    Here $D^\delta$ is the same **chemical diffusivity** used in Modules 05 and
    06. The dc and ac views are the same scaling idea:

    $$\ell_{D,\rm dc}\sim\sqrt{D^\delta t},\qquad
    \ell_{D,\rm ac}\sim\sqrt{D^\delta/\omega}.$$

    High frequency samples only a thin region; low frequency can reach the far
    boundary. For a semi-infinite sample,

    $$
    Z_W=\frac{W}{\sqrt{\mathrm{i}\omega}}
       =\frac{W(1-\mathrm{i})}{\sqrt{2\omega}},
    $$

    so $Z'=-Z''$: the $45^\circ$ line follows from diffusion rather than being
    inserted as an equivalent-circuit element.

    With $\widetilde\omega=\omega L^2/D^\delta$, the general resistance scale is

    $$R_D=\frac{L}{zFS D^\delta}\left|\frac{\partial E}{\partial c}\right|.$$

    Here $c$ is molar concentration. In the dilute ideal limit,
    $|\partial E/\partial c|=RT/(zFc_0)$, giving
    $R_D=RTL/(z^2F^2c_0SD^\delta)$. The far-boundary cases are

    $$
    \underbrace{\widehat{\Delta c}(L)=0}_{\text{open: fixed composition}},
    \quad
    \widetilde Z_{\rm open}
      =\frac{\tanh\sqrt{\mathrm{i}\widetilde\omega}}{\sqrt{\mathrm{i}\widetilde\omega}},
    $$

    $$
    \underbrace{d\widehat{\Delta c}/dx|_L=0}_{\text{blocked / FSW: zero flux}},
    \quad
    \widetilde Z_{\rm blocked}
      =\frac{\coth\sqrt{\mathrm{i}\widetilde\omega}}{\sqrt{\mathrm{i}\widetilde\omega}}.
    $$

    **Open boundary does not mean electrical open circuit.** It means a
    composition reservoir is maintained at the far face. The blocked case is
    also called finite-space Warburg (FSW).
    """)
    return


@app.cell
def _(mo):
    mo.md(r"""
    **Predict before moving the controls.** At high frequency the concentration
    wave should penetrate only a short distance and should not feel the far
    boundary. At low frequency, decide whether a fixed-composition reservoir or
    a zero-flux wall should produce the larger accumulation.
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
    dc_log_time_07 = mo.ui.slider(
        start=-2.0, stop=6.0, step=0.25, value=2.0,
        label=r"dc: $\log_{10}(t/\mathrm{s})$",
    )
    ac_log_frequency_07 = mo.ui.slider(
        start=-4.0, stop=5.0, step=0.25, value=-1.0,
        label=r"ac: $\log_{10}(f/\mathrm{Hz})$",
    )
    mo.vstack(
        [
            mo.hstack(
                [warburg_boundary_07, warburg_log_omega_07],
                justify="start", gap=1.4,
            ),
            mo.hstack(
                [warburg_log_diffusivity_07, warburg_length_07, dc_log_time_07, ac_log_frequency_07],
                justify="start", gap=1.4,
            ),
        ]
    )
    return (
        ac_log_frequency_07,
        dc_log_time_07,
        warburg_boundary_07,
        warburg_length_07,
        warburg_log_diffusivity_07,
        warburg_log_omega_07,
    )


@app.cell
def _(
    ac_log_frequency_07,
    dc_log_time_07,
    mo,
    np,
    plt,
    warburg_length_07,
    warburg_log_diffusivity_07,
):
    _diffusivity_m2_s = 10.0 ** warburg_log_diffusivity_07.value * 1.0e-4
    _length_m = warburg_length_07.value * 1.0e-6
    _time_s = 10.0 ** dc_log_time_07.value
    _frequency_hz = 10.0 ** ac_log_frequency_07.value
    _omega = 2.0 * np.pi * _frequency_hz
    _dc_length = np.sqrt(_diffusivity_m2_s * _time_s)
    _ac_length = np.sqrt(_diffusivity_m2_s / _omega)
    _ratios = np.array([_dc_length, _ac_length]) / _length_m
    _visible = np.minimum(_ratios, 1.0)

    _figure, _axis = plt.subplots(figsize=(11.8, 3.7), constrained_layout=True)
    _axis.barh(
        [1.0, 0.0], [1.0, 1.0], height=0.34,
        color="#E7EAED", edgecolor="#B8C0C8", label="specimen thickness",
    )
    _axis.barh(
        [1.0, 0.0], _visible, height=0.34,
        color=["#6B86A5", "#B8734A"], alpha=0.82,
    )
    for _y, _ratio in zip([1.0, 0.0], _ratios):
        _axis.text(
            min(_ratio, 1.0) + 0.025, _y,
            rf"$\ell_D/L={_ratio:.3g}$" + (" (far face reached)" if _ratio >= 1.0 else ""),
            va="center", color="#526173",
        )
    _axis.set(
        xlim=(0.0, 1.43),
        yticks=[0.0, 1.0],
        yticklabels=[r"ac: $\sqrt{D^\delta/\omega}$", r"dc: $\sqrt{D^\delta t}$"],
        xlabel=r"Penetration length / specimen thickness, $\ell_D/L$",
        title="DC time and AC frequency ask how far composition can respond",
    )
    _axis.grid(axis="x", alpha=0.22)
    mo.vstack([
        _figure,
        mo.md(r"""
        **Figure takeaway.** Longer dc time and lower ac
        frequency both let the composition disturbance travel farther. The far
        boundary can affect the response only when $\ell_D$ becomes comparable
        with $L$.
        """),
    ])
    return


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
    format_nyquist_axis_07,
    np,
    plt,
    warburg_phases_07,
    warburg_position_07,
    warburg_reduced_omega_07,
    warburg_selected_omega_07,
    warburg_spectra_07,
    warburg_profile_07,
):
    _semi_profile = warburg_profile_07(
        warburg_selected_omega_07, warburg_position_07, "semi-infinite"
    )
    _semi_impedance = warburg_spectra_07["semi-infinite"]
    _figure, _axes = plt.subplots(
        1, 3, figsize=(15.0, 4.8), constrained_layout=True
    )
    _profile_axis, _nyquist_axis, _bode_axis = _axes

    _profile_styles = ("-", "--", "-.", ":")
    _profile_colors = ("#4C7C86", "#B8734A", "#7C6A91", "#5F8A6B")
    for _phase, _style, _color in zip(
        warburg_phases_07, _profile_styles, _profile_colors
    ):
        _profile_axis.plot(
            warburg_position_07,
            np.real(_semi_profile * np.exp(1j * _phase)),
            color=_color,
            ls=_style,
            lw=1.8,
            label=rf"$\omega t={_phase / np.pi:.1f}\pi$",
        )
    _profile_axis.axhline(0.0, color="#73808C", lw=0.9, ls=":")
    _profile_axis.set(
        xlabel=r"Normalized position, $x/L$ (dimensionless)",
        ylabel=r"Normalized concentration, $\Delta c/|\widehat{\Delta c}(0)|$",
        title="Concentration wave decays into the solid",
    )
    _profile_axis.grid(alpha=0.24)
    _profile_axis.legend(frameon=False, fontsize=10.5, ncols=2)

    _nyquist_axis.plot(
        _semi_impedance.real,
        -_semi_impedance.imag,
        color="#4C7C86",
        lw=1.9,
        marker="o",
        markevery=55,
        ms=4,
    )
    _nyquist_axis.set(
        xlim=(-0.02, 2.2),
        ylim=(-0.02, 2.2),
        xlabel=r"Normalized real impedance, $\widetilde Z'$",
        ylabel=r"Normalized negative imaginary impedance, $-\widetilde Z''$",
        title=r"Diffusion gives $Z'=-Z''$",
    )
    format_nyquist_axis_07(_nyquist_axis)

    _bode_axis.loglog(
        warburg_reduced_omega_07,
        np.abs(_semi_impedance),
        color="#4C7C86",
        lw=1.9,
        label=r"$|\widetilde Z|$",
    )
    _phase_axis = _bode_axis.twinx()
    _phase_axis.semilogx(
        warburg_reduced_omega_07,
        np.angle(_semi_impedance, deg=True),
        color="#B8734A",
        ls="--",
        lw=1.7,
    )
    _bode_axis.set(
        xlabel=r"Reduced angular frequency, $\widetilde\omega$ (dimensionless)",
        ylabel=r"Magnitude, $|\widetilde Z|$ (dimensionless)",
        title=r"Magnitude falls as $\widetilde\omega^{-1/2}$",
    )
    _phase_axis.set_ylabel(
        r"Phase of $\widetilde Z$ (degrees)", color="#B8734A"
    )
    _phase_axis.set_ylim(-55.0, -35.0)
    _phase_axis.tick_params(axis="y", colors="#B8734A")
    _bode_axis.grid(alpha=0.24, which="both")
    plt.close(_figure)
    _figure
    return


@app.cell
def _(
    format_nyquist_axis_07,
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
    for _boundary in ("open", "blocked"):
        _impedance = warburg_spectra_07[_boundary]
        _plot_mask = (-_impedance.imag <= 4.0) & (_impedance.real <= 2.2)
        _axis_warburg.plot(
            _impedance.real[_plot_mask],
            -_impedance.imag[_plot_mask],
            lw=1.8,
            ls="-" if _boundary == "open" else "--",
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
    _axis_warburg.set(
        xlim=(-0.04, 1.55),
        ylim=(-0.04, 2.0),
        xlabel=r"$\widetilde Z'$",
        ylabel=r"$-\widetilde Z''$",
        title="Finite boundaries separate at low frequency",
    )
    format_nyquist_axis_07(_axis_warburg)
    _axis_warburg.legend(loc="upper left")

    _profile_colors = ["#4C7C86", "#B8734A", "#7C6A91", "#5F8A6B"]
    for _phase, _color in zip(warburg_phases_07, _profile_colors):
        _snapshot = np.real(warburg_selected_profile_07 * np.exp(1j * _phase))
        _axis_profile.plot(
            warburg_position_07,
            _snapshot,
            lw=1.8,
            color=_color,
            label=rf"$\omega t={_phase / np.pi:.1f}\pi$",
        )
    _axis_profile.axhline(0.0, color="#73808C", lw=0.9)
    _axis_profile.set(
        xlabel=r"Position $x/L$",
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

    **Figure takeaway.** The $45^\circ$ segment is the high-frequency
    semi-infinite asymptote. Only the low-frequency end reveals whether the far
    boundary is a fixed-composition reservoir or a zero-flux wall.

    $$
    \begin{aligned}
    \widetilde Z_{\rm open}&\longrightarrow 1-
      \frac{\mathrm{i}\widetilde\omega}{3},\\[2mm]
    \widetilde Z_{\rm blocked}&\longrightarrow
      \frac{1}{3}+\frac{1}{\mathrm{i}\widetilde\omega},
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
    and $z=1$; those fixed values affect $R_D$ but not the dimensionless curve.
    """)
    return


@app.cell
def _(mo):
    mo.accordion({
        "Advanced interpretation — why DRT needs caution": mo.md(r"""
        A finite-length diffusion response is one physical process, yet its
        $45^\circ$ region requires a distribution of relaxation times when
        represented by many $R\parallel C$ elements. A DRT calculation may
        therefore show several peaks or shoulders for this single diffusion
        problem. Conversely, a peak does not by itself identify one microscopic
        mechanism.

        Use DRT as a representation aid only after checking geometry, scaling,
        boundary conditions, and whether the proposed process reproduces both
        Nyquist and Bode behavior.
        """)
    })
    return



@app.cell
def _(mo):
    mo.md(r"""
    ## 3. Why a mixed conductor needs a transmission line


    A single lumped resistor or capacitor can reproduce a spectral shape,
    but it cannot show **where** ionic and electronic current flow, where current
    transfers between carriers, or how selective contacts change that transfer.
    A transmission line keeps this spatial information.

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
    \frac{dI_e}{dx}=-\mathrm{i}\omega c_{\rm chem}(u_e-u_i),\qquad
    \frac{dI_i}{dx}=+\mathrm{i}\omega c_{\rm chem}(u_e-u_i).
    $$

    $u_e$ and $u_i$ are **voltage-equivalent electrochemical potentials** in
    volts. Their difference drives local chemical storage. Adding the last two
    equations gives $d(I_e+I_i)/dx=0$: total current is conserved even though
    its division between carriers changes with position.


    For monovalent carriers the rail voltages are defined by

    $$u_e=-\widetilde\mu_e/F,\qquad
    u_i=+\widetilde\mu_i/F,\qquad
    u_e-u_i=-\mu_M/F.$$

    These signs make the rail difference the voltage equivalent of the neutral
    composition chemical potential. They also keep conventional electronic and
    ionic currents in the same circuit-current orientation.

    The equation-to-circuit map is direct:

    | transport statement | circuit representation |
    |---|---|
    | $du_e/dx=-r_e I_e$ | electronic rail resistance $r_e\,dx$ |
    | $du_i/dx=-r_i I_i$ | ionic rail resistance $r_i\,dx$ |
    | $\mathrm{i}\omega c_{\rm chem}(u_e-u_i)$ | local chemical storage $c_{\rm chem}\,dx$ between rails |
    | terminal carrier fluxes | contact impedances or ideal passing/blocking conditions |

    The general TLM places four contact impedances at its ends:
    $Z_A$ and $Z_B$ connect the left electrode to the electronic and ionic
    rails; $Z_C$ and $Z_D$ do the same on the right. Changing these impedances
    changes the boundary conditions without changing the bulk rails. That is
    why a TLM is a representation of transport equations and contacts rather
    than an arbitrary collection of fit elements.
    A useful time scale is

    $$
    \tau_{\rm chem}=(R_e+R_i)C_{\rm chem},\qquad
    f_{\rm chem}=\frac{1}{2\pi\tau_{\rm chem}},
    $$

    Here $r_e,r_i$ have units $\Omega$ m$^{-1}$ and $c_{\rm chem}$ has units
    F m$^{-1}$. Total quantities are $R_e=r_eL$, $R_i=r_iL$, and
    $C_{\rm chem}=c_{\rm chem}L$. The area $S$ is already included in the rail
    resistances. The resulting $f_{\rm chem}$ is a scaling frequency, not a
    promise that every boundary condition has a peak exactly there.

    A dielectric capacitance $C_{\rm diel}$ can be added for an electrolyte,
    but it is shown only as an application below and is not in this teaching model.
    """)
    return


@app.cell
def _(mo, np, plt):
    _figure, _axis = plt.subplots(figsize=(12.6, 3.2), constrained_layout=True)
    _x_nodes = np.linspace(0.12, 0.88, 6)
    _y_e, _y_i = 0.72, 0.25
    _axis.plot([0.05, 0.95], [_y_e, _y_e], color="#B8734A", lw=2.0)
    _axis.plot([0.05, 0.95], [_y_i, _y_i], color="#4C7C86", lw=2.0)
    for _node in _x_nodes:
        _axis.plot([_node, _node], [_y_i + 0.06, _y_e - 0.06], color="#7C6A91", lw=1.8)
        _axis.plot([_node - 0.018, _node + 0.018], [0.50, 0.50], color="#7C6A91", lw=2.0)
        _axis.plot([_node - 0.018, _node + 0.018], [0.46, 0.46], color="#7C6A91", lw=2.0)
    _axis.text(0.50, 0.82, r"electronic rail: $r_e\,dx$", ha="center", color="#8F5638", fontsize=15)
    _axis.text(0.50, 0.10, r"ionic rail: $r_i\,dx$", ha="center", color="#3D6972", fontsize=15)
    _axis.text(0.50, 0.53, r"distributed $c_{\rm chem}\,dx$", ha="center", color="#665777", fontsize=14)
    _axis.text(
        0.015, 0.48, "Left contacts", ha="center", va="center",
        rotation=90, color="#526173",
    )
    _axis.text(
        0.985, 0.48, "Right contacts", ha="center", va="center",
        rotation=90, color="#526173",
    )
    _axis.set(xlim=(0.0, 1.0), ylim=(0.0, 1.0), title="A continuous two-rail transmission line")
    _axis.text(0.06, 0.78, r"$Z_A$", color="#8F5638", weight="bold")
    _axis.text(0.06, 0.18, r"$Z_B$", color="#3D6972", weight="bold")
    _axis.text(0.94, 0.78, r"$Z_C$", color="#8F5638", weight="bold", ha="right")
    _axis.text(0.94, 0.18, r"$Z_D$", color="#3D6972", weight="bold", ha="right")
    _axis.text(
        0.50, 0.94, r"optional $C_{\rm diel}$ is not included here",
        ha="center", color="#73808C", fontsize=11,
    )
    _axis.axis("off")
    mo.vstack([
        _figure,
        mo.md("**Figure takeaway.** The two rails are transport equations made visible: rail resistances carry electronic and ionic current, distributed chemical capacitance stores neutral composition, and the four terminal impedances impose boundary conditions."),
    ])
    return


@app.cell
def _(mo):
    mo.md(r"""
    ### Contact presets used by the interactive line

    | preset | left face | right face | main consequence |
    |---|---|---|---|
    | electron-reversible, ions blocked | $u_e$ fixed; $I_i=0$ | $u_e$ fixed; $I_i=0$ | electronic dc path plus chemical polarization |
    | cross-selective | $u_e$ fixed; $I_i=0$ | $I_e=0$; $u_i$ fixed | no single carrier spans the specimen; blocking low-frequency response |
    | both carriers reversible | $u_e=u_i$ fixed | $u_e=u_i$ fixed | no chemical rail difference; $R_e\parallel R_i$ |

    **Prediction.** Before changing the contact preset, decide whether a dc
    path remains and whether the contacts can drive $u_e-u_i$. Those two
    questions predict the low-frequency response more reliably than naming the
    shape afterward.
    """)
    return


@app.cell
def _(mo):
    tlm_contact_case_07 = mo.ui.dropdown(
        options={
            "Electron-reversible contacts; ions blocked at both faces": "electron contacts, ions blocked",
            "Cross-selective: electron contact left; ion contact right": "cross-selective contacts",
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
    mo,
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
        lw=1.9,
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
    _axis_nyquist.set(
        ylim=(-0.03 * _nyquist_cap, _nyquist_cap),
        xlabel=r"$Z'/R_\parallel$",
        ylabel=r"$-Z''/R_\parallel$",
        title=tlm_contact_case_07.value.capitalize()
        + "\nNyquist response selected by the contacts",
    )
    format_nyquist_axis_07(_axis_nyquist)
    if _visible[_selected_index]:
        _axis_nyquist.legend(loc="best")
    _axis_bode.loglog(
        tlm_frequency_ratio_07,
        np.abs(_normalized_impedance),
        lw=1.9,
        color="#4C7C86",
    )
    _axis_phase = _axis_bode.twinx()
    _axis_phase.semilogx(
        tlm_frequency_ratio_07,
        np.angle(_normalized_impedance, deg=True),
        lw=1.7,
        color="#B8734A",
    )
    _axis_bode.axvline(1.0, color="#9AA3AB", lw=1.2, ls="--")
    _axis_bode.set(
        xlabel=r"$f/f_{\rm chem}=\omega\tau_{\rm chem}$",
        ylabel=r"$|Z|/R_\parallel$",
        title="Frequency keeps the time-scale information",
    )
    _axis_phase.set_ylabel(r"Phase of $Z$ (degrees)", color="#B8734A")
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
    mo.vstack([
        _figure,
        mo.md("**Figure takeaway.** The same MIEC interior produces conducting, polarized, or blocking low-frequency behavior when only the contact boundary conditions change."),
    ])
    return


@app.cell
def _(
    mo,
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
    _axis_potential.plot(tlm_position_07, _u_e_snapshot, lw=1.9, color="#B8734A", label=r"$u_e$")
    _axis_potential.plot(tlm_position_07, _u_i_snapshot, lw=1.9, color="#4C7C86", label=r"$u_i$")
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
        xlabel=r"Position $x/L$",
        ylabel="Instantaneous voltage-equivalent potential (V)",
        title="Potentials show where chemical storage is driven",
    )
    _axis_potential.grid(True)
    _axis_potential.legend(loc="best")

    _axis_current.plot(tlm_position_07, _i_e_snapshot, lw=1.9, color="#B8734A", label=r"$I_e$")
    _axis_current.plot(tlm_position_07, _i_i_snapshot, lw=1.9, color="#4C7C86", label=r"$I_i$")
    _axis_current.plot(
        tlm_position_07,
        _i_total_snapshot,
        lw=1.7,
        ls="--",
        color="#5F8A6B",
        label=r"$I_e+I_i$",
    )
    _axis_current.axhline(0.0, color="#73808C", lw=0.9)
    _axis_current.set(
        xlabel=r"Position $x/L$",
        ylabel="Instantaneous current / common scale",
        title="Carrier currents exchange, but their sum stays constant",
    )
    _axis_current.grid(True)
    _axis_current.legend(loc="best")
    _figure.suptitle(
        rf"Internal TLM state at $f/f_{{\rm chem}}={tlm_selected_omega_07:.3g}$ and $\omega t={tlm_profile_phase_07.value}^\circ$",
        fontsize=17,
    )
    mo.vstack([
        _figure,
        mo.md(r"**Figure takeaway.** Chemical storage is driven by $u_e-u_i=-\mu_M/F$; current can transfer between rails with position even though $I_e+I_i$ remains constant."),
    ])
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
            "and capacitive. It is the mirror-coordinate version of the selective-contact "
            "geometry displayed in Module 06 and in many chemical-diffusion measurements."
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
def _(mo):
    _scope = mo.callout(
        mo.md(r"""
        **Scope of the interactive TLM model.** It includes uniform,
        one-dimensional electronic and ionic rails, distributed chemical
        capacitance, and three ideal contact presets. It does **not** fit finite
        $Z_A$–$Z_D$, dielectric/stray capacitance, surface reaction impedance,
        constant-phase elements, microstructural distributions, or nonlinear
        large-signal response. Those belong in a model only when the experiment
        supplies evidence for them.
        """),
        kind="warn",
    )
    _applications = mo.md(r"""
    ### Three application cards for the same TLM anatomy

    **1. Warburg measurement in an electron-dominant MIEC**

    Let $r_e\ll r_i$ and use one ion-selective face. Chemical capacitance and
    ionic transport set the diffusion response; changing the far contact selects
    the fixed-composition or blocked finite-length limit.

    **2. Solid electrolyte between ion-blocking metal electrodes**

    Let $r_e\gg r_i$. Chemical capacitance can become small enough that the
    dielectric specimen capacitance $C_{\rm diel}$ and stray capacitance matter.
    The present teaching model does not add them, but the general anatomy shows where
    they enter.

    **3. Thin MIEC electrode with surface exchange**

    If both rail resistances are small compared with a surface reaction
    resistance, the body equilibrates nearly uniformly. Surface reaction
    resistance can then appear in parallel with the total
    $C_{\rm chem}=c_{\rm chem}L$, explaining why a surface process can be paired
    with a volume-scaling capacitance.

    **Takeaway.** These are controlled simplifications of one transport model,
    not three unrelated equivalent circuits.
    """)
    mo.vstack([_scope, _applications])
    return



@app.cell
def _(
    capacitor_impedance_07,
    format_nyquist_axis_07,
    np,
    plt,
    rc_impedance_07,
    resistor_impedance_07,
    tlm_parameters_07,
    tlm_profile_07,
    tlm_solution_07,
    tlm_spectrum_07,
    warburg_impedance_07,
    warburg_scales_07,
    waveform_data_07,
):
    _check_frequency = np.logspace(-5.0, 5.0, 800)
    _check_resistor = resistor_impedance_07(_check_frequency, 37.0)
    resistor_limit_error_07 = max(
        np.max(np.abs(_check_resistor.real - 37.0)),
        np.max(np.abs(_check_resistor.imag)),
    )
    _check_capacitor = capacitor_impedance_07(_check_frequency, 2.5e-6)
    capacitor_limit_error_07 = max(
        np.max(np.abs(_check_capacitor.real)),
        np.max(np.abs(
            (-_check_capacitor.imag)
            * (2.0 * np.pi * _check_frequency * 2.5e-6)
            - 1.0
        )),
    )
    _nyquist_test_figure, _nyquist_test_axis = plt.subplots()
    format_nyquist_axis_07(_nyquist_test_axis)
    nyquist_equal_axis_error_07 = abs(float(_nyquist_test_axis.get_aspect()) - 1.0)
    plt.close(_nyquist_test_figure)
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
    tlm_voltage_mapping_error_07 = np.max(np.abs(
        _tlm_check_profile["u_e"]
        - _tlm_check_profile["u_i"]
        - _tlm_check_profile["u_chemical"]
    ))
    _length_check_m = 2.5e-4
    _r_e_per_m = _tlm_check_parameters["R_e_ohm"] / _length_check_m
    _c_per_m = _tlm_check_parameters["C_chemical_f"] / _length_check_m
    tlm_total_distributed_error_07 = max(
        abs(_r_e_per_m * _length_check_m / _tlm_check_parameters["R_e_ohm"] - 1.0),
        abs(_c_per_m * _length_check_m / _tlm_check_parameters["C_chemical_f"] - 1.0),
    )
    _scale_check = warburg_scales_07(120.0, 2.0e-8, 250.0, 0.8, 800.0, 1.0)
    warburg_resistance_scale_error_07 = abs(
        _scale_check["resistance_general_ohm"]
        / _scale_check["resistance_diffusion_ohm"]
        - 1.0
    )
    _dc_time = 3.7
    _dc_ac_diffusivity = 4.2e-11
    dc_ac_length_error_07 = abs(
        np.sqrt(_dc_ac_diffusivity * _dc_time)
        / np.sqrt(_dc_ac_diffusivity / (1.0 / _dc_time))
        - 1.0
    )
    return (
        capacitor_limit_error_07,
        dc_ac_length_error_07,
        nyquist_equal_axis_error_07,
        resistor_limit_error_07,
        phasor_sign_error_07,
        rc_circle_error_07,
        rc_peak_error_07,
        tlm_boundary_error_07,
        tlm_current_conservation_error_07,
        tlm_finiteness_07,
        tlm_passivity_margin_07,
        tlm_reversible_error_07,
        tlm_total_distributed_error_07,
        tlm_voltage_mapping_error_07,
        warburg_high_frequency_error_07,
        warburg_low_frequency_error_07,
        warburg_passivity_margin_07,
        warburg_resistance_scale_error_07,
    )


@app.cell
def _(
    capacitor_limit_error_07,
    dc_ac_length_error_07,
    mo,
    nyquist_equal_axis_error_07,
    resistor_limit_error_07,
    phasor_sign_error_07,
    rc_circle_error_07,
    rc_peak_error_07,
    tlm_boundary_error_07,
    tlm_current_conservation_error_07,
    tlm_finiteness_07,
    tlm_passivity_margin_07,
    tlm_reversible_error_07,
    tlm_total_distributed_error_07,
    tlm_voltage_mapping_error_07,
    warburg_high_frequency_error_07,
    warburg_low_frequency_error_07,
    warburg_passivity_margin_07,
    warburg_resistance_scale_error_07,
):
    _checks = [
        (
            "Ideal resistor and capacitor limits",
            max(resistor_limit_error_07, capacitor_limit_error_07) < 1.0e-12,
            r"The resistor point and capacitive line must follow directly from $Z_R=R$ and $Z_C=1/(\mathrm{i}\omega C)$.",
        ),
        (
            "Equal Nyquist scaling",
            nyquist_equal_axis_error_07 < 1.0e-14,
            "Equal data scale keeps angles, semicircles, and diffusion slopes geometrically honest.",
        ),
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
            "DC/AC diffusion-length scaling",
            dc_ac_length_error_07 < 1.0e-14,
            r"Choosing $\omega=1/t$ must give the same penetration length in time and frequency views.",
        ),
        (
            "General and dilute Warburg resistance scales",
            warburg_resistance_scale_error_07 < 1.0e-14,
            "The general thermodynamic slope must reduce to the ideal dilute expression under the stated assumption.",
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
            "TLM voltage and unit mappings",
            max(tlm_voltage_mapping_error_07, tlm_total_distributed_error_07) < 1.0e-12,
            "$u_e-u_i$ must be the stored chemical voltage, and distributed quantities must recover their totals after multiplication by $L$.",
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
    _checks_table = mo.md("\n".join(_heading + _rows))
    mo.accordion({"Physical consistency checks": _checks_table})
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## 5. Three messages to keep

    $$
    \boxed{
    \text{small sinusoid}
    \rightarrow Z(\omega)
    \rightarrow \text{time scale}
    \rightarrow \text{transport + storage + boundaries}
    }
    $$

    1. **Frequency selects a time and length scale.** A semicircle marks a
       relaxation; Warburg behavior appears when chemical diffusion is the
       distributed response.
    2. **The far boundary controls the low-frequency end.** Fixed composition
       gives a finite resistance; zero flux gives capacitive accumulation. Both
       share the high-frequency $45^\circ$ limit.
    3. **A TLM connects physics, storage, and contacts.** The same MIEC rails can
       produce different spectra when $Z_A$–$Z_D$ change, so geometry, sign
       convention, units, and boundary conditions must precede feature labels.

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


@app.cell
def _(mo):
    mo.md(r"""
    **Reference:** [shared notation and sign conventions](https://github.com/qiyanglu/Solid-State-Ionics-Interactive/blob/main/NOTATION.md)
    """)
    return



if __name__ == "__main__":
    app.run()
