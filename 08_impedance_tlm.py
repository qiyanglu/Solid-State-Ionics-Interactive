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
    from matplotlib.patches import Rectangle
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
            "lines.solid_capstyle": "round",
            "figure.dpi": 115,
        }
    )

    return Rectangle, mo, np, plt


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
    # Module 08 · Impedance spectroscopy: from a sine wave to a transport model

    **What can respond at a chosen frequency?**

    Electrochemical impedance spectroscopy (EIS) applies a small sinusoidal
    voltage and measures the sinusoidal current. Fast processes can follow a
    rapid oscillation; slower transport appears only as the frequency is
    lowered. An impedance spectrum is therefore a map of **dynamics**, not just
    a collection of arcs to label.

    We will first learn the complex-number language, then watch one-dimensional
    chemical diffusion create Warburg impedance, and finally connect that
    picture to the two-rail transmission line of a mixed conductor.


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

    def waveform_data_07(frequency_hz, current_lead_deg, sample_count=1200):
        """Return normalized voltage and current over a fixed two-second window."""
        frequency = _positive_07("frequency_hz", frequency_hz)
        lead_rad = np.deg2rad(float(current_lead_deg))
        time_s = np.linspace(0.0, 2.0, int(sample_count))
        phase = 2.0 * np.pi * frequency * time_s
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

    def series_rc_impedance_07(frequency_hz, resistance_ohm, capacitance_f):
        """Return the impedance of an ideal resistor and capacitor in series."""
        return resistor_impedance_07(
            frequency_hz, resistance_ohm
        ) + capacitor_impedance_07(frequency_hz, capacitance_f)

    def parallel_rc_impedance_07(frequency_hz, resistance_ohm, capacitance_f):
        """Return the impedance of an ideal resistor and capacitor in parallel."""
        frequencies = np.asarray(frequency_hz, dtype=float)
        resistance = _positive_07("resistance_ohm", resistance_ohm)
        capacitance = _positive_07("capacitance_f", capacitance_f)
        return resistance / (
            1.0 + 1j * 2.0 * np.pi * frequencies * resistance * capacitance
        )

    def set_equal_nyquist_limits(axis, real_values, minus_imaginary_values, margin=0.08):
        """Set equal numerical spans with data margins on every Nyquist axis."""
        real = np.asarray(real_values, dtype=float).ravel()
        minus_imaginary = np.asarray(minus_imaginary_values, dtype=float).ravel()
        finite = np.isfinite(real) & np.isfinite(minus_imaginary)
        if not np.any(finite):
            raise ValueError("Nyquist data must contain at least one finite point")
        real = real[finite]
        minus_imaginary = minus_imaginary[finite]
        x_min, x_max = float(real.min()), float(real.max())
        y_min, y_max = float(minus_imaginary.min()), float(minus_imaginary.max())
        absolute_scale = max(
            float(np.max(np.abs(real))),
            float(np.max(np.abs(minus_imaginary))),
            1.0,
        )
        raw_span = max(x_max - x_min, y_max - y_min)
        data_span = (
            max(0.25 * absolute_scale, 1.0e-6)
            if raw_span <= 1.0e-9 * absolute_scale
            else raw_span
        )
        full_span = (1.0 + 2.0 * float(margin)) * data_span
        near_zero_real = abs(x_min) <= max(0.02 * data_span, 1.0e-10 * absolute_scale)
        if near_zero_real:
            x_lower = min(x_min - margin * data_span, -0.04 * data_span)
            x_limits = (x_lower, x_lower + full_span)
        else:
            x_center = 0.5 * (x_min + x_max)
            x_limits = (x_center - 0.5 * full_span, x_center + 0.5 * full_span)
        y_center = 0.5 * (y_min + y_max)
        y_limits = (y_center - 0.5 * full_span, y_center + 0.5 * full_span)
        axis.set_xlim(*x_limits)
        axis.set_ylim(*y_limits)
        axis.set_aspect("equal", adjustable="box")
        return x_limits, y_limits

    def format_nyquist_axis_07(axis):
        """Apply labels-independent Nyquist geometry and grid styling."""
        axis.set_aspect("equal", adjustable="box")
        axis.ticklabel_format(axis="both", style="plain", useOffset=False)
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
            if q_value.real > 20.0:
                return (
                    np.exp(-q_value * xi)
                    * (1.0 - np.exp(-2.0 * q_value * (1.0 - xi)))
                    / (1.0 - np.exp(-2.0 * q_value))
                )
            return np.sinh(q_value * (1.0 - xi)) / np.sinh(q_value)
        if boundary == "blocked":
            if q_value.real > 20.0:
                return (
                    np.exp(-q_value * xi)
                    * (1.0 + np.exp(-2.0 * q_value * (1.0 - xi)))
                    / (1.0 + np.exp(-2.0 * q_value))
                )
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
        if abs(k_value) < 1.0e-1:
            k_squared = k_value**2
            mean_factor = k_squared / 2.0 - k_squared**2 / 24.0
            difference_factor = 1.0 + k_squared / 12.0 - k_squared**2 / 720.0
        else:
            tanh_half = _stable_tanh_positive_real_07(
                np.array([k_value / 2.0])
            )[0]
            mean_factor = k_value * tanh_half
            difference_factor = (k_value / 2.0) / tanh_half
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

    def tlm_distributed_parameters_07(
        parallel_resistance_ohm_per_m,
        chemical_capacitance_f_per_m,
        length_um,
        conductivity_ratio,
    ):
        """Convert the distributed teaching controls to whole-sample values."""
        length_m = _positive_07("length_um", length_um) * 1.0e-6
        r_parallel = _positive_07(
            "parallel_resistance_ohm_per_m", parallel_resistance_ohm_per_m
        )
        c_chemical = _positive_07(
            "chemical_capacitance_f_per_m", chemical_capacitance_f_per_m
        )
        parameters = tlm_parameters_07(
            r_parallel * length_m,
            c_chemical * length_m,
            conductivity_ratio,
        )
        r_e_per_m = parameters["R_e_ohm"] / length_m
        r_i_per_m = parameters["R_i_ohm"] / length_m
        diffusivity = 1.0 / ((r_e_per_m + r_i_per_m) * c_chemical)
        parameters.update(
            {
                "length_m": length_m,
                "r_parallel_ohm_per_m": r_parallel,
                "c_chemical_f_per_m": c_chemical,
                "r_e_ohm_per_m": r_e_per_m,
                "r_i_ohm_per_m": r_i_per_m,
                "D_chemical_m2_per_s": diffusivity,
            }
        )
        return parameters

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
        row_norms = np.linalg.norm(boundary_matrix, axis=1)
        row_scaled_matrix = boundary_matrix / row_norms[:, np.newaxis]
        row_scaled_target = target / row_norms
        column_norms = np.linalg.norm(row_scaled_matrix, axis=0)
        balanced_matrix = row_scaled_matrix / column_norms[np.newaxis, :]
        balanced_coefficients = np.linalg.solve(balanced_matrix, row_scaled_target)
        coefficients = balanced_coefficients / column_norms
        total_current = coefficients[1] / parameters["R_parallel_ohm"]
        impedance = 1.0 / total_current
        balanced_residual = balanced_matrix @ balanced_coefficients - row_scaled_target
        residual_scale = (
            np.linalg.norm(balanced_matrix, ord=np.inf)
            * np.linalg.norm(balanced_coefficients, ord=np.inf)
            + np.linalg.norm(row_scaled_target, ord=np.inf)
        )
        residual = np.linalg.norm(balanced_residual, ord=np.inf) / max(
            residual_scale, np.finfo(float).tiny
        )
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
        if k_value.real > 20.0:
            exp_from_right = np.exp(k_value * (zeta - 0.5))
            exp_from_left = np.exp(-k_value * (zeta + 0.5))
            exp_across = np.exp(-k_value)
            cosh_over_cosh = (
                exp_from_right + exp_from_left
            ) / (1.0 + exp_across)
            sinh_over_cosh = (
                exp_from_right - exp_from_left
            ) / (1.0 + exp_across)
            sinh_over_two_sinh = 0.5 * (
                exp_from_right - exp_from_left
            ) / (1.0 - exp_across)
            cosh_over_sinh = (
                exp_from_right + exp_from_left
            ) / (1.0 - exp_across)
            chemical_voltage = (
                mean_value * cosh_over_cosh
                - difference_value * sinh_over_two_sinh
            )
            weighted_current = (
                -mean_value * k_value * sinh_over_cosh
                + difference_value * (k_value / 2.0) * cosh_over_sinh
            ) / r_sum
        else:
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
        parallel_rc_impedance_07,
        rc_impedance_07,
        resistor_impedance_07,
        series_rc_impedance_07,
        set_equal_nyquist_limits,
        tlm_distributed_parameters_07,
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

    **Prediction for Figure 1.** Raising $f$ shortens the period in seconds,
    but it does not change the phase angle between the two phasors. Which part
    of the figure should therefore move when you change frequency?
    """)
    return


@app.cell
def _(mo):
    waveform_frequency_07 = mo.ui.slider(
        start=-0.3, stop=1.0, step=0.1, value=0.0, label="Wave-frequency exponent"
    )
    waveform_lead_07 = mo.ui.slider(
        start=0, stop=90, step=5, value=45, label="Current lead (degrees)"
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
def _(mo, np, plt, waveform_frequency_07, waveform_result_07):
    _time_s, _voltage, _current, _z_phase = waveform_result_07
    _frequency_hz = 10.0 ** waveform_frequency_07.value
    _period_s = 1.0 / _frequency_hz
    _figure, (_axis_wave, _axis_phasor) = plt.subplots(
        1, 2, figsize=(12.5, 4.4), constrained_layout=True
    )
    _axis_wave.plot(_time_s, _voltage, lw=1.9, color="#4C7C86", label=r"$\Delta V/V_a$")
    _axis_wave.plot(_time_s, _current, lw=1.9, color="#B8734A", label=r"$\Delta I/I_a$")
    _axis_wave.axhline(0.0, color="#73808C", lw=0.9)
    _axis_wave.set(
        xlabel=r"Time, $t$ (s)",
        xlim=(0.0, 2.0),
        ylabel="Normalized signal (dimensionless)",
        title=rf"$f={_frequency_hz:.3g}$ Hz, period $={_period_s:.3g}$ s",
    )
    _axis_wave.ticklabel_format(axis="x", style="sci", scilimits=(-2, 3))
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
        mo.md(r"The two-second window stays fixed, so frequency now changes the visible number of cycles. The selected phase offset remains the same, and the phasors store that relation without carrying a time axis."),
    ])
    return


@app.cell
def _(mo):
    mo.md(r"""
    ### The same resistor and capacitor can make two different Nyquist shapes

    | connection | impedance | Nyquist signature |
    |---|---:|---|
    | resistor | $Z_R=R$ | point on the real axis |
    | capacitor | $Z_C=1/(\mathrm{i}\omega C)$ | vertical capacitive line |
    | series RC | $Z=R+1/(\mathrm{i}\omega C)$ | vertical line at $Z'=R$ |
    | parallel $R\parallel C$ | $Z=R/(1+\mathrm{i}\omega RC)$ | semicircle |

    In series, the same current passes through both elements and the capacitor
    adds a vertical contribution to the fixed real resistance. In parallel,
    the current can divide between the two paths and the spectrum bends into a
    semicircle. In both cases, $\tau=RC$ and $\omega\tau=1$ marks the crossover.

    **Predict for Figure 2.** Which connection should retain the real part
    $Z'=R$ at every frequency, and which should approach $R$ only at low
    frequency?
    """)
    return


@app.cell
def _(
    capacitor_impedance_07,
    format_nyquist_axis_07,
    mo,
    np,
    plt,
    resistor_impedance_07,
    set_equal_nyquist_limits,
):
    _frequency_hz = np.logspace(-1.0, 5.0, 500)
    resistor = resistor_impedance_07(_frequency_hz, 500.0)
    capacitor = capacitor_impedance_07(_frequency_hz, 10.0e-6)
    figure, (resistor_axis, capacitor_axis) = plt.subplots(
        1, 2, figsize=(12.4, 5.0), constrained_layout=True
    )
    resistor_axis.scatter(
        resistor.real[:1], -resistor.imag[:1], s=95,
        color="#4C7C86", edgecolor="white", zorder=4,
    )
    resistor_axis.set(
        xlabel=r"$Z'$ ($\Omega$)", ylabel=r"$-Z''$ ($\Omega$)",
        title="Ideal resistor: one point",
    )
    set_equal_nyquist_limits(
        resistor_axis,
        np.array([0.0, resistor.real[0]]),
        np.zeros(2),
    )
    format_nyquist_axis_07(resistor_axis)

    capacitor_axis.plot(
        capacitor.real, -capacitor.imag, color="#B8734A", lw=2.0
    )
    capacitor_axis.set(
        xlabel=r"$Z'$ ($\Omega$)", ylabel=r"$-Z''$ ($\Omega$)",
        title="Ideal capacitor: vertical line at $Z'=0$",
    )
    set_equal_nyquist_limits(capacitor_axis, capacitor.real, -capacitor.imag)
    format_nyquist_axis_07(capacitor_axis)
    plt.close(figure)
    mo.vstack([
        mo.md(r"""
        ### Start with the two ideal elements

        A resistor has no phase lag, so all frequencies collapse to one real-axis
        point. A capacitor has zero real impedance and a frequency-dependent
        negative imaginary part, so it traces a vertical Nyquist line. The small
        negative real-axis margin keeps that line visible instead of hiding it
        under the plot spine.
        """),
        figure,
    ])
    return (figure,)


@app.cell
def _(
    format_nyquist_axis_07,
    mo,
    np,
    parallel_rc_impedance_07,
    plt,
    series_rc_impedance_07,
    set_equal_nyquist_limits,
):
    _frequency_hz = np.logspace(-1.0, 5.0, 600)
    resistance_ohm = 500.0
    capacitance_f = 10.0e-6
    tau_s = resistance_ohm * capacitance_f
    z_parallel = parallel_rc_impedance_07(_frequency_hz, resistance_ohm, capacitance_f)
    z_series = series_rc_impedance_07(_frequency_hz, resistance_ohm, capacitance_f)
    apex_frequency_hz = 1.0 / (2.0 * np.pi * tau_s)
    apex_index = int(np.argmin(np.abs(np.log(_frequency_hz / apex_frequency_hz))))

    parallel_figure, parallel_axis = plt.subplots(figsize=(8.2, 5.4), constrained_layout=True)
    parallel_axis.plot(z_parallel.real, -z_parallel.imag, color="#4C7C86", lw=2.0)
    parallel_axis.scatter(
        z_parallel.real[apex_index], -z_parallel.imag[apex_index],
        s=90, color="#C49345", edgecolor="white", zorder=4,
        label=r"$\omega RC=1$",
    )
    parallel_axis.set(
        xlabel=r"Real impedance, $Z'$ ($\Omega$)",
        ylabel=r"Negative imaginary impedance, $-Z''$ ($\Omega$)",
        title=r"A parallel $R\parallel C$ branch gives a semicircle",
    )
    set_equal_nyquist_limits(parallel_axis, z_parallel.real, -z_parallel.imag)
    format_nyquist_axis_07(parallel_axis)
    parallel_axis.legend(loc="best")
    plt.close(parallel_figure)

    series_visible = -z_series.imag <= np.percentile(-z_series.imag, 75.0)
    series_figure, series_axis = plt.subplots(figsize=(6.8, 5.0), constrained_layout=True)
    series_axis.plot(
        z_series.real[series_visible], -z_series.imag[series_visible],
        color="#B8734A", lw=2.0,
    )
    series_axis.set(
        xlabel=r"$Z'$ ($\Omega$)", ylabel=r"$-Z''$ ($\Omega$)",
        title=r"Series RC keeps $Z'=R$",
    )
    set_equal_nyquist_limits(
        series_axis, z_series.real[series_visible], -z_series.imag[series_visible]
    )
    format_nyquist_axis_07(series_axis)
    plt.close(series_figure)

    mo.vstack([
        parallel_figure,
        mo.md(r"""
        In parallel, current divides between the resistive and capacitive paths.
        The arc top occurs at $\omega RC=1$, where the real and capacitive
        responses are equally important.
        """),
        mo.accordion({
            "Explore further — the same resistor and capacitor in series":
            mo.vstack([
                series_figure,
                mo.md(r"In series, the capacitor adds a vertical contribution while $Z'=R$ remains fixed."),
            ])
        }),
    ])
    return


@app.cell
def _(mo):
    _section_heading = mo.md(r"""
    ### When two parallel-RC relaxations overlap

    One ideal $R\parallel C$ branch gives one relaxation time. Add a second
    branch in series and predict whether the Nyquist arcs will remain distinct
    as their time constants move closer together. The Bode view retains the
    frequency information even when the arcs overlap.
    """)
    rc_series_07 = mo.ui.slider(start=0, stop=200, step=5, value=25, label=r"$R_s$ ($\Omega$)")
    rc_resistance_1_07 = mo.ui.slider(
        start=100, stop=1500, step=50, value=600, label=r"$R_1$ ($\Omega$)"
    )
    rc_log_tau_1_07 = mo.ui.slider(
        start=-4.0, stop=1.0, step=0.1, value=-1.7, label="First relaxation-time exponent"
    )
    rc_show_second_07 = mo.ui.checkbox(value=True, label="Show a second relaxation")
    rc_resistance_ratio_07 = mo.ui.slider(
        start=0.25, stop=2.0, step=0.05, value=1.4, label=r"$R_2/R_1$"
    )
    rc_log_separation_07 = mo.ui.slider(
        start=0.0,
        stop=4.0,
        step=0.1,
        value=2.0,
        label="Relaxation-time separation exponent",
    )
    rc_section_heading_07 = _section_heading
    return (
        rc_log_separation_07,
        rc_log_tau_1_07,
        rc_resistance_1_07,
        rc_resistance_ratio_07,
        rc_section_heading_07,
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
    rc_log_separation_07,
    rc_log_tau_1_07,
    rc_resistance_1_07,
    rc_resistance_ratio_07,
    rc_section_heading_07,
    rc_series_07,
    rc_show_second_07,
    rc_spectrum_07,
    rc_tau_1_07,
    rc_tau_2_07,
    set_equal_nyquist_limits,
):
    rc_primary_controls_07 = mo.hstack(
        [rc_series_07, rc_resistance_1_07, rc_log_tau_1_07],
        justify="start",
        gap=1.5,
    )
    _second_control_items = [rc_show_second_07]
    if rc_show_second_07.value:
        _second_control_items.extend(
            [rc_resistance_ratio_07, rc_log_separation_07]
        )
    rc_secondary_controls_07 = mo.hstack(
        _second_control_items,
        justify="start",
        gap=1.5,
    )
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
    set_equal_nyquist_limits(_axis_nyquist, rc_spectrum_07.real, -rc_spectrum_07.imag)
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
    mo.accordion({"Explore further — overlapping parallel-RC relaxations": mo.vstack([
        rc_section_heading_07,
        rc_primary_controls_07,
        rc_secondary_controls_07,
        _figure,
        mo.md(r"Each ideal $R\parallel C$ branch contributes one relaxation time; Nyquist shows shape, while Bode plots preserve the frequency location needed to interpret it."),
    ])})
    return


@app.cell
def _(mo):
    interpretation = mo.md(r"""
    ### A ceramic interpretation — with a caution

    In a brick-layer picture, a high-frequency arc is often
    associated with grain interiors and a lower-frequency arc with grain
    boundaries. The assignment is plausible when the two regions have distinct
    $R C$ times and their fitted capacitances scale sensibly with geometry.
    When the times overlap, the Nyquist curve need not reveal two complete
    semicircles. Frequency range, capacitance, thickness scaling, and independent
    microstructural knowledge should support the assignment; arc order alone is
    not enough.
    """)
    mo.accordion({"Explore further — assigning ceramic arcs": interpretation})
    return


@app.cell
def _(mo):
    warburg_details = mo.md(r"""
    ## 2. Warburg impedance: diffusion written in the frequency domain

    Consider a neutral composition variable $c$ in a one-dimensional slab of
    length $L$ and cross-sectional area $S$:

    $$
    \frac{\partial c}{\partial t}=D^\delta\frac{\partial^2c}{\partial x^2}
    \quad\Longrightarrow\quad
    \mathrm{i}\omega\widehat{\Delta c}=D^\delta\frac{d^2\widehat{\Delta c}}{dx^2}.
    $$

    Here $D^\delta$ is the same **chemical diffusivity** used in Modules 05 and
    06. A harmonic perturbation introduces the propagation constant
    $q=\sqrt{\mathrm{i}\omega/D^\delta}$. A larger frequency therefore makes
    the concentration wave decay more rapidly with position. At sufficiently
    low frequency the wave reaches the far face, so its boundary condition
    becomes part of the measured response. For a semi-infinite sample,

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
    mo.vstack([
        mo.md(r"""
        ## 2. Warburg impedance: when diffusion cannot keep up

        Use the same one-dimensional chemical diffusivity $D^\delta$ as in
        Modules 05 and 06. A sinusoidal surface perturbation creates a
        concentration wave inside the specimen. High frequency probes only a
        shallow region; lower frequency lets the wave travel farther.

        For a semi-infinite specimen,

        $$Z_W=\frac{W}{\sqrt{\mathrm{i}\omega}}
        =\frac{W(1-\mathrm{i})}{\sqrt{2\omega}}.$$

        Therefore $Z'=-Z''$: the familiar $45^\circ$ line is a consequence of
        diffusion, not an assumed line drawn on the Nyquist plot.
        """),
        mo.accordion({
            "Model details — diffusion equation, scales, and finite boundaries": warburg_details
        }),
    ])
    return


@app.cell
def _(mo):
    mo.md(r"""
    ### Finite length: make the far boundary visible

    **Predict before moving the controls.** Lowering $f$, shortening $L$, or
    raising $D^\delta$ all reduce $f/f_D$ and give the concentration wave more
    opportunity to reach the far face. Decide whether a fixed-composition
    reservoir or a zero-flux wall should then produce the larger accumulation.
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
        label="Far boundary condition",
    )
    warburg_log_frequency_07 = mo.ui.slider(
        start=-9.0,
        stop=3.0,
        step=0.25,
        value=-4.75,
        label="Warburg-frequency exponent",
    )
    warburg_phase_07 = mo.ui.dropdown(
        options={"0 cycle": 0.0, "1/4 cycle": 90.0, "1/2 cycle": 180.0, "3/4 cycle": 270.0},
        value="0 cycle",
        label="Snapshot phase in cycle",
    )
    warburg_compare_boundaries_07 = mo.ui.checkbox(
        value=False, label="Compare both far boundaries"
    )
    warburg_log_diffusivity_07 = mo.ui.slider(
        start=-12.0,
        stop=-5.0,
        step=0.25,
        value=-8.0,
        label="Chemical-diffusivity exponent",
    )
    warburg_length_07 = mo.ui.slider(
        start=10, stop=500, step=10, value=100, label="Diffusion length (micrometers)"
    )
    warburg_temperature_07 = mo.ui.slider(
        start=400, stop=1400, step=50, value=800, label=r"$T$ (K)"
    )
    warburg_log_concentration_07 = mo.ui.slider(
        start=0.0,
        stop=5.0,
        step=0.25,
        value=2.25,
        label=r"$\log_{10}(c_0/\mathrm{mol\,m^{-3}})$",
    )
    warburg_area_07 = mo.ui.slider(
        start=0.1, stop=2.0, step=0.1, value=0.5, label=r"$S$ (cm$^2$)"
    )
    core_controls_07 = mo.hstack(
        [warburg_boundary_07, warburg_log_frequency_07, warburg_log_diffusivity_07],
        justify="start", align="center", wrap=True, gap=1.4,
    )
    advanced_controls_07 = mo.vstack([
        mo.hstack(
            [warburg_phase_07, warburg_length_07, warburg_compare_boundaries_07],
            justify="start", align="center", wrap=True, gap=1.4,
        ),
        mo.hstack(
            [warburg_temperature_07, warburg_log_concentration_07, warburg_area_07],
            justify="start", align="center", wrap=True, gap=1.4,
        ),
    ])
    mo.vstack([
        core_controls_07,
        mo.accordion({"Explore further — phase, thickness, and impedance scale": advanced_controls_07}),
    ])
    return (
        warburg_area_07,
        warburg_boundary_07,
        warburg_compare_boundaries_07,
        warburg_length_07,
        warburg_log_concentration_07,
        warburg_log_diffusivity_07,
        warburg_log_frequency_07,
        warburg_phase_07,
        warburg_temperature_07,
    )


@app.cell
def _(
    np,
    warburg_area_07,
    warburg_boundary_07,
    warburg_impedance_07,
    warburg_length_07,
    warburg_log_concentration_07,
    warburg_log_diffusivity_07,
    warburg_log_frequency_07,
    warburg_phase_07,
    warburg_profile_07,
    warburg_scales_07,
    warburg_temperature_07,
):
    warburg_scale_data_07 = warburg_scales_07(
        warburg_length_07.value,
        10.0 ** warburg_log_diffusivity_07.value,
        concentration_mol_per_m3=10.0 ** warburg_log_concentration_07.value,
        area_cm2=warburg_area_07.value,
        temperature_k=warburg_temperature_07.value,
        charge_number=1.0,
    )
    warburg_selected_frequency_hz_07 = 10.0 ** warburg_log_frequency_07.value
    warburg_selected_omega_07 = (
        warburg_selected_frequency_hz_07
        / warburg_scale_data_07["frequency_diffusion_hz"]
    )
    _log_f_d = np.log10(warburg_scale_data_07["frequency_diffusion_hz"])
    _log_f_min = min(_log_f_d - 4.0, warburg_log_frequency_07.value - 0.25)
    _log_f_max = max(_log_f_d + 4.0, warburg_log_frequency_07.value + 0.25)
    warburg_frequency_hz_07 = np.logspace(_log_f_min, _log_f_max, 600)
    warburg_reduced_omega_07 = (
        warburg_frequency_hz_07
        / warburg_scale_data_07["frequency_diffusion_hz"]
    )
    _resistance_scale = warburg_scale_data_07["resistance_diffusion_ohm"]
    warburg_spectra_07 = {
        _boundary: _resistance_scale
        * warburg_impedance_07(warburg_reduced_omega_07, _boundary)
        for _boundary in ("open", "blocked")
    }
    warburg_selected_impedance_07 = _resistance_scale * warburg_impedance_07(
        np.array([warburg_selected_omega_07]), warburg_boundary_07.value
    )[0]
    warburg_position_07 = np.linspace(0.0, 1.0, 300)
    warburg_position_um_07 = warburg_length_07.value * warburg_position_07
    warburg_selected_profile_07 = warburg_profile_07(
        warburg_selected_omega_07, warburg_position_07, warburg_boundary_07.value
    )
    warburg_phase_radians_07 = np.deg2rad(warburg_phase_07.value)
    return (
        warburg_frequency_hz_07,
        warburg_phase_radians_07,
        warburg_position_07,
        warburg_position_um_07,
        warburg_reduced_omega_07,
        warburg_scale_data_07,
        warburg_selected_frequency_hz_07,
        warburg_selected_impedance_07,
        warburg_selected_omega_07,
        warburg_selected_profile_07,
        warburg_spectra_07,
    )


@app.cell
def _(
    format_nyquist_axis_07,
    mo,
    np,
    plt,
    set_equal_nyquist_limits,
    warburg_impedance_07,
):
    _semi_reduced_omega = np.logspace(-4.0, 4.0, 500)
    _semi_impedance = warburg_impedance_07(
        _semi_reduced_omega, "semi-infinite"
    )
    _reduced_depth = np.linspace(0.0, 4.0, 260)
    _semi_profile = np.exp(-(1.0 + 1j) * _reduced_depth)
    _phases = np.array([0.0, 0.5, 1.0, 1.5]) * np.pi
    _figure, _axes = plt.subplots(
        1, 2, figsize=(12.8, 4.8), constrained_layout=True
    )
    _profile_axis, _nyquist_axis = _axes

    _profile_styles = ("-", "--", "-.", ":")
    _profile_colors = ("#4C7C86", "#B8734A", "#7C6A91", "#5F8A6B")
    for _phase, _style, _color in zip(
        _phases, _profile_styles, _profile_colors
    ):
        _profile_axis.plot(
            _reduced_depth,
            np.real(_semi_profile * np.exp(1j * _phase)),
            color=_color,
            ls=_style,
            lw=1.8,
            label=rf"$\omega t={_phase / np.pi:.1f}\pi$",
        )
    _profile_axis.axhline(0.0, color="#73808C", lw=0.9, ls=":")
    _profile_axis.set(
        xlabel=r"Reduced depth, $x\sqrt{\omega/(2D^\delta)}$ (dimensionless)",
        ylabel=r"Normalized concentration, $\Delta c/|\widehat{\Delta c}(0)|$",
        title="The oscillation decays with depth",
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
        xlabel=r"Normalized real impedance, $\widetilde Z'$",
        ylabel=r"Normalized negative imaginary impedance, $-\widetilde Z''$",
        title=r"Diffusion gives $Z'=-Z''$",
    )
    set_equal_nyquist_limits(_nyquist_axis, _semi_impedance.real, -_semi_impedance.imag)
    format_nyquist_axis_07(_nyquist_axis)


    plt.close(_figure)
    mo.vstack([
        _figure,
        mo.md(r"""
        A semi-infinite concentration wave decays and
        accumulates phase with depth. That same solution gives both the
        $45^\circ$ Nyquist line. The phase lag accumulates with depth.
        """),
    ])
    return


@app.cell
def _(
    format_nyquist_axis_07,
    mo,
    np,
    plt,
    warburg_boundary_07,
    warburg_compare_boundaries_07,
    warburg_frequency_hz_07,
    warburg_phase_radians_07,
    warburg_position_um_07,
    warburg_reduced_omega_07,
    warburg_scale_data_07,
    warburg_selected_frequency_hz_07,
    warburg_selected_impedance_07,
    warburg_selected_omega_07,
    warburg_selected_profile_07,
    warburg_spectra_07,
    set_equal_nyquist_limits,
):
    _colors = {
        "open": "#4C7C86",
        "blocked": "#B8734A",
    }
    _labels = {
        "open": r"Fixed composition: $\widehat{\Delta c}(L)=0$",
        "blocked": r"Zero flux: $d\widehat{\Delta c}/dx|_L=0$",
    }
    _styles = {"open": "-", "blocked": "--"}
    _boundaries = (
        ("open", "blocked")
        if warburg_compare_boundaries_07.value
        else (warburg_boundary_07.value,)
    )
    _selected_index = int(
        np.argmin(
            np.abs(
                np.log(warburg_frequency_hz_07 / warburg_selected_frequency_hz_07)
            )
        )
    )
    _resistance_scale = warburg_scale_data_07["resistance_diffusion_ohm"]

    _figure, (_axis_profile, _axis_warburg) = plt.subplots(
        1, 2, figsize=(12.8, 5.1), constrained_layout=True
    )

    _profile_envelope = np.abs(warburg_selected_profile_07)
    _profile_snapshot = np.real(
        warburg_selected_profile_07 * np.exp(1j * warburg_phase_radians_07)
    )
    _axis_profile.fill_between(
        warburg_position_um_07,
        -_profile_envelope,
        _profile_envelope,
        color="#B8734A",
        alpha=0.18,
        label=r"envelope $\pm|\widehat{\Delta c}|$",
    )
    _axis_profile.plot(
        warburg_position_um_07,
        _profile_snapshot,
        color="#4C7C86",
        lw=1.9,
        label="selected phase",
    )
    _axis_profile.plot(
        warburg_position_um_07,
        _profile_envelope,
        color="#B8734A",
        ls="--",
        lw=1.5,
    )
    _axis_profile.plot(
        warburg_position_um_07,
        -_profile_envelope,
        color="#B8734A",
        ls="--",
        lw=1.5,
    )
    _axis_profile.axhline(0.0, color="#73808C", lw=0.9)
    _axis_profile.axvline(warburg_position_um_07[-1], color="#73808C", lw=1.0, ls=":")
    _axis_profile.set(
        xlabel=r"Position, $x$ ($\mu$m)",
        ylabel=r"Normalized $\Delta c$ (dimensionless)",
        title=_labels[warburg_boundary_07.value],
    )
    _axis_profile.grid(True, alpha=0.24)
    _axis_profile.legend(loc="best", fontsize=10.5)

    _selected_extent = max(
        abs(warburg_selected_impedance_07.real),
        abs(warburg_selected_impedance_07.imag),
    )
    _nyquist_limit = max(
        2.0 * _resistance_scale,
        min(8.0 * _resistance_scale, 1.18 * _selected_extent),
    )
    for _boundary in _boundaries:
        _impedance = warburg_spectra_07[_boundary]
        _plot_mask = (
            (-_impedance.imag >= -1.0e-12)
            & (-_impedance.imag <= _nyquist_limit)
            & (_impedance.real <= _nyquist_limit)
        )
        _axis_warburg.plot(
            _impedance.real[_plot_mask],
            -_impedance.imag[_plot_mask],
            lw=1.8,
            ls=_styles[_boundary],
            color=_colors[_boundary],
            label=_labels[_boundary],
        )
    if _selected_extent <= _nyquist_limit:
        _axis_warburg.scatter(
            warburg_selected_impedance_07.real,
            -warburg_selected_impedance_07.imag,
            s=95,
            facecolor="none",
            edgecolor="#30343B",
            linewidth=1.8,
            zorder=5,
            label="selected frequency",
        )
    else:
        _axis_warburg.text(
            0.22,
            0.94,
            "Selected low-frequency point\nis above this viewport",
            transform=_axis_warburg.transAxes,
            va="top",
            color="#526173",
        )
    plotted_warburg = np.concatenate([
        warburg_spectra_07[boundary][
            (-warburg_spectra_07[boundary].imag >= -1.0e-12)
            & (-warburg_spectra_07[boundary].imag <= _nyquist_limit)
            & (warburg_spectra_07[boundary].real <= _nyquist_limit)
        ]
        for boundary in _boundaries
    ])
    _axis_warburg.set(
        xlabel=r"Real impedance, $Z'$ ($\Omega$)",
        ylabel=r"Negative imaginary impedance, $-Z''$ ($\Omega$)",
        title="Far boundary controls the low-frequency end",
    )
    set_equal_nyquist_limits(
        _axis_warburg, plotted_warburg.real, -plotted_warburg.imag
    )
    format_nyquist_axis_07(_axis_warburg)
    _axis_warburg.legend(loc="best", fontsize=10.0)


    plt.close(_figure)

    _boundary_sentence = {
        "open": "The reservoir pins the far-face composition.",
        "blocked": "The zero-flux wall stores material and becomes capacitive at low frequency.",
    }[warburg_boundary_07.value]
    mo.vstack([
        _figure,
        mo.md(
            rf"""
            {_boundary_sentence} The selected point has
            $f/f_D=\widetilde\omega={warburg_selected_omega_07:.3g}$,
            with $f_D={warburg_scale_data_07['frequency_diffusion_hz']:.3g}$ Hz
            and $R_D={warburg_scale_data_07['resistance_diffusion_ohm']:.3g}\ \Omega$.
            Frequency, $D^\delta$, and $L$ change whether the profile reaches
            the far face; $T$, $c_0$, and $S$ change the impedance scale.
            """
        ),
    ])
    return


@app.cell
def _(mo):
    finite_limits = mo.md(r"""
    ### Read the two finite-length limits

    Both finite boundaries recover the semi-infinite $45^\circ$ response at
    high frequency. Their low-frequency limits explain why the two curves then
    separate:

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

    The dimensionless shape is controlled by $f/f_D$, with
    $f_D=D^\delta/(2\pi L^2)$. The dilute resistance scale
    $R_D=RTL/(F^2c_0SD^\delta)$ converts that shape to ohms. This is why the
    profile controls ($f$, $D^\delta$, $L$) and amplitude controls ($T$, $c_0$,
    $S$) teach different parts of the same response.
    """)
    mo.accordion({"Explore further — finite-length limits and impedance scaling": finite_limits})
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
    tlm_details = mo.md(r"""
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
    D^\delta=\frac{1}{(r_e+r_i)c_{\rm chem}},\qquad
    \boxed{(R_e+R_i)C_{\rm chem}=\frac{L^2}{D^\delta}},
    \qquad
    f_{\rm chem}=\frac{1}{2\pi(R_e+R_i)C_{\rm chem}}.
    $$

    Module 07 derived the storage term. For the ideal monovalent pair used in
    Modules 03, 05, and 06,

    $$
    C_{\rm chem}=c_{\rm chem}^VSL
    =\frac{F^2SLc}{2RT},\qquad c_{\rm chem}=c_{\rm chem}^V S.
    $$

    The TLM therefore recovers the same chemical diffusivity and the same
    finite-slab clock rather than introducing a second transport coefficient.

    Here $r_e,r_i$ have units $\Omega$ m$^{-1}$ and $c_{\rm chem}$ has units
    F m$^{-1}$. Total quantities are $R_e=r_eL$, $R_i=r_iL$, and
    $C_{\rm chem}=c_{\rm chem}L$. The area $S$ is already included in the rail
    resistances. The resulting $f_{\rm chem}$ is a scaling frequency, not a
    promise that every boundary condition has a peak exactly there.

    A dielectric capacitance $C_{\rm diel}$ can be added for an electrolyte,
    but it is shown only as an application below and is not in this teaching model.
    """)
    mo.vstack([
        mo.md(r"""
        ## 3. A mixed conductor needs two connected pathways

        A MIEC carries electronic current on one rail and ionic current on the
        other. Distributed chemical capacitance connects the rails because a
        local change in stoichiometry stores chemical free energy. Current may
        transfer between the rails, but total current is conserved:

        $$\frac{d(I_e+I_i)}{dx}=0.$$

        The end contacts decide which carrier may enter or leave. The same bulk
        material can therefore show a different spectrum when its boundary
        conditions change. In an appropriate selective-contact limit, this
        distributed model becomes the finite-diffusion response seen above.
        """),
        mo.accordion({"Model details — TLM equations, signs, units, and time scale": tlm_details}),
    ])
    return


@app.cell
def _(Rectangle, mo, np, plt):
    _figure, _axis = plt.subplots(figsize=(12.8, 4.2), constrained_layout=True)
    _axis.add_patch(
        Rectangle(
            (-0.14, 0.16), 0.14, 0.62,
            facecolor="#D9DDE0", edgecolor="#66717B", linewidth=1.0,
        )
    )
    _axis.add_patch(
        Rectangle(
            (0.0, 0.16), 1.0, 0.62,
            facecolor="#F3EEDB", edgecolor="#66717B", linewidth=1.0,
        )
    )
    _axis.add_patch(
        Rectangle(
            (1.0, 0.16), 0.14, 0.62,
            facecolor="#D9DDE0", edgecolor="#66717B", linewidth=1.0,
        )
    )

    _electronic_y, _ionic_y = 0.60, 0.32
    _axis.plot([0.0, 1.0], [_electronic_y, _electronic_y], color="#4C7C86", lw=2.0)
    _axis.plot([0.0, 1.0], [_ionic_y, _ionic_y], color="#B8734A", lw=2.0)
    for _capacitor_x in np.linspace(0.14, 0.86, 6):
        _axis.plot(
            [_capacitor_x, _capacitor_x], [_ionic_y + 0.025, 0.43],
            color="#665777", lw=1.2,
        )
        _axis.plot(
            [_capacitor_x - 0.022, _capacitor_x + 0.022], [0.43, 0.43],
            color="#665777", lw=1.8,
        )
        _axis.plot(
            [_capacitor_x - 0.022, _capacitor_x + 0.022], [0.49, 0.49],
            color="#665777", lw=1.8,
        )
        _axis.plot(
            [_capacitor_x, _capacitor_x], [0.49, _electronic_y - 0.025],
            color="#665777", lw=1.2,
        )

    _terminal_style = {
        "boxstyle": "round,pad=0.18",
        "facecolor": "white",
        "edgecolor": "#9AA3AB",
        "alpha": 0.96,
    }
    for _x, _y, _label, _side, _color in (
        (0.0, _electronic_y, r"$Z_A$", "left", "#4C7C86"),
        (0.0, _ionic_y, r"$Z_B$", "left", "#B8734A"),
        (1.0, _electronic_y, r"$Z_C$", "right", "#4C7C86"),
        (1.0, _ionic_y, r"$Z_D$", "right", "#B8734A"),
    ):
        _direction = -1.0 if _side == "left" else 1.0
        _axis.plot([_x, _x + 0.14 * _direction], [_y, _y], color=_color, lw=1.8)
        _axis.text(
            _x + 0.21 * _direction,
            _y,
            _label,
            ha="right" if _side == "left" else "left",
            va="center",
            color=_color,
            bbox=_terminal_style,
        )

    _axis.plot(
        [-0.07, -0.07, 0.485], [0.78, 0.96, 0.96],
        color="#9AA3AB", ls="--", lw=1.3,
    )
    _axis.plot(
        [0.515, 1.07, 1.07], [0.96, 0.96, 0.78],
        color="#9AA3AB", ls="--", lw=1.3,
    )
    _axis.plot([0.485, 0.485], [0.91, 1.01], color="#9AA3AB", lw=1.7)
    _axis.plot([0.515, 0.515], [0.91, 1.01], color="#9AA3AB", lw=1.7)

    _axis.text(
        0.50, 0.735, "Uniform mixed ionic-electronic conductor",
        ha="center", va="center", fontsize=13, weight="bold",
    )
    _axis.text(
        0.50, 0.645, r"Electronic rail  $r_e$",
        ha="center", va="center", color="#3D6972",
    )
    _axis.text(
        0.50, 0.275, r"Ionic rail  $r_i$",
        ha="center", va="center", color="#8F5638",
    )
    _axis.text(
        0.50, 0.46, r"Distributed $c_{\rm chem}$",
        ha="center", va="center", color="#665777",
        bbox={"boxstyle": "round,pad=0.16", "facecolor": "#F3EEDB", "edgecolor": "none"},
    )
    _axis.text(-0.07, 0.07, "LEFT ELECTRODE", ha="center", va="center", fontsize=10)
    _axis.text(1.07, 0.07, "RIGHT ELECTRODE", ha="center", va="center", fontsize=10)
    _axis.text(
        0.50, 1.075, r"$C_{\rm diel}$ is not included in this simplified model",
        ha="center", va="center", color="#73808C", fontsize=11,
    )
    _axis.annotate(
        "Positive $x$", xy=(0.68, 0.09), xytext=(0.42, 0.09),
        ha="center", va="center",
        arrowprops={"arrowstyle": "->", "color": "#526173", "linewidth": 1.1},
    )
    _axis.set(xlim=(-0.38, 1.38), ylim=(0.0, 1.15))
    _axis.axis("off")
    plt.close(_figure)
    mo.vstack([
        _figure,
        mo.md(r"""
        The schematic makes the distributed model explicit:
        electrodes surround a uniform MIEC, the two rails carry electronic and
        ionic current, and each properly drawn $c_{\rm chem}$ element stores
        composition between—not along—the rails. $Z_A$–$Z_D$ set the four
        rail-end boundary conditions.
        """),
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
        },
        value="Electron-reversible contacts; ions blocked at both faces",
        label="Contact boundary conditions",
    )
    tlm_log_ratio_07 = mo.ui.slider(
        start=-3.0,
        stop=3.0,
        step=0.25,
        value=0.0,
        label="Conductivity-ratio exponent",
    )
    tlm_log_r_parallel_per_m_07 = mo.ui.slider(
        start=3.0,
        stop=9.0,
        step=0.25,
        value=6.0,
        label=r"$\log_{10}(r_\parallel/\Omega\,\mathrm{m}^{-1})$",
    )
    tlm_log_c_chemical_per_m_07 = mo.ui.slider(
        start=-4.0,
        stop=4.0,
        step=0.25,
        value=0.0,
        label=r"$\log_{10}(c_{\rm chem}/\mathrm{F\,m^{-1}})$",
    )
    tlm_length_07 = mo.ui.slider(
        start=10, stop=500, step=10, value=100, label="TLM length (micrometers)"
    )
    tlm_log_selected_frequency_07 = mo.ui.slider(
        start=-8.0,
        stop=6.0,
        step=0.25,
        value=0.5,
        label="TLM-frequency exponent",
    )
    tlm_profile_phase_07 = mo.ui.slider(
        start=0, stop=330, step=30, value=60, label="Snapshot phase (degrees)"
    )
    tlm_internal_view_07 = mo.ui.dropdown(
        options=["Composition response", "Voltage-equivalent potentials", "Rail currents"],
        value="Composition response",
        label="Internal view",
    )
    core_tlm_controls_07 = mo.hstack(
        [tlm_contact_case_07, tlm_log_ratio_07, tlm_log_selected_frequency_07],
        justify="start", align="center", wrap=True, gap=1.4,
    )
    advanced_tlm_controls_07 = mo.hstack(
        [tlm_log_r_parallel_per_m_07, tlm_log_c_chemical_per_m_07, tlm_length_07],
        justify="start", align="center", wrap=True, gap=1.4,
    )
    mo.vstack([
        core_tlm_controls_07,
        mo.accordion({"Explore further — distributed transport scale": advanced_tlm_controls_07}),
    ])
    return (
        tlm_contact_case_07,
        tlm_internal_view_07,
        tlm_length_07,
        tlm_log_c_chemical_per_m_07,
        tlm_log_ratio_07,
        tlm_log_r_parallel_per_m_07,
        tlm_log_selected_frequency_07,
        tlm_profile_phase_07,
    )


@app.cell
def _(
    np,
    tlm_contact_case_07,
    tlm_distributed_parameters_07,
    tlm_length_07,
    tlm_log_c_chemical_per_m_07,
    tlm_log_ratio_07,
    tlm_log_r_parallel_per_m_07,
    tlm_log_selected_frequency_07,
    tlm_profile_07,
    tlm_spectrum_07,
):
    tlm_parameter_data_07 = tlm_distributed_parameters_07(
        10.0 ** tlm_log_r_parallel_per_m_07.value,
        10.0 ** tlm_log_c_chemical_per_m_07.value,
        tlm_length_07.value,
        10.0 ** tlm_log_ratio_07.value,
    )
    _log_f_chemical = np.log10(tlm_parameter_data_07["frequency_chemical_hz"])
    _log_f_min = min(_log_f_chemical - 4.0, tlm_log_selected_frequency_07.value - 0.25)
    _log_f_max = max(_log_f_chemical + 4.0, tlm_log_selected_frequency_07.value + 0.25)
    tlm_frequency_hz_07 = np.logspace(
        _log_f_min, _log_f_max, 520
    )
    tlm_frequency_ratio_07 = (
        tlm_frequency_hz_07 / tlm_parameter_data_07["frequency_chemical_hz"]
    )
    tlm_selected_frequency_hz_07 = 10.0 ** tlm_log_selected_frequency_07.value
    tlm_selected_omega_07 = (
        tlm_selected_frequency_hz_07
        / tlm_parameter_data_07["frequency_chemical_hz"]
    )
    tlm_spectrum_data_07 = tlm_spectrum_07(
        tlm_frequency_ratio_07, tlm_parameter_data_07, tlm_contact_case_07.value
    )
    tlm_position_07 = np.linspace(0.0, 1.0, 300)
    tlm_position_um_07 = tlm_length_07.value * tlm_position_07
    tlm_profile_data_07 = tlm_profile_07(
        tlm_position_07,
        tlm_selected_omega_07,
        tlm_parameter_data_07,
        tlm_contact_case_07.value,
    )
    return (
        tlm_frequency_hz_07,
        tlm_frequency_ratio_07,
        tlm_parameter_data_07,
        tlm_position_07,
        tlm_position_um_07,
        tlm_profile_data_07,
        tlm_selected_frequency_hz_07,
        tlm_selected_omega_07,
        tlm_spectrum_data_07,
    )


@app.cell
def _(mo, tlm_parameter_data_07):
    mo.callout(
        mo.md(
            rf"""
            The distributed controls imply
            $D^\delta=\mathbf{{{tlm_parameter_data_07['D_chemical_m2_per_s'] * 1.0e4:.3e}}}$
            cm² s$^{{-1}}$ and
            $(R_e+R_i)C_{{\rm chem}}=\mathbf{{{tlm_parameter_data_07['tau_chemical_s']:.3e}}}$ s.
            These are the same storage–transport identities developed in Module 07.
            """
        ),
        kind="info",
    )
    return


@app.cell
def _(
    format_nyquist_axis_07,
    mo,
    np,
    plt,
    tlm_contact_case_07,
    tlm_frequency_hz_07,
    tlm_parameter_data_07,
    tlm_selected_frequency_hz_07,
    tlm_spectrum_data_07,
    set_equal_nyquist_limits,
):
    _impedance = tlm_spectrum_data_07
    _real = _impedance.real
    _minus_imaginary = -_impedance.imag
    _selected_index = int(
        np.argmin(np.abs(np.log(tlm_frequency_hz_07 / tlm_selected_frequency_hz_07)))
    )
    _representative_frequencies = np.logspace(
        np.log10(tlm_frequency_hz_07[0]), np.log10(tlm_frequency_hz_07[-1]), 7
    )
    _representative_indices = np.array([
        int(np.argmin(np.abs(tlm_frequency_hz_07 - frequency)))
        for frequency in _representative_frequencies
    ])

    _nyquist_figure, _nyquist_axis = plt.subplots(
        figsize=(8.4, 6.0), constrained_layout=True
    )
    _nyquist_axis.plot(_real, _minus_imaginary, color="#4C7C86", lw=1.9)
    _nyquist_axis.scatter(
        _real[_representative_indices],
        _minus_imaginary[_representative_indices],
        s=42,
        color="#B8734A",
        edgecolors="white",
        linewidths=0.6,
        zorder=3,
        label="representative frequencies",
    )
    _nyquist_axis.scatter(
        _real[_selected_index],
        _minus_imaginary[_selected_index],
        s=125,
        facecolors="none",
        edgecolors="#30343B",
        linewidths=1.8,
        zorder=4,
        label="selected frequency",
    )
    _nyquist_axis.set(
        xlabel=r"Real impedance, $Z'$ ($\Omega$)",
        ylabel=r"Negative imaginary impedance, $-Z''$ ($\Omega$)",
        title="Distributed transport and storage create one spectrum",
    )
    set_equal_nyquist_limits(_nyquist_axis, _real, _minus_imaginary)
    format_nyquist_axis_07(_nyquist_axis)
    _nyquist_axis.legend(frameon=False, loc="best")

    _bode_figure, (_magnitude_axis, _phase_axis) = plt.subplots(
        1, 2, figsize=(12.8, 4.4), constrained_layout=True
    )
    _magnitude_axis.loglog(
        tlm_frequency_hz_07, np.abs(_impedance), color="#4C7C86", lw=1.8
    )
    _phase_axis.semilogx(
        tlm_frequency_hz_07, np.angle(_impedance, deg=True), color="#B8734A", lw=1.8
    )
    _magnitude_axis.set(
        xlabel=r"Frequency, $f$ (Hz)", ylabel=r"Magnitude, $|Z|$ ($\Omega$)",
        title="Magnitude",
    )
    _phase_axis.set(
        xlabel=r"Frequency, $f$ (Hz)", ylabel=r"Phase of $Z$ (degrees)",
        title="Phase",
    )
    for _axis in (_magnitude_axis, _phase_axis):
        _axis.axvline(
            tlm_parameter_data_07["frequency_chemical_hz"],
            color="#73808C", ls=":", lw=1.2,
        )
        _axis.grid(True, which="both", alpha=0.24)
    plt.close(_nyquist_figure)
    plt.close(_bode_figure)

    mo.vstack([
        mo.md("### Read one contact-dependent TLM spectrum"),
        _nyquist_figure,
        mo.md(
            rf"""
            The selected contact preset is **{tlm_contact_case_07.value}**.
            The line comes from distributed rail transport and chemical
            storage; no Warburg segment is inserted. The black ring selects the
            frequency used for the spatial view below. Here
            $f_{{\rm chem}}={tlm_parameter_data_07['frequency_chemical_hz']:.3g}$ Hz.
            """
        ),
        mo.accordion({"Explore further — TLM Bode magnitude and phase": _bode_figure}),
    ])
    return

@app.cell
def _(
    mo,
    np,
    plt,
    tlm_internal_view_07,
    tlm_parameter_data_07,
    tlm_position_um_07,
    tlm_profile_data_07,
    tlm_profile_phase_07,
    tlm_selected_frequency_hz_07,
    tlm_selected_omega_07,
):
    _phase_factor = np.exp(1j * np.deg2rad(tlm_profile_phase_07.value))
    _chemical_hat = tlm_profile_data_07["u_chemical"]
    _chemical_scale = max(float(np.max(np.abs(_chemical_hat))), 1.0e-30)
    _concentration_hat = _chemical_hat / _chemical_scale
    _concentration_snapshot = np.real(_concentration_hat * _phase_factor)
    _concentration_envelope = np.abs(_concentration_hat)
    _u_e_snapshot = np.real(tlm_profile_data_07["u_e"] * _phase_factor)
    _u_i_snapshot = np.real(tlm_profile_data_07["u_i"] * _phase_factor)
    _chemical_snapshot = np.real(_chemical_hat * _phase_factor)
    _current_scale = max(
        float(np.max(np.abs(tlm_profile_data_07["I_total"]))),
        1.0 / tlm_parameter_data_07["R_parallel_ohm"],
    )
    _i_e_snapshot = np.real(tlm_profile_data_07["I_e"] * _phase_factor) / _current_scale
    _i_i_snapshot = np.real(tlm_profile_data_07["I_i"] * _phase_factor) / _current_scale
    _i_total_snapshot = np.real(tlm_profile_data_07["I_total"] * _phase_factor) / _current_scale

    _figure, _axis = plt.subplots(figsize=(10.8, 4.8), constrained_layout=True)
    _view = tlm_internal_view_07.value
    if _view == "Composition response":
        _length_um = tlm_position_um_07[-1]
        _electrode_width_um = 0.035 * _length_um
        _axis.axvspan(-_electrode_width_um, 0.0, color="#D9DDE0", alpha=0.9)
        _axis.axvspan(_length_um, _length_um + _electrode_width_um, color="#D9DDE0", alpha=0.9)
        _axis.fill_between(
            tlm_position_um_07, -_concentration_envelope, _concentration_envelope,
            color="#B8734A", alpha=0.16, label="full-cycle envelope",
        )
        _axis.plot(
            tlm_position_um_07, _concentration_snapshot,
            color="#4C7C86", lw=1.9, label="selected-phase snapshot",
        )
        _axis.set(
            xlim=(-_electrode_width_um, _length_um + _electrode_width_um),
            ylim=(-1.08, 1.08),
            ylabel="Normalized composition response",
            title="Chemical storage profile",
        )
        _explanation = (
            "The stored rail difference is the voltage equivalent of the neutral "
            "composition chemical potential. This spatial wave connects directly to Warburg diffusion."
        )
    elif _view == "Voltage-equivalent potentials":
        _axis.plot(tlm_position_um_07, _u_e_snapshot, color="#4C7C86", lw=1.9, label=r"$u_e$")
        _axis.plot(tlm_position_um_07, _u_i_snapshot, color="#B8734A", lw=1.9, label=r"$u_i$")
        _axis.plot(
            tlm_position_um_07, _chemical_snapshot,
            color="#7C6A91", lw=1.6, ls="--", label=r"$u_e-u_i$",
        )
        _axis.set(ylabel="Voltage-equivalent potential (V)", title="Potentials along the two rails")
        _explanation = (
            r"The difference $u_e-u_i=-\mu_M/F$ drives local chemical storage; "
            "the two rail potentials also determine their separate current gradients."
        )
    else:
        _axis.plot(tlm_position_um_07, _i_e_snapshot, color="#4C7C86", lw=1.9, label=r"$I_e$")
        _axis.plot(tlm_position_um_07, _i_i_snapshot, color="#B8734A", lw=1.9, label=r"$I_i$")
        _axis.plot(
            tlm_position_um_07, _i_total_snapshot,
            color="#5F8A6B", lw=1.6, ls="--", label=r"$I_e+I_i$",
        )
        _axis.set(ylabel="Current / common scale", title="Current transfers between rails")
        _explanation = (
            "Electronic and ionic currents can exchange with position, but their "
            "sum remains constant everywhere in the line."
        )
    _axis.axhline(0.0, color="#73808C", lw=0.9)
    _axis.set_xlabel(r"Position, $x$ ($\mu$m)")
    _axis.grid(True, alpha=0.24)
    _axis.legend(frameon=False, loc="best")
    _figure.suptitle(
        rf"$f={tlm_selected_frequency_hz_07:.3g}$ Hz, "
        + rf"$f/f_{{\rm chem}}={tlm_selected_omega_07:.3g}$, "
        + rf"phase $={tlm_profile_phase_07.value}^\circ$",
        fontsize=14,
    )
    plt.close(_figure)

    mo.vstack([
        mo.md("### Look inside the line at one frequency"),
        mo.hstack(
            [tlm_internal_view_07, tlm_profile_phase_07],
            justify="start",
            align="center",
            wrap=True,
            gap=1.4,
        ),
        _figure,
        mo.md(_explanation),
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
        capacitance, and two ideal contact presets. It does **not** fit finite
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

    These are controlled simplifications of one transport model,
    not three unrelated equivalent circuits.
    """)
    mo.accordion({"Explore further — TLM model scope and applications": mo.vstack([_scope, _applications])})
    return



@app.cell
def _(
    capacitor_impedance_07,
    format_nyquist_axis_07,
    np,
    parallel_rc_impedance_07,
    plt,
    resistor_impedance_07,
    series_rc_impedance_07,
    set_equal_nyquist_limits,
    tlm_distributed_parameters_07,
    tlm_parameters_07,
    tlm_profile_07,
    tlm_solution_07,
    tlm_spectrum_07,
    warburg_impedance_07,
    warburg_profile_07,
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
    _capacitor_y = -_check_capacitor.imag
    _capacitor_limits = set_equal_nyquist_limits(
        _nyquist_test_axis, _check_capacitor.real, _capacitor_y
    )
    format_nyquist_axis_07(_nyquist_test_axis)
    _x_limits, _y_limits = _capacitor_limits
    nyquist_equal_axis_error_07 = abs(
        (_x_limits[1] - _x_limits[0]) - (_y_limits[1] - _y_limits[0])
    ) / (_y_limits[1] - _y_limits[0])
    nyquist_capacitor_margin_pass_07 = bool(
        _x_limits[0] < 0.0 < _x_limits[1]
        and np.all(_capacitor_y > _y_limits[0])
        and np.all(_capacitor_y < _y_limits[1])
    )
    plt.close(_nyquist_test_figure)
    _check_series_rc = series_rc_impedance_07(_check_frequency, 37.0, 2.5e-6)
    series_rc_limit_error_07 = max(
        np.max(np.abs(_check_series_rc.real - 37.0)),
        np.max(np.abs(_check_series_rc.imag - _check_capacitor.imag)),
    )
    _check_rc = parallel_rc_impedance_07(
        _check_frequency, 1.0, 1.0
    )
    rc_circle_error_07 = np.max(
        np.abs((_check_rc.real - 0.5) ** 2 + _check_rc.imag**2 - 0.25)
    )
    _peak_index = int(np.argmax(-_check_rc.imag))
    rc_peak_error_07 = abs(2.0 * np.pi * _check_frequency[_peak_index] - 1.0)

    _waveform_check = waveform_data_07(3.0, 37.0)
    phasor_sign_error_07 = abs(_waveform_check[3] + 37.0)
    _waveform_slow = waveform_data_07(0.2, 20.0)
    _waveform_fast = waveform_data_07(20.0, 20.0)
    waveform_frequency_error_07 = max(
        abs((_waveform_slow[0][-1] - _waveform_slow[0][0]) - 2.0),
        abs((_waveform_fast[0][-1] - _waveform_fast[0][0]) - 2.0),
        float(np.array_equal(_waveform_slow[1], _waveform_fast[1])),
    )

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
    _warburg_profile_positions = np.linspace(0.0, 1.0, 301)
    _warburg_extreme_profiles = {
        _boundary: warburg_profile_07(
            1.0e13, _warburg_profile_positions, _boundary
        )
        for _boundary in ("semi-infinite", "open", "blocked")
    }
    warburg_profile_finiteness_07 = all(
        np.all(np.isfinite(_profile))
        for _profile in _warburg_extreme_profiles.values()
    )
    warburg_profile_boundary_error_07 = max(
        max(abs(_profile[0] - 1.0) for _profile in _warburg_extreme_profiles.values()),
        abs(_warburg_extreme_profiles["open"][-1]),
    )

    _tlm_check_parameters = tlm_parameters_07(100.0, 1.0e-4, 100.0)
    _tlm_check_frequency = np.logspace(-12.0, 12.0, 260)
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
    _tlm_capacitive_line = tlm_spectrum_07(
        _tlm_check_frequency, _tlm_check_parameters, "cross-selective contacts"
    )
    _tlm_margin_figure, _tlm_margin_axis = plt.subplots()
    _tlm_x_limits, _tlm_y_limits = set_equal_nyquist_limits(
        _tlm_margin_axis, _tlm_capacitive_line.real, -_tlm_capacitive_line.imag
    )
    nyquist_tlm_margin_pass_07 = bool(
        _tlm_x_limits[0] < float(np.min(_tlm_capacitive_line.real))
        and np.all(_tlm_capacitive_line.real > _tlm_x_limits[0])
        and np.all(_tlm_capacitive_line.real < _tlm_x_limits[1])
        and np.all(-_tlm_capacitive_line.imag > _tlm_y_limits[0])
        and np.all(-_tlm_capacitive_line.imag < _tlm_y_limits[1])
        and abs(
            (_tlm_x_limits[1] - _tlm_x_limits[0])
            - (_tlm_y_limits[1] - _tlm_y_limits[0])
        ) < 1.0e-12 * (_tlm_y_limits[1] - _tlm_y_limits[0])
    )
    plt.close(_tlm_margin_figure)
    tlm_passivity_margin_07 = min(
        np.min(tlm_spectrum_07(_tlm_check_frequency, _tlm_check_parameters, _case).real)
        for _case in _tlm_cases
    )
    tlm_boundary_error_07 = max(
        tlm_solution_07(10.0 ** _frequency, _tlm_check_parameters, _case)["boundary_residual"]
        for _case in _tlm_cases
        for _frequency in (-12.0, -3.0, 0.0, 3.0, 12.0)
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
    _tlm_extreme_profile = tlm_profile_07(
        np.linspace(0.0, 1.0, 301),
        1.0e13,
        _tlm_check_parameters,
        "cross-selective contacts",
    )
    tlm_extreme_profile_finiteness_07 = all(
        np.all(np.isfinite(_tlm_extreme_profile[_key]))
        for _key in ("u_e", "u_i", "u_chemical", "I_e", "I_i", "I_total")
    )
    _tlm_extreme_current_scale = max(
        abs(_tlm_extreme_profile["I_total"][0]), 1.0e-30
    )
    tlm_extreme_profile_error_07 = max(
        abs(_tlm_extreme_profile["u_e"][0] - 1.0),
        abs(_tlm_extreme_profile["I_i"][0]),
        abs(_tlm_extreme_profile["I_e"][-1]),
        abs(_tlm_extreme_profile["u_i"][-1]),
        np.ptp(_tlm_extreme_profile["I_total"]) / _tlm_extreme_current_scale,
        np.max(np.abs(
            _tlm_extreme_profile["u_e"]
            - _tlm_extreme_profile["u_i"]
            - _tlm_extreme_profile["u_chemical"]
        )),
    )
    _tlm_distributed_check = tlm_distributed_parameters_07(
        4.0e5, 2.5, 250.0, 100.0
    )
    _length_check_m = _tlm_distributed_check["length_m"]
    tlm_total_distributed_error_07 = max(
        abs(_tlm_distributed_check["r_parallel_ohm_per_m"] * _length_check_m / _tlm_distributed_check["R_parallel_ohm"] - 1.0),
        abs(_tlm_distributed_check["c_chemical_f_per_m"] * _length_check_m / _tlm_distributed_check["C_chemical_f"] - 1.0),
        abs(_tlm_distributed_check["r_e_ohm_per_m"] * _length_check_m / _tlm_distributed_check["R_e_ohm"] - 1.0),
        abs(_tlm_distributed_check["r_i_ohm_per_m"] * _length_check_m / _tlm_distributed_check["R_i_ohm"] - 1.0),
    )
    _distributed_d_from_line = 1.0 / (
        (_tlm_distributed_check["r_e_ohm_per_m"] + _tlm_distributed_check["r_i_ohm_per_m"])
        * _tlm_distributed_check["c_chemical_f_per_m"]
    )
    _distributed_d_from_time = _length_check_m**2 / _tlm_distributed_check["tau_chemical_s"]
    tlm_diffusivity_mapping_error_07 = max(
        abs(_tlm_distributed_check["D_chemical_m2_per_s"] / _distributed_d_from_line - 1.0),
        abs(_tlm_distributed_check["D_chemical_m2_per_s"] / _distributed_d_from_time - 1.0),
    )
    _scale_check = warburg_scales_07(120.0, 2.0e-8, 250.0, 0.8, 800.0, 1.0)
    warburg_resistance_scale_error_07 = abs(
        _scale_check["resistance_general_ohm"]
        / _scale_check["resistance_diffusion_ohm"]
        - 1.0
    )
    return (
        capacitor_limit_error_07,
        nyquist_capacitor_margin_pass_07,
        nyquist_equal_axis_error_07,
        nyquist_tlm_margin_pass_07,
        resistor_limit_error_07,
        phasor_sign_error_07,
        rc_circle_error_07,
        rc_peak_error_07,
        series_rc_limit_error_07,
        tlm_boundary_error_07,
        tlm_current_conservation_error_07,
        tlm_diffusivity_mapping_error_07,
        tlm_extreme_profile_error_07,
        tlm_extreme_profile_finiteness_07,
        tlm_finiteness_07,
        tlm_passivity_margin_07,
        tlm_reversible_error_07,
        tlm_total_distributed_error_07,
        tlm_voltage_mapping_error_07,
        warburg_high_frequency_error_07,
        warburg_low_frequency_error_07,
        warburg_passivity_margin_07,
        warburg_profile_boundary_error_07,
        warburg_profile_finiteness_07,
        warburg_resistance_scale_error_07,
        waveform_frequency_error_07,
    )


@app.cell
def _(
    capacitor_limit_error_07,
    series_rc_limit_error_07,
    mo,
    nyquist_capacitor_margin_pass_07,
    nyquist_equal_axis_error_07,
    nyquist_tlm_margin_pass_07,
    resistor_limit_error_07,
    phasor_sign_error_07,
    rc_circle_error_07,
    rc_peak_error_07,
    tlm_boundary_error_07,
    tlm_current_conservation_error_07,
    tlm_diffusivity_mapping_error_07,
    tlm_extreme_profile_error_07,
    tlm_extreme_profile_finiteness_07,
    tlm_finiteness_07,
    tlm_passivity_margin_07,
    tlm_reversible_error_07,
    tlm_total_distributed_error_07,
    tlm_voltage_mapping_error_07,
    warburg_high_frequency_error_07,
    warburg_low_frequency_error_07,
    warburg_passivity_margin_07,
    warburg_profile_boundary_error_07,
    warburg_profile_finiteness_07,
    warburg_resistance_scale_error_07,
    waveform_frequency_error_07,
):
    _checks = [
        (
            "Ideal elements and series-RC limit",
            max(resistor_limit_error_07, capacitor_limit_error_07, series_rc_limit_error_07) < 1.0e-12,
            r"The series example must equal the same resistor plus the same capacitive impedance at every frequency.",
        ),
        (
            "Equal Nyquist scaling and visible capacitive axis",
            nyquist_equal_axis_error_07 < 1.0e-14 and nyquist_capacitor_margin_pass_07 and nyquist_tlm_margin_pass_07,
            "Equal data scale preserves geometry, and a small negative real-axis margin keeps a vertical capacitive line inside the axes.",
        ),
        (
            "Waveform frequency and phasor sign",
            max(phasor_sign_error_07, waveform_frequency_error_07) < 1.0e-12,
            "The observation window stays at two seconds while frequency changes the visible cycle count; a leading current gives negative impedance phase for $e^{i\\omega t}$.",
        ),
        (
            "Ideal $R\\parallel C$ semicircle",
            max(rc_circle_error_07, rc_peak_error_07) < 2.0e-2,
            "The geometric arc and its $\\omega RC=1$ apex must come from the same circuit equation.",
        ),
        (
            "General and dilute Warburg resistance scales",
            warburg_resistance_scale_error_07 < 1.0e-14,
            "The general thermodynamic slope must reduce to the ideal dilute expression under the stated assumption.",
        ),
        (
            "Finite Warburg limits",
            warburg_profile_finiteness_07
            and max(
                warburg_high_frequency_error_07,
                warburg_low_frequency_error_07,
                warburg_profile_boundary_error_07,
            ) < 2.0e-6,
            "Both boundaries keep finite profiles, share the high-frequency 45 degree limit, and approach their own low-frequency forms.",
        ),
        (
            "Passive Warburg response",
            warburg_passivity_margin_07 > 0.0,
            "Diffusion dissipates energy, so the real part of its impedance cannot be negative.",
        ),
        (
            "TLM boundary conditions",
            tlm_extreme_profile_finiteness_07
            and max(tlm_boundary_error_07, tlm_extreme_profile_error_07) < 2.0e-10,
            "Every spectrum and profile must stay finite and satisfy the selected passing or blocking conditions at both faces.",
        ),
        (
            "TLM total-current conservation",
            abs(tlm_current_conservation_error_07) < 1.0e-12,
            "Current may transfer between rails, but $I_e+I_i$ must be independent of position.",
        ),
        (
            "TLM voltage and unit mappings",
            max(tlm_voltage_mapping_error_07, tlm_total_distributed_error_07, tlm_diffusivity_mapping_error_07) < 1.0e-12,
            r"$u_e-u_i$ must be the stored chemical voltage; distributed quantities must recover their totals and $D^\delta=1/[(r_e+r_i)c_{\rm chem}]$.",
        ),
        (
            "Collapsed reversible-contact regression",
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
    ## What to carry forward

    $$
    \boxed{
    \text{small sinusoid}
    \rightarrow Z(\omega)
    \rightarrow \text{time scale}
    \rightarrow \text{transport + storage + boundaries}
    }
    $$

    1. **Frequency separates response regimes.** A semicircle marks a
       relaxation, while a distributed chemical-diffusion response produces
       the Warburg frequency dependence.
    2. **The far boundary controls the low-frequency end.** Fixed composition
       gives a finite resistance; zero flux gives capacitive accumulation. Both
       share the high-frequency $45^\circ$ limit.
    3. **A TLM connects physics, storage, and contacts.** The same MIEC rails can
       produce different spectra when $Z_A$–$Z_D$ change, so geometry, sign
       convention, units, and boundary conditions must precede feature labels.

    This module deliberately uses ideal capacitors, uniform one-dimensional
    transport, linear response, and ideal contacts. Constant-phase elements,
    electrode kinetics, microstructural distributions, and nonlinear large-signal
    effects belong in later model extensions, not in the first interpretation.

    ### Continue with the full TLM teaching tool

    This module closely adapts the original tool's schematic, distributed
    parameter controls, frequency-colored spectra, and spatial teaching views,
    while limiting the boundary choices to two transparent ideal contact
    cases. The separate
    [TLM teaching tool](https://qiyanglu.github.io/TLM-teaching-tool/) exposes all
    four terminals, dielectric storage, and more general boundary impedances;
    its
    [source repository](https://github.com/qiyanglu/TLM-teaching-tool) documents
    the terminal signs and boundary conventions.

    ### Sources and further reading

    - Q. Lu, *Solid State Ionics, Lecture 8: Impedance Spectroscopy* (course
      slides). The phasor, Nyquist, $R\parallel C$, time-constant separation,
      and brick-layer narrative use the notation defined above.
    - Q. Lu, [Warburg impedance: more than a 45-degree line](https://mp.weixin.qq.com/s/CyCXnWWEoX586lzMGl0A9Q).
      This article motivates the finite-length discussion; the equations above
      state both boundary conditions explicitly.
    - Q. Lu, [Transmission lines for mixed ionic-electronic conductors](https://mp.weixin.qq.com/s/zR9QI0GnGUvPz2VRc9bHOA).
    - A. E. Bumberger, A. Nenning, and J. Fleig, “Transmission line revisited —
      the impedance of mixed ionic and electronic conductors,” *PCCP* **26** (2024),
      [doi:10.1039/D4CP00975D](https://doi.org/10.1039/D4CP00975D).
    - Course context and related tutorials:
      [Solid State Ionics teaching page](https://ssi-westlake.com/teaching/) and
      [tutorial index](https://ssi-westlake.com/tutorial/).

    Together, these sources motivate the examples above. The equations and
    limiting cases use the conventions stated here and are checked in the final
    accordion.
    """)
    return


@app.cell
def _(mo):
    mo.md(r"""
    **Previous:** [Module 07 — Chemical Capacitance](https://qiyanglu.github.io/Solid-State-Ionics-Interactive/07-chemical-capacitance/)

    **Reference:** [shared notation and sign conventions](https://github.com/qiyanglu/Solid-State-Ionics-Interactive/blob/main/docs/NOTATION.md)
    """)
    return



if __name__ == "__main__":
    app.run()
