"""Cross-check the ideal-pair chemical-capacitance identities used in Modules 03, 05–08.

This script intentionally imports no notebook or marimo runtime.  It evaluates
the same formulas on one shared SI-unit parameter set and exits nonzero if an
identity drifts during future edits.
"""

from __future__ import annotations

import math


R_J_PER_MOL_K = 8.314462618
F_C_PER_MOL = 96485.33212


def relative_error(value: float, reference: float) -> float:
    return abs(value - reference) / max(abs(reference), 1.0e-300)


def require_identity(name: str, value: float, reference: float, tolerance: float = 2.0e-14) -> None:
    error = relative_error(value, reference)
    if not math.isfinite(error) or error > tolerance:
        raise AssertionError(f"{name}: relative error {error:.3e} exceeds {tolerance:.1e}")
    print(f"PASS  {name}: relative error {error:.3e}")


def main() -> None:
    temperature_k = 800.0
    pair_concentration_mol_per_m3 = 350.0
    area_m2 = 1.7e-4
    length_m = 80.0e-6
    sigma_i_s_per_m = 1.2e-2
    sigma_e_s_per_m = 8.4e-1

    c_chem_v_f_per_m3 = (
        F_C_PER_MOL**2
        * pair_concentration_mol_per_m3
        / (2.0 * R_J_PER_MOL_K * temperature_k)
    )
    c_chem_total_f = c_chem_v_f_per_m3 * area_m2 * length_m
    c_chem_pair_formula_f = (
        F_C_PER_MOL**2
        * area_m2
        * length_m
        * pair_concentration_mol_per_m3
        / (2.0 * R_J_PER_MOL_K * temperature_k)
    )
    require_identity(
        "Module 05 total capacitance equals Module 07 volumetric storage times volume",
        c_chem_total_f,
        c_chem_pair_formula_f,
    )

    diffusivity_i_m2_per_s = (
        R_J_PER_MOL_K
        * temperature_k
        * sigma_i_s_per_m
        / (F_C_PER_MOL**2 * pair_concentration_mol_per_m3)
    )
    diffusivity_e_m2_per_s = (
        R_J_PER_MOL_K
        * temperature_k
        * sigma_e_s_per_m
        / (F_C_PER_MOL**2 * pair_concentration_mol_per_m3)
    )
    diffusivity_pair_m2_per_s = (
        2.0
        * diffusivity_i_m2_per_s
        * diffusivity_e_m2_per_s
        / (diffusivity_i_m2_per_s + diffusivity_e_m2_per_s)
    )
    sigma_amb_s_per_m = (
        sigma_i_s_per_m
        * sigma_e_s_per_m
        / (sigma_i_s_per_m + sigma_e_s_per_m)
    )
    diffusivity_from_storage_m2_per_s = sigma_amb_s_per_m / c_chem_v_f_per_m3
    require_identity(
        "Modules 03/05/06 ideal-pair diffusivity equals conductivity divided by storage",
        diffusivity_pair_m2_per_s,
        diffusivity_from_storage_m2_per_s,
    )

    r_i_ohm_per_m = 1.0 / (sigma_i_s_per_m * area_m2)
    r_e_ohm_per_m = 1.0 / (sigma_e_s_per_m * area_m2)
    c_chem_f_per_m = c_chem_v_f_per_m3 * area_m2
    diffusivity_tlm_m2_per_s = 1.0 / (
        (r_e_ohm_per_m + r_i_ohm_per_m) * c_chem_f_per_m
    )
    require_identity(
        "Module 08 distributed TLM diffusivity equals the Module 07 transport identity",
        diffusivity_tlm_m2_per_s,
        diffusivity_from_storage_m2_per_s,
    )

    resistance_sum_ohm = (r_e_ohm_per_m + r_i_ohm_per_m) * length_m
    tlm_time_s = resistance_sum_ohm * c_chem_total_f
    diffusion_time_s = length_m**2 / diffusivity_pair_m2_per_s
    require_identity(
        "TLM RC time equals the finite-slab chemical-diffusion time",
        tlm_time_s,
        diffusion_time_s,
    )


if __name__ == "__main__":
    main()
