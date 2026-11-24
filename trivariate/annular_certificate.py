"""Continuous bounds for the canonical binary64 annular pair."""

from fractions import Fraction
from math import sqrt

from .annular_audit import audit as audit_coefficients


def certify(model):
    """Audit every coefficient before using the canonical map factorization."""
    result = audit_coefficients(model)
    if result["status"] != "valid":
        return result
    heights = model["heights_mm"]
    q = Fraction.from_float(sqrt(0.5))
    q_square_error = abs(q * q - Fraction(1, 2))
    # D is a Bernstein combination of (1,q,1), hence q<=D<=1.
    # X^2+Y^2-D^2=(4q^2-2)s^2t^2 and st<=1/4.
    radial_error = 2 * abs(4 * q * q - 2) / (16 * q * q)
    # At plane n, only the rounded middle homogeneous product differs from q*n.
    # D>=q and 2st<=1/2, hence |2st/D|<=1/(2q).
    levels = (0, heights[0], heights[0] + heights[1])
    plane_error = max(
        abs(Fraction.from_float(float(q * n)) - q * n) / (2 * q)
        for n in levels
    )
    axial_step = min(heights) - 2 * plane_error
    angular_numerator_min = q  # 2q(s^2+t^2)+2st >= q.
    determinant_min = angular_numerator_min * axial_step

    def fraction_json(value):
        return {"numerator": value.numerator, "denominator": value.denominator}

    if radial_error > Fraction(1, 10**12) or 2 * plane_error > Fraction(1, 10**12):
        return {"status": "invalid", "reason": "fidelity_threshold_exceeded"}
    if determinant_min <= 0 or axial_step <= 0:
        return {"status": "invalid", "reason": "nonpositive_jacobian_bound"}
    result["certificate"] = {
        "status": "valid",
        "method": "exact_rational_factor_bounds_with_algebraic_sqrt_enclosure",
        "domain": "two closed unit cubes; 1..1000 mm per block",
        "bounds": {
            "sqrt_half_weight_error": fraction_json(q_square_error / Fraction(7, 5)),
            "radius_error_mm": fraction_json(radial_error),
            "height_error_mm": fraction_json(2 * plane_error),
            "minimum_weight": fraction_json(q),
            "determinant_lower_bound": fraction_json(determinant_min),
            "angular_monotonicity_margin": fraction_json(angular_numerator_min / 2),
            # X+Y>=D, so sqrt(X^2+Y^2)/D >= 1/sqrt(2) > 7/10.
            "radial_monotonicity_margin": fraction_json(Fraction(7, 10)),
            "axial_monotonicity_margin": fraction_json(axial_step),
        },
        "claims": {
            "representation_fidelity": "valid",
            "strictly_positive_jacobian": "valid",
            "block_injectivity": "valid_by_sequential_polar_and_axial_inversion",
            "shared_trace": "valid_exact_coefficient_equality",
            "partition": "valid_two_axial_intervals_with_common_face",
        },
        "assumptions": ["all stored coefficients exactly match canonical binary64 reconstruction",
                        "standard correctly rounded binary64 multiplication"],
    }
    return result
