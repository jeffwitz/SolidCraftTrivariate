"""Tensor-product rational Bezier maps for small trivariate patches.

Controls are homogeneous ``(wx, wy, wz, w)`` coefficients in millimetres,
with the third (``w``) index varying fastest.  Derivative columns are returned
in parameter-axis order ``(u, v, w)``; face restrictions therefore retain
the order of their remaining parent axes.
"""

from math import comb, isfinite


_SUPPORTED_DEGREES = {(1, 1, 1), (2, 1, 1)}


def _bernstein(degree, parameter):
    """Return the Bernstein basis of ``degree`` at ``parameter``."""
    return tuple(
        comb(degree, index)
        * parameter**index
        * (1.0 - parameter) ** (degree - index)
        for index in range(degree + 1)
    )


def _bernstein_derivative(degree, parameter):
    """Return derivatives of the Bernstein basis of ``degree``."""
    if degree == 0:
        return (0.0,)
    lower = _bernstein(degree - 1, parameter)
    return tuple(
        degree * ((lower[index - 1] if index else 0.0) - (lower[index] if index < degree else 0.0))
        for index in range(degree + 1)
    )


def _number(value):
    if not isinstance(value, (int, float)) or isinstance(value, bool):
        return False
    try:
        return isfinite(float(value))
    except (OverflowError, ValueError):
        return False


class RationalBezierMap:
    """A rational tensor-product Bezier map on the closed unit cube.

    Args:
        degrees: Nonnegative integer degrees in ``(u, v, w)`` order.
        controls: Flat homogeneous coefficients ``(wx, wy, wz, weight)``;
            the ``w`` control index varies fastest, then ``v``, then ``u``.
            Coordinates are in millimetres.

    A malformed degree tuple or coefficient array raises :class:`ValueError`.
    A well-formed degree family outside the currently supported families
    raises :class:`NotImplementedError`.
    """

    def __init__(self, degrees, controls):
        try:
            degree_values = tuple(degrees)
        except TypeError as error:
            raise ValueError("degrees must be three nonnegative integers") from error
        if (
            len(degree_values) != 3
            or any(type(degree) is not int or degree < 0 for degree in degree_values)
        ):
            raise ValueError("degrees must be three nonnegative integers")

        try:
            coefficient_values = tuple(tuple(point) for point in controls)
        except (TypeError, ValueError) as error:
            raise ValueError("controls must be flat homogeneous quadruples") from error
        expected_count = (degree_values[0] + 1) * (degree_values[1] + 1) * (degree_values[2] + 1)
        if len(coefficient_values) != expected_count or any(
            len(point) != 4 for point in coefficient_values
        ):
            raise ValueError("control count or homogeneous coefficient dimension is invalid")
        if any(not _number(value) for point in coefficient_values for value in point):
            raise ValueError("homogeneous coefficients must be finite real numbers")
        coefficient_values = tuple(tuple(float(value) for value in point) for point in coefficient_values)
        if any(point[3] <= 0 for point in coefficient_values):
            raise ValueError("homogeneous weights must be strictly positive")
        if degree_values not in _SUPPORTED_DEGREES:
            raise NotImplementedError(f"unsupported degree family {degree_values!r}")

        self.degrees = degree_values
        self.controls = coefficient_values

    def evaluate(self, u, v, w):
        """Evaluate the map and its analytic Jacobian columns at ``(u,v,w)``.

        Returns ``(xyz, (du, dv, dw))``.  All three parameters must be finite
        and lie in the closed unit interval.
        """
        parameters = (u, v, w)
        if any(not _number(value) or not 0.0 <= value <= 1.0 for value in parameters):
            raise ValueError("parameters must be finite and lie in the unit cube")
        parameters = tuple(float(value) for value in parameters)
        bases = tuple(_bernstein(degree, value) for degree, value in zip(self.degrees, parameters))
        derivative_bases = tuple(
            _bernstein_derivative(degree, value)
            for degree, value in zip(self.degrees, parameters)
        )
        homogeneous = [0.0, 0.0, 0.0, 0.0]
        derivatives = [[0.0, 0.0, 0.0, 0.0] for _ in range(3)]
        du, dv, dw = self.degrees

        for i in range(du + 1):
            for j in range(dv + 1):
                for k in range(dw + 1):
                    index = (i * (dv + 1) + j) * (dw + 1) + k
                    point = self.controls[index]
                    basis = bases[0][i] * bases[1][j] * bases[2][k]
                    for component in range(4):
                        homogeneous[component] += basis * point[component]

                    factors = (
                        derivative_bases[0][i] * bases[1][j] * bases[2][k],
                        bases[0][i] * derivative_bases[1][j] * bases[2][k],
                        bases[0][i] * bases[1][j] * derivative_bases[2][k],
                    )
                    for column, factor in zip(derivatives, factors):
                        for component in range(4):
                            column[component] += factor * point[component]

        if any(not isfinite(value) for value in homogeneous) or any(
            not isfinite(value) for derivative in derivatives for value in derivative
        ):
            raise ValueError("evaluation overflowed finite double precision")
        weight = homogeneous[3]
        if weight <= 0.0:
            raise ValueError("evaluated homogeneous weight is not positive in double precision")
        xyz = tuple(value / weight for value in homogeneous[:3])
        columns = tuple(
            tuple(
                (derivative[component] - xyz[component] * derivative[3]) / weight
                for component in range(3)
            )
            for derivative in derivatives
        )
        if any(not isfinite(value) for value in xyz) or any(
            not isfinite(value) for column in columns for value in column
        ):
            raise ValueError("evaluation produced a nonfinite result")
        return xyz, columns
