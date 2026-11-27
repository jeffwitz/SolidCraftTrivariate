import unittest

import numpy as np

from pocket_lateral_c1_sweep import (
    Patch,
    UnionFind,
    c1_constraints,
    edge_basis_indices,
    geometry,
    knots_for,
    spline_pair,
)


class SplineRefinementTests(unittest.TestCase):
    def test_c0_and_c1_dimension_counts(self):
        for p in (3, 4):
            for n in (1, 2, 4, 8):
                c0 = len(knots_for(p, n, 0))-p-1
                c1 = len(knots_for(p, n, 1))-p-1
                self.assertEqual(c0, p*n+1)
                self.assertEqual(c1, p+1+(n-1)*(p-1))

    def test_basis_partitions_unity_and_derivatives_sum_to_zero(self):
        for p in (3, 4):
            for c in (0, 1):
                knots = knots_for(p, 4, c)
                basis, deriv = spline_pair(knots, p)
                for x in np.linspace(0.0, 1.0, 33):
                    self.assertAlmostEqual(float(np.sum(basis(x))), 1.0, places=12)
                    self.assertAlmostEqual(float(np.sum(deriv(x))), 0.0, places=11)

    def test_c1_derivative_is_continuous_at_each_internal_knot(self):
        eps = 1e-8
        for p in (3, 4):
            knots = knots_for(p, 4, 1)
            _, deriv = spline_pair(knots, p)
            for x in np.unique(knots)[1:-1]:
                left = np.asarray(deriv(float(x-eps)))
                right = np.asarray(deriv(float(x+eps)))
                self.assertLess(float(np.max(np.abs(left-right))), 2e-6)

    def test_c1_interface_rows_cover_each_tangential_span(self):
        p, n = 3, 2
        ku = knots_for(p, n, 1)
        kv = knots_for(p, n, 1)
        nb = len(ku)-p-1
        left = Patch("left", np.array([[0., 0.], [1., 0.], [1., 1.], [0., 1.]]),
                     1., p, n, 1, ku, kv, nb, nb, 0)
        right = Patch("right", np.array([[1., 0.], [2., 0.], [2., 1.], [1., 1.]]),
                      1., p, n, 1, ku, kv, nb, nb, nb*nb)
        patches = {"left": left, "right": right}
        interfaces = [("left", "u1", "right", "u0")]
        uf = UnionFind(2*nb*nb)
        for a, b in zip(edge_basis_indices(left, "u1"), edge_basis_indices(right, "u0")):
            uf.join(a, b)
        roots = sorted({uf.find(i) for i in range(2*nb*nb)})
        root_index = {r: i for i, r in enumerate(roots)}
        mapping = np.array([root_index[uf.find(i)] for i in range(2*nb*nb)])
        constraints = c1_constraints(patches, interfaces, mapping, len(roots), set(), p)
        # p+3 collocation nodes on each of the two tangential spans yields a
        # nonempty reduced system; a represented physical affine field is in
        # its exact nullspace.
        self.assertGreater(constraints.shape[0], 0)
        coeff = np.zeros(len(roots))
        for patch in (left, right):
            gu = np.array([np.mean(patch.knots_u[i+1:i+p+1]) for i in range(nb)])
            gv = np.array([np.mean(patch.knots_v[j+1:j+p+1]) for j in range(nb)])
            values = np.empty(nb*nb)
            for i, u in enumerate(gu):
                for j, v in enumerate(gv):
                    x, y = geometry(patch, float(u), float(v))[0]
                    values[i*nb+j] = x+y
            for local, value in enumerate(values):
                coeff[mapping[patch.offset+local]] = value
        self.assertLess(float(np.max(np.abs(constraints@coeff))), 1e-10)


if __name__ == "__main__":
    unittest.main()
