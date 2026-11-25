"""Independent volume integration of the multi-span JSON control nets."""
import itertools
import json
import math
import sys
from pathlib import Path

import numpy as np


def basis(degree, t):
    if degree == 1:
        return np.stack([1-t, t], axis=1), np.tile([-1, 1], (len(t), 1))
    assert degree == 2
    return (np.stack([(1-t)**2, 2*t*(1-t), t*t], axis=1),
            np.stack([2*t-2, 2-4*t, 2*t], axis=1))


def volume(chart, order):
    nodes, weights = np.polynomial.legendre.leggauss(order)
    basis_pairs = [basis(p, (nodes+1)/2) for p in chart['degrees']]
    coeff = np.asarray(chart['xyzw']).reshape(*chart['shape'], 4)
    knots = [np.asarray(k) for k in chart['knots']]
    active = [[i for i in range(len(k)-1) if k[i] < k[i+1]] for k in knots]
    total = 0.
    for ijk in itertools.product(*active):
        controls = coeff[tuple(slice(i-p, i+1) for i, p in zip(ijk, chart['degrees']))]
        values = []
        for axis in [-1, 0, 1, 2]:
            matrices = [pair[int(j == axis)] for j, pair in enumerate(basis_pairs)]
            values.append(np.einsum('ai,bj,ck,ijkh->abch', *matrices, controls, optimize=True))
        h = values[0]
        xyz = h[..., :3]/h[..., 3:]
        jac = np.stack([(v[..., :3]-v[..., 3:]*xyz)/h[..., 3:] for v in values[1:]], axis=-1)
        total += np.einsum('i,j,k,ijk', weights/2, weights/2, weights/2, np.linalg.det(jac))
    return total


if __name__ == '__main__':
    # Order of the three files is sharp, chamfer, fillet.
    for name, path, target in zip(('sharp', 'chamfer', 'fillet'), sys.argv[1:],
                                   (96*math.pi, 299.2367002544223, 300.5816074146435), strict=True):
        model = json.loads(Path(path).read_text())
        v12, v16 = [sum(volume(c, q) for c in model['charts']) for q in (12, 16)]
        assert max(abs(v16-target), abs(v16-v12)) < 1e-8, (name, v12, v16, target)
        print(f'{name}: {len(model["charts"])} charts, V12={v12:.12f}, V16={v16:.12f}, reference={target:.12f}')
