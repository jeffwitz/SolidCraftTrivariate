"""Exact rational pointwise checks on an independently loaded narrow-row model."""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import json
import sys


def determinant(block, u, v, w):
    values = [[(1-u)**2, 2*u*(1-u), u*u], [1-v,v], [1-w,w]]
    derivatives = [[-2*(1-u), 2-4*u, 2*u], [-1,1], [-1,1]]
    homogeneous = [[F(0)]*4 for _ in range(4)]
    for coord, control in zip(product(range(3),range(2),range(2)),block["controls"],strict=True):
        for axis in range(4):
            weight = F(1)
            for i in range(3):
                weight *= (derivatives if axis == i else values)[i][coord[i]]
            for c in range(4):
                homogeneous[axis][c] += weight*control[c]
    pos = homogeneous[3]
    x,y,z = [[(row[i]*pos[3]-row[3]*pos[i])/pos[3]**2 for i in range(3)] for row in homogeneous[:3]]
    return sum(x[i]*(y[(i+1)%3]*z[(i+2)%3]-y[(i+2)%3]*z[(i+1)%3]) for i in range(3))


if __name__ == "__main__":
    model = json.loads(Path(sys.argv[1]).read_text())
    bound = F(sys.argv[2])
    if bound <= 0:
        raise ValueError("strictly positive certificate required")
    nodes = [F(i,4) for i in range(5)]
    grid = [determinant(b,u,v,F(1,2)) for b in model["blocks"] for u,v in product(nodes,repeat=2)]
    assert min(grid) >= bound
    print(json.dumps({"samples":len(grid),"minimum_exact":str(min(grid)),"minimum_mm3":float(min(grid))}))
