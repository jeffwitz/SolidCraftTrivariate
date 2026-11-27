"""Integrate the gallery pocket's continuum energy and surface work by boundary moments."""

import argparse
import json
from fractions import Fraction as F
from pathlib import Path

from sympy import Poly, QQ, symbols

ROOT = Path(__file__).resolve().parents[1]
MODEL = ROOT / "docs/site/data/plate_blind_pocket/plate_blind_pocket.model.json"
T = symbols("t")
HEIGHT = 20
LONGITUDINAL = F(7, 2)  # lambda + 2 mu, with shear modulus one


def polynomial(x):
    return Poly(x, T, domain=QQ)


def integral(poly):
    primitive = poly.integrate()
    return F(primitive.eval(1) - primitive.eval(0))


def boundary(model, feature, layer):
    """CCW orientation of each material section includes clockwise hole edges."""
    blocks = {b["id"]: b for b in model["blocks"]}
    arcs = []
    for face in model["design_faces"]:
        if face["feature"] != feature or not face["block"].endswith(layer):
            continue
        block = blocks[face["block"]]
        assert block["degrees"] == [2, 2, 1]
        net = block["controls"]
        assert len({p[-1] for p in net}) == 1 and net[0][-1] > 0
        triples = {"v0": [0, 3, 6], "u1": [6, 7, 8],
                   "v1": [8, 5, 2], "u0": [2, 1, 0]}
        knots = [net[2*i] for i in triples[face["side"]]]
        coord = []
        for d in range(2):
            a, b, c = (F(p[d], p[3]) for p in knots)
            coord.append(polynomial(a*(1-T)**2 + 2*b*T*(1-T) + c*T**2))
        arcs.append(coord)
    assert len(arcs) == 4
    starts = sorted((x.eval(0), y.eval(0)) for x, y in arcs)
    assert starts == sorted((x.eval(1), y.eval(1)) for x, y in arcs)
    return arcs


def measure(model, power, offset, bx, by, scale):
    """u_z = (z/20)*t^power where t = (offset+bx*x+by*y)/scale."""
    volumes = dict(volume=F(0), squared_displacement=F(0), twice_strain_energy=F(0))
    work = dict(body=F(0), outer_wall=F(0), pocket_wall=F(0), floor=F(0), top=F(0))
    pocket_area = None
    for layer, lo, hi in [("below", 0, 14), ("beside", 14, HEIGHT)]:
        outer = boundary(model, "outer/body", layer)
        inner = boundary(model, "cut/pocket/wall", "beside")
        material = outer if lo == 0 else outer + inner

        def planar(arcs, degree, flux=False):
            # Green's theorem with explicit x primitive; outward n ds=(dy,-dx).
            value = F(0)
            for x, y in arcs:
                t = (polynomial(offset) + bx*x + by*y).mul_ground(F(1, scale))
                if flux:
                    value += integral(t**degree * (bx*y.diff()-by*x.diff())) / scale
                else:
                    value += integral(t**(degree+1) * y.diff()) * F(scale, bx*(degree+1))
            return value

        height_moment = F(hi**3-lo**3, 3*HEIGHT**2)
        norm = F(bx**2+by**2, scale**2)
        moment = planar(material, 2*power)
        lower = planar(material, 2*power-2)
        volumes["volume"] += (hi-lo)*planar(material, 0)
        volumes["squared_displacement"] += height_moment*moment
        volumes["twice_strain_energy"] += (height_moment*power**2*norm*lower
                                            + LONGITUDINAL*moment*F(hi-lo, HEIGHT**2))
        work["body"] -= height_moment*power*(power-1)*norm*lower
        work["outer_wall"] += height_moment*power*planar(outer, 2*power-1, True)
        if lo:
            work["pocket_wall"] += height_moment*power*planar(inner, 2*power-1, True)
            work["top"] += LONGITUDINAL*moment*F(1, HEIGHT)
            work["floor"] -= LONGITUDINAL*planar(inner, 2*power)*F(lo, HEIGHT**2)
            pocket_area = -planar(inner, 0)
    assert sum(work.values()) == volumes["twice_strain_energy"]
    assert volumes["volume"] == 48976 and pocket_area > 0
    return {"integrals": {k: str(v) for k, v in volumes.items()},
            "load_work_by_boundary_N_mm": {k: str(v) for k, v in work.items()},
            "pocket_floor_area_mm2": str(pocket_area)}


def force_controls(model, power, offset, bx, by, scale):
    """Integrate oriented stresses and the first moment of the body load exactly."""
    from functools import lru_cache
    from sympy import Matrix, expand

    x, y, z = symbols("x y z")
    phi = ((offset + bx*x + by*y) / scale)**power
    uz = phi*z/HEIGHT
    grad = Matrix([uz.diff(x), uz.diff(y), uz.diff(z)])
    stress = Matrix.diag(*[F(3, 2)*grad[2]]*3)
    for i in range(3):
        stress[2,i] += grad[i]
        stress[i,2] += grad[i]
    body = -Matrix([sum(stress[i,j].diff(t) for j,t in enumerate((x,y,z))) for i in range(3)])
    moment = Matrix([x,y,z]).cross(body)
    arcs = {"outer":boundary(model,"outer/body","below"), "hole":boundary(model,"cut/pocket/wall","beside")}
    assert arcs["outer"] == boundary(model,"outer/body","beside")

    @lru_cache(maxsize=None)
    def line(group, a, b, direction):
        return sum(integral(xx**a * yy**b * (yy.diff() if direction==0 else -xx.diff()))
                   for xx,yy in arcs[group])

    def integrate(expr, group, direction, axial):
        total = F(0)
        for (i,j,k),v in Poly(expand(expr),x,y,z,domain=QQ).terms():
            section = (line(group,i+1,j,0)/(i+1) if direction is None
                       else line(group,i,j,direction))
            total += F(v)*section*axial(k)
        return total

    vectors = {k:[F(0)]*3 for k in ("body","body_moment","bottom","top","floor","wall")}
    for lo,hi in ((0,14),(14,20)):
        for group in (("outer",) if lo==0 else ("outer","hole")):
            def axial(k):
                return F(hi**(k+1)-lo**(k+1),k+1)
            for i in range(3):
                for name,field in (("body",body),("body_moment",moment)):
                    vectors[name][i] += integrate(field[i],group,None,axial)
                for axis in range(2):
                    vectors["wall"][i] += integrate(stress[i,axis],group,axis,axial)
                if lo==0:
                    vectors["bottom"][i] -= integrate(stress[i,2],group,None,lambda k: 0**k)
                else:
                    vectors["top"][i] += integrate(stress[i,2],group,None,lambda k:hi**k)
                    if group=="hole":
                        vectors["floor"][i] -= integrate(stress[i,2],group,None,lambda k:lo**k)
    assert all(sum(v[i] for k,v in vectors.items() if k!="body_moment")==0 for i in range(3))
    return {k:list(map(str,v)) for k,v in vectors.items()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    model = json.loads(MODEL.read_text())
    states = {"affine": (1, 0, 1, 2, 140), "quintic": (5, 10, 1, 2, 150),
              "affine_shifted": (1, 10, 1, 2, 150), "quadratic_mixed": (2, 10, 1, 2, 150),
              "historical_quadratic": (2, -30, 1, 0, 60)}
    result = {name: {**measure(model, *s), "force_and_body_moment": force_controls(model, *s)}
              for name, s in states.items()}
    data = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(data)
    else:
        print(data, end="")


if __name__ == "__main__":
    main()
