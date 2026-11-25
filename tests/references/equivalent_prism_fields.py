"""Compare same-domain elasticity fields at independently located physical points."""
import argparse
import json
import time
from pathlib import Path

import numpy as np
from scipy.stats import qmc, t as student

from bracket_field_distance import q2
from equivalent_prism_mechanics import read, save


class Geometry:
    def __init__(self, folder):
        self.blocks = {}
        for block in read(folder / "model.json")["blocks"]:
            h = np.asarray(block["controls"], dtype=float).reshape(2, 2, 2, 4)
            assert block["degrees"] == [1, 1, 1] and np.all(h[..., 3] == h.flat[3])
            points = h[..., :3] / h[..., 3:]
            self.blocks[block["id"]] = points

    def locate(self, x):
        owners = np.full(len(x), -1, dtype=int)
        coords = np.zeros_like(x)
        hits = np.zeros(len(x), dtype=int)
        for index, corners in enumerate(self.blocks.values()):
            bottom = corners[:, :, 0]
            shift = corners[0, 0, 1] - bottom[0, 0]
            assert np.allclose(corners[:, :, 1] - bottom, shift) and shift[2] > 0
            w = (x[:, 2] - bottom[0, 0, 2]) / shift[2]
            p = x[:, :2] - w[:, None] * shift[:2]
            xy = bottom[..., :2]
            subset = np.flatnonzero(np.all((p >= xy.min((0, 1))-1e-10) & (p <= xy.max((0, 1))+1e-10), axis=1))
            if not len(subset):
                continue
            a = xy[0, 0]
            b, c = xy[1, 0]-a, xy[0, 1]-a
            d = xy[1, 1] - a - b - c
            u, v = np.full((2, len(subset)), .5)
            def residual():
                return a + u[:, None]*b + v[:, None]*c + (u*v)[:, None]*d - p[subset]
            for _ in range(20):
                error = residual()
                if np.max(np.abs(error)) < 5e-12:
                    break
                du, dv = b + v[:, None]*d, c + u[:, None]*d
                det = du[:, 0]*dv[:, 1] - du[:, 1]*dv[:, 0]
                det = np.where(np.abs(det) > 1e-14, det, np.inf)
                u -= np.clip((dv[:, 1]*error[:, 0] - dv[:, 0]*error[:, 1])/det, -1, 1)
                v -= np.clip((du[:, 0]*error[:, 1] - du[:, 1]*error[:, 0])/det, -1, 1)
            local = np.column_stack((u, v, w[subset]))
            valid = (np.max(np.abs(residual()), axis=1) < 1e-9) & np.all((local >= -1e-11) & (local <= 1+1e-11), axis=1)
            chosen = subset[valid]
            coords[chosen] = np.clip(local[valid], 0, 1)
            owners[chosen] = index
            hits[chosen] += 1
        if np.any(hits != 1):
            raise ValueError(f"Physical inverse: {sum(hits == 0)} missing, {sum(hits > 1)} multiple")
        return owners, coords


class Field:
    def __init__(self, folder, geometry, level):
        root_index = {k:i for i,k in enumerate(geometry.blocks)}
        self.keys = []
        with np.load(folder / "field.npz") as data:
            self.geom = data["geom"]
            self.coeff = data["field"]
            self.names = data["names"].tolist()
        for name in self.names:
            root = name
            bits = []
            while root not in root_index:
                root,sep,child = root.rpartition("/")
                if not sep or len(child)!=3 or set(child)-set("01"):
                    raise ValueError(f"Unknown source for {name}")
                bits.append(child)
            index = np.zeros(3,dtype=int)
            for child in reversed(bits):
                index = 2*index + np.array(list(child),dtype=int)
            self.keys.append((root,len(bits),index))
        levels = np.array([key[1] for key in self.keys])
        if level is not None and np.any(levels != level):
            raise ValueError("Uniform field level disagrees with cell name")
        self.n = 1 << int(levels.max())
        self.scale = 1 << levels
        self.indices = np.array([key[2] for key in self.keys])
        self.lookup = np.full((len(root_index),self.n,self.n,self.n),-1,dtype=np.int32)
        for i,(root,_,index) in enumerate(self.keys):
            step = self.n//int(self.scale[i])
            slices = tuple(slice(k*step,(k+1)*step) for k in index)
            region = self.lookup[(root_index[root],*slices)]
            if np.any(region != -1):
                raise ValueError("Overlapping adaptive cells")
            region.fill(i)
        if np.any(self.lookup < 0):
            raise ValueError("Incomplete adaptive cover")
        self.physical = read(folder / "physical.json")

    def cell_ids(self, location):
        owners, uvw = location
        indices = np.minimum((uvw*self.n).astype(int), self.n-1)
        return self.lookup[(owners, *indices.T)]

    def evaluate(self, location, x):
        ids = self.cell_ids(location)
        t = location[1]*self.scale[ids,None]-self.indices[ids]
        geom, coeff = self.geom[ids], self.coeff[ids]
        b, db = q2(t)
        weights = geom[:, :, 3]
        cart = geom[:, :, :3]/weights[:, :, None]
        wb = weights*b
        total = wb.sum(axis=1)
        basis = wb/total[:, None]
        dwb = weights[:, :, None]*db
        deriv = (dwb - basis[:, :, None]*dwb.sum(axis=1)[:, None, :])/total[:, None, None]
        pos = np.einsum("qi,qij->qj", basis, cart)
        jac = np.einsum("qic,qia->qca", cart, deriv)
        if np.max(np.abs(pos-x)) > 1e-9 or np.any(np.linalg.det(jac) <= 0):
            raise ValueError("Physical evaluation disagrees with root inverse")
        u = np.einsum("qi,qij->qj", basis, coeff)
        grad = np.einsum("qic,qia->qca", coeff, deriv) @ np.linalg.inv(jac)
        return u, grad


def energy(g):
    e = (g + g.transpose(0, 2, 1))*.5
    # E = 100000 MPa, nu = 1/4 give lambda = mu = 40000 MPa.
    return 40000*np.sum(e*e, axis=(1, 2)) + 20000*np.trace(e, axis1=1, axis2=2)**2


def clouds(geometries, seed=20261125, log2=11):
    batches = []
    for i, corners in enumerate(geometries[0].blocks.values()):
        lo, hi = corners.min((0, 1, 2)), corners.max((0, 1, 2))
        assert len(np.unique(corners.reshape(-1,3), axis=0)) == 8
        assert np.all((corners == lo) | (corners == hi))
        for rep in range(8):
            x = lo + (hi-lo)*qmc.Sobol(3, scramble=True, seed=seed+1000*i+rep).random_base2(log2)
            batches.append((rep, np.prod(hi-lo), x, [g.locate(x) for g in geometries]))
    return batches


def compare(fields, batches, affine, counts=(128,512,2048)):
    total = np.zeros((len(counts), 8, 5))
    manufactured = np.array([[1,1,0],[0,-2,1],[2,0,1]])/1000
    exact_error = 0.
    for rep, volume, x, locations in batches:
        (u,g),(v,h) = [f.evaluate(t, x) for f,t in zip(fields,locations)]
        if affine:
            exact = x@manufactured.T+np.array([3,-2,1])/1000
            exact_error = max(exact_error, float(np.abs(u-exact).max()),float(np.abs(v-exact).max()))
        samples = np.column_stack((energy(g),energy(h),energy(g-h),np.sum((u-v)**2,axis=1),np.sum(v*v,axis=1)))
        for i,n in enumerate(counts):
            total[i,rep] += volume*samples[:n].mean(axis=0)
    metrics = []
    for n, values in zip(counts,total):
        ratios = np.sqrt(values[:,[2,3]]/values[:,[1,4]])
        mean = ratios.mean(axis=0)
        width = student.ppf(.975,7)*ratios.std(axis=0,ddof=1)/np.sqrt(8)
        recovery = values[:,:2].mean(axis=0)/[f.physical["strain_energy_n_mm"] for f in fields] - 1.
        metrics.append(dict(n_per_box_scramble=n,energy_distance=float(mean[0]),l2_distance=float(mean[1]),
                            empirical_half_width_95=width.tolist(),relative_energy_recovery=recovery.tolist()))
    return dict(metrics=metrics,affine_max_mm=exact_error if affine else None)


def run(root):
    start = time.monotonic()
    dirs = [root/tag for tag in ("boxes","contour")]
    geometry = [Geometry(d) for d in dirs]
    batches = clouds(geometry)
    solves = read(root/"solves.json")
    rows = {}
    for name,level in [("affine",0),("h0-q4",0),("h1-q4",1),("h2-q4",2),("h2-q6",2),("h3-q4",3)]:
        folders = [d/name for d in dirs]
        if any(not (d/"field.npz").exists() for d in folders):
            rows[name] = {"status":"unavailable"}
            continue
        fields = [Field(d,g,level) for d,g in zip(folders,geometry)]
        rows[name] = dict(status="completed",**compare(fields,batches,name=="affine"),
                          work=[f.physical["external_work_n_mm"] for f in fields])
        del fields
        print(root.name,name,rows[name]["metrics"][-1],flush=True)
    checks = dict(loading_balance=all(r.get("qualified_equilibrium_and_loading",False) for r in solves) and len(solves)==12,
                  coverage=all(r["status"]=="completed" for r in rows.values()))
    if checks["coverage"]:
        last = [rows[f"h{h}-q4"]["metrics"][-1] for h in (1,2,3)]
        work = rows["h3-q4"]["work"]
        checks.update(affine=rows["affine"]["affine_max_mm"]<=1e-9,
                      recovery=all(max(map(abs,r["metrics"][-1]["relative_energy_recovery"]))<=.015 for r in rows.values()),
                      energy_decay=last[0]["energy_distance"]>last[1]["energy_distance"]>last[2]["energy_distance"],
                      l2_decay=last[0]["l2_distance"]>last[1]["l2_distance"]>last[2]["l2_distance"],
                      final_energy=last[2]["energy_distance"]<=.05,final_work=abs(work[0]/work[1]-1)<=.01,
                      quadrature=all(abs(a/b-1)<=.001 for a,b in zip(rows["h2-q4"]["work"],rows["h2-q6"]["work"])))
    return dict(case=root.name,solves=solves,fields=rows,checks=checks,qualified=all(checks.values()),
                points_located_once=len(batches)*2048,elapsed_s=time.monotonic()-start)


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root",type=Path)
    parser.add_argument("output",type=Path)
    args=parser.parse_args()
    save(args.output,run(args.root))
