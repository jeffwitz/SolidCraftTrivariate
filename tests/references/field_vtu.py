"""Read V2 documents with VTK and independently evaluate each one-sided field."""
import argparse
from itertools import product
import json
from pathlib import Path
import subprocess

import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy

from mechanical_job import shape
from projected_bulk import element_matrices, q1_basis

HERE = Path(__file__).resolve().parent


def invoke(*args):
    return subprocess.run(list(map(str, args)), text=True, capture_output=True, check=True)


def read_chart_mesh(path):
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(str(path))
    reader.Update()
    grid = reader.GetOutput()
    if grid.GetNumberOfCells() == 0 or reader.GetErrorCode():
        raise ValueError(f"Unreadable VTU at {path}")
    return grid


def compare(grid, result, n, analytic):
    blocks = result["analysis_blocks"]
    points = vtk_to_numpy(grid.GetPoints().GetData())
    data = grid.GetPointData()
    fields = {k: vtk_to_numpy(data.GetArray(k)) for k in ("uvw", "stress_mpa", "displacement_mm", "von_mises_mpa")}
    block_indices = vtk_to_numpy(grid.GetCellData().GetArray("analysis_block_index"))
    if (len(points), len(block_indices)) != (len(blocks) * (n+1)**3, len(blocks)*n**3):
        raise ValueError("VTK piece counts differ")
    unknowns = np.array(result["displacement_mm_xyz"]).reshape(-1, 3)
    job = result["job"]
    young, nu = (job["material"][k] for k in ("young_mpa", "poisson"))
    mu = young/(2*(1+nu))
    bulk = young/(3*(1-2*nu))
    max_x = max_u = max_s = max_vm = max_a = max_jump = max_stress_jump = 0.0
    min_linear_det = float("inf")
    common_pairs = 0
    coincident = {}
    samples = {}
    for ci, bi in enumerate(block_indices):
        cell = grid.GetCell(ci)
        if cell.GetCellType() != vtk.VTK_HEXAHEDRON:
            raise ValueError("Unexpected linear cell type")
        ids = [cell.GetPointId(i) for i in range(8)]
        if len(set(ids)) != 8:
            raise ValueError("Repeated cell node")
        low = fields["uvw"][ids[0]]
        expected = np.array([(0,0,0),(1,0,0),(1,1,0),(0,1,0),(0,0,1),(1,0,1),(1,1,1),(0,1,1)])/n + low
        if np.max(np.abs(fields["uvw"][ids]-expected)) > 1e-14:
            raise ValueError("Incorrect VTK node orientation")
        p = points[ids]
        for a,b,c in product((0,1),repeat=3):
            # Trilinear Jacobian in physical coordinates at each of the eight corners.
            lookup = np.array([[[0,4],[3,7]],[[1,5],[2,6]]])
            jac = np.array([p[lookup[1,b,c]]-p[lookup[0,b,c]],p[lookup[a,1,c]]-p[lookup[a,0,c]],p[lookup[a,b,1]]-p[lookup[a,b,0]]])
            min_linear_det = min(min_linear_det, float(np.linalg.det(jac)))
        for pi in ids:
            if pi in samples and samples[pi] != bi:
                raise ValueError("Chart boundary was welded")
            samples[pi] = bi
    if min_linear_det <= 0:
        raise ValueError("Inverted linear sample cell")
    pressures = {}
    if job["formulation"] == "projected_bulk_q1":
        for i, block in enumerate(blocks):
            _, _, mass, coupling = element_matrices(block["controls"], young, nu, job["order"])
            pressures[i] = bulk*np.linalg.solve(mass, coupling@unknowns[block["dofs"]].ravel())
    for pi, bi in samples.items():
        block = blocks[bi]
        uvw = fields["uvw"][pi]
        basis, xyz, _, grad, _ = shape(block["controls"], uvw)
        coeff = unknowns[block["dofs"]]
        u = basis @ coeff
        gradient = coeff.T @ grad
        div = np.trace(gradient)
        mean = q1_basis(uvw) @ pressures[bi] if bi in pressures else bulk*div
        stress = mu*(gradient+gradient.T) + np.eye(3)*(mean-2*mu*div/3)
        dev = stress - np.eye(3)*np.trace(stress)/3
        vm = np.sqrt(1.5*np.sum(dev*dev))
        stream_s = fields["stress_mpa"][pi].reshape(3,3)
        stream_u = fields["displacement_mm"][pi]
        max_x = max(max_x, np.linalg.norm(points[pi]-xyz))
        max_u = max(max_u, np.linalg.norm(stream_u-u))
        max_s = max(max_s, np.linalg.norm(stream_s-stress)/max(1.,np.linalg.norm(stress)))
        max_vm = max(max_vm, abs(fields["von_mises_mpa"][pi]-vm)/max(1.,vm))
        if analytic is not None:
            exact = np.array(analytic["gradient"]) @ xyz + np.array(analytic["offset"])
            max_a = max(max_a, np.linalg.norm(stream_u-exact))
        key = tuple(np.rint(xyz*1e9).astype(np.int64))
        for prev_b, prev_x, prev_u, prev_s in coincident.get(key, []):
            if prev_b != bi and np.linalg.norm(prev_x-xyz)<1e-10:
                max_jump = max(max_jump, np.linalg.norm(stream_u-prev_u))
                max_stress_jump = max(max_stress_jump, np.linalg.norm(stream_s-prev_s))
                common_pairs += 1
        coincident.setdefault(key, []).append((bi,xyz,stream_u,stream_s))
    sizes = vtk.vtkCellSizeFilter()
    sizes.SetInputData(grid)
    sizes.ComputeVolumeOn()
    sizes.Update()
    volume = float(vtk_to_numpy(sizes.GetOutput().GetCellData().GetArray("Volume")).sum())
    if max(max_x,max_u,max_jump,max_a) > 1e-9 or max(max_s,max_vm) > 1e-8:
        raise ValueError(f"Physical field discrepancy: {max_x,max_u,max_s,max_vm,max_a,max_jump}")
    return dict(points=len(points), cells=len(block_indices), volume_mm3=volume, common_chart_pairs=common_pairs,
                max_position_mm=max_x, max_displacement_mm=max_u, max_affine_error_mm=max_a,
                relative_stress_error=max_s, relative_von_mises_error=max_vm, max_interface_jump_mm=max_jump,
                max_one_sided_stress_difference_mpa=max_stress_jump, min_linear_corner_determinant_mm3=min_linear_det)


def run(cli, work):
    work.mkdir(parents=True, exist_ok=True)
    rows = []
    for src in sorted((HERE/"mechanical_jobs").glob("*.json")):
        fixture = json.loads(src.read_text())
        if fixture["id"].startswith("refus") or src.name.startswith("refusal"):
            continue
        forms = [fixture["job"]["formulation"]]
        if src.name in ("cube.json","open-bend-affine.json","bracket.json"):
            forms.append("projected_bulk_q1")
        for form in dict.fromkeys(forms):
            case = fixture["id"]+"-"+form
            folder = work/case
            folder.mkdir(exist_ok=True)
            recipe, job, model, result, vtu = (folder/k for k in ("recipe.json","job.json","model.json","result.json","sample.vtu"))
            recipe.write_text(json.dumps(fixture["recipe"]))
            job.write_text(json.dumps(dict(fixture["job"],formulation=form)))
            invoke(cli,"build",recipe,model)
            invoke(cli,"solve-job",model,job,result)
            solution = json.loads(result.read_text())
            analytic = fixture["oracle"] if fixture["oracle"]["kind"] == "affine_displacement" else None
            levels = (1,2,4,8) if src.name == "open-bend-affine.json" and form=="displacement" else (2,)
            samples = []
            for n in levels:
                invoke(cli,"export-field-vtu",result,n,vtu)
                m = compare(read_chart_mesh(vtu),solution,n,analytic)
                samples.append(dict(subdivisions=n,**m))
            if len(samples)>1:
                # Exact toroidal sweep volume: cross section area times centroid arc length.
                r = fixture["recipe"]["major_radius_mm"]
                number = fixture["recipe"]["quadrants"]
                exact = 400*r*number*np.pi/2
                for m in samples:
                    m["linear_volume_relative_error"] = abs(m["volume_mm3"]/exact-1.)
                if not all(b["linear_volume_relative_error"]<a["linear_volume_relative_error"] for a,b in zip(samples,samples[1:])):
                    raise ValueError("Chordal volume refinement failed")
            before = vtu.read_bytes()
            refused = subprocess.run([cli,"export-field-vtu",result,"0",vtu],capture_output=True)
            if not refused.returncode or before!=vtu.read_bytes():
                raise ValueError("Invalid sampling replaced the previous VTU")
            rows.append(dict(case=case,samples=samples,preserved_after_refusal=True))
            print(case,flush=True)
    return {"schema":"mechanical_vtu_field_oracle_v1","vtk_version":vtk.vtkVersion.GetVTKVersion(),"cases":rows}


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    for k in ("cli","work","output"):
        parser.add_argument(k,type=Path)
    a=parser.parse_args()
    a.output.write_text(json.dumps(run(a.cli.resolve(),a.work),indent=2,allow_nan=False)+"\n")
