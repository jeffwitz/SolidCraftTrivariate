"""Read native chart diagnostics in VTK and verify rational derivatives with NumPy SVD."""
import argparse
import copy
from itertools import product
import json
import math
from pathlib import Path
import subprocess

import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy

HERE=Path(__file__).resolve().parent


def charts():
    def archive(file,id):
        return next(r["recipe"] for r in json.loads((HERE/file).read_text()) if r["id"]==id)
    base=json.loads((HERE/"affine_pair.json").read_text())
    yield "affine_pair",base
    yield "annular_pair",json.loads((HERE/"annular_pair.json").read_text())
    yield "bracket",json.loads((HERE/"filleted_bracket.json").read_text())
    perforated=archive("perforated_prism_corpus.json","two_staggered")
    yield "two_tunnels",perforated
    yield "anisotropic",dict(family="affine_image_v1",source=dict(family="reference_rectangle_v1"),
                             transform=dict(matrix=[[32,7,0],[0,1,0],[0,0,1]],translation=[1,2,3],denominator=4))
    yield "torus",dict(family="circular_bend_v1",section="full_annulus_v1",major_radius_mm=16,quadrants=4)
    yield "inclined_bores",archive("circular_polygon_corpus.json","extreme_shear")


def basis(degree,t):
    return np.column_stack([math.comb(degree,j)*t**j*(1-t)**(degree-j) for j in range(degree+1)])


def analytic(b,uvw):
    degree=b.get("degrees",[1,1,1])
    H=np.array(b["controls"],dtype=float).reshape(*(np.array(degree)+1),4)
    val=[basis(d,t) for d,t in zip(degree,uvw.T)]
    def tensor(v,h):return np.einsum("pi,pj,pk,ijkc->pc",*v,h,optimize=True)
    h=tensor(val,H)
    position=h[:,:3]/h[:,3,None]
    jac=[]
    for k,d in enumerate(degree):
        pol=val.copy();pol[k]=basis(d-1,uvw[:,k])
        diff=tensor(pol,d*np.diff(H,axis=k))
        jac.append((diff[:,:3]-position*diff[:,3,None])/h[:,3,None])
    return position,np.stack(jac,axis=-1)


def run(cli,folder):
    folder.mkdir(parents=True,exist_ok=True)
    rec=[]
    for name,recipe in charts():
        rp=folder/(name+".recipe");mp=folder/(name+".model")
        rp.write_text(json.dumps(recipe))
        subprocess.run([cli,"build",rp,mp],check=True,capture_output=True)
        model=json.loads(mp.read_text())
        for n in (1,3,6):
            vt=folder/f"{name}-{n}.vtu"
            subprocess.run([cli,"export-geometry-vtu",mp,str(n),vt],check=True,capture_output=True)
            reader=vtk.vtkXMLUnstructuredGridReader();reader.SetFileName(str(vt));reader.Update()
            grid=reader.GetOutput()
            assert not reader.GetErrorCode() and grid.GetNumberOfCells()==len(model["blocks"])*n**3
            actual=vtk_to_numpy(grid.GetPoints().GetData())
            data=grid.GetPointData()
            uvw=vtk_to_numpy(data.GetArray("uvw"));sidx=vtk_to_numpy(grid.GetCellData().GetArray("geometry_block_index"))
            quality={k:vtk_to_numpy(data.GetArray(k)) for k in ("singular_values_mm","condition_number","scaled_jacobian","jacobian_determinant_mm3")}
            assert np.array_equal(sidx,np.repeat(np.arange(len(model["blocks"])),n**3))
            assert len(uvw)==len(model["blocks"])*(n+1)**3
            errors=dict(position_mm=0.,singular_values_mm=0.,condition_number=0.,scaled_jacobian=0.,jacobian_determinant_mm3=0.)
            expected_uvw=np.array(list(product(np.linspace(0,1,n+1),repeat=3)))
            for j,b in enumerate(model["blocks"]):
                sub=slice(j*(n+1)**3,(j+1)*(n+1)**3)
                assert np.allclose(uvw[sub],expected_uvw,rtol=0,atol=2e-16)
                x,J=analytic(b,uvw[sub]);sv=np.linalg.svd(J,compute_uv=False);det=np.linalg.det(J)
                vals=dict(singular_values_mm=sv,condition_number=sv[:,0]/sv[:,-1],
                          scaled_jacobian=det/np.prod(np.linalg.norm(J,axis=1),axis=1),jacobian_determinant_mm3=det)
                errors["position_mm"]=max(errors["position_mm"],float(abs(x-actual[sub]).max()))
                for key,value in vals.items():
                    errors[key]=max(errors[key],float(abs(quality[key][sub]/value-1).max()))
                for cell_id in (j*n**3,(j+1)*n**3-1):
                    c=grid.GetCell(cell_id)
                    assert c.GetCellType()==vtk.VTK_HEXAHEDRON
                    indices=[c.GetPointId(k) for k in range(8)]
                    U=uvw[indices];basis_edges=U[[1,3,4]]-U[0]
                    assert np.allclose(basis_edges,np.eye(3)/n)
                    assert len(set(indices))==8
            rec.append(dict(model=name,subdivisions=n,blocks=len(model["blocks"]),points=len(actual),cells=len(sidx),
                            maximum_errors=errors,condition_range=[float(quality["condition_number"].min()),float(quality["condition_number"].max())],
                            passed=errors["position_mm"]<1e-10 and max(v for k,v in errors.items() if k!="position_mm")<1e-9))
    mp=folder/"inclined_bores.model";bad=copy.deepcopy(json.loads(mp.read_text()))
    bad["blocks"][0]["controls"][0][0]+=1
    bp=folder/"altered.model";bp.write_text(json.dumps(bad))
    checks=[]
    for event,path,n in [("altered",bp,3),("budget",mp,32),("out_of_range",mp,0)]:
        out=folder/"protected.vtu";out.write_bytes(b"accepted")
        cmd=subprocess.run([cli,"export-geometry-vtu",path,str(n),out],capture_output=True,text=True)
        checks.append(dict(case=event,refused=cmd.returncode!=0,preserved=out.read_bytes()==b"accepted",diagnostic=cmd.stdout.strip()))
    return dict(schema="geometry_vtu_independent_oracle_v1",vtk=vtk.vtkVersion.GetVTKVersion(),cases=rec,refusals=checks,
                qualified=all(r["passed"] for r in rec) and all(r["refused"] and r["preserved"] for r in checks))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("cli","work","output"):p.add_argument(k,type=Path)
    a=p.parse_args();d=run(a.cli.resolve(),a.work);a.output.write_text(json.dumps(d,indent=2,allow_nan=False)+"\n")
    print("Geometry VTU",d["qualified"])
    raise SystemExit(int(not d["qualified"]))
