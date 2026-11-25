"""Physical, CAD and rotation-covariance oracles for rational affine images."""
import FreeCAD as App
import Part

import argparse
import json
from pathlib import Path
import subprocess

import numpy as np

from circular_bend import geometry, reference_solid
from coaxial_program import compile_mass_oracle
from perforated_prism import occ_solid

HERE = Path(__file__).resolve().parent


def read(p):
    return json.loads(p.read_text())


def write(p, data):
    p.write_text(json.dumps(data, indent=2, allow_nan=False)+"\n")


def call(cli,*args):
    p=subprocess.run([str(cli),*map(str,args)],capture_output=True,text=True)
    return p.returncode,json.loads(p.stdout)


def corpus():
    box=dict(family="orthogonal_feature_program_v1",base_mm=[[0,12],[0,10],[0,8]],features=[])
    cavity=dict(box,features=[dict(id="void",operation="cut",bounds_mm=[[3,7],[2,9],[2,5]])])
    profile=next(r["recipe"] for r in read(HERE/"perforated_prism_corpus.json") if r["id"]=="two_staggered")
    torus=dict(family="circular_bend_v1",section="full_annulus_v1",major_radius_mm=16,quadrants=4)
    identity=np.eye(3,dtype=int).tolist()
    rot=[[3,-4,0],[4,3,0],[0,0,5]]
    shear=[[2,1,-1],[0,2,1],[0,0,2]]
    yield "translation",box,identity,[7,-3,2],1
    yield "rotation",profile,rot,[13,-9,6],5
    yield "sheared_cavity",cavity,shear,[2,1,-3],2
    yield "anisotropy",profile,[[3,0,0],[0,2,0],[0,0,4]],[-7,5,2],2
    yield "rotated_sheared_torus",torus,[[0,-3,1],[4,1,2],[0,0,5]],[13,-7,4],3
    yield "rotated_sheared_torus_i64",torus,[[0,-3,1],[4,1,2],[0,0,5]],[13,-7,4],3
    yield "torus_overflow_i64",torus,[[64,0,0],[0,64,0],[0,0,64]],[4096,4096,4096],1
    yield "cyclic_axes",cavity,[[0,0,1],[1,0,0],[0,1,0]],[-4,7,2],1


def cad_source(recipe):
    if recipe["family"]=="circular_bend_v1":
        return reference_solid(recipe["section"],recipe["major_radius_mm"],recipe["quadrants"])
    if recipe["family"]=="polygonal_perforated_prism_v1":
        loops=[recipe["outer"],*recipe["holes"]]
        return occ_solid([np.array([v["xy_mm"] for v in l["vertices"]]) for l in loops],recipe)
    def box(b):
        lo,hi=np.array(b).T
        return Part.makeBox(*map(float,hi-lo),App.Vector(*map(float,lo)))
    body=box(recipe["base_mm"])
    for feature in recipe["features"]:
        assert feature["operation"]=="cut"
        body=body.cut(box(feature["bounds_mm"]))
    return body.removeSplitter()


def transformed_torus(recipe,a,shift):
    """Build the two analytic toroidal shells from independent biquadratic patches."""
    quarter=np.array([[1.,0.],[1.,1.],[0.,1.]])
    quarter_turn=np.array([[0.,-1.],[1.,0.]])
    curves=[quarter@np.linalg.matrix_power(quarter_turn,i).T for i in range(4)]
    w=np.array([1.,np.sqrt(.5),1.])
    weights=np.outer(w,w).tolist()
    def torus(minor):
        faces=[]
        for outer in curves:
            for inner in curves:
                pts=[]
                for cs,sn in outer:
                    row=[]
                    for cq,sq in inner:
                        radius=recipe["major_radius_mm"]+minor*cq
                        x=np.array([cs*radius,sn*radius,minor*sq])
                        row.append(App.Vector(*(a@x+shift)))
                    pts.append(row)
                surface=Part.BSplineSurface()
                surface.buildFromPolesMultsKnots(pts,[3,3],[3,3],[0.,1.],[0.,1.],False,False,2,2,weights)
                faces.append(surface.toShape())
        shell=Part.makeShell(faces)
        solid=Part.makeSolid(shell)
        if solid.Volume<0:solid.reverse()
        if not solid.isValid() or not shell.isClosed():raise ValueError("Open or invalid analytic torus shell")
        return solid
    return torus(4.).cut(torus(2.))


def mechanics(cli,folder,model,rotation,offset):
    source=model["source_model"]
    tags=[f["feature"] for f in source["design_faces"]]
    assert "start" in tags and "end" in tags
    base=dict(schema="solidcraft_mechanical_job_v2",material=dict(young_mpa=100000,poisson=.25),
              formulation="displacement",order=5,plan=dict(schema="adaptive_refinement_plan_v1",
              max_depth=0,max_cells=len(source["blocks"]),steps=[]))
    data=[]
    for i,(m,Q,t) in enumerate(((source,np.eye(3),np.zeros(3)),(model,rotation,offset))):
        mp=folder/f"mech-{i}.model.json";write(mp,m)
        def solve(job,name):
            jp=folder/f"{name}.job.json";out=folder/f"{name}.result.json";write(jp,job)
            rc,diag=call(cli,"solve-job",mp,jp,out)
            if rc:raise ValueError(diag)
            return read(out)
        affine=dict(denominator=1000,numerators=[[3,1,1,0],[-2,0,-2,1],[1,2,0,1]])
        patch=solve(dict(base,supports=[dict(select=dict(feature="*"),components=[0,1,2],affine_mm=affine)],loads=[]),f"patch-{i}")
        traction=Q@np.array([1.,2.,-3.])
        loaded=solve(dict(base,supports=[dict(select=dict(feature="start"),components=[0,1,2])],
                         loads=[dict(select=dict(feature="end"),vector_mpa=traction.tolist())]),f"traction-{i}")
        controls=np.array(patch["displacement_mm_xyz"]).reshape(-1,3);error=0.
        for b in patch["analysis_blocks"]:
            h=np.array(b["controls"]);x=h[:,:3]/h[:,3,None]
            exact=(x@np.array(affine["numerators"])[:,1:].T+np.array(affine["numerators"])[:,0])/1000
            error=max(error,float(np.abs(controls[b["dofs"]]-exact).max()))
        data.append((loaded,error))
    a,b=data[0][0],data[1][0]
    def rel(x,y):return float(np.linalg.norm(x-y)/max(np.linalg.norm(y),1e-30))
    r0=np.array(a["applied_resultant_n"]);rotated=rotation@r0
    u=np.array(a["displacement_mm_xyz"]).reshape(-1,3)
    v=np.array(b["displacement_mm_xyz"]).reshape(-1,3)
    errs=dict(affine_control_max_mm=max(d[1] for d in data),work_relative=abs(b["external_work_n_mm"]/a["external_work_n_mm"]-1),
              displacement_relative=rel(v,u@rotation.T),force_relative=rel(np.array(b["applied_resultant_n"]),rotated),
              moment_relative=rel(np.array(b["applied_moment_n_mm"]),rotation@np.array(a["applied_moment_n_mm"])+np.cross(offset,rotated)))
    return dict(errors=errs,qualified=max(errs.values())<=1e-8 and errs["affine_control_max_mm"]<=1e-9)


def run(cli,work):
    work.mkdir(parents=True,exist_ok=True);helper=compile_mass_oracle(work)
    rows=[]
    for name,recipe,ints,t,d in corpus():
        folder=work/name;folder.mkdir(exist_ok=True)
        transform=dict(matrix=ints,translation=t,denominator=d)
        if name.endswith("_i64"):transform["homogeneous_range"]="i64"
        rp=folder/"recipe.json";mp=folder/"model.json";sp=folder/"native.step"
        write(rp,dict(family="affine_image_v1",source=recipe,transform=transform))
        expected_refusal=name=="torus_overflow_i64"
        if expected_refusal:
            previous=(work/"translation/model.json").read_bytes()
            mp.write_bytes(previous)
        rc,audit=call(cli,"build",rp,mp)
        row=dict(id=name,recipe=read(rp),audit=audit,qualified=False);rows.append(row)
        if expected_refusal:
            row.update(expected_admission=False,preserved=mp.read_bytes()==previous)
            row["qualified"]=bool(rc and audit.get("reason")=="affine_image_i64_overflow" and row["preserved"])
            continue
        row["expected_admission"]=True
        if rc:continue
        m=read(mp);a=np.array(ints,dtype=float)/d;shift=np.array(t)/d
        orig=cad_source(recipe);det=np.linalg.det(a)
        expected_volume=orig.Volume*det
        expected_centroid=a@np.array(list(orig.Solids[0].CenterOfMass))+shift
        matrix=np.eye(4);matrix[:3,:3]=a;matrix[:3,3]=shift
        cad_matrix=App.Matrix(*matrix.flatten().tolist())
        try:
            target=orig.transformGeometry(cad_matrix)
            row["occ_transform"]="transformGeometry"
        except Part.OCCError as exc:
            row["initial_occ_transform_error"]=str(exc)
            if recipe["family"]!="circular_bend_v1":raise
            # OCCT's located analytic primitive can integrate mass yet misclassify
            # an affine non-similarity. Construct the analytic image in Cartesian coordinates.
            misplaced=orig.transformed(cad_matrix,True)
            witnesses=[]
            for witness in ([16,0,3],[0,13,0],[0,-13,0]):
                y=a@witness+shift
                witnesses.append(dict(source_mm=witness,transformed_mm=y.tolist(),analytic_source=orig.isInside(App.Vector(*witness),1e-8,False),
                                      located_primitive=misplaced.isInside(App.Vector(*y),1e-8,False)))
            row["located_primitive_classification_witnesses"]=witnesses
            target=transformed_torus(recipe,a,shift)
            row["occ_transform"]="independent toroidal tensor patches with transformed poles"
        rc,step=call(cli,"export-step",mp,sp);row["step"]=step
        if rc:continue
        physical=Part.read(str(sp))
        valid=all(s.isValid() and len(s.Solids)==1 and len(s.Shells)==len(orig.Shells) for s in (target,physical))
        masses=[]
        for i,s in enumerate((target,physical)):
            bp=folder/f"reference-{i}.brep";s.exportBrep(str(bp))
            masses.extend(json.loads(subprocess.check_output([str(helper),str(bp),"--centroid"],text=True)))
        vol=max(abs(q["volume_mm3"]/expected_volume-1) for q in masses)
        cen=max(np.linalg.norm(q["centroid_mm"]-expected_centroid) for q in masses)
        bbox=orig.BoundBox;lo=np.array([bbox.XMin,bbox.YMin,bbox.ZMin]);hi=np.array([bbox.XMax,bbox.YMax,bbox.ZMax])
        sampled=np.random.default_rng(2026112519).uniform(lo-1,hi+1,size=(256,3))
        errors=[0,0];witness=[]
        for x in sampled:
            p=App.Vector(*x);expected=orig.isInside(p,1e-8,False);y=App.Vector(*(a@x+shift))
            values=[s.isInside(y,1e-8,False) for s in (target,physical)]
            for k,v in enumerate(values):errors[k]+=int(v!=expected)
            if values!=[expected,expected] and len(witness)<10:
                witness.append(dict(x=list(x),expected=expected,cad=values))
        mismatches=sum(errors)
        xe=je=0.;min_det=float("inf");uvw=[np.linspace(0.,1.,4)]*3
        for old,new in zip(m["source_model"]["blocks"],m["blocks"]):
            x,j=geometry(old,uvw);y,k=geometry(new,uvw)
            xe=max(xe,float(np.abs(y-(x@a.T+shift)).max()))
            je=max(je,float(np.abs(k-a@j).max()));min_det=min(min_det,float(np.linalg.det(k).min()))
        row.update(valid=valid,masses=masses,volume_relative=vol,centroid_error_mm=cen,
                   points=len(sampled),classification_mismatches=mismatches,classification_by_oracle=errors,witness=witness,max_map_error_mm=xe,
                   max_jacobian_error_mm=je,min_sampled_det_mm3=min_det,
                   qualified=bool(valid and vol<1e-8 and cen<1e-6 and not mismatches and max(xe,je)<1e-10 and
                                  audit["determinant_lower_bound_mm3"]<=min_det*(1+1e-12)))
        if name=="rotation":
            row["mechanics"]=mechanics(cli,folder,m,a,shift)
            row["qualified"] &= row["mechanics"]["qualified"]
        print(name,row["qualified"],vol,cen,flush=True)
    return dict(schema="affine_image_independent_cad_and_covariance_v1",occ=Part.OCC_VERSION,cases=rows)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("cli","work","output"):p.add_argument(k,type=Path)
    args=p.parse_args();report=run(args.cli.resolve(),args.work)
    write(args.output,report)
    if not all(row["qualified"] for row in report["cases"]):raise SystemExit(1)
