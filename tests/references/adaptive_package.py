"""Reassemble the locally refined Q2 problem and solve its exported constraints."""
import argparse
from fractions import Fraction
import json
from pathlib import Path
import resource
import subprocess
import sys
import time

# This import fixes BLAS to one thread before importing NumPy and SciPy.
from solver_package import assemble_operators, face_rule, shape
import numpy as np
from scipy import linalg, sparse
from scipy.sparse.linalg import spsolve


def scalar_lift(data, ndof):
    roots = np.asarray(data["free_roots"])
    rows, cols, values = [], [], []
    assert np.array_equal(np.unique(roots),roots)
    assert len(data["lines"]) == ndof
    for i, line in enumerate(data["lines"]):
        assert line["dof"] == i and Fraction(*line["offset"]) == 0
        terms = [(col, Fraction(num,den)) for col,num,den in line["terms"]]
        assert sum(x for _,x in terms) == 1
        assert len(set(col for col,_ in terms))==len(terms)
        for col, x in terms:
            rows.append(i)
            cols.append(col)
            values.append(float(x))
    matrix = sparse.csr_matrix((values,(rows,cols)),shape=(ndof,len(roots)))
    assert (matrix[roots] != sparse.eye(len(roots))).nnz == 0
    return matrix, roots


def leaf_coverage(mesh,block_ids):
    occupied={}
    denom=mesh["grid_denominator"]
    assert {leaf["block_id"] for leaf in mesh["leaves"]}==block_ids
    for leaf in mesh["leaves"]:
        d=2**leaf["level"]
        ijk=np.asarray(leaf["index"])
        bounds=np.asarray(leaf["domain"]["numerator"])/leaf["domain"]["denominator"]
        assert np.array_equal(bounds,np.stack((ijk,ijk+1),axis=1)/d)
        at=occupied.setdefault(leaf["root"],set())
        key=(leaf["level"],*map(int,ijk))
        assert key not in at
        at.add(key)
    for cells in occupied.values():
        volume=Fraction(0)
        for level,x,y,z in cells:
            assert 2**level<=denom
            assert all((level-k,x>>k,y>>k,z>>k) not in cells for k in range(1,level+1))
            volume+=Fraction(1,8**level)
        assert volume==1


def trace_constraints(package, mesh, lift):
    """Recover scalar trace equations by independently sampling paired block faces."""
    blocks={b["id"]:b for b in package["analysis_space"]["blocks"]}
    ndof=lift.shape[0]
    grid=np.asarray([(a,b) for a in (0.,0.5,1.) for b in (0.,0.5,1.)])
    equations=[]
    max_position=0.
    for patch in mesh["interfaces"]:
        rect=np.asarray(patch["first_subrect"],dtype=float)/mesh["grid_denominator"]
        first=rect[:,0]+grid*(rect[:,1]-rect[:,0])
        orient=patch["orientation"]
        second=np.where([orient["flip_first"],orient["flip_second"]],1-first,first)
        if orient["swap"]:
            second=second[:,::-1]
        rect2=np.asarray(patch["second_subrect"])/mesh["grid_denominator"]
        assert np.array_equal(second.min(axis=0),rect2[:,0]) and np.array_equal(second.max(axis=0),rect2[:,1])
        constraints=np.zeros((len(grid),ndof))
        xyz=[]
        for entry,points,sign in [(patch["first"],first,1),(patch["second"],second,-1)]:
            b=blocks[entry["block"]]
            axis="uvw".index(entry["side"][0])
            tangents=[a for a in range(3) if a!=axis]
            leaf=entry["leaf"]
            local=points*2**leaf["level"]-np.asarray(leaf["index"])[tangents]
            assert np.min(local)>=0 and np.max(local)<=1
            uvw,_=face_rule(entry["side"],local)
            values,pos,_,_=shape(uvw,b["controls"])
            constraints[:,b["dofs"]]+=sign*values
            xyz.append(pos)
        max_position=max(max_position,float(np.max(np.abs(xyz[0]-xyz[1]))))
        equations.extend(row for row in constraints if np.max(np.abs(row))>1e-12)
    max_kernel_error=0.
    rank=0
    if equations:
        equations=np.array(equations)
        singular=linalg.svdvals(equations)
        rank=int(np.sum(singular>1e-9))
        max_kernel_error=float(np.max(np.abs(equations@lift)))
    assert max_position<1e-10 and max_kernel_error<1e-10 and rank==ndof-lift.shape[1]
    return dict(rank=rank,unisolvent_rows=len(equations),position_max_mm=max_position,lift_kernel_max=max_kernel_error)


def component_major(lift):
    # Full rows alternate x,y,z, while reduced columns are grouped by component.
    order=np.arange(3*lift.shape[0]).reshape(3,-1).T.ravel()
    return sparse.block_diag([lift]*3,format="csr")[order]


def compare(package,native):
    start=time.monotonic()
    assert package["schema"]=="solidcraft_adaptive_solver_package_v1"
    k,f,coords,tips,vol,area,affine,ferr=assemble_operators(package)
    ndof=len(coords)
    mesh=package["adaptive_mesh"]
    leaf_coverage(mesh,{b["id"] for b in package["analysis_space"]["blocks"]})
    lift,roots=scalar_lift(mesh["exact_lift"],ndof)
    interface=trace_constraints(package,mesh,lift)
    reconstructed=lift@coords[roots]
    geom_error=np.max(np.abs(reconstructed-coords))
    assert geom_error<1e-10
    clamped=np.flatnonzero(np.abs(coords[:,1]-50)<1e-10)
    fixed=package["reduction"]["fixed_vector_dofs"]
    assert np.array_equal(np.sort(fixed),(3*clamped[:,None]+np.arange(3)).ravel())
    live=np.flatnonzero(~np.isin(roots,clamped))
    mounted=component_major(lift[:,live])
    c=package["reduction"]["C"]
    entries=c["unit_entries_row_column_numerator_denominator"]
    r,j,num,den=zip(*entries)
    provided=sparse.coo_matrix(([float(Fraction(n,d)) for n,d in zip(num,den)],(r,j)),shape=c["shape"]).tocsr()
    assert (provided-mounted).nnz==0
    full=component_major(lift)
    body=(full.T@k@full).tocsc()
    modes=np.column_stack([np.tile(a,ndof) for a in np.eye(3)]+[np.cross(a,coords).ravel() for a in np.eye(3)])
    rigid_action=np.linalg.norm(k@modes)/(np.linalg.norm(k.data)*np.linalg.norm(modes))
    # K before constraint elimination need not annihilate interface forces for
    # a strained state, but a rigid motion has zero element strain everywhere.
    assert rigid_action<1e-12
    eig=linalg.eigvalsh(body.toarray(),driver="evd",check_finite=False)
    threshold=eig[-1]*1e-10
    assert np.sum(abs(eig)<threshold)==6 and eig[0]>-threshold and eig[6]>threshold
    lhs=(mounted.T@k@mounted).tocsc()
    rhs=mounted.T@f
    u=mounted@spsolve(lhs,rhs)
    residual=np.linalg.norm(mounted.T@(k@u-f))/np.linalg.norm(rhs)
    tip_values=[float(weights@u[3*ids+1]) for ids,weights in tips]
    assert np.ptp(tip_values)<1e-9
    energy=float(u@(k@u)/2)
    oracle=np.asarray(native["full_interleaved_displacement_mm"])
    uerr=np.linalg.norm(u-oracle)/np.linalg.norm(oracle)
    out=dict(volume_mm3=vol,area_mm2=area,point_tip_y_mm=tip_values[0],strain_energy_n_mm=energy)
    native_values=dict(volume_mm3=native["volume_mm3"],area_mm2=native["x60_area_mm2"],
                       point_tip_y_mm=native["tip_displacement_mm"][1],strain_energy_n_mm=native["strain_energy_n_mm"])
    relative={name:abs(value/native_values[name]-1) for name,value in out.items()}
    assert uerr<1e-7 and max(relative.values())<1e-7 and residual<1e-10
    return dict(**out,trace_checks=interface,full_scalar_dofs=ndof,scalar_masters=len(roots),free_vector_dofs=mounted.shape[1],
                number_of_rigid_modes=6,free_eigenvalue_7=float(eig[6]),spectral_threshold=float(threshold),
                affine_gradient_error=float(affine),affine_lift_error=float(geom_error),
                rigid_action_relative=float(rigid_action),rhs_max_error_n=float(ferr),
                field_relative_error=float(uerr),relative_error=relative,projected_residual_relative=float(residual),
                elapsed_s=time.monotonic()-start,process_peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)


def campaign(root,work,output):
    work.mkdir(parents=True,exist_ok=True)
    cli=root/"target/release/solidcraft-trivariate"
    solve=root/"target/release/examples/solve_adaptive_plan"
    cases=[]
    checks=[]
    def invoke(cmd,ok=True):
        start=time.monotonic()
        run=subprocess.run(list(map(str,cmd)),capture_output=True,text=True)
        if (run.returncode==0)!=ok:
            raise RuntimeError(f"{cmd}: {run.stdout} {run.stderr}")
        return time.monotonic()-start
    for rx in (64,56):
        recipe=work/f"x{rx}.recipe"
        model=work/f"x{rx}.model"
        recipe.write_text(json.dumps(dict(family="parameterized_bracket_v1",radii_sixteenths_mm=[rx,64])))
        invoke([cli,"build",recipe,model])
        for depth in range(3):
            leaf="hole_x_sector_east_z10_30"
            steps=[dict(refine=[dict(root=leaf,level=i,index=([0,0,0] if i==0 else [1,0,0]))]) for i in range(depth)]
            plan=dict(schema="adaptive_refinement_plan_v1",max_depth=4,max_cells=1500,steps=steps)
            name=f"x{rx}-local{depth}"
            pp=work/f"{name}.plan"
            pp.write_text(json.dumps(plan))
            pkg=work/f"{name}.pkg"
            result=work/f"{name}.native"
            report=work/f"{name}.json"
            export_s=invoke([cli,"export-adaptive",model,pp,8,pkg])
            native_s=invoke([solve,model,pp,8,result])
            reader_s=invoke([sys.executable,Path(__file__).resolve(),pkg,result,report])
            row=json.loads(report.read_text())
            cases.append(dict(radii_ticks=[rx,64],depth=depth,export_s=export_s,native_s=native_s,
                              reader_process_wall_s=reader_s,**row))
        prior=pkg.read_bytes()
        invalid_plan=dict(plan,max_cells=36)
        pp.write_text(json.dumps(invalid_plan))
        invoke([cli,"export-adaptive",model,pp,8,pkg],False)
        assert pkg.read_bytes()==prior
        pp.write_text(json.dumps(plan))
        bad_model=json.loads(model.read_text());bad_model["blocks"][0]["controls"][0][0]+=1
        bad=work/f"x{rx}-corrupted.model";bad.write_text(json.dumps(bad_model))
        invoke([cli,"export-adaptive",bad,pp,8,pkg],False)
        assert pkg.read_bytes()==prior
        checks.append(dict(radius_ticks=rx,budget_refusal_preserved_previous_package=True,
                           altered_geometry_preserved_previous_package=True))
    output.write_text(json.dumps(dict(schema="x2-independent-assembly-v1",q=8,cases=cases,
                                      transaction_tests=checks),indent=2,allow_nan=False)+"\n")


if __name__=="__main__":
    if len(sys.argv)>1 and sys.argv[1]=="campaign":
        parser=argparse.ArgumentParser(description="Run all six adaptive package oracles")
        parser.add_argument("work",type=Path)
        parser.add_argument("output",type=Path)
        a=parser.parse_args(sys.argv[2:])
        campaign(Path(__file__).resolve().parents[2],a.work.resolve(),a.output)
    else:
        p=argparse.ArgumentParser(description=__doc__)
        p.add_argument("package",type=Path)
        p.add_argument("native",type=Path)
        p.add_argument("output",type=Path)
        a=p.parse_args()
        report=compare(json.loads(a.package.read_text()),json.loads(a.native.read_text()))
        a.output.write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
