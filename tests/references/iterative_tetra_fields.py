"""Calibrate a CalculiX iterative solver and replay independent quadratic tetrahedral decks."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import sys
import time

import numpy as np
from bracket_field_distance import TetraField
from refine_tetra_fields import HERE, measure


def displacements(directory):
    output=(directory/"bracket.dat").read_text()
    block=output.split("for set NALL and time")[-1]
    return np.array([[float(x) for x in line.split()] for line in block.splitlines()
                     if len(line.split())==4 and line.split()[0].isdigit()])


def physics(directory):
    text=(directory/"bracket.dat").read_text()
    forces=text.split("for set CLAMP and time")[1].split("displacements")[0]
    rows=[list(map(float,line.split()[1:])) for line in forces.splitlines()
          if len(line.split())==4 and line.split()[0].isdigit()]
    reaction=np.sum(rows,axis=0)
    u=displacements(directory)
    uy={int(row[0]):row[2] for row in u}
    deck=(directory/"bracket.inp").read_text()
    loads=np.fromstring(deck.split("*CLOAD\n")[1].split("*")[0].replace(","," "),sep=" ").reshape(-1,3)
    assert np.all(loads[:,1]==2) and abs(loads[:,2].sum()+1000)<1e-8
    work=sum(uy[int(n)]*f for n,_,f in loads)
    tetra=TetraField(directory)
    tip,_,inside,_,_=tetra.evaluate(np.array([[60.,4.,20.]]))
    if not inside[0]:
        raise ValueError("Tip not covered by tetrahedral mesh")
    return {"mean_tip_y_mm":-work/1000,"point_tip_y_mm":float(tip[0,1]),
            "strain_energy_n_mm":work/2,"reaction_n":reaction.tolist(),
            "reaction_relative_error":float(np.linalg.norm(reaction-[0,1000,0])/1000),
            "nodes":len(u),"tet10":len(tetra.x)}


def solver(source,work,ccx,env,iterative,seconds):
    work.mkdir(parents=True,exist_ok=True)
    deck=(source/"bracket.inp").read_text()
    assert deck.count("*STATIC\n")==1
    if "*NODE PRINT,NSET=Nall" not in deck:
        deck=deck.replace("*END STEP","*NODE PRINT,NSET=Nall\nU\n*END STEP")
    if iterative:
        deck=deck.replace("*STATIC\n","*STATIC,SOLVER=ITERATIVE CHOLESKY\n")
    (work/"bracket.inp").write_text(deck)
    run=measure([str(ccx),"-i","bracket"],env,work/"ccx.log",cwd=work,seconds=seconds)
    run["input_sha256"]=hashlib.sha256(deck.encode()).hexdigest()
    log=(work/"ccx.log").read_text(errors="replace")
    run["solver_log"]=[l.strip() for l in log.splitlines() if any(s in l.lower() for s in
                       ("system of equations using","preconditioning","iteration","error","tolerance"))][-18:]
    if run["status"]=="completed":
        if iterative and "Cholesky preconditioning" not in log:
            run["status"]="requested_preconditioner_unverified"
        else:
            result=physics(work)
            (work/"result.json").write_text(json.dumps(result,indent=2)+"\n")
            run["result"]=result
    return run


def report_rows(text, quantity, group):
    import io
    import re
    header=re.compile(r"^\s*"+quantity+r"[^\n]*for set "+group+r" and time[^\n]*$",re.M|re.I)
    hit=header.search(text)
    if not hit:
        raise ValueError(f"Missing {quantity} on {group}")
    lines=text[hit.end():].lstrip().split("\n\n",1)[0]
    return np.loadtxt(io.StringIO(lines),ndmin=2)


def defect_corrected(source,work,ccx,env,seconds,max_corrections,retained_first=None):
    start=time.monotonic()
    work.mkdir(parents=True,exist_ok=True)
    original=(source/"bracket.inp").read_text()
    start_load=original.index("*CLOAD\n")+len("*CLOAD\n")
    end_load=original.index("*",start_load)
    loads=np.fromstring(original[start_load:end_load].replace(","," "),sep=" ").reshape(-1,3)
    clamp_text=original.split("*NSET,NSET=CLAMP\n")[1].split("*")[0].strip()
    clamp=np.fromstring(clamp_text.replace(","," "),sep=" ",dtype=int)
    deck=original[:end_load]+"*NODE PRINT,NSET=Nall\nU,RF\n*END STEP\n"
    deck=deck.replace("*STATIC\n","*STATIC,SOLVER=ITERATIVE CHOLESKY\n")
    passes=[]
    accumulated_u=accumulated_rf=None
    residuals=[]
    for iteration in range(1+max_corrections):
        path=work/f"pass-{iteration+1}"
        path.mkdir(exist_ok=True)
        if iteration==0 and retained_first is not None:
            first,record=retained_first
            if record["status"]!="completed" or (first/"bracket.inp").read_text()!=deck:
                raise ValueError("retained first correction deck or status does not match")
            path=first
            run=dict(record,reused=True,displacement_force_sha256=hashlib.sha256((first/"bracket.dat").read_bytes()).hexdigest())
        else:
            if iteration==1 and retained_first is not None:
                earlier=retained_first[0].parent/"pass-2/bracket.inp"
                if earlier.read_text()!=deck:raise ValueError("regenerated correction load differs from the timed-out input")
            (path/"bracket.inp").write_text(deck)
            run=measure([str(ccx),"-i","bracket"],env,path/"ccx.log",cwd=path,
                        seconds=max(0,seconds-(time.monotonic()-start)))
        log=(path/"ccx.log").read_text(errors="replace")
        run["solver_log"]=[l.strip() for l in log.splitlines() if any(x in l for x in
                           ("Cholesky preconditioning","error condition","# of iterations"))]
        passes.append(run)
        if run["status"]!="completed" or "Cholesky preconditioning" not in log or "error condition (0=good, 1=bad) = 0" not in log:
            return {"status":"solver_error","passes":passes}
        dat=(path/"bracket.dat").read_text()
        u=report_rows(dat,"displacements","NALL")
        rf=report_rows(dat,"forces","NALL")
        assert np.array_equal(u[:,0],rf[:,0])
        if iteration==0:
            ids=u[:,0].astype(int)
            assert np.all(np.diff(ids)>0)
            free=~np.isin(ids,clamp)
            applied=np.zeros((len(ids),3))
            load_dof=np.searchsorted(ids,loads[:,0].astype(int))
            assert np.array_equal(ids[load_dof],loads[:,0].astype(int))
            np.add.at(applied,(load_dof,loads[:,1].astype(int)-1),loads[:,2])
            accumulated_u=u[:,1:].copy()
            accumulated_rf=rf[:,1:].copy()
        else:
            assert np.array_equal(ids,u[:,0])
            accumulated_u+=u[:,1:]
            accumulated_rf+=rf[:,1:]
        defect=applied-accumulated_rf
        residuals.append(float(np.linalg.norm(defect[free])/np.linalg.norm(applied[free])))
        if residuals[-1]<=1e-6:
            break
        if iteration<max_corrections:
            # The prescribed row reactions must never enter the correction RHS.
            node,axis=np.nonzero((defect!=0)&free[:,None])
            correction="".join(f"{ids[i]},{d+1},{defect[i,d]:.10e}\n" for i,d in zip(node,axis))
            deck=deck[:deck.index("*CLOAD\n")+7]+correction+deck[deck.index("*NODE PRINT"):]
    # Write the physically superposed results as a separate derived data file.
    # The two original solver output files and both input decks stay intact.
    (work/"bracket.inp").write_text(original)
    with (work/"bracket.dat").open("w") as out:
        out.write(" forces (fx,fy,fz) for set CLAMP and time 1.0\n\n")
        np.savetxt(out,np.column_stack((ids[~free],accumulated_rf[~free])),fmt=["%d","%.16e","%.16e","%.16e"])
        out.write("\n displacements (vx,vy,vz) for set NALL and time 1.0\n\n")
        np.savetxt(out,np.column_stack((ids,accumulated_u)),fmt=["%d","%.16e","%.16e","%.16e"])
    data=physics(work)
    (work/"result.json").write_text(json.dumps(data,indent=2)+"\n")
    return {"status":"completed" if residuals[-1]<=1e-6 else "residual_above_limit",
            "passes":passes,"free_residual_relative":residuals,"result":data,
            "elapsed_s":time.monotonic()-start,"input_sha256":hashlib.sha256(original.encode()).hexdigest()}


def exercise(root,prefix,corrected=False,max_corrections=1):
    root=root.resolve()
    build=root.parent.parent
    ccx=prefix/"bin/ccx"
    env=dict(os.environ,LD_LIBRARY_PATH=str(prefix/"lib"),OMP_NUM_THREADS="2",OPENBLAS_NUM_THREADS="1")
    root.mkdir(parents=True,exist_ok=True)
    metrics=[]
    def save(row):
        metrics.append(row)
        dest=root/"metrics.tmp"
        dest.write_text(json.dumps(metrics,indent=2)+"\n")
        dest.replace(root/"metrics.json")
        print(row["stage"],row.get("size_mm"),row.get("qualified"),flush=True)
    for size in (2.,1.25):
        source=build/f"ef1-h{size:g}"
        direct=source
        earlier=build/f"field-comparison/tet-iterative/direct-h{size:g}"
        if corrected and (earlier/"result.json").exists():
            direct=earlier
        row={"stage":"calibration","size_mm":size}
        if "for set NALL and time" not in (direct/"bracket.dat").read_text():
            direct=root/f"direct-h{size:g}"
            row["direct"]=solver(source,direct,ccx,env,False,1800)
            if row["direct"]["status"]!="completed":
                row["qualified"]=False;save(row);return
        work=root/f"cal-h{size:g}"
        row["iterative"]=defect_corrected(source,work,ccx,env,1800,max_corrections) if corrected else solver(source,work,ccx,env,True,1800)
        row["qualified"]=False
        if row["iterative"]["status"]=="completed":
            a,b=displacements(direct),displacements(work)
            assert np.array_equal(a[:,0],b[:,0])
            nodal=float(np.linalg.norm(a[:,1:]-b[:,1:])/np.linalg.norm(a[:,1:]))
            dp=physics(direct);ip=row["iterative"]["result"]
            delta={k:abs(dp[k]-ip[k])/abs(dp[k]) for k in ("mean_tip_y_mm","point_tip_y_mm","strain_energy_n_mm")}
            probe=[]
            for loc in (direct,work):
                out=root/(loc.name+"-calibration-probe")
                post=measure([sys.executable,str(HERE/"fillet_probe.py"),str(loc),str(out),"--ccx",str(ccx)],
                             dict(os.environ,OPENBLAS_NUM_THREADS="1"),root/(out.name+".log"))
                if post["status"]!="completed":
                    row["probe_error"]=post;save(row);return
                probe.append(json.loads((out/"result.json").read_text())["von_mises_mpa"])
            stress=abs(probe[1]/probe[0]-1.)
            row.update(nodal_relative=nodal,scalar_relative=delta,fillet_stress_mpa=probe,stress_relative=stress)
            row["qualified"]=bool(max(nodal,*delta.values(),ip["reaction_relative_error"])<=1e-5 and stress<=1e-3)
        save(row)
        if not row["qualified"]:return
    for ticks in (64,56):
        for size in (1.,.9):
            start=time.monotonic()
            source=build/f"field-comparison/tet-fine/x{ticks}-h{size:g}"
            row={"stage":"fine","radius_x_ticks":ticks,"size_mm":size}
            work=root/f"x{ticks}-h{size:g}"
            if not (source/"bracket.inp").exists():
                source=root/f"deck-x{ticks}-h{size:g}"
                row["mesh"]=measure([str(prefix/"bin/python"),str(HERE/"calculix_bracket.py"),
                    "--ccx",str(ccx),"--deck-only","--store-displacements","--radius-x",str(ticks/16),
                    "--size",str(size),str(source)],dict(env,PYTHONPATH=str(prefix/"lib")),root/(source.name+".log"))
                if row["mesh"]["status"]!="completed":
                    row["qualified"]=False;save(row);break
            seconds=max(0,1800-(time.monotonic()-start))
            row["iterative"]=defect_corrected(source,work,ccx,env,seconds,max_corrections) if corrected else solver(source,work,ccx,env,True,seconds)
            row["qualified"]=row["iterative"]["status"]=="completed"
            if row["qualified"]:
                out=work/"probe"
                row["post"]=measure([sys.executable,str(HERE/"fillet_probe.py"),str(work),str(out),"--ccx",str(ccx)],
                                      dict(os.environ,OPENBLAS_NUM_THREADS="1"),work/"probe.log",seconds=max(0,1800-(time.monotonic()-start)))
                row["qualified"]=bool(row["post"]["status"]=="completed" and row["iterative"]["result"]["reaction_relative_error"]<=1e-5)
                if row["post"]["status"]=="completed":row["probe"]=json.loads((out/"result.json").read_text())
            row["elapsed_s"]=time.monotonic()-start
            save(row)
            if not row["qualified"]:break


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory",type=Path)
    parser.add_argument("--freecad-prefix",type=Path,required=True)
    parser.add_argument("--defect-correction",action="store_true")
    parser.add_argument("--max-corrections",type=int,choices=(1,2,3),default=1)
    args=parser.parse_args()
    exercise(args.directory,args.freecad_prefix.resolve(),args.defect_correction,args.max_corrections)
