"""Apply the recorded EF3 thresholds to the direct and multigrid reference runs."""
import argparse
import json
import math
from pathlib import Path


def summarize(directory):
    loads=lambda path:json.loads(path.read_text())
    direct=loads(directory/"full/result.json")
    calibration=loads(directory/"mg-calibration/result.json")
    fine=loads(directory/"mg-fine/result.json")
    ef=loads(Path(__file__).resolve().parents[2]/"docs/results/bracket-stress-20261125.json")
    reference=min(ef["independent_tetra"],key=lambda item:item["displacement"]["size_mm"])
    tip_ref=reference["displacement"]["point_tip_y_mm"]
    stress_ref=reference["stress"]["von_mises_mpa"]
    records={}
    for run in direct["runs"]+fine["runs"]:
        if run["status"]=="completed":
            records[run["h_mm"],run["order"]]=run["result"]
    checks=[]
    for i,h in enumerate([5.,2.5,1.25,1.]):
        for order in (4,6,8):
            r=records.get((h,order))
            if r is None:
                continue
            previous=records.get(([5.,2.5,1.25,1.][i-1],order)) if i else None
            other=records.get((h,6 if order==8 else 8))
            def rel(x,y):return abs(x-y)/abs(y)
            tip,stress=r["point_tip_y_mm"],r["von_mises_mpa"]
            change_u=rel(tip,previous["point_tip_y_mm"]) if previous else None
            change_s=rel(stress,previous["von_mises_mpa"]) if previous else None
            delta_q=max(rel(tip,other["point_tip_y_mm"]),rel(stress,other["von_mises_mpa"])) if other else None
            volume=rel(r["volume_mm3"],34080.-616.*math.pi)
            force=max(abs(x-y) for x,y in zip(r["force_n"],[0.,-1000.,0.]))/1000.
            balance=max(abs(x+y) for x,y in zip(r["force_n"],r["reaction_n"]))/1000.
            passed=(volume<=1e-8 and max(force,balance)<=1e-6 and r["relative_residual"]<=1e-8
                    and rel(tip,tip_ref)<=.01 and rel(stress,stress_ref)<=.05
                    and change_u is not None and change_u<=.01 and change_s<=.03
                    and delta_q is not None and delta_q<=.005)
            checks.append({"h_mm":h,"order":order,"qualified":passed,"tip_error_relative":rel(tip,tip_ref),
                           "stress_error_relative":rel(stress,stress_ref),"volume_relative":volume,
                           "force_relative":force,"balance_relative":balance,"tip_h_change":change_u,
                           "stress_h_change":change_s,"max_quadrature_qoi_change":delta_q})
    replay=[]
    for run in calibration["runs"]:
        assert run["status"]=="completed"
        r=run["result"]
        baseline=records[(r["h_mm"],r["gauss_order"])]
        errs={key:abs(r[key]/baseline[key]-1.) for key in ("strain_energy_n_mm","point_tip_y_mm","von_mises_mpa")}
        assert max(errs.values())<1e-7
        replay.append({"h_mm":r["h_mm"],"agreement":errs})
    assert any(c["qualified"] for c in checks)
    return {"target":{"reference":"independent OCC/Netgen C3D10 h=1.25 mm","tip_mm":tip_ref,"fillet_vm_mpa":stress_ref},
            "checks":checks,"multigrid_calibration":replay,"transfer_checks":loads(directory/"transfer.json"),
            "direct_series":direct,"multigrid_coarse_series":calibration,"multigrid_fine_series":fine}


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("directory",type=Path)
    print(json.dumps(summarize(p.parse_args().directory),indent=2,allow_nan=False))
