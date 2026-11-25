"""Collect raw Kirsch pilot trials and the first stable refinement for each method."""

import argparse
import json
from pathlib import Path
from kirsch_campaign import MANIFEST

ORDER = {"native": (4,8), "tetra": (6,), "immersed": (6,12)}
LEVELS = {"native": ("l1","l2","l3","l4","l5"), "tetra": ("h2","h1","h0.7","h0.5","h0.4"),
          "immersed": ("n4","n8","n16","n32","n64")}


def interval(trial):
    m = trial["metrics"]
    if "hoop_mpa" in m:
        return (m["hoop_mpa"],)*2
    return m["hoop_min_mpa"],m["hoop_max_mpa"]


def summarize(directory):
    manifest=json.loads(MANIFEST.read_text())
    limits=manifest["analysis"]
    goal=limits["reference_hoop_stress"]
    output={"corpus":manifest["version"],"quarter_symmetry":True,"methods":{}}
    for method,levels in LEVELS.items():
        cases={}
        for case in manifest["mandatory"]:
            trials=[]
            refinements=[]
            previous=None
            elapsed=0.
            for level in levels:
                group=[]
                for q in ORDER[method]:
                    path=directory/method/case["id"]/f"{level}-q{q}"/"trial.json"
                    t=json.loads(path.read_text())
                    group.append(t)
                    trials.append(t)
                    elapsed+=t["elapsed_s"]
                if not all(t["status"]=="measured" for t in group):
                    refinements.append(dict(level=level,passes=False,cumulative_s=elapsed))
                    previous=None
                    continue
                stress=interval(group[-1])
                # Compare matching extrema of all one-sided limits.
                variation=(max(abs(x-y)/abs(x) for x,y in zip(stress,previous))
                           if previous is not None else None)
                quadrature=max(abs(x-y)/goal for x,y in zip(stress,interval(group[0])))
                error=max(abs(x/goal-1) for x in stress)
                balance=max(t["metrics"]["force_balance_relative"] for t in group)
                passes=(variation is not None and variation<limits["successive_qoi_relative_change"]
                        and error<=limits["qoi_relative_tolerance"]
                        and balance<=limits["boundary_resultant_relative_tolerance"]
                        and quadrature<.005 and elapsed<limits["wall_clock_budget_seconds_per_method_case"])
                refinements.append(dict(level=level,passes=passes,hoop_range_mpa=stress,
                                        qoi_error_relative=error,successive_change_relative=variation,
                                        quadrature_change_relative=quadrature if len(group)>1 else None,
                                        process_peak_rss_kib=max(t["process_peak_rss_kib"] for t in group),
                                        cumulative_s=elapsed,configuration_s=group[-1]["elapsed_s"]))
                previous=stress
            first=next((r["level"] for r in refinements if r["passes"]),None)
            stable=next((r["level"] for i,r in enumerate(refinements) if all(v["passes"] for v in refinements[i:])),None)
            cases[case["id"]]=dict(status="qualified" if stable else "unqualified",first_passing_level=first,
                                    first_level_with_no_later_failure=stable,refinements=refinements,trials=trials)
            if method=="tetra":
                check=json.loads((directory/method/case["id"]/"h0.4-q12/trial.json").read_text())
                c=cases[case["id"]]
                c["final_quadrature_trial"]=check
                if check["status"]=="measured":
                    delta=max(abs(a-b)/goal for a,b in zip(interval(check),interval(trials[-1])))
                    c["final_quadrature_relative_change"]=delta
                    if delta>=.005:
                        c["status"]="unqualified"
                else:
                    c["status"]="unqualified"
        output["methods"][method]=cases
    return output


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("directory",type=Path)
    p.add_argument("output",type=Path)
    args=p.parse_args()
    result=summarize(args.directory)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2)+"\n")
    for method,cases in result["methods"].items():
        for name,c in cases.items():
            print(method,name,c["first_passing_level"],c["first_level_with_no_later_failure"])
