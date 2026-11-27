"""Independently replay both single cuts and their overlap on the stock SolidCraft host."""
import argparse
import json
from pathlib import Path
import subprocess


ROOT = Path(__file__).resolve().parents[1]
# Below A's mouth, below B's mouth, common void, gap between mouths,
# below both floors, untouched far corner, A below B's floor.
PROBES = [[12,20,39], [28,20,39], [20,20,25], [20,20,39], [22,20,17], [35,35,30], [22,20,20]]
EXPECTED_VOID = {"A": {0,2,6}, "B": {1,2}, "both": {0,1,2,6}}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--work", type=Path, default=ROOT/"build"/"solidcraft-oblique-reference")
    p.add_argument("--host", default="solidcraft-cli")
    args = p.parse_args()
    args.work.mkdir(parents=True, exist_ok=True)
    cmds = json.loads((ROOT/"tests/references/solidcraft-plugin/oblique-cube-host.json").read_text())["commands"]
    first = next(i+1 for i,c in enumerate(cmds) if c["command"] == "solid.loft")
    base = next(i+1 for i,c in enumerate(cmds) if c["command"] == "solid.extrude")
    probe = args.work/"probes.json"
    probe.write_text(json.dumps([{"command":"inspect.classify_points", "params":{"body":"Cube","points":PROBES}}]))
    volumes = {}
    for case, program in [("A",cmds[:first]), ("B",cmds[:base]+cmds[first:]), ("both",cmds)]:
        path = args.work/(case+".json")
        path.write_text(json.dumps(program))
        design = path.with_suffix(".sc")
        subprocess.run([args.host,"run",path,"--save",design,"--quiet"],check=True)
        measured = json.loads(subprocess.check_output([args.host,"eval",design]))
        body, = measured["bodies"]
        if not body["closed"]:
            raise ValueError(f"{case}: open body")
        volumes[case] = body["volume_mm3"]
        result = json.loads(subprocess.check_output([args.host,"run",probe,"--in",design,"--results"]))
        classes = [r["class"] for r in result[0]["results"]]
        expected = ["outside" if i in EXPECTED_VOID[case] else "inside" for i in range(len(PROBES))]
        if classes != expected:
            raise ValueError(f"{case}: material classification {classes} differs from {expected}")
    overlap = 40**3 - volumes["A"] - volumes["B"] + volumes["both"]
    if not 0 < overlap < min(40**3-volumes["A"],40**3-volumes["B"]):
        raise ValueError(f"expected a partial overlap, found {overlap}")
    screenshot = args.work/"both.png"
    subprocess.run([args.host,"snapshot","--in",args.work/"both.sc","--out",screenshot,
                    "--no-axes","--no-grid","--width","900","--height","700"],check=True)
    subprocess.run([args.host,"snapshot","--in",args.work/"both.sc","--out",args.work/"section.png",
                    "--section","y=20,flip","--view","front","--no-axes","--no-grid",
                    "--width","740","--height","600"],check=True)
    print(json.dumps({"material_volumes_mm3":volumes,"overlap_mm3":overlap,"probes":3*len(PROBES)}))


if __name__ == "__main__":
    main()
