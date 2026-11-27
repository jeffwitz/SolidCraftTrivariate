"""Exercise the installed SolidCraft guest, native solver and FieldCraft results reader."""
import argparse
import json
import math
import os
from pathlib import Path
import subprocess
import sys

from check_fieldcraft_roundtrip import check


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("kit", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    kit, output = args.kit.resolve(), args.output.resolve()
    (package,) = kit.glob("solidcraft-trivariate-*.zip")
    profile = output / "profile"
    companion = profile / "companion"
    launch = [sys.executable, kit / "scripts/solidcraft_compact_desktop.py",
              package, "--profile", profile]
    subprocess.run([*launch, "--prepare-only", "--external-consumer"], cwd=kit, check=True)
    env = {**os.environ, "XDG_DATA_HOME": str(profile / "data"),
           "SOLIDCRAFT_SETTINGS_DIR": str(profile / "settings")}

    def host(command, params):
        return json.loads(subprocess.check_output(["solidcraft-cli", "--keep-plugin-approvals",
            "exec", command, json.dumps(params)], env=env))

    grants = host("plugin.grants", {})["grants"]
    grant, = [g for g in grants if g["id"] == "solidcraft-trivariate" and g["status"] == "current"]
    params = {"path": str(profile/"data/solidcraft/plugins/solidcraft-trivariate"),
              "approval": grant["package_hash"]}
    host("plugin.run", {**params, "entry": "fieldcraft-iga-solve-display",
                        "params": {"case_name": "loaded-bracket"}})
    subprocess.run([*launch, "--process-once"], cwd=kit, check=True)
    status = json.loads((companion/"fieldcraft-analysis-status.json").read_text())
    if status.get("state") != "complete":
        raise RuntimeError(status)
    archive = Path(status["archive"])
    case = archive.parent
    result = json.loads((case/"loaded-bracket.result.json").read_text())
    if (result["scalar_dofs"] != 171 or result["relative_force_balance"] > 1e-7
            or not math.isclose(result["support_reaction_n"][0], -2000, abs_tol=1e-4)
            or not math.isclose(result["strain_energy_n_mm"], 135.44475386, rel_tol=1e-7)):
        raise RuntimeError("unexpected native bracket response")
    check(os.environ.get("FIELDCRAFT_CLI", "fieldcraft-cli"), archive,
          case/"loaded-bracket.field.vtu", output/"readback")


if __name__ == "__main__":
    main()
