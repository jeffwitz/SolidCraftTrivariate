"""Open a document in the isolated SolidCraft trivariate desktop profile."""
import argparse
import json
import os
from pathlib import Path
import subprocess

from solidcraft_compact_profile import seed


ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("package", type=Path, help="Installable trivariate plugin ZIP")
    parser.add_argument("document", type=Path, nargs="?",
                        default=ROOT / "plugins/solidcraft-trivariate/tests/ribbon-bracket.sc")
    parser.add_argument("--profile", type=Path, default=ROOT / "build/compact-solidcraft")
    parser.add_argument("--prepare-only", action="store_true", help="Install and verify without opening a window")
    args = parser.parse_args()
    profile = args.profile.resolve()
    document = args.document.resolve(strict=True)
    package = args.package.resolve(strict=True)
    settings = profile / "settings"
    seed(settings)
    companion = profile / "companion"
    companion.mkdir(exist_ok=True)
    env = {**os.environ, "XDG_DATA_HOME": str(profile / "data"), "SOLIDCRAFT_SETTINGS_DIR": str(settings)}

    def command(name, params):
        raw = subprocess.check_output(["solidcraft-cli", "--keep-plugin-approvals", "exec",
                                       name, json.dumps(params)], env=env)
        return json.loads(raw)

    review = command("plugin.install", {"path": str(package), "dry_run": True})
    if review["id"] != "solidcraft-trivariate":
        parser.error("expected the SolidCraft Trivariate plugin package")
    grants = command("plugin.grants", {})["grants"]
    if not any(g["id"] == review["id"] and g["status"] == "current" and
               g["package_hash"] == review["approval"] and "file" in g["command_groups"] and
               any(f["name"] == "companion" and f["folder"] == str(companion) and f["read"] and f["write"]
                   for f in g["folders"]) for g in grants):
        command("plugin.install", {"path": str(package), "approval": review["approval"],
                                   "archive": review["archive"], "replace": True,
                                   "approve": {"command_groups": ["file"], "folders": [
                                       {"name": "companion", "path": str(companion), "read": True, "write": True}]}})
    inventory = command("plugin.list", {})["plugins"]
    if not any(p["manifest"]["id"] == review["id"] and not p["disabled"] for p in inventory):
        parser.error("the trivariate plugin is disabled in this profile")
    if not args.prepare_only:
        os.execvpe("solidcraft", ["solidcraft", "--new-instance", "--view", "model", str(document)], env)


if __name__ == "__main__":
    main()
