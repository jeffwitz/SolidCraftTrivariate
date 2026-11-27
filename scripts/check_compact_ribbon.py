"""Verify the buttons actually drawn in the stock SolidCraft window under Xvfb."""
import argparse
import json
import os
from pathlib import Path
import secrets
import socket
import subprocess
import time

from solidcraft_desktop_smoke import screenshot


def control(port, token, method, params=None):
    request = {"id": 1, "method": method, "params": params or {}, "token": token}
    with socket.create_connection(("127.0.0.1", port), timeout=30) as client:
        client.sendall((json.dumps(request) + "\n").encode())
        reply = json.loads(client.makefile("rb").readline())
    if not reply.get("ok"):
        raise RuntimeError(f"{method}: {reply.get('error')}")
    return reply["result"]


def click_handle(port, token, name):
    for _ in range(40):
        rows = control(port, token, "ui.handles", {"prefix": name})
        hit = next((r for r in rows if r["name"] == name), None)
        if hit:
            control(port, token, "ui.click", {"x": hit["x"], "y": hit["y"]})
            return
        time.sleep(0.15)
    raise RuntimeError(f"unavailable desktop control: {name}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("package", type=Path)
    parser.add_argument("out", type=Path)
    args = parser.parse_args()
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=True)
    token = secrets.token_hex(32)
    env = dict(os.environ, XDG_DATA_HOME=str(out / "data"), SOLIDCRAFT_SETTINGS_DIR=str(out / "settings"),
               SOLIDCRAFT_CONTROL_DIR=str(out / "control"), SOLIDCRAFT_CONTROL_TOKEN=token)
    subprocess.run(["python3", "scripts/solidcraft_compact_desktop.py", str(args.package.resolve()),
                    "--profile", str(out), "--prepare-only"], check=True)
    port = 27893
    with (out / "desktop.log").open("w") as log:
        app = subprocess.Popen(["solidcraft", "--new-instance", "--control", str(port),
                                "plugins/solidcraft-trivariate/tests/ribbon-bracket.sc"], env=env, stdout=log, stderr=log)
        try:
            expected = {"solidcraft-trivariate." + name for name in
                        ("drill", "boolean", "fillet", "chamfer", "extrude")}
            names = set()
            for _ in range(90):
                if app.poll() is not None:
                    raise RuntimeError("SolidCraft exited before exposing its ribbon")
                try:
                    names = {row["name"] for row in control(port, token, "ui.handles")}
                except (OSError, ValueError, RuntimeError):
                    pass
                visible = {s.removeprefix("button:") for s in names if s.startswith("button:")}
                if expected <= visible:
                    break
                time.sleep(0.5)
            windows = subprocess.check_output(["xdotool","search","--onlyvisible","--name","SolidCraft"],text=True).split()
            screenshot(windows[-1], out / "ribbon.png")
            (out / "handles.json").write_text(json.dumps(sorted(names),indent=2))
            print("Visible commands:", sorted(visible), flush=True)
            if not expected <= visible or "tab:SKETCH" not in names:
                raise RuntimeError("missing compact ribbon buttons: " + str(sorted(expected-visible)))
            unwanted = {v for v in visible if not v.startswith("sketch.") and v not in expected}
            tabs = names & {"tab:"+s for s in "MODEL SURFACE ASSEMBLE INSERT MANAGE VIEW INSPECT WELDMENTS".split()+["SHEET METAL"]}
            if unwanted or tabs:
                raise RuntimeError(f"unexpected ribbon commands or tabs: {unwanted}, {tabs}")
            source = control(port, token, "engine.execute", {"command": "timeline.export_script", "params": {}})
            original = source["script"]["commands"][0]["params"]["feature"]["step"]
            click_handle(port, token, "button:solidcraft-trivariate.fillet")
            click_handle(port, token, "plugin_dialog:run")
            # The GUI applies plugin edits asynchronously. Observe the change on
            # the original feature before saving and independently reopening it.
            for _ in range(120):
                export = control(port, token, "engine.execute", {"command": "timeline.export_script", "params": {}})
                steps = export["script"]["commands"]
                if len(steps) == 1 and steps[0]["params"]["feature"]["step"] != original:
                    break
                time.sleep(0.25)
            else:
                raise RuntimeError("the Fillet ribbon click did not update the bracket import")
            document = out / "clicked-fillet.sc"
            control(port, token, "engine.execute", {"command": "file.save_as", "params": {"path": str(document)}})
            measurement = json.loads(subprocess.check_output(["solidcraft-cli", "eval", str(document)],env=env))
            if measurement["body_count"] != 1 or abs(measurement["total"]["volume_mm3"] - 32143.411718185) > 1e-6:
                raise RuntimeError("reopened GUI fillet differs from the qualified trivariate result")
        finally:
            app.terminate()
            app.wait(timeout=10)


if __name__ == "__main__":
    main()
