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


def execute(port, token, command, params=None):
    return control(port, token, "engine.execute", {"command": command, "params": params or {}})


def open_stock(port, token, fixture):
    execute(port, token, "doc.open", {"path": str(fixture)})
    time.sleep(0.75)


def wait_volume(port, token, expected, label):
    for _ in range(120):
        measure = execute(port, token, "inspect.measure")
        if abs(measure["total"]["volume_mm3"] - expected) < 1e-5:
            return measure
        time.sleep(0.25)
    raise RuntimeError(f"{label} volume did not reach its reference: {measure['total']}")


def save_and_reopen(port, token, out, label, expected, env):
    document = out / f"clicked-{label}.sc"
    execute(port, token, "file.save_as", {"path": str(document)})
    reopened = json.loads(subprocess.check_output(["solidcraft-cli", "eval", str(document)],
                                                   env=env))
    if reopened["body_count"] != 1 or abs(reopened["total"]["volume_mm3"] - expected) > 1e-5:
        raise RuntimeError(f"saved {label} differs from its qualified reference: {reopened}")


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
        fixture = Path("plugins/solidcraft-trivariate/tests/face-box-add.sc").resolve()
        app = subprocess.Popen(["solidcraft", "--new-instance", "--control", str(port),
                                str(fixture)], env=env, stdout=log, stderr=log)
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
            # Select a real vertical box edge through the host picker contract.
            # The plugin must turn this into the host's sparse edge reference.
            execute(port, token, "select.set", {"items": [{"type": "edge", "body": "FaceIntegrationBox",
                                                            "index": 8, "point": [0, 0, 6]}]})
            click_handle(port, token, "button:solidcraft-trivariate.fillet")
            time.sleep(1)
            screenshot(windows[-1], out / "dialog.png")
            # The dialog reads the selection and defaults radius_mm to 0.25.
            click_handle(port, token, "plugin_dialog:run")
            fillet_volume = 35999.83904862254
            wait_volume(port, token, fillet_volume, "edge fillet")
            save_and_reopen(port, token, out, "fillet", fillet_volume, env)
            screenshot(windows[-1], out / "fillet-result.png")

            # Reload the unmodified stock, select the top face and chamfer its
            # boundary with the default 0.25 mm distance.
            open_stock(port, token, fixture)
            execute(port, token, "select.set", {"items": [{"type": "face", "body": "FaceIntegrationBox",
                                                            "index": 0,
                                                            "point": [15.94875162046673, 15.94875162046673, 12]}]})
            click_handle(port, token, "button:solidcraft-trivariate.chamfer")
            time.sleep(1)
            screenshot(windows[-1], out / "chamfer-dialog.png")
            click_handle(port, token, "plugin_dialog:run")
            chamfer_volume = 35993.14583333333
            wait_volume(port, token, chamfer_volume, "face chamfer")
            save_and_reopen(port, token, out, "chamfer", chamfer_volume, env)
            screenshot(windows[-1], out / "chamfer-result.png")

            # Finally exercise the face-attached sketch extrusion with its
            # default direction and length, then remeasure after reopening.
            open_stock(port, token, fixture)
            click_handle(port, token, "button:solidcraft-trivariate.extrude")
            time.sleep(1)
            screenshot(windows[-1], out / "face-extrude-dialog.png")
            click_handle(port, token, "plugin_dialog:run")
            # The fixture oracle measures +21.68888888889 mm³ at 2 mm.
            expected_volume = 36000 + (36021.68888888889 - 36000) * (12 / 2)
            try:
                wait_volume(port, token, expected_volume, "face extrusion")
            except RuntimeError:
                screenshot(windows[-1], out / "face-extrude-failure.png")
                raise
            save_and_reopen(port, token, out, "face-extrude", expected_volume, env)
            screenshot(windows[-1], out / "face-extrude-result.png")
        finally:
            app.terminate()
            app.wait(timeout=10)


if __name__ == "__main__":
    main()
