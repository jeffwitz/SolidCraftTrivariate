#!/usr/bin/env python3
"""Exercise a native SolidCraft desktop sketch overlay and capture its window."""
import argparse
import json
import os
from pathlib import Path
import select
import socket
import subprocess
import sys
import time


class McpClient:
    """Use the installed SolidCraft CLI bridge, not a headless model session."""
    def __init__(self, port: int, env: dict):
        self.process = subprocess.Popen(
            ["solidcraft-cli", "mcp", "--connect", f"127.0.0.1:{port}"],
            stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            text=True, bufsize=1, env=env,
        )
        self.next_id = 1
        self.call("initialize", {"protocolVersion": "2025-06-18", "capabilities": {}, "clientInfo": {"name": "desktop-smoke", "version": "1"}})
        self.notify("notifications/initialized")

    def notify(self, method: str):
        self.process.stdin.write(json.dumps({"jsonrpc": "2.0", "method": method}) + "\n")
        self.process.stdin.flush()

    def call(self, method: str, params: dict) -> dict:
        req_id = self.next_id
        self.next_id += 1
        self.process.stdin.write(json.dumps({"jsonrpc": "2.0", "id": req_id, "method": method, "params": params}) + "\n")
        self.process.stdin.flush()
        while True:
            ready, _, _ = select.select([self.process.stdout], [], [], 40)
            if not ready:
                raise TimeoutError(f"MCP {method} timed out")
            line = self.process.stdout.readline()
            if not line:
                raise RuntimeError(f"MCP bridge exited ({self.process.poll()}): {self.process.stderr.read()[-2000:]}")
            msg = json.loads(line)
            if msg.get("id") == req_id:
                if "error" in msg:
                    raise RuntimeError(f"MCP {method} error: {msg['error']}")
                return msg.get("result", {})

    def execute(self, command: str, params: dict) -> dict:
        reply = self.call("tools/call", {"name": "execute", "arguments": {"command": command, "params": params}})
        texts = [item.get("text", "") for item in reply.get("content", []) if item.get("type") == "text"]
        parsed = [json.loads(text) for text in texts if text.startswith("{")]
        if not parsed or parsed[-1].get("ok") is False:
            raise RuntimeError(f"execute {command} failed: {texts}")
        return parsed[-1]

    def close(self):
        self.process.terminate()
        try:
            self.process.wait(timeout=5)
        except subprocess.TimeoutExpired:
            self.process.kill()
            self.process.wait(timeout=5)


def screenshot(window_id: str, output: Path) -> None:
    raw = output.with_suffix(".xwd")
    with raw.open("wb") as stream:
        subprocess.run(["xwd", "-silent", "-id", window_id, "-out", str(raw)], check=True, stdout=subprocess.DEVNULL)
    subprocess.run(["convert", str(raw), str(output)], check=True)
    raw.unlink()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--timeline", type=Path, default=Path("tests/references/solidcraft_sketch_bracket.timeline.json"))
    parser.add_argument("--kernel", default="target/release/solidcraft-trivariate")
    parser.add_argument("--port", type=int, default=27891)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    design = args.out / "preview-bracket.sc"
    subprocess.run([sys.executable, "scripts/solidcraft_preview_design.py", str(args.timeline), "1", str(design), "--kernel", args.kernel], check=True)

    control_dir = args.out / "control"
    control_dir.mkdir(exist_ok=True)
    token = "a" * 64
    env = os.environ.copy()
    env["SOLIDCRAFT_CONTROL_DIR"] = str(control_dir)
    env["SOLIDCRAFT_CONTROL_TOKEN"] = token
    log = (args.out / "solidcraft.log").open("wb")
    app = subprocess.Popen(["solidcraft", "--new-instance", "--control", str(args.port), str(design)], env=env, stdout=log, stderr=subprocess.STDOUT)
    try:
        ready = False
        deadline = time.monotonic() + 45
        client = None
        while time.monotonic() < deadline and app.poll() is None:
            try:
                client = McpClient(args.port, env)
                client.execute("document.inspect", {})
                ready = True
                break
            except (OSError, RuntimeError, json.JSONDecodeError, TimeoutError):
                if client is not None:
                    client.close()
                    client = None
                time.sleep(0.25)
        if not ready:
            raise RuntimeError(f"desktop control did not become ready (exit={app.poll()})")

        windows = subprocess.run(["xdotool", "search", "--onlyvisible", "--class", "solidcraft"], text=True, capture_output=True).stdout.split()
        if not windows:
            windows = subprocess.run(["xdotool", "search", "--onlyvisible", "--name", "SolidCraft"], text=True, capture_output=True).stdout.split()
        if not windows:
            raise RuntimeError("SolidCraft control is ready but no visible desktop window was found")
        wid = windows[-1]
        activate = subprocess.run(["xdotool", "windowactivate", "--sync", wid], text=True, capture_output=True)
        if activate.returncode:
            # Xvfb has no window manager, so _NET_ACTIVE_WINDOW is commonly unsupported.
            # XSetInputFocus still targets the real SolidCraft X11 window for capture.
            subprocess.run(["xdotool", "windowfocus", "--sync", wid], check=True)
        execute = client.execute
        execute("view.apply", {"view": "top"})
        # Keep the test view unambiguous: hide model bodies and source profile, leaving only
        # the generated native partition sketch to contribute viewport lines.
        execute("browser.visibility", {"folder": "bodies", "visible": False})
        execute("browser.visibility", {"items": [1], "visible": False})
        execute("browser.visibility", {"items": [6], "visible": False})
        time.sleep(1.0)
        screenshot(wid, args.out / "hidden.png")
        execute("browser.visibility", {"items": [6], "visible": True})
        time.sleep(1.0)
        screenshot(wid, args.out / "shown.png")
        execute("browser.visibility", {"items": [6], "visible": False})
        time.sleep(0.5)
        screenshot(wid, args.out / "hidden-again.png")

        from PIL import Image, ImageChops
        hidden = Image.open(args.out / "hidden.png").convert("RGB")
        shown = Image.open(args.out / "shown.png").convert("RGB")
        hidden_again = Image.open(args.out / "hidden-again.png").convert("RGB")
        diff = ImageChops.difference(hidden, shown)
        changed = sum(1 for px in diff.getdata() if max(px) >= 40)
        repeat_diff = ImageChops.difference(hidden, hidden_again)
        repeat_changed = sum(1 for px in repeat_diff.getdata() if max(px) >= 12)
        report = {
            "desktop_window": wid,
            "ready_via_control": True,
            "native_host": "SolidCraft desktop",
            "renderer": "Xvfb + Mesa software rendering",
            "design": str(design),
            "camera": "top",
            "only_visible_geometry": "native sketch feature 6; bodies and source sketch hidden",
            "screenshots": ["hidden.png", "shown.png", "hidden-again.png"],
            "changed_pixels_threshold_40": changed,
            "repeat_changed_pixels_threshold_12": repeat_changed,
        }
        (args.out / "desktop-smoke.json").write_text(json.dumps(report, indent=2) + "\n")
        # Require a broad, clear line contribution and a reversible hide, not a handful of UI pixels.
        if changed < 250 or repeat_changed > 100:
            raise RuntimeError(f"weak/non-reversible native sketch rendering: {report}")
        print(json.dumps(report, indent=2))
        return 0
    finally:
        if "client" in locals() and client is not None:
            client.close()
        app.terminate()
        try:
            app.wait(timeout=10)
        except subprocess.TimeoutExpired:
            app.kill()
            app.wait(timeout=5)
        log.close()


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"desktop smoke failed: {exc}", file=sys.stderr)
        raise
