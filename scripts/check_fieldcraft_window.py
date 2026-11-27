"""Open the sampled IGA archive in a separate real FieldCraft desktop and capture its contour."""
import argparse
import json
import socket
import subprocess
import time
from pathlib import Path


def call(port, method, params):
    with socket.create_connection(("127.0.0.1", port), timeout=15) as connection:
        payload = {"id": 1, "method": method, "params": params}
        connection.sendall((json.dumps(payload) + "\n").encode())
        reply = json.loads(connection.makefile("rb").readline())
    if not reply["ok"]:
        raise RuntimeError(reply["error"])
    return reply["result"]


def capture(archive, output):
    with socket.socket() as endpoint:
        endpoint.bind(("127.0.0.1", 0))
        port = endpoint.getsockname()[1]
    app = subprocess.Popen(["fieldcraft", "--control", str(port), str(archive)])
    try:
        for _ in range(120):
            if app.poll() is not None:
                raise RuntimeError("FieldCraft exited before opening the archive")
            try:
                info = call(port, "engine.execute", {"command": "model.info", "params": {}})
                if info.get("nodes", 0):
                    break
            except (ConnectionError, TimeoutError):
                pass
            time.sleep(0.5)
        else:
            raise TimeoutError("no model in the FieldCraft desktop")
        call(port, "engine.execute", {"command": "results.contour",
                                     "params": {"step": "IGA", "field": "displacement"}})
        call(port, "ui.capture", {"path": str(output.resolve()), "width": 1100, "height": 750})
        if output.read_bytes()[:8] != b"\x89PNG\r\n\x1a\n":
            raise RuntimeError("desktop capture is missing")
    finally:
        app.terminate()
        try:
            app.wait(timeout=10)
        except subprocess.TimeoutExpired:
            app.kill()
            app.wait()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    capture(args.archive.resolve(strict=True), args.output)
