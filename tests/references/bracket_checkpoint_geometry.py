"""Reconstruct every homogeneous leaf net exactly from the archived root package."""
import argparse
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import struct
import subprocess
import time

import numpy as np


def subdivide(net, bits):
    result = net
    for axis, bit in enumerate(bits):
        a, b, c = np.moveaxis(result, axis, 0)
        middle = a + 2*b + c
        halves = (4*a, 2*(a+b), middle) if bit == "0" else (middle, 2*(b+c), 4*c)
        result = np.moveaxis(np.stack(halves), 0, axis)
    return result


def verify_geometry(package, binary):
    roots = {b["id"]: np.array(b["controls"], dtype=object).reshape(3, 3, 3, 4)
             for b in package["analysis_space"]["blocks"]}

    @lru_cache(None)
    def expected(identifier):
        if "/" not in identifier:
            return roots[identifier]
        parent, bits = identifier.rsplit("/", 1)
        assert len(bits) == 3 and set(bits) <= {"0", "1"}
        return subdivide(expected(parent), bits)

    seen = set()
    first_bad = None
    with binary.open("rb") as f:
        assert f.read(8) == b"SCQ2GEO1"
        count, = struct.unpack("<I", f.read(4))
        for _ in range(count):
            n, = struct.unpack("<H", f.read(2))
            name = f.read(n).decode()
            assert name not in seen
            seen.add(name)
            data = f.read(27 * 4 * 16)
            assert len(data) == 27 * 4 * 16
            actual = np.array([int.from_bytes(data[j:j+16], "little", signed=True)
                               for j in range(0, len(data), 16)], dtype=object)
            target = expected(name).ravel()
            if min(actual[3::4]) <= 0 or np.any(actual * target[3] != target * actual[3]):
                first_bad = first_bad or name
        assert f.read() == b""
    return seen, first_bad


def checkpoint(path):
    """Read the archive in physical leaf order, preserving all binary64 values."""
    rows = {}
    with path.open("rb") as f:
        assert f.read(8) == b"SCQ2ADP1"
        header = struct.unpack("<BBBBI", f.read(8))
        for _ in range(header[-1]):
            n, = struct.unpack("<H", f.read(2))
            key = f.read(n).decode()
            assert key not in rows
            rows[key] = np.frombuffer(f.read(27*4*8), dtype="<f8").reshape(27, 4)
        assert not f.read()
    return header, rows


def run(root, work):
    rows = []
    for rx, ry in ((40,80),(80,40),(50,74),(51,74)):
        case = f"x{rx}-y{ry}"
        folder = work / case
        folder.mkdir(parents=True, exist_ok=True)
        field = folder / "field.bin"
        cmd = ["/usr/bin/taskset", "-c", "2,3", str(root/"target/release/examples/local_fillet_cost"),
               "--field-tier", str(rx), str(ry), "2", "8", str(folder/"report.json"), str(field)]
        start = time.monotonic()
        subprocess.run(cmd, check=True, timeout=1800)
        elapsed = time.monotonic() - start
        start = time.monotonic()
        package = json.loads((root/"build/j4e"/case/"native.package.json").read_text())
        ids, bad = verify_geometry(package, field.with_suffix(".geometry.bin"))
        header, data = checkpoint(field)
        reference = root/"build/j4f"/case/"s2.bin"
        same = hashlib.sha256(field.read_bytes()).digest() == hashlib.sha256(reference.read_bytes()).digest()
        passed = (header[1:4] == (rx,ry,8) and same and bad is None and ids == data.keys() and len(ids) == 2528)
        rows.append(dict(radii_ticks=[rx,ry], leaves=len(ids), homogeneous_coefficients=108*len(ids),
                         first_geometry_mismatch=bad, archived_checkpoint_equal=same,
                         new_run_s=elapsed, independent_geometry_check_s=time.monotonic()-start, passed=passed))
        print(case, passed, flush=True)
    return dict(schema="native_exact_checkpoint_geometry_v1", cases=rows, qualified=all(r["passed"] for r in rows))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("work", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    result = run(Path(__file__).resolve().parents[2], args.work.resolve())
    args.output.write_text(json.dumps(result, indent=2)+"\n")
    raise SystemExit(int(not result["qualified"]))
