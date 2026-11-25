#!/usr/bin/env python3
"""Build, monitor, and compare the A6 enriched-field campaigns."""

import argparse
import json
import os
from pathlib import Path
import platform
import selectors
import subprocess
import sys
import time


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT = ROOT / "build/enriched-bracket"
DEFAULT_REPORT = ROOT / "docs/results/enriched-bracket-20261125.json"
FIELD_REFERENCE = ROOT / "build/field-comparison"
DISTANCE_ORACLE = ROOT / "tests/references/adaptive_field_distance.py"
MIB = 1024 * 1024
GIB = 1024 * MIB


def process_rss_bytes(pid):
    try:
        status = Path(f"/proc/{pid}/status").read_text()
    except OSError:
        return None
    for line in status.splitlines():
        if line.startswith("VmRSS:"):
            return int(line.split()[1]) * 1024
    return None


def run_monitored(binary, radius, directory, rounds, memory_limit, solve_limit):
    process = subprocess.Popen(
        [str(binary), str(radius), str(directory), str(rounds)],
        cwd=ROOT,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        bufsize=0,
    )
    selector = selectors.DefaultSelector()
    selector.register(process.stdout, selectors.EVENT_READ)
    os.set_blocking(process.stdout.fileno(), False)
    buffer = b""
    stage = None
    stage_start = None
    stage_records = []
    peak_rss = 0
    started = time.monotonic()
    stop_reason = None
    log_path = directory / "run.log"
    directory.mkdir(parents=True, exist_ok=True)
    with log_path.open("wb") as log:
        while process.poll() is None:
            for key, _ in selector.select(timeout=0.2):
                chunk = os.read(key.fd, 65536)
                if chunk:
                    buffer += chunk
                    while b"\n" in buffer:
                        line, buffer = buffer.split(b"\n", 1)
                        log.write(line + b"\n")
                        sys.stdout.buffer.write(line + b"\n")
                        sys.stdout.buffer.flush()
                        text = line.decode(errors="replace")
                        if text.startswith("A6_RESOLUTION_START|"):
                            stage, stage_start = text, time.monotonic()
                        elif text.startswith("A6_RESOLUTION_END|"):
                            elapsed = time.monotonic() - stage_start if stage_start else None
                            stage_records.append({"start": stage, "end": text,
                                                  "monitored_wall_s": elapsed})
                            stage, stage_start = None, None
            rss = process_rss_bytes(process.pid)
            if rss is not None:
                peak_rss = max(peak_rss, rss)
                if rss > memory_limit:
                    stop_reason = "rss_limit_exceeded"
                    process.kill()
                    break
            if stage_start is not None and time.monotonic() - stage_start > solve_limit:
                stop_reason = "per_resolution_time_limit_exceeded"
                process.kill()
                break
        if buffer:
            log.write(buffer)
            sys.stdout.buffer.write(buffer)
            sys.stdout.buffer.flush()
        process.wait()
    return {
        "returncode": process.returncode,
        "stop_reason": stop_reason,
        "campaign_wall_s": time.monotonic() - started,
        "sampled_peak_rss_bytes": peak_rss,
        "memory_limit_bytes": memory_limit,
        "per_resolution_time_limit_s": solve_limit,
        "resolution_stages": stage_records,
        "log": str(log_path.relative_to(ROOT)),
    }


def run_distance(radius, adaptive_dir, log2):
    command = [
        sys.executable,
        str(DISTANCE_ORACLE),
        str(FIELD_REFERENCE),
        str(adaptive_dir / "fields"),
        "--radius-ticks",
        str(radius),
        "--log2",
        str(log2),
    ]
    result = subprocess.run(command, cwd=ROOT, capture_output=True, text=True, check=False,
                            env={**os.environ, "OPENBLAS_NUM_THREADS": "1"})
    if result.returncode:
        return {"status": "failed", "returncode": result.returncode,
                "stdout": result.stdout[-4000:], "stderr": result.stderr[-4000:]}
    return {"status": "complete", **json.loads(result.stdout)}


def should_extend_sampling(result):
    levels = result.get("levels", [])
    if len(levels) < 2:
        return False
    previous, last = levels[-2]["methods"], levels[-1]["methods"]
    for method, current in last.items():
        prior = previous[method]
        for key in ("mean_l2_relative", "mean_energy_relative"):
            value = current[key]
            delta = abs(value - prior[key])
            if delta > max(0.001, 0.05 * abs(value)):
                return True
    return False


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--rounds", type=int, default=6)
    parser.add_argument("--memory-limit-gib", type=float, default=6.0)
    parser.add_argument("--resolution-limit-seconds", type=float, default=1800.0)
    args = parser.parse_args()
    output = args.output if args.output.is_absolute() else ROOT / args.output
    report_path = args.report if args.report.is_absolute() else ROOT / args.report
    build = subprocess.run(
        ["rtk", "cargo", "build", "--release", "--example", "enriched_bracket"],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    if build.returncode:
        raise SystemExit(build.stderr or build.stdout)
    binary = ROOT / "target/release/examples/enriched_bracket"
    revision = subprocess.run(["rtk", "git", "rev-parse", "HEAD"], cwd=ROOT,
                              capture_output=True, text=True, check=True).stdout.strip()
    rustc = subprocess.run(["rtk", "rustc", "--version"], cwd=ROOT,
                           capture_output=True, text=True, check=True).stdout.strip()
    report = {
        "contract": "A6 eb9b568 + A2/A4/A5",
        "contract_base_revision": "eb9b568",
        "workspace_head_revision": revision,
        "worktree_has_local_changes": True,
        "environment": {
            "platform": platform.platform(),
            "python": sys.version,
            "rustc": rustc,
            "distance_oracle": str(DISTANCE_ORACLE.relative_to(ROOT)),
            "distance_reference": "native h3/q8 field; eight Sobol scrambles with fixed seeds",
        },
        "build_command": "rtk cargo build --release --example enriched_bracket",
        "build_stdout": build.stdout,
        "build_stderr": build.stderr,
        "limits": {
            "rounds": args.rounds,
            "max_accepted_leaves": 4000,
            "max_trial_leaves": 16000,
            "rss_bytes": int(args.memory_limit_gib * GIB),
            "per_resolution_wall_seconds": args.resolution_limit_seconds,
        },
        "radii": {},
    }
    for radius in (64, 56):
        directory = output / f"x{radius}"
        run = run_monitored(binary, radius, directory, args.rounds,
                            int(args.memory_limit_gib * GIB), args.resolution_limit_seconds)
        record = {"run": run}
        complete_path = directory / "completed.json"
        if complete_path.exists():
            record["campaign"] = json.loads(complete_path.read_text())
        elif (directory / "progress.json").exists():
            record["campaign"] = json.loads((directory / "progress.json").read_text())
        if run["returncode"] == 0 and run["stop_reason"] is None:
            distance = run_distance(radius, directory, 11)
            if distance.get("status") == "complete" and should_extend_sampling(distance):
                distance["extended_sampling_triggered"] = True
                distance = run_distance(radius, directory, 13)
            else:
                distance["extended_sampling_triggered"] = False
            record["field_distance"] = distance
        else:
            record["field_distance"] = {"status": "not_run_due_to_campaign_stop"}
        report["radii"][str(radius)] = record
        report_path.parent.mkdir(parents=True, exist_ok=True)
        report_path.write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n")
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print(json.dumps({
        radius: {
            "stop_reason": value["run"]["stop_reason"],
            "returncode": value["run"]["returncode"],
            "peak_rss_bytes": value["run"]["sampled_peak_rss_bytes"],
            "distance_status": value["field_distance"].get("status"),
        }
        for radius, value in report["radii"].items()
    }, indent=2))


if __name__ == "__main__":
    main()
