"""Create a SolidCraft document with a native, initially hidden construction sketch."""
import argparse
import json
from pathlib import Path
import subprocess
import tempfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("timeline", type=Path)
    parser.add_argument("sketch_id", type=int)
    parser.add_argument("output", type=Path)
    parser.add_argument("--kernel", default="target/release/solidcraft-trivariate")
    parser.add_argument("--host", default="solidcraft-cli")
    args = parser.parse_args()
    if args.output.suffix != ".sc":
        parser.error("output must be a SolidCraft .sc document")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=args.output.parent) as folder:
        folder = Path(folder)
        descriptor = folder / "preview.json"
        subprocess.run([args.kernel, "sketch-preview", args.timeline, str(args.sketch_id), descriptor], check=True)
        preview = json.loads(descriptor.read_text())
        commands = list(preview["source_script"]["commands"])
        feature = preview["feature"]
        commands.extend([
            {"command": "timeline.add", "params": {"feature": feature}},
            {"command": "browser.visibility", "params": {"items": [feature["id"]], "visible": False}},
        ])
        script = folder / "design.json"
        script.write_text(json.dumps({"commands": commands}))
        design = folder / "design.sc"
        subprocess.run([args.host, "run", script, "--save", design, "--quiet"], check=True)
        design.replace(args.output)
        print(f"{args.output} — toggle {feature['name']} in the SolidCraft browser")


if __name__ == "__main__":
    main()
