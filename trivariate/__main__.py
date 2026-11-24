import argparse
import json
import os
from pathlib import Path
import sys
import tempfile

from .audit import audit
from .build import build


def main():
    parser = argparse.ArgumentParser()
    commands = parser.add_subparsers(dest="command", required=True)
    commands.add_parser("capabilities")
    create = commands.add_parser("build")
    create.add_argument("recipe", type=Path)
    create.add_argument("output", type=Path)
    commands.add_parser("audit").add_argument("model", type=Path)
    args = parser.parse_args()
    try:
        if args.command == "capabilities":
            result = {"families": ["affine_pair_v1"], "commands": ["build", "audit", "capabilities"]}
        elif args.command == "audit":
            result = audit(json.loads(args.model.read_text()))
        else:
            model = build(json.loads(args.recipe.read_text()))
            result = audit(model)
            if result["status"] == "valid":
                args.output.parent.mkdir(parents=True, exist_ok=True)
                with tempfile.NamedTemporaryFile(dir=args.output.parent, delete=False) as stream:
                    temporary = Path(stream.name)
                try:
                    temporary.write_text(json.dumps(model, indent=2) + "\n")
                    os.replace(temporary, args.output)
                finally:
                    temporary.unlink(missing_ok=True)
    except (ValueError, OSError) as error:
        result = {"status": "error", "reason": str(error)}
    print(json.dumps(result, sort_keys=True))
    return int(result.get("status", "valid") != "valid")


if __name__ == "__main__":
    sys.exit(main())
