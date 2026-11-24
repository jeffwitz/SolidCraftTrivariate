import argparse
import json
import os
from pathlib import Path
import sys
import tempfile

from .audit import audit
from .build import build
from .spline import RationalBezierMap


def main():
    parser = argparse.ArgumentParser()
    commands = parser.add_subparsers(dest="command", required=True)
    commands.add_parser("capabilities")
    create = commands.add_parser("build")
    create.add_argument("recipe", type=Path)
    create.add_argument("output", type=Path)
    commands.add_parser("audit").add_argument("model", type=Path)
    evaluate = commands.add_parser("evaluate")
    evaluate.add_argument("mapping", type=Path)
    evaluate.add_argument("u")
    evaluate.add_argument("v")
    evaluate.add_argument("w")
    args = parser.parse_args()
    try:
        if args.command == "capabilities":
            result = {
                "families": ["affine_pair_v1", "annular_pair_v1", "filleted_bracket_v1"],
                "commands": ["build", "audit", "evaluate", "capabilities"],
                "evaluators": {
                    "rational_bezier_tensor": {
                        "degrees": [[1, 1, 1], [2, 1, 1]],
                        "domain": [[0, 1], [0, 1], [0, 1]],
                    }
                },
            }
        elif args.command == "audit":
            result = audit(json.loads(args.model.read_text()))
        elif args.command == "evaluate":
            mapping = json.loads(args.mapping.read_text())
            if not isinstance(mapping, dict) or set(mapping) != {"degrees", "controls"}:
                raise ValueError("mapping JSON must contain exactly degrees and controls")
            position, columns = RationalBezierMap(
                mapping["degrees"], mapping["controls"]
            ).evaluate(float(args.u), float(args.v), float(args.w))
            result = {
                "status": "evaluated",
                "position": position,
                "jacobian_columns": columns,
            }
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
    except NotImplementedError as error:
        result = {"status": "unsupported", "reason": str(error)}
    except (ValueError, OSError) as error:
        result = {"status": "error", "reason": str(error)}
    print(json.dumps(result, sort_keys=True))
    return int(result.get("status", "valid") not in {"valid", "evaluated"})


if __name__ == "__main__":
    sys.exit(main())
