"""Bundle the Linux desktop launcher with its two native executables and sample design."""
import argparse
import json
from pathlib import Path
import tarfile
from zipfile import ZipFile


ROOT = Path(__file__).resolve().parents[1]


def bundle(package, kernel, helper, notices, dependency_licenses, output):
    with ZipFile(package) as archive:
        manifest = json.loads(archive.read("plugin.json"))
    stem = f"solidcraft-trivariate-{manifest['version']}"
    files = {
        package: f"{stem}.zip",
        kernel: "target/release/solidcraft-trivariate",
        helper: "target/release/fieldcraft-iga-archive-bridge",
        ROOT / "plugins/solidcraft-trivariate/README.md": "README.md",
        dependency_licenses: "notices/fieldcraft/THIRD-PARTY-LICENSES.txt",
    }
    for rel in ("scripts/solidcraft_compact_desktop.py", "scripts/solidcraft_compact_profile.py",
                "scripts/fieldcraft_iga_bridge.py", "plugins/solidcraft-trivariate/tests/ribbon-bracket.sc",
                "docs/site/data/fieldcraft_x64_y64/fieldcraft_x64_y64.recipe.json",
                "plugins/solidcraft-trivariate/examples/fieldcraft_bracket_hex_pa5.json"):
        files[ROOT / rel] = rel
    for name in ("LICENSE-MIT", "LICENSE-APACHE", "ATTRIBUTION.md", "NOTICE"):
        files[notices / name] = f"notices/fieldcraft/{name}"
    output.mkdir(parents=True, exist_ok=True)
    target = output / f"{stem}-linux-x86_64.tar.gz"
    with tarfile.open(target, "w:gz") as archive:
        for source, destination in files.items():
            archive.add(source.resolve(strict=True), f"{stem}/{destination}", recursive=False)
    return target


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("package", "kernel", "helper", "notices", "dependency_licenses", "output"):
        parser.add_argument(name, type=Path)
    args = parser.parse_args()
    print(bundle(args.package, args.kernel, args.helper, args.notices, args.dependency_licenses, args.output))
