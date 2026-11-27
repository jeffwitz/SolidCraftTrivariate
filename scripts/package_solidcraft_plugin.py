"""Build the installable ZIP from the same manifest and component checked by the host."""
import argparse
import json
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile, ZipInfo


def package(staged, destination):
    manifest = json.loads((staged / "plugin.json").read_text())
    module = manifest["runtime"]["module"]
    if Path(module).name != module:
        raise ValueError("the release component must be at the package root")
    destination.mkdir(parents=True, exist_ok=True)
    target = destination / f"{manifest['id']}-{manifest['version']}.zip"
    files = [staged / "plugin.json", staged / module, *sorted((staged / "icons").glob("*.svg"))]
    with ZipFile(target, "w") as archive:
        for file in files:
            entry = ZipInfo(file.relative_to(staged).as_posix(), (1980, 1, 1, 0, 0, 0))
            entry.create_system = 3
            entry.external_attr = 0o100644 << 16
            archive.writestr(entry, file.read_bytes(), compress_type=ZIP_DEFLATED, compresslevel=9)
    return target


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("staged", type=Path)
    parser.add_argument("destination", type=Path)
    args = parser.parse_args()
    print(package(args.staged, args.destination))
