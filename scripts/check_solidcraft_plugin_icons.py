#!/usr/bin/env python3
"""Check the SolidCraft Trivariate plugin's 0.48.1 icon contract."""

import json
from pathlib import Path, PurePosixPath
import re
import sys
import xml.etree.ElementTree as ET


EXPECTED = {
    "extrude-active-sketch": "icons/extrude.svg",
    "drill": "icons/drill.svg",
    "boolean": "icons/boolean.svg",
    "chamfer": "icons/chamfer.svg",
    "fillet": "icons/fillet.svg",
    "extrude": "icons/extrude.svg",
    "fieldcraft-iga-solve-display": "icons/compute.svg",
}
SVG_NS = "http://www.w3.org/2000/svg"
ROOT_TAG = f"{{{SVG_NS}}}svg"
ALLOWED_TAGS = {
    ROOT_TAG,
    f"{{{SVG_NS}}}g",
    f"{{{SVG_NS}}}path",
    f"{{{SVG_NS}}}circle",
    f"{{{SVG_NS}}}defs",
    f"{{{SVG_NS}}}mask",
    f"{{{SVG_NS}}}clipPath",
}
SAFE_NAME = re.compile(r"[A-Za-z0-9_.-]+\.svg\Z")


def fail(message: str) -> None:
    raise ValueError(message)


def safe_icon_path(relative: str) -> PurePosixPath:
    path = PurePosixPath(relative)
    if path.is_absolute() or len(path.parts) != 2 or path.parts[0] != "icons":
        fail(f"icon path must be a package-local icons/<name>.svg path: {relative}")
    if not SAFE_NAME.fullmatch(path.name):
        fail(f"icon path contains an unsupported name: {relative}")
    return path


def check_icon(root: Path, relative: str) -> None:
    path = safe_icon_path(relative)
    target = root.joinpath(*path.parts)
    if target.is_symlink() or not target.resolve().is_relative_to(root):
        fail(f"icon must be a regular file inside the plugin package: {relative}")
    if not target.is_file():
        fail(f"manifest references a missing icon: {relative}")
    raw = target.read_bytes()
    if len(raw) > 64 * 1024:
        fail(f"icon exceeds SolidCraft's 64 KiB limit: {relative}")
    try:
        svg = ET.fromstring(raw)
    except ET.ParseError as error:
        fail(f"invalid SVG in {relative}: {error}")
    if svg.tag != ROOT_TAG or svg.get("viewBox") != "0 0 64 64":
        fail(f"{relative} must be an SVG with viewBox='0 0 64 64'")
    elements = list(svg.iter())
    if any(element.tag not in ALLOWED_TAGS for element in elements):
        bad = next(element.tag for element in elements if element.tag not in ALLOWED_TAGS)
        fail(f"{relative} uses an unsupported SVG element: {bad}")
    if not any(element.tag == f"{{{SVG_NS}}}path" and element.get("d", "").strip()
               for element in elements):
        fail(f"{relative} contains no non-empty vector path")
    for element in elements:
        for attribute, value in element.attrib.items():
            local_attribute = attribute.rsplit("}", 1)[-1].lower()
            if local_attribute.startswith("on") or local_attribute in {"href", "src", "style"}:
                fail(f"{relative} contains a link or event handler")
            if attribute in {"fill", "stroke", "color"} and value not in {"none", "currentColor"}:
                fail(f"{relative} uses a non-portable color: {value}")


def main() -> int:
    if len(sys.argv) != 2:
        print(f"usage: {Path(sys.argv[0]).name} PACKAGE_DIR", file=sys.stderr)
        return 2
    root = Path(sys.argv[1]).resolve()
    try:
        manifest = json.loads((root / "plugin.json").read_text(encoding="utf-8"))
        entries = {entry["name"]: entry for entry in manifest["runtime"]["entries"]}
        missing = EXPECTED.keys() - entries.keys()
        if missing:
            fail(f"manifest is missing icon-bearing commands: {', '.join(sorted(missing))}")
        referenced = set()
        for name, expected in EXPECTED.items():
            schema = entries[name].get("input_schema", {})
            actual = schema.get("x-solidcraft-icon")
            if actual != expected:
                fail(f"{name}: expected icon {expected!r}, got {actual!r}")
            check_icon(root, actual)
            referenced.add(actual)
        # Check every additional package image declared by another plugin entry too. This
        # keeps this audit useful when optional integration commands are added later.
        def collect_icons(value):
            if isinstance(value, dict):
                for key, item in value.items():
                    if key == "x-solidcraft-icon":
                        if not isinstance(item, str):
                            fail("x-solidcraft-icon must be a relative path string")
                        yield item
                    else:
                        yield from collect_icons(item)
            elif isinstance(value, list):
                for item in value:
                    yield from collect_icons(item)

        all_declared = set(collect_icons(manifest))
        for image in all_declared:
            safe_icon_path(image)
            check_icon(root, image)
            referenced.add(image)
        icons = {path.relative_to(root).as_posix() for path in (root / "icons").glob("*.svg")}
        if len(icons) > 64:
            fail("the package exceeds SolidCraft's 64-image limit")
        if icons != referenced:
            fail(f"package icon set differs from referenced assets; unused={sorted(icons-referenced)}, "
                 f"missing={sorted(referenced-icons)}")
    except (OSError, KeyError, TypeError, json.JSONDecodeError, ValueError) as error:
        print(f"plugin icon check failed: {error}", file=sys.stderr)
        return 1
    print(f"SolidCraft 0.48.1 plugin icons: {len(icons)} SVGs, {len(EXPECTED)} exact command assignments; paths and portable SVG subset OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
