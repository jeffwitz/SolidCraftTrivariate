"""Seed an isolated SolidCraft desktop profile with the trivariate sketch ribbon."""
import argparse
import json
from pathlib import Path


def seed(folder):
    hidden = ["MODEL", "SURFACE", "ASSEMBLE", "INSERT", "INSPECT", "MANAGE",
              "VIEW", "SHEET METAL", "WELDMENTS"]
    buttons = {"CREATE": ["sketch.create", "sketch.line", "sketch.rectangle.two_point"],
               "MODIFY": [], "PATTERN": [], "REFERENCE": [], "INSPECT": [], "CONFIGURE": [],
               "CONSTRAIN": ["sketch.constraint.coincident", "sketch.constraint.horizontal_vertical"],
               "FINISH SKETCH": ["sketch.finish", "sketch.cancel"]}
    profile = {"preferences": {"hidden_toolbar_tabs": hidden, "sketch_tab_always": True},
               "toolbar": {"version": 4, "panels": {"DESIGN/SKETCH/" + k: v for k, v in buttons.items()}}}
    folder.mkdir(parents=True, exist_ok=True)
    target = folder / "app.ron"
    # eframe stores a map of string keys to string values in RON.
    # This map uses only ASCII JSON strings, which are also valid RON strings.
    if not target.exists():
        with target.open("x") as stream:
            json.dump({"solidcraft.prefs": json.dumps(profile)}, stream)
    return target


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("folder", type=Path)
    print(seed(parser.parse_args().folder))
