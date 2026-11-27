"""HTML manual for the native kernel and its scientific qualification."""
from pathlib import Path
import os

project = "SolidCraft Trivariate"
author = "SolidCraft Trivariate contributors"
root_doc = "index"
language = "en"
extensions = ["sphinx_rtd_theme", "sphinx.ext.mathjax"]
html_theme = "sphinx_rtd_theme"
html_title = "SolidCraft Trivariate"
html_theme_options = {"navigation_depth": 4, "titles_only": False}
exclude_patterns = ["requirements.txt", "_generated/**"]
html_extra_path = [str(Path(__file__).resolve().parents[2] / "build/docs-extra")]
html_baseurl = os.environ.get("READTHEDOCS_CANONICAL_URL", "")
nitpicky = True
