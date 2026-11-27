"""Check the built Sphinx manual as a reader navigating a project subdirectory."""
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import unquote, urlsplit

ROOT = Path(__file__).resolve().parents[1] / "build/manual"


class Links(HTMLParser):
    def __init__(self):
        super().__init__()
        self.links = []

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if tag in {"a", "link", "img", "script"}:
            link = attrs.get("href") or attrs.get("src")
            if link:
                self.links.append(link)


def check():
    pages = sorted(p for p in ROOT.rglob("*.html")
                   if p.relative_to(ROOT).parts[0] not in {"api", "_static"})
    if not pages:
        raise ValueError("No built manual found")
    targets = set()
    for page in pages:
        parser = Links()
        parser.feed(page.read_text())
        for href in parser.links:
            url = urlsplit(href)
            if url.netloc or url.scheme or not url.path:
                continue
            dest = (page.parent / unquote(url.path)).resolve()
            if not dest.is_relative_to(ROOT) or not dest.exists():
                raise ValueError(f"{page.relative_to(ROOT)} links to missing {href}")
            targets.add(dest)
    for required in ["tutorials/index.html", "how-to/index.html",
                     "reference/index.html", "explanation/index.html",
                     "tutorials/solidcraft-plugin.html", "how-to/solidcraft-fieldcraft.html",
                     "reference/fieldcraft-plugin.html", "explanation/plugin-boundaries.html",
                     "api/solidcraft_trivariate/index.html", "search.html"]:
        if ROOT / required not in targets:
            raise ValueError(f"Unreachable manual entry: {required}")
    print(f"Verified {len(pages)} manual pages and {len(targets)} local link targets.")


if __name__ == "__main__":
    check()
