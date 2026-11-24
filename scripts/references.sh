#!/usr/bin/env bash
set -eu
cd "$(dirname "$0")/.."
mode=${1:-check}
case "$mode" in check|clone) ;; *) exit 2 ;; esac
while read -r repo rev dirs; do
    [[ "$repo" == "#"* || -z "$repo" ]] && continue
    dir="upstream/${repo##*/}"
    if [[ ! -e "$dir" && "$mode" == clone ]]; then
        git clone --filter=blob:none --no-checkout "https://github.com/$repo.git" "$dir"
        if [[ "$dirs" != "-" ]]; then git -C "$dir" sparse-checkout set $dirs; fi
        git -C "$dir" checkout --detach "$rev"
    fi
    [[ -d "$dir/.git" && "$(git -C "$dir" rev-parse HEAD)" == "$rev" && -z "$(git -C "$dir" status --porcelain)" ]] || {
        echo "Missing, modified or unexpected checkout: $dir (expected $rev)" >&2
        exit 1
    }
    echo "$dir $rev"
done < references.lock
test -f upstream/solidcraft/solidcore/docs/coordinating.md
test -f upstream/fieldcraft/tests/corpus/filleted-bracket-blocks/recipe.json
