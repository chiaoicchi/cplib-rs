#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -ne 1 ]; then
  echo "usage: bd <problem>" >&2
  exit 1
fi

p="$1"

if ! root="$(git rev-parse --show-toplevel 2>/dev/null)"; then
  echo "error: run inside the git repository" >&2
  exit 1
fi
if [ "$(dirname "$PWD")" != "$root" ] || [ "$(basename "$PWD")" != "verify" ]; then
  echo "error: run inside verify/" >&2
  exit 1
fi

if [ ! -f "src/bin/$p.rs" ]; then
  echo "error: src/bin/$p.rs not found" >&2
  exit 1
fi

lib="$root"
edition="$(sed -n 's/^edition *= *"\([0-9]*\)".*/\1/p' "$lib/Cargo.toml" | head -n 1)"

out="$(bundle-rs --lib "$lib" "src/bin/$p.rs" | rustfmt --edition "$edition")"
printf '%s\n' "$out" | wl-copy
bytes="$(printf '%s\n' "$out" | wc -c)"
echo "copied: src/bin/$p.rs -> clipboard ($bytes bytes)"
