#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -ne 1 ]; then
  echo "usage: ck <problem>" >&2
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

cache="${XDG_CACHE_HOME:-$HOME/.cache}/cplib-verify"
lcp="$cache/library-checker-problems"
mkdir -p "$cache"
if [ ! -d "$lcp" ]; then
  git clone --depth 1 https://github.com/yosupo06/library-checker-problems.git "$lcp"
else
  git -C "$lcp" pull --ff-only || echo "warning: cannot update library-checker-problems" >&2
fi

echo "generate test cases"
(cd "$lcp" && ./generate.py -p "$p")

dir="$(find "$lcp" -mindepth 2 -maxdepth 2 -type d -name "$p" | head -n 1)"
if [ -z "$dir" ]; then
  echo "error: problem $p not found in library-checker-problems" >&2
  exit 1
fi

edition="$(sed -n 's/^edition *= *"\([0-9]*\)".*/\1/p' Cargo.toml | head -n 1)"
bundled="$(mktemp --suffix=.rs)"
bin="$(mktemp)"
tmp="$(mktemp)"
trap 'rm -f "$bundled" "$bin" "$tmp"' EXIT

echo "bundle"
bundle-rs --lib "$root" "src/bin/$p.rs" | rustfmt --edition "$edition" >"$bundled"
echo "build"
rustc --edition "$edition" -O --crate-name "$p" -o "$bin" "$bundled"

echo "test"
ulimit -s unlimited
passed=0
failed=0
for input in "$dir"/in/*.in; do
  [ -e "$input" ] || continue
  base="$(basename "$input" .in)"
  expected="$dir/out/$base.out"
  if [ ! -f "$expected" ]; then
    echo "warning: $expected not found, skipping" >&2
    continue
  fi
  if ! "$bin" <"$input" >"$tmp" 2>/dev/null; then
    echo "FAILED (runtime error): $base"
    failed=$((failed + 1))
    continue
  fi
  if [ -x "$dir/checker" ]; then
    if "$dir/checker" "$input" "$tmp" "$expected" >/dev/null 2>&1; then
      echo "PASSED: $base"
      passed=$((passed + 1))
    else
      echo "FAILED: $base"
      failed=$((failed + 1))
    fi
  else
    if diff -q "$tmp" "$expected" >/dev/null 2>&1; then
      echo "PASSED: $base"
      passed=$((passed + 1))
    else
      echo "FAILED: $base"
      failed=$((failed + 1))
    fi
  fi
done

echo
if [ "$passed" -eq 0 ] && [ "$failed" -eq 0 ]; then
  echo "no test cases found"
  exit 1
elif [ "$failed" -eq 0 ]; then
  echo "all $passed tests passed"
else
  echo "results: $passed passed, $failed failed"
  exit 1
fi
