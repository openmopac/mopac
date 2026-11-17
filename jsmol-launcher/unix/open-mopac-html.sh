#!/usr/bin/env bash
# open-mopac-html.sh - open a MOPAC-generated HTML file with JSmol (macOS / Linux).
#
# Same logic as the Windows launcher:
#   1. walk up from the HTML file's folder until a folder containing
#      jsmol/JSmol.min.js is found  (the "workspace root")
#   2. start  python3 -m http.server  rooted there, unless one is already running
#   3. open  http://localhost:<port>/<relative path>  in the default browser
#
# Usage:  open-mopac-html.sh [--no-open] [--port N] file.html
#         open-mopac-html.sh --stop
set -euo pipefail

PORT=8000; OPEN=1
STATE_DIR="${XDG_STATE_HOME:-$HOME/.local/state}/mopac-jsmol"; mkdir -p "$STATE_DIR"
STATE="$STATE_DIR/server"

while [ $# -gt 0 ]; do
  case "$1" in
    --no-open) OPEN=0; shift ;;
    --port) PORT="$2"; shift 2 ;;
    --stop)
      if [ -f "$STATE" ]; then read -r pid root port < "$STATE"; kill "$pid" 2>/dev/null && echo "Stopped server for $root"; rm -f "$STATE"
      else echo "No server is running."; fi; exit 0 ;;
    *) FILE="$1"; shift ;;
  esac
done
[ -n "${FILE:-}" ] || { echo "usage: $0 [--no-open] [--port N] file.html | --stop" >&2; exit 1; }
[ -f "$FILE" ] || { echo "File not found: $FILE" >&2; exit 1; }

# 1. absolute path and workspace root
FILE="$(cd "$(dirname "$FILE")" && pwd)/$(basename "$FILE")"
dir="$(dirname "$FILE")"; root=""
probe="$dir"
while :; do
  if [ -f "$probe/jsmol/JSmol.min.js" ]; then root="$probe"; break; fi
  [ "$probe" = "/" ] && break
  probe="$(dirname "$probe")"
done
if [ -z "$root" ] && [ -n "${MOPAC_JSMOL_ROOT:-}" ] && [ -f "$MOPAC_JSMOL_ROOT/jsmol/JSmol.min.js" ]; then root="$MOPAC_JSMOL_ROOT"; fi
if [ -z "$root" ]; then
  cat >&2 <<EOF
No 'jsmol' folder was found in $dir or any folder above it.
Put one copy of the JSmol 'jsmol' folder (the one containing JSmol.min.js) in a
folder above your MOPAC working folders, e.g. ~/MOPAC/jsmol, or set
MOPAC_JSMOL_ROOT to that parent folder.
EOF
  exit 1
fi

# 2. server
ok=0
if [ -f "$STATE" ]; then
  read -r pid sroot sport < "$STATE"
  if kill -0 "$pid" 2>/dev/null; then
    if [ "$sroot" = "$root" ]; then ok=1; PORT="$sport"; else kill "$pid" 2>/dev/null; fi
  fi
fi
if [ "$ok" = 0 ]; then
  while lsof -iTCP:"$PORT" -sTCP:LISTEN >/dev/null 2>&1; do PORT=$((PORT+1)); done
  (cd "$root" && nohup python3 -m http.server "$PORT" --bind 127.0.0.1 >"$STATE_DIR/server.log" 2>&1 & echo $! "$root" "$PORT" > "$STATE")
  for _ in $(seq 1 50); do
    curl -sI "http://localhost:$PORT/jsmol/JSmol.min.js" | grep -q '200' && { ok=1; break; }
    sleep 0.2
  done
  [ "$ok" = 1 ] || { echo "The local web server did not start (port $PORT)." >&2; exit 1; }
fi

# 3. URL
rel="${FILE#"$root"/}"
enc="$(printf '%s' "$rel" | python3 -c 'import sys,urllib.parse; print("/".join(urllib.parse.quote(p) for p in sys.stdin.read().split("/")))')"
url="http://localhost:$PORT/$enc"
echo "JSmol root : $root"
echo "Opening    : $url"
if [ "$OPEN" = 1 ]; then
  if command -v open >/dev/null; then open "$url"; else xdg-open "$url"; fi
fi
