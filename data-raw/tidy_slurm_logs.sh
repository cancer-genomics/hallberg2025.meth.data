#!/usr/bin/env bash
# Move all but the most recent slurm-*.out to logs/archive/
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ARCHIVE="$SCRIPT_DIR/logs/archive"
mkdir -p "$ARCHIVE"

# Find all slurm logs sorted oldest-first; skip the last one
mapfile -t logs < <(ls -t "$SCRIPT_DIR"/slurm-*.out 2>/dev/null)

if [ ${#logs[@]} -le 1 ]; then
    echo "Nothing to archive (${#logs[@]} log(s) found)."
    exit 0
fi

to_archive=("${logs[@]:1}")  # all but the first (newest)
for f in "${to_archive[@]}"; do
    echo "Archiving $(basename "$f")"
    mv "$f" "$ARCHIVE/"
done

echo "Done. Kept: $(basename "${logs[0]}")"
