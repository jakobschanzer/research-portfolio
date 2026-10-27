#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 2 || $# -gt 3 ]]; then
  echo "Usage: $0 TRAJECTORY_ROOT CONFIG_DIR [OUTPUT_DIR]" >&2
  exit 2
fi

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
project_root=$(cd "$script_dir/.." && pwd)
trajectory_root=$(cd "$1" && pwd)
config_dir=$(cd "$2" && pwd)
output_dir=${3:-"$project_root/data/generated/hbond"}
hbond_executable="$project_root/build/hbb"
temperatures=(5 15 25 35 45 55 65 75 85 95)
max_jobs=${MAX_JOBS:-4}

command -v catdcd >/dev/null || {
  echo "ERROR: catdcd is required but was not found on PATH." >&2
  exit 1
}
[[ -x "$hbond_executable" ]] || {
  echo "ERROR: $hbond_executable is missing. Run 'make build' first." >&2
  exit 1
}
for required in drive top; do
  [[ -f "$config_dir/$required" ]] || {
    echo "ERROR: missing $config_dir/$required" >&2
    exit 1
  }
done

mkdir -p "$output_dir"

wait_for_slot() {
  while (( $(jobs -r | wc -l) >= max_jobs )); do
    sleep 1
  done
}

for temperature in "${temperatures[@]}"; do
  trajectory="$trajectory_root/$temperature/prd.xtc"
  [[ -f "$trajectory" ]] || {
    echo "ERROR: missing trajectory $trajectory" >&2
    exit 1
  }
  wait_for_slot
  (
    work_dir=$(mktemp -d "${TMPDIR:-/tmp}/hbond-${temperature}.XXXXXX")
    trap 'rm -rf "$work_dir"' EXIT
    cp "$config_dir/drive" "$config_dir/top" "$work_dir/"
    catdcd -o "$work_dir/prd.crdbox" -otype crdbox -xtc "$trajectory"
    (
      cd "$work_dir"
      "$hbond_executable" < prd.crdbox > hbb.log 2>&1
    )
    cp "$work_dir/dist.out" "$output_dir/dist.$temperature"
    cp "$work_dir/hbb.log" "$output_dir/hbb.$temperature.log"
    echo "Completed ${temperature} °C"
  ) &
done

wait
python3 "$project_root/src/python/compile_hbond_outputs.py" \
  --input-dir "$output_dir" \
  --output-dir "$output_dir"
echo "Hydrogen-bond workflow complete: $output_dir"
