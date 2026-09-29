#!/usr/bin/env bash
set -euo pipefail
if [[ $# -ne 3 ]]; then
    echo 'Usage: bash run_mscart.sh /absolute/path/to/MSCART ground_low|ground_high|satellite /absolute/new/output'
    exit 2
fi
task_executable=$(realpath -- "$1")
task_root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
task_case=$2
case "$task_case" in
    ground_low|ground_high)
        task_level=${task_case#ground_}
        task_nml="$task_root/ground/${task_level}_input.nml"
        task_optics="$task_root/ground/${task_level}_input.nc"
        task_field=optics.nc;task_batches=40;task_photons=50000;;
    satellite)
        task_nml="$task_root/satellite/input.nml"
        task_optics="$task_root/satellite/input.nc"
        task_field=Lidar_Hu2001_OP.nc;task_batches=100;task_photons=100000;;
    *) echo 'Unknown case.' >&2;exit 2;;
esac
[[ -x "$task_executable" && ! -e "$3" ]] || { echo 'Use a valid executable and a new output directory.' >&2;exit 2; }
mkdir -p -- "$3"
task_output=$(realpath -- "$3")
cp -- "$task_nml" "$task_output/input.nml"
cp -- "$task_optics" "$task_output/$task_field"
cd -- "$task_output"
"$task_executable" "$task_batches" "$task_photons" 1 input.nml result.nc >run.log 2>&1
echo "Finished: $task_output/result.nc"
