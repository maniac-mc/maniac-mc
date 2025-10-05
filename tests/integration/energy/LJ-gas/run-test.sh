#!/bin/bash
set -e  # Exit on error

topology_path="../../../../mc-topology/testcase-energy/LJ-gas/"
build_path="../../../../build/maniac"

input=$topology_path"input.maniac"
data=$topology_path"topology.data"
inc=$topology_path"parameters.inc"
outputs="outputs/" # optional

$build_path -i $input -d $data -p $inc -o $outputs >> log.maniac

maniac_log=$outputs"log.maniac"
lammps_log=$topology_path"log.lammps"

# Extract reference TotEng from LAMMPS log (last numeric value in that column)
ref_energy=$(grep "TotEng" "$lammps_log" -A 1 | tail -n 1 | awk '{print $2}')

tolerance=0.01 # allowed deviation

# Get the last TotEng value from the log
last_toteng=$(awk '/TotEng/ {getline; gsub(/^[| ]+/,""); print $2}' "$maniac_log")

if [ -z "$last_toteng" ]; then
    echo "❌ Test failed: 'TotEng' not found in $maniac_log"
    exit 1
fi

# Compute absolute difference
diff=$(awk -v a="$last_toteng" -v b="$ref_energy" 'BEGIN { print (a > b ? a - b : b - a) }')

# Check tolerance
if (( $(echo "$diff > $tolerance" | bc -l) )); then
    echo "❌ Test failed: TotEng=$last_toteng differs from ref=$ref_energy by $diff"
    exit 1
else
    echo "✅ Test passed: TotEng=$last_toteng (diff=$diff within tolerance $tolerance)"
fi
