#!/bin/bash
set -euo pipefail

topology_path="../../../mc-topology/testcase-adsorption/ZIF8-CH4O-H2O/"
build_path="../../../build/maniac"

data=$topology_path"topology.data"
inc=$topology_path"parameters.inc"
outputs="../../../outputs/"

# Loop over all good-input files
for input in good-input-*.maniac; do

    # Clean outputs before running
    rm -rf "$outputs"
    mkdir -p "$outputs"

    # Run silently (everything goes into log.maniac)
    $build_path -i "$input" -d "$data" -p "$inc" -o "$outputs" > log.maniac 2>&1

    # === Test: program termination
    if grep -q "Simulation Completed" log.maniac; then
        echo "✅ [PASS] $input : Simulation terminated normally"
    else
        echo "❌ [FAIL] $input : No termination message found"
        exit 1
    fi
done

# Run bad inputs (expect failure)
for input in bad-input-*.maniac; do
    rm -rf "$outputs"
    mkdir -p "$outputs"

    $build_path -i "$input" -d "$data" -p "$inc" -o "$outputs" > log.maniac 2>&1 || true  # don't exit on error

    if grep -E -q "(Error|STOP)" log.maniac; then
        echo "✅ [PASS] $input : Correctly failed with error"
    else
        echo "❌ [FAIL] $input : Unexpectedly succeeded (no error in log)"
        exit 1
    fi
done
