#!/bin/bash
set -euo pipefail

topology_path="../../../mc-topology/molecule-reservoir/CH4O-H2O/"
build_path="../../../build/maniac"

inc=$topology_path"parameters.inc"
input="input.maniac"
outputs="../../../outputs/"

# Loop over all good-input files
for data in good-*.data; do
    # Clean outputs before running
    rm -rf "$outputs"
    mkdir -p "$outputs"

    # Run silently (everything goes into log.maniac)
    $build_path -i "$input" -d "$data" -p "$inc" -o "$outputs" > log.maniac 2>&1

    # === Test: program termination
    if grep -q "Simulation Completed" log.maniac; then
        echo "✅ [PASS] $data : Simulation terminated normally"
    else
        echo "❌ [FAIL] $data : No termination message found"
        exit 1
    fi
done

# Run bad inputs (expect failure)
for data in bad-*.data; do
    rm -rf "$outputs"
    mkdir -p "$outputs"

    $build_path -i "$input" -d "$data" -p "$inc" -o "$outputs" > log.maniac 2>&1 || true  # don't exit on error

    if grep -E -q "(Error|STOP)" log.maniac; then
        echo "✅ [PASS] $data : Correctly failed with error"
    else
        echo "❌ [FAIL] $data : Unexpectedly succeeded (no error in log)"
        exit 1
    fi
done
