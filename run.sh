#!/bin/bash
set -e  # Exit on error

case="H2O-gas"

base_energy="mc-topology/testcase-energy"
base_adsorption="mc-topology/testcase-adsorption"
base_reservoir="mc-topology/molecule-reservoir"

case "$case" in
  "SLIT")
    folder="$base_adsorption/SLIT"
    input="$folder/input.maniac"
    data="$folder/topology.data"
    inc="$folder/parameters.inc"
    ;;
  "CH4O-H2O")
    folder="$base_reservoir/CH4O-H2O"
    input="$folder/input.maniac"
    data="$folder/topology.data"
    inc="$folder/parameters.inc"
    ;;
  "ZIF8-CH4O-H2O")
    folder="$base_adsorption/ZIF8-CH4O-H2O"
    input="$folder/input.maniac"
    data="$folder/topology.data"
    inc="$folder/parameters.inc"
    res="$base_reservoir/CH4O-H2O/topology.data"
    ;;
  "ZIF8-H2O-1")
    folder="$base_energy/ZIF8-H2O"
    input="$folder/zif8-water.maniac"
    data="$folder/zif8-water.data"
    inc="$folder/zif8-water.inc"
    res="$base_reservoir/tip4p-water/water.data"
    ;;
  "ZIF8-H2O-2")
    folder="$base_adsorption/ZIF8-H2O"
    input="$folder/zif8-water.maniac"
    data="$folder/zif8-water.data"
    inc="$folder/zif8-water.inc"
    res="$base_reservoir/tip4p-water/water.data"
    ;;
  "MFI-CO2")
    folder="$base_adsorption/MFI-CO2"
    input="$folder/MFI-CO2.maniac"
    data="$folder/MFI-CO2.data"
    inc="$folder/MFI-CO2.inc"
    res="$base_reservoir/co2/co2.data"
    ;;
  "LJ-gas")
    folder="$base_energy/LJ-gas"
    input="$folder/input.maniac"
    data="$folder/topology.data"
    inc="$folder/parameters.inc"
    res=""
    ;;
  "H2O-gas")
    folder="$base_adsorption/H2O"
    input="$folder/input.maniac"
    data="$folder/topology.data"
    inc="$folder/parameters.inc"
    res=""
    ;;
  "DIPOLE-orthorhombic")
    folder="$base_energy/DIPOLE-orthorhombic"
    input="$folder/dipole.maniac"
    data="$folder/dipole.data"
    inc="$folder/dipole.inc"
    res=""
    ;;
  *)
    echo "Unknown case: $case"
    exit 1
    ;;
esac

outputs="outputs/"

echo "Running maniac with system=$case ..."
./build/maniac -i "$input" -d "$data" -p "$inc" -o "$outputs"  ${res:+-r "$res"}
