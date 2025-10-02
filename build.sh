#!/bin/bash
set -e  # Exit immediately if any command fails

#!/bin/bash
set -e

# Read version and git commit ID
VERSION=$(cat version.txt | xargs)
GIT_COMMIT=$(git rev-parse --short HEAD)

# Generate Fortran module from template
sed -e "s/@VERSION@/$VERSION/" \
    -e "s/@GIT_COMMIT@/$GIT_COMMIT/" \
    src/version_module.f90.in > src/version_module.f90

echo "Cleaning with make..."
make distclean # Use distclean to remove dependencies.d and include/ as well

echo "Generating dependencies..."
make depend # Ensure dependencies.d is regenerated

echo "Building project..."
make