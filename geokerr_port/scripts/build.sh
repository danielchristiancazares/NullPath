#!/bin/bash

# GEOKERR Build Script
# Compiles FORTRAN reference implementations and batch testers

set -e  # Exit on any error

echo "GEOKERR REFERENCE BUILD SYSTEM"
echo "==============================="

# Check for FORTRAN compiler
if command -v gfortran >/dev/null 2>&1; then
    FC=gfortran
    FFLAGS="-O3 -fdefault-real-8 -fdefault-double-8"
elif command -v f77 >/dev/null 2>&1; then
    FC=f77
    FFLAGS="-O3"
elif command -v ifort >/dev/null 2>&1; then
    FC=ifort
    FFLAGS="-O3 -r8"
else
    echo "Error: No FORTRAN compiler found (gfortran, f77, or ifort required)"
    exit 1
fi

echo "Using compiler: $FC"
echo "Compiler flags: $FFLAGS"

# Create build directory
cd "$(dirname "$0")/.."  # Go to geokerr_port root
BUILD_DIR="build"
mkdir -p "$BUILD_DIR"

echo "Building in: $(pwd)/$BUILD_DIR"

# Build elliptic integral batch tester
echo "Building elliptic integral batch tester..."
cd reference
$FC $FFLAGS -o "../$BUILD_DIR/elliptic_batch" geokerr_batch.f
echo "✓ elliptic_batch executable created"

# Build geodesic batch tester (placeholder - needs full geokerr integration)
echo "Building geodesic batch tester..."
$FC $FFLAGS -o "../$BUILD_DIR/geodesic_batch" geodesic_batch.f
echo "✓ geodesic_batch executable created"

# Make executables runnable
chmod +x "../$BUILD_DIR/elliptic_batch"
chmod +x "../$BUILD_DIR/geodesic_batch"

echo ""
echo "Build complete!"
echo "Executables created in $BUILD_DIR/:"
ls -la "../$BUILD_DIR/"

echo ""
echo "Next steps:"
echo "1. cd validation && ../build/elliptic_batch"
echo "2. cd validation && ../build/geodesic_batch" 
echo "3. Compare outputs with CUDA implementation"