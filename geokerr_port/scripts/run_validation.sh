#!/bin/bash

# GEOKERR Validation Pipeline
# Generates reference data and validates implementations

set -e

echo "GEOKERR VALIDATION PIPELINE"
echo "==========================="

cd "$(dirname "$0")/.."  # Go to geokerr_port root

# Ensure build exists
if [ ! -d "build" ] || [ ! -f "build/elliptic_batch" ]; then
    echo "Building FORTRAN reference implementations..."
    ./scripts/build.sh
fi

# Generate test cases if needed
if [ ! -f "validation/test_data/elliptic_inputs.dat" ]; then
    echo "Generating test cases..."
    cd validation
    source venv/bin/activate
    python3 generate_test_cases.py
    cd ..
fi

# Run elliptic integral validation
echo "Generating elliptic integral reference data..."
cd validation
cp test_data/elliptic_inputs.dat elliptic_inputs.dat
../build/elliptic_batch
mv elliptic_outputs.dat test_data/elliptic_reference.dat

# Run geodesic validation  
echo "Generating geodesic reference data..."
cp test_data/geodesic_inputs.dat geodesic_inputs.dat
../build/geodesic_batch
mv geodesic_outputs.dat test_data/geodesic_reference.dat

# Clean up temporary files
rm elliptic_inputs.dat geodesic_inputs.dat

echo ""
echo "Validation complete!"
echo "Reference data generated:"
echo "  - validation/test_data/elliptic_reference.dat (5000 cases)"
echo "  - validation/test_data/geodesic_reference.dat (24 cases)"
echo ""
echo "Next steps:"
echo "1. Implement CUDA Carlson elliptic integrals in cuda/carlson_elliptic.cu"
echo "2. Implement CUDA geodesic solver in cuda/geokerr_cuda.cu"
echo "3. Create validation suite in cuda/geokerr_validation.cu"
echo "4. Run performance benchmarks"