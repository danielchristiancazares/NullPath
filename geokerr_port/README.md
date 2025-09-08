# GEOKERR CUDA PORT PROJECT

Semi-analytic Kerr geodesic computation ported from FORTRAN to CUDA for high-precision black hole ray tracing.

## Project Structure

```
geokerr_port/
├── README.md                    # This file
├── reference/                   # Original FORTRAN geokerr code
│   ├── geokerr.f               # Original Dexter & Agol (2009) code
│   ├── geokerr_batch.f         # Modified for batch testing
│   └── carlson_test.f          # Standalone elliptic integral tester
├── cuda/                       # CUDA implementation
│   ├── carlson_elliptic.cu     # Carlson elliptic integrals
│   ├── geokerr_cuda.cu         # Semi-analytic geodesic solver
│   └── geokerr_validation.cu   # End-to-end validation suite
├── validation/                 # Reference datasets and validation
│   ├── generate_test_cases.py  # Parameter sweep generator
│   ├── reference_data/         # Golden reference outputs
│   └── validation_reports/     # Test results and comparisons
├── tests/                      # Unit tests and benchmarks
│   ├── test_elliptic.cu        # Carlson function unit tests
│   ├── test_geodesics.cu       # Geodesic solver tests
│   └── benchmark.cu            # Performance comparisons
└── scripts/                    # Build and automation scripts
    ├── build.sh               # Compilation scripts
    ├── run_validation.sh      # Validation pipeline
    └── generate_reference.sh  # Reference data generation
```

## Implementation Timeline

### Phase 1: Reference Implementation (Weeks 1-2)
- [x] Set up directory structure
- [ ] Download and compile original geokerr
- [ ] Create FORTRAN batch testing harness
- [ ] Generate comprehensive parameter sweep
- [ ] Create golden reference datasets

### Phase 2: Elliptic Integrals (Week 3)
- [ ] Implement Carlson RF, RD, RJ, RC functions in CUDA
- [ ] Unit test against FORTRAN reference
- [ ] Performance benchmarking

### Phase 3: Semi-Analytic Geodesics (Weeks 4-6)
- [ ] Port core geodesic solver (GEOMU, GEOR, GEOPHITIME)
- [ ] Handle edge cases and orbit classification
- [ ] End-to-end validation against reference

### Phase 4: Integration (Week 7)
- [ ] Integrate with existing black hole ray tracer
- [ ] Replace numerical integration with semi-analytic
- [ ] Validate conservation monitoring

### Phase 5: Optimization (Week 8)
- [ ] GPU performance optimization
- [ ] Memory access optimization
- [ ] Final benchmarking and documentation

## References

- Dexter & Agol (2009): "A Fast New Public Code for Computing Photon Orbits in a Kerr Spacetime"
- Carlson (1995): "Numerical computation of real or complex elliptic integrals"
- Original geokerr: https://faculty.washington.edu/agol/geokerr/

## Target Accuracy

- Elliptic integrals: < 1e-12 relative error
- Geodesic coordinates: < 1e-10 relative error  
- Performance: >10x speedup over CPU geokerr
- Conservation: Machine precision Hamiltonian conservation