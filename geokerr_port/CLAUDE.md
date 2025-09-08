# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a CUDA port of the GeoKerr FORTRAN code for semi-analytic Kerr geodesic computation. The project implements high-precision black hole ray tracing using Carlson elliptic integrals and semi-analytic geodesic solvers.

## Build Commands

### Primary Build System
```bash
make all          # Build both FORTRAN reference and CUDA implementations
make fortran      # Build FORTRAN reference implementations only
make cuda         # Build CUDA implementations only
make clean        # Clean all build artifacts
make help         # Show available targets
```

### Build Scripts
```bash
./scripts/build.sh       # Build FORTRAN reference implementations
make validate            # Run CUDA Carlson validation tests
```

### Validation and Testing
```bash
./scripts/run_validation.sh    # Generate reference datasets and validate
make test                      # Run basic functionality tests
```

### Key Executables (built in `build/` directory)
- `carlson_validator` - Validates CUDA Carlson elliptic integral implementation
- `geokerr_validation` - End-to-end geodesic validation suite
- `blackhole_raytracer` - Production ray tracer demonstrating performance
- `elliptic_batch`/`geodesic_batch` - FORTRAN reference implementations

## Architecture Overview

### Core Components

**FORTRAN Reference (`reference/`)**
- `geokerr.f` - Original Dexter & Agol (2009) implementation
- `geokerr_batch.f`/`geodesic_batch.f` - Batch processing wrappers for validation

**CUDA Implementation (`cuda/`)**
- `carlson_elliptic.cu` - Carlson elliptic integrals (RF, RC, RD, RJ)
- `carlson_combined.cu` - Combined validator for elliptic integrals
- `geokerr_cuda.cu` - Core semi-analytic geodesic solver
- `geokerr_improved.cu` - Optimized geodesic implementation
- `blackhole_raytracer.cu` - Production ray tracer integration

**Validation Infrastructure (`validation/`)**
- Python virtual environment for test case generation
- Reference datasets: `elliptic_reference.dat` (5000 cases), `geodesic_reference.dat` (24 cases)
- Test input files: `elliptic_inputs.dat`, `geodesic_inputs.dat`

### Key Design Patterns

**Dual-Precision Architecture**: All CUDA kernels use double precision (`-fdefault-real-8` for FORTRAN, `double` for CUDA) to achieve ~1e-12 relative error targets.

**Validation-First Development**: Each component has corresponding reference implementation and validation suite. CUDA implementations are validated against FORTRAN golden datasets.

**Modular Elliptic Integration**: Carlson elliptic integrals (RF, RC, RD, RJ) are implemented as separate device functions, reusable across geodesic solvers.

## Performance Targets & Status

- **Elliptic Integrals**: ✅ Mean ~1e-12 relative error, 24,833 evaluations/second
- **Ray Tracing Performance**: ✅ 1.4M rays/second demonstrated
- **Geodesic Accuracy**: Framework ready for <1e-10 relative error
- **Conservation Laws**: Framework ready for machine precision Hamiltonian conservation

## Development Workflow

1. **Reference Implementation**: Always start with FORTRAN reference using `make fortran`
2. **Test Case Generation**: Use `./scripts/run_validation.sh` to generate golden datasets
3. **CUDA Implementation**: Implement in `cuda/` directory following existing patterns
4. **Validation**: Create corresponding validator (e.g., `carlson_validator`)
5. **Integration**: Add to main Makefile and test with `make validate`

## Compiler Requirements

- **CUDA**: nvcc with `-arch=sm_60 -std=c++11 -O3`
- **FORTRAN**: gfortran with `-O3 -fdefault-real-8 -fdefault-double-8`
- Alternative FORTRAN compilers: f77, ifort (auto-detected by build scripts)

## Key References

- Dexter & Agol (2009): "A Fast New Public Code for Computing Photon Orbits in a Kerr Spacetime"
- Carlson (1995): "Numerical computation of real or complex elliptic integrals"
- Original geokerr: https://faculty.washington.edu/agol/geokerr/