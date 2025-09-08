# GEOKERR VALIDATION FRAMEWORK

Comprehensive cryptographic validation system for comparing GEOKERR FORTRAN reference against CUDA port implementations.

## Features

- **Cryptographic Validation**: BLAKE3 hashing with Merkle tree proofs
- **ULP-Based Quantization**: Deterministic floating-point comparison
- **Crash-Safe Resume**: Automatic recovery from interruptions  
- **Fail Window Analysis**: ±32 sample context around failures
- **Minimal Storage**: ~1-2KB per passing case, cryptographic audit trails
- **Production Ready**: Handles millions of geodesic cases

## Quick Start

```bash
# Build the framework
make

# Initialize validation run
./build/gkval init --run-id TEST --size 1024 --atol 1e-12

# Run validation (CUDA vs FORTRAN)  
./build/gkval validate TEST ../build/geodesic_batch

# Verify integrity
./build/gkval audit TEST

# Analyze failures
./build/gkval diff TEST 12345 phi
```

## Architecture

### Core Components
- `gkval_core.{c,h}` - ULP quantization, BLAKE3 hashing, Merkle trees
- `gkval_validator.{c,h}` - Streaming validation with crash-safe resume
- `gkval_adapters.{c,h}` - Interfaces to CUDA and FORTRAN implementations
- `gkval_io.c` - File I/O for manifests, ledgers, fail windows
- `gkval_cli.c` - Command-line interface

### Data Model
```
runs/<run_id>/
  manifest.json     # Configuration and environment
  ledger.ndjson     # Per-case results with Merkle roots
  cursor.json       # Crash-safe resume state
  fail_windows/     # Binary diff data for mismatches
```

## Validation Process

1. **Initialize**: Set tolerances, grid parameters, quantization levels
2. **Stream**: Process geodesic segments with real-time comparison
3. **Hash**: Quantize and hash reference data for cryptographic proofs
4. **Merkle**: Build integrity trees from segments to case level
5. **Resume**: Automatic crash recovery from any point
6. **Audit**: Verify all stored hashes against recomputed values

## Configuration

Default tolerances:
- `atol = 1e-12` (absolute tolerance)
- `rtol = 1e-12` (relative tolerance)  
- `max_ulps = 2` (maximum ULP distance)
- `quant_lsb = 8` (quantization bits)

## Integration

The framework connects to existing implementations:
- **CUDA**: Calls `compute_geokerr_complete_cuda()` 
- **FORTRAN**: Executes `geodesic_batch` with proper file I/O

## Storage Efficiency

- **Passing cases**: ~1KB each (Merkle roots only)
- **Failing cases**: ~32KB each (with context window)
- **Total for 1M cases**: ~10GB maximum storage

## Build Requirements

- GCC with C99 support
- NVIDIA CUDA toolkit (nvcc)
- Existing GEOKERR CUDA implementation
- FORTRAN geodesic batch processor

## Files

### Framework Source
- Core validation library (66KB)
- Build system and executables (280KB)  
- Documentation and examples

### Test Data  
- Reference datasets in `test_data/` (2.8MB)
- FORTRAN validation inputs/outputs
- Test case generators

**Total Framework Size: 3.1MB**