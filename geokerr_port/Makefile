# GEOKERR CUDA PORT MAKEFILE
# Builds FORTRAN reference, CUDA implementations, and validation tools

.PHONY: all clean fortran cuda test validate help ref-lib

# Compiler settings
FC = gfortran
# Support long fixed-form lines in legacy .f sources and default real kinds
FFLAGS = -O3 -fdefault-real-8 -fdefault-double-8 -ffixed-line-length-132

NVCC = nvcc
NVFLAGS = -O3 -arch=sm_60 -std=c++11
NVLIBS = -lcudart

# Directories
BUILD_DIR = build
CUDA_DIR = cuda
REF_DIR = reference
SCRIPTS_DIR = scripts
VALIDATION_DIR = validation
REF_BUILD_DIR = $(REF_DIR)/build
REFLIB = $(REF_BUILD_DIR)/libgeokerr_ref.so

# Default target
all: fortran cuda

help:
	@echo "GEOKERR Build System"
	@echo "==================="
	@echo ""
	@echo "Targets:"
	@echo "  all        - Build FORTRAN and CUDA implementations"
	@echo "  fortran    - Build FORTRAN reference implementations"
	@echo "  ref-lib    - Build FORTRAN shared library (libgeokerr_ref.so)"
	@echo "  cuda       - Build CUDA implementations"
	@echo "  validate   - Run validation tests"
	@echo "  test       - Run quick functionality tests"
	@echo "  clean      - Clean build artifacts"
	@echo ""
	@echo "Files built in $(BUILD_DIR)/ and $(REF_BUILD_DIR)/"

# Create build directories
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(REF_BUILD_DIR):
	mkdir -p $(REF_BUILD_DIR)

# FORTRAN reference implementations
fortran: $(BUILD_DIR) $(BUILD_DIR)/elliptic_batch $(BUILD_DIR)/geodesic_batch

$(BUILD_DIR)/elliptic_batch: $(REF_DIR)/geokerr_batch.f
	cd $(REF_DIR) && $(FC) $(FFLAGS) -o "../$(BUILD_DIR)/elliptic_batch" geokerr_batch.f

$(BUILD_DIR)/geodesic_batch: $(REF_DIR)/geodesic_batch.f  
	cd $(REF_DIR) && $(FC) $(FFLAGS) -o "../$(BUILD_DIR)/geodesic_batch" geodesic_batch.f

# FORTRAN shared library (ISO_C_BINDING API)
ref-lib: $(REFLIB)


# Build shared lib with analytic GEOKERR and helpers
$(REFLIB): $(REF_BUILD_DIR) $(REF_DIR)/geokerr_eval_lambda.f90 $(REF_DIR)/geokerr_original.f
	$(FC) -shared -fPIC $(FFLAGS) -o $@ \
		$(REF_DIR)/geokerr_original.f \
		$(REF_DIR)/geokerr_eval_lambda.f90
	@echo "Built $@"

# CUDA implementations
cuda: $(BUILD_DIR) $(BUILD_DIR)/carlson_validator $(BUILD_DIR)/geokerr_test $(BUILD_DIR)/geokerr_validation $(BUILD_DIR)/geokerr_improved $(BUILD_DIR)/geokerr_final_validation $(BUILD_DIR)/blackhole_raytracer $(BUILD_DIR)/geokerr_compare

$(BUILD_DIR)/carlson_validator: $(CUDA_DIR)/carlson_combined.cu
	$(NVCC) $(NVFLAGS) -o $@ $< $(NVLIBS)

$(BUILD_DIR)/geokerr_test: $(CUDA_DIR)/geokerr_cuda.cu
	$(NVCC) $(NVFLAGS) -DTEST_MAIN -o $@ $< $(NVLIBS)

$(BUILD_DIR)/geokerr_validation: $(CUDA_DIR)/geokerr_validation.cu $(CUDA_DIR)/geokerr_cuda.cu
	$(NVCC) $(NVFLAGS) -o $@ $^ $(NVLIBS)

$(BUILD_DIR)/geokerr_improved: $(CUDA_DIR)/geokerr_improved.cu
	$(NVCC) $(NVFLAGS) -DTEST_IMPROVED -o $@ $< $(NVLIBS)

$(BUILD_DIR)/geokerr_final_validation: $(CUDA_DIR)/geokerr_final_validation.cu $(CUDA_DIR)/geokerr_improved.cu
	$(NVCC) $(NVFLAGS) -o $@ $^ $(NVLIBS)

$(BUILD_DIR)/blackhole_raytracer: $(CUDA_DIR)/blackhole_raytracer.cu $(CUDA_DIR)/geokerr_improved.cu
	$(NVCC) $(NVFLAGS) -o $@ $^ $(NVLIBS)

$(BUILD_DIR)/geokerr_compare: $(CUDA_DIR)/geokerr_compare.cu $(CUDA_DIR)/geokerr_complete.cu $(CUDA_DIR)/carlson_elliptic.cu
	$(NVCC) $(NVFLAGS) -dc $(CUDA_DIR)/carlson_elliptic.cu -o $(BUILD_DIR)/carlson_elliptic.o
	$(NVCC) $(NVFLAGS) -dc $(CUDA_DIR)/geokerr_complete.cu -o $(BUILD_DIR)/geokerr_complete.o
	$(NVCC) $(NVFLAGS) $(CUDA_DIR)/geokerr_compare.cu $(BUILD_DIR)/carlson_elliptic.o $(BUILD_DIR)/geokerr_complete.o -o $@ $(NVLIBS)

# Validation and testing
validate: cuda fortran
	@echo "Running CUDA Carlson validation..."
	cd $(BUILD_DIR) && ./carlson_validator

test: $(BUILD_DIR)/carlson_validator
	@echo "Running basic CUDA functionality test..."
	cd $(BUILD_DIR) && echo "Test functionality placeholder"

# Reference data generation (requires Python environment)
reference-data: fortran
	cd $(VALIDATION_DIR) && ./scripts/run_validation.sh

# Clean build artifacts
clean:
	rm -rf $(BUILD_DIR)
	rm -rf $(REF_BUILD_DIR)
	find . -name "*.o" -delete
	find . -name "*.mod" -delete

# Performance benchmark (future)
benchmark: cuda
	@echo "Performance benchmark placeholder"

# Documentation generation (future)
docs:
	@echo "Documentation generation placeholder"