# GEOKERR VALIDATION FRAMEWORK MAKEFILE
# Builds the comprehensive validation system with cryptographic proofs

.PHONY: all clean test install

# Compiler settings
CC = gcc
CFLAGS = -std=c99 -O3 -Wall -Wextra -fPIC
LDFLAGS = -lm

# For CUDA integration (optional)
NVCC = nvcc
# CUDA compilation flags: disable FMAD contraction and FTZ for determinism; keep precise div/sqrt
NVFLAGS = -O3 -arch=sm_60 -std=c++11 --fmad=false --prec-div=true --prec-sqrt=true --ftz=false

# Directories
BUILD_DIR = build
SRC_DIR = .
INCLUDE_DIR = .
REFLIB_DIR = ../reference/build

# Source files
CORE_SRCS = gkval_core.c gkval_io.c
VALIDATOR_SRCS = gkval_validator.c gkval_adapters.c
CLI_SRCS = gkval_cli.c
ALL_SRCS = $(CORE_SRCS) $(VALIDATOR_SRCS) $(CLI_SRCS)
TEST_SRCS = gkval_fail_window_test.c

# Object files
CORE_OBJS = $(CORE_SRCS:%.c=$(BUILD_DIR)/%.o)
VALIDATOR_OBJS = $(VALIDATOR_SRCS:%.c=$(BUILD_DIR)/%.o)
CLI_OBJS = $(CLI_SRCS:%.c=$(BUILD_DIR)/%.o)
ALL_OBJS = $(ALL_SRCS:%.c=$(BUILD_DIR)/%.o)
TEST_OBJS = $(TEST_SRCS:%.c=$(BUILD_DIR)/%.o)

# Targets
GKVAL_CLI = $(BUILD_DIR)/gkval
LIBGKVAL_CORE = $(BUILD_DIR)/libgkval_core.a
LIBGKVAL_VALIDATOR = $(BUILD_DIR)/libgkval_validator.a
FAILWIN_TEST = $(BUILD_DIR)/gkval_fail_window_test

# Default target
all: $(BUILD_DIR) $(GKVAL_CLI) $(LIBGKVAL_CORE) $(LIBGKVAL_VALIDATOR)

# Create build directory
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

# Core library (hash, quantization, merkle)
$(LIBGKVAL_CORE): $(CORE_OBJS)
	ar rcs $@ $^

# Validator library (streaming validation, adapters)
$(LIBGKVAL_VALIDATOR): $(VALIDATOR_OBJS) $(LIBGKVAL_CORE)
	ar rcs $@ $(VALIDATOR_OBJS)

# Main CLI tool (with CUDA linking and rpath to reference lib)
$(GKVAL_CLI): $(CLI_OBJS) $(LIBGKVAL_VALIDATOR) $(LIBGKVAL_CORE) $(BUILD_DIR)/geokerr_complete.o $(BUILD_DIR)/carlson_elliptic.o $(BUILD_DIR)/geokerr_stream.o
	$(NVCC) $(NVFLAGS) -o $@ $(CLI_OBJS) $(LIBGKVAL_VALIDATOR) $(LIBGKVAL_CORE) $(BUILD_DIR)/geokerr_complete.o $(BUILD_DIR)/carlson_elliptic.o $(BUILD_DIR)/geokerr_stream.o $(LDFLAGS) \
		-L$(REFLIB_DIR) -lgeokerr_ref \
		-Xlinker -rpath -Xlinker $(abspath $(REFLIB_DIR))

# CUDA object compilation
$(BUILD_DIR)/geokerr_complete.o: ../cuda/geokerr_complete.cu
	$(NVCC) $(NVFLAGS) -dc $< -o $@

$(BUILD_DIR)/carlson_elliptic.o: ../cuda/carlson_elliptic.cu
	$(NVCC) $(NVFLAGS) -dc $< -o $@

$(BUILD_DIR)/geokerr_stream.o: ../cuda/geokerr_stream.cu
	$(NVCC) $(NVFLAGS) -dc $< -o $@

# Object file compilation
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c
	$(CC) $(CFLAGS) -I$(INCLUDE_DIR) -c $< -o $@

# Installation
install: all
	install -D $(GKVAL_CLI) /usr/local/bin/gkval
	install -D $(LIBGKVAL_CORE) /usr/local/lib/libgkval_core.a
	install -D $(LIBGKVAL_VALIDATOR) /usr/local/lib/libgkval_validator.a
	install -D gkval_core.h /usr/local/include/gkval_core.h
	install -D gkval_validator.h /usr/local/include/gkval_validator.h
	install -D gkval_adapters.h /usr/local/include/gkval_adapters.h

# Testing
test: all $(FAILWIN_TEST)
	@echo "Running validation framework tests..."
	./$(GKVAL_CLI) init --run-id TEST001 --size 256 --delta 1e-2
	@echo "✅ Test initialization complete"
	@echo "Running fail window test..."
	./$(FAILWIN_TEST)
	@echo "✅ Fail window test passed"
	@echo ""
	@echo "To run a full validation test:"
	@echo "  ./$(GKVAL_CLI) validate TEST001 ../build/geodesic_batch"
	@echo ""
	@echo "To audit results:"
	@echo "  ./$(GKVAL_CLI) audit TEST001"

# Demos and examples  
demo: all
	@echo "GEOKERR VALIDATION FRAMEWORK DEMO"
	@echo "================================="
	@echo ""
	@echo "1. Initialize validation run..."
	./$(GKVAL_CLI) init --run-id DEMO --size 128 --delta 1e-2 --atol 1e-10 --rtol 1e-10
	@echo ""
	@echo "2. Show generated manifest..."
	@cat runs/DEMO/manifest.json
	@echo ""
	@echo "3. Run validation (demo mode with placeholder solvers)..."
	./$(GKVAL_CLI) validate DEMO /dev/null || true
	@echo ""
	@echo "✅ Demo complete! Check runs/DEMO/ for generated files"

# Performance testing
benchmark: all
	@echo "GEOKERR VALIDATION PERFORMANCE BENCHMARK"
	@echo "======================================="
	@echo "Testing ULP quantization performance..."
	@time -p bash -c 'for i in {1..10000}; do echo "1.2345678901234567e-12" | ./$(GKVAL_CLI) init --run-id BENCH$$i --size 1; done' 2>/dev/null || true
	@echo "✅ Benchmark complete"

# Development utilities
format:
	clang-format -i *.c *.h

lint:
	cppcheck --enable=all --std=c99 *.c *.h

# Documentation generation
docs:
	doxygen Doxyfile || echo "Install doxygen to generate documentation"

# Clean build artifacts
clean:
	rm -rf $(BUILD_DIR)
	rm -rf runs/
	find . -name "*.o" -delete
	find . -name "*.a" -delete
	find . -name "*~" -delete

# Integration with parent project
integration: all
	@echo "Integrating with GEOKERR CUDA port..."
	cp $(GKVAL_CLI) ../build/
	cp $(LIBGKVAL_CORE) ../build/
	cp $(LIBGKVAL_VALIDATOR) ../build/
	@echo "✅ Integration complete"

# Help
help:
	@echo "GEOKERR Validation Framework Build System"
	@echo "========================================"
	@echo ""
	@echo "Targets:"
	@echo "  all         - Build all components"
	@echo "  test        - Run basic functionality tests"
	@echo "  demo        - Interactive demonstration"
	@echo "  benchmark   - Performance testing"
	@echo "  install     - Install to system directories"
	@echo "  integration - Copy to parent project build/"
	@echo "  clean       - Remove build artifacts"
	@echo "  help        - Show this help"
	@echo ""
	@echo "The main executable is: $(GKVAL_CLI)"
	@echo ""
	@echo "Usage examples:"
	@echo "  make test                    # Quick functionality test"
	@echo "  make demo                    # Interactive demo"
	@echo "  make benchmark              # Performance testing"
	@echo "  make integration            # Install to parent project"

# Dependencies
$(BUILD_DIR)/gkval_core.o: gkval_core.c gkval_core.h
$(BUILD_DIR)/gkval_io.o: gkval_io.c gkval_core.h
$(BUILD_DIR)/gkval_validator.o: gkval_validator.c gkval_validator.h gkval_core.h
$(BUILD_DIR)/gkval_adapters.o: gkval_adapters.c gkval_adapters.h gkval_validator.h gkval_core.h
$(BUILD_DIR)/gkval_cli.o: gkval_cli.c gkval_core.h gkval_validator.h gkval_adapters.h
$(FAILWIN_TEST): $(TEST_OBJS) $(CORE_OBJS) $(BUILD_DIR)/gkval_validator.o
	$(CC) $(CFLAGS) -o $@ $(TEST_OBJS) $(CORE_OBJS) $(BUILD_DIR)/gkval_validator.o $(LDFLAGS)
