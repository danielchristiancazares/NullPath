# Simple CUDA Makefile for NullPath (blackhole project)

NVCC ?= nvcc
STD ?= c++17
SM ?= 89
ARCH := sm_$(SM)
CXXFLAGS ?= -O3 -std=$(STD) -Xcompiler -Wall
DBGFLAGS ?= -G -lineinfo -std=$(STD) -Xcompiler "-g -O0 -Wall -Wextra"

BIN := bin
APP := $(BIN)/nullpath
TEST := $(BIN)/bh_sanity_test
TEST_DERIVS := $(BIN)/bh_kerr_derivs
RENDER := $(BIN)/render

.PHONY: all app test run-test clean debug

all: app test

$(BIN):
	@mkdir -p $(BIN)

app: $(BIN)
	$(NVCC) $(CXXFLAGS) -arch=$(ARCH) -Iinclude schwarzschild_blackhole.cu -o $(APP)

debug: $(BIN)
	$(NVCC) $(DBGFLAGS) -arch=$(ARCH) -Iinclude schwarzschild_blackhole.cu -o $(APP)

test: $(BIN)
	$(NVCC) -O2 -std=$(STD) -arch=$(ARCH) -Iinclude -DBLACKHOLE_NO_MAIN \
		tests/bh_sanity_test.cu schwarzschild_blackhole.cu -o $(TEST)

run-test: test
	./$(TEST)

quick-test: test
	BH_NUM_RAYS=$${BH_NUM_RAYS:-4096} ./$(TEST)

# Analytic derivative self-test (host-only)
.PHONY: test-derivs run-derivs
test-derivs: $(BIN)
	$(NVCC) -O2 -std=$(STD) -arch=$(ARCH) -Iinclude tests/test_kerr_derivs.cu -o $(TEST_DERIVS)

run-derivs: test-derivs
	./$(TEST_DERIVS)

clean:
	rm -rf $(BIN)

render: $(BIN)
	$(NVCC) -O3 -std=$(STD) -arch=$(ARCH) -Iinclude render.cu -o $(RENDER)

render-run: render
	./$(RENDER) --out render.ppm --w 1280 --h 720 --fov 50 --mass 10 --rcam 50 --samples 1

render-fast: $(BIN)
	$(NVCC) -O3 -use_fast_math -std=$(STD) -arch=$(ARCH) -Iinclude render.cu -o $(RENDER)

render-4k: render-fast
	./$(RENDER) --out render_4k.ppm --w 3840 --h 2160 --fov 55 --mass 10 --rcam 50 --samples 2 --tile 256

# Auto-detect SM via nvidia-smi (first GPU). Use: make auto
SM_DETECTED := $(shell nvidia-smi --query-gpu=compute_cap --format=csv,noheader 2>/dev/null | head -n1 | tr -d '.' )

.PHONY: auto auto-render auto-app auto-test
auto:
	@if [ -n "$(SM_DETECTED)" ]; then \
		echo "[auto] Detected SM=$(SM_DETECTED)"; \
		$(MAKE) SM=$(SM_DETECTED) all; \
	else \
		echo "[auto] Could not detect SM; building with SM=$(SM)"; \
		$(MAKE) all; \
	fi

auto-render:
	@if [ -n "$(SM_DETECTED)" ]; then \
		echo "[auto] Detected SM=$(SM_DETECTED)"; \
		$(MAKE) SM=$(SM_DETECTED) render; \
	else \
		$(MAKE) render; \
	fi

auto-app:
	@if [ -n "$(SM_DETECTED)" ]; then \
		echo "[auto] Detected SM=$(SM_DETECTED)"; \
		$(MAKE) SM=$(SM_DETECTED) app; \
	else \
		$(MAKE) app; \
	fi

auto-test:
	@if [ -n "$(SM_DETECTED)" ]; then \
		echo "[auto] Detected SM=$(SM_DETECTED)"; \
		$(MAKE) SM=$(SM_DETECTED) test; \
	else \
		$(MAKE) test; \
	fi

# Encode frames to MP4 using ffmpeg
.PHONY: encode
encode:
	@if ! command -v ffmpeg >/dev/null 2>&1; then \
		echo "ffmpeg not found. Install it or run: sudo apt install ffmpeg"; \
		exit 1; \
	fi
	@if [ -z "$(PAT)" ]; then \
		echo "Usage: make encode PAT=frame_%05d.ppm OUT=out.mp4 FPS=24"; \
		echo "       PAT uses printf-style numbering (e.g., frame_%05d.ppm)"; \
		exit 2; \
	fi
	@if [ -z "$(OUT)" ]; then OUT=out.mp4; else OUT=$(OUT); fi; \
	FPS=$${FPS:-24}; \
	ffmpeg -y -framerate $$FPS -i $(PAT) -c:v libx264 -preset slow -crf 18 -pix_fmt yuv420p $$OUT

.PHONY: render-long
render-long: render-fast
	@echo "[info] Starting long render: 4K, spin=0.9, spp-total=64 (4 spp x 16 passes)"
	./$(RENDER) \
	  --out kerr_4k.ppm --w 3840 --h 2160 --fov 55 \
	  --mass 10 --rcam 50 --theta 80 --spin 0.9 \
	  --spp-total 64 --samples 4 --checkpoint-every 4 --tile 256
