# Simple CUDA Makefile for NullPath (blackhole project)

NVCC ?= nvcc
STD ?= c++17
SM ?= 89
ARCH := sm_$(SM)
COMPAT := -allow-unsupported-compiler
CXXFLAGS ?= -O3 -std=$(STD) $(COMPAT) -Xcompiler -Wall
DBGFLAGS ?= -G -lineinfo -std=$(STD) $(COMPAT) -Xcompiler "-g -O0 -Wall -Wextra"

# Defaults for run targets (override like: make quick-test BH_NUM_RAYS=8192)
BH_NUM_RAYS ?= 4096
FPS ?= 24
OUT ?= out.mp4

ifeq ($(OS),Windows_NT)
EXEEXT := .exe
else
EXEEXT :=
endif

BIN := bin
APP := $(BIN)/nullpath
TEST := $(BIN)/bh_sanity_test
TEST_DERIVS := $(BIN)/bh_kerr_derivs
RENDER := $(BIN)/render

APP_EXE := $(APP)$(EXEEXT)
TEST_EXE := $(TEST)$(EXEEXT)
TEST_DERIVS_EXE := $(TEST_DERIVS)$(EXEEXT)
RENDER_EXE := $(RENDER)$(EXEEXT)

.PHONY: all app test run-test clean debug

all: app test

$(BIN):
ifeq ($(OS),Windows_NT)
	@pwsh -NoProfile -Command "New-Item -ItemType Directory -Force -Path '$(BIN)' | Out-Null"
else
	@mkdir -p $(BIN)
endif

app: $(BIN)
	$(NVCC) $(CXXFLAGS) -arch=$(ARCH) -Iinclude schwarzschild_blackhole.cu -o $(APP)

debug: $(BIN)
	$(NVCC) $(DBGFLAGS) -arch=$(ARCH) -Iinclude schwarzschild_blackhole.cu -o $(APP)

test: $(BIN)
	$(NVCC) -O2 -std=$(STD) $(COMPAT) -arch=$(ARCH) -Iinclude -DBLACKHOLE_NO_MAIN \
		tests/bh_sanity_test.cu schwarzschild_blackhole.cu -o $(TEST)

run-test: test
ifeq ($(OS),Windows_NT)
	@pwsh -NoProfile -Command "$$env:BH_NUM_RAYS='$(BH_NUM_RAYS)'; & '$(TEST_EXE)'"
else
	BH_NUM_RAYS=$(BH_NUM_RAYS) ./$(TEST_EXE)
endif

quick-test: test
ifeq ($(OS),Windows_NT)
	@pwsh -NoProfile -Command "$$env:BH_NUM_RAYS='$(BH_NUM_RAYS)'; & '$(TEST_EXE)'"
else
	BH_NUM_RAYS=$(BH_NUM_RAYS) ./$(TEST_EXE)
endif

# Analytic derivative self-test (host-only)
.PHONY: test-derivs run-derivs
test-derivs: $(BIN)
	$(NVCC) -O2 -std=$(STD) $(COMPAT) -arch=$(ARCH) -Iinclude tests/test_kerr_derivs.cu -o $(TEST_DERIVS)

run-derivs: test-derivs
ifeq ($(OS),Windows_NT)
	@pwsh -NoProfile -Command "& '$(TEST_DERIVS_EXE)'"
else
	./$(TEST_DERIVS_EXE)
endif

clean:
ifeq ($(OS),Windows_NT)
	@pwsh -NoProfile -Command "if (Test-Path '$(BIN)') { Remove-Item -Recurse -Force '$(BIN)' }"
else
	rm -rf $(BIN)
endif

render: $(BIN)
	$(NVCC) -O3 -std=$(STD) $(COMPAT) -arch=$(ARCH) -Iinclude render.cu -o $(RENDER)

render-run: render
ifeq ($(OS),Windows_NT)
	@pwsh -NoProfile -Command "& '$(RENDER_EXE)' --out render.ppm --w 1280 --h 720 --fov 50 --mass 10 --rcam 50 --samples 1"
else
	./$(RENDER_EXE) --out render.ppm --w 1280 --h 720 --fov 50 --mass 10 --rcam 50 --samples 1
endif

render-fast: $(BIN)
	$(NVCC) -O3 -use_fast_math -std=$(STD) $(COMPAT) -arch=$(ARCH) -Iinclude render.cu -o $(RENDER)

render-4k: render-fast
ifeq ($(OS),Windows_NT)
	@pwsh -NoProfile -Command "& '$(RENDER_EXE)' --out render_4k.ppm --w 3840 --h 2160 --fov 55 --mass 10 --rcam 50 --samples 2 --tile 256"
else
	./$(RENDER_EXE) --out render_4k.ppm --w 3840 --h 2160 --fov 55 --mass 10 --rcam 50 --samples 2 --tile 256
endif

# Auto-detect SM via nvidia-smi (first GPU). Use: make auto
ifeq ($(OS),Windows_NT)
SM_DETECTED := $(shell pwsh -NoProfile -Command "$$cap = (nvidia-smi --query-gpu=compute_cap --format=csv,noheader 2>$$null | Select-Object -First 1); if ($$cap) { ($$cap.Trim()).Replace('.', '') }")
else
SM_DETECTED := $(shell nvidia-smi --query-gpu=compute_cap --format=csv,noheader 2>/dev/null | head -n1 | tr -d '.' )
endif

.PHONY: auto auto-render auto-app auto-test
auto:
ifneq ($(strip $(SM_DETECTED)),)
	@echo "[auto] Detected SM=$(SM_DETECTED)"
	@$(MAKE) SM=$(SM_DETECTED) all
else
	@echo "[auto] Could not detect SM; building with SM=$(SM)"
	@$(MAKE) all
endif

auto-render:
ifneq ($(strip $(SM_DETECTED)),)
	@echo "[auto] Detected SM=$(SM_DETECTED)"
	@$(MAKE) SM=$(SM_DETECTED) render
else
	@$(MAKE) render
endif

auto-app:
ifneq ($(strip $(SM_DETECTED)),)
	@echo "[auto] Detected SM=$(SM_DETECTED)"
	@$(MAKE) SM=$(SM_DETECTED) app
else
	@$(MAKE) app
endif

auto-test:
ifneq ($(strip $(SM_DETECTED)),)
	@echo "[auto] Detected SM=$(SM_DETECTED)"
	@$(MAKE) SM=$(SM_DETECTED) test
else
	@$(MAKE) test
endif

# Encode frames to MP4 using ffmpeg
.PHONY: encode
encode:
ifeq ($(strip $(PAT)),)
	@echo "Usage: make encode PAT=frame_%05d.ppm OUT=out.mp4 FPS=24"
	@echo "       PAT uses printf-style numbering (e.g., frame_%05d.ppm)"
	@exit 2
endif
ifeq ($(OS),Windows_NT)
	@pwsh -NoProfile -Command "if (-not (Get-Command ffmpeg -ErrorAction SilentlyContinue)) { Write-Error 'ffmpeg not found (not in PATH). Install it and try again.'; exit 1 }"
else
	@command -v ffmpeg >/dev/null 2>&1 || { echo "ffmpeg not found. Install it and try again."; exit 1; }
endif
	ffmpeg -y -framerate $(FPS) -i $(PAT) -c:v libx264 -preset slow -crf 18 -pix_fmt yuv420p $(OUT)

.PHONY: render-long
render-long: render-fast
	@echo "[info] Starting long render: 4K, spin=0.9, spp-total=64 (4 spp x 16 passes)"
ifeq ($(OS),Windows_NT)
	@pwsh -NoProfile -Command "& '$(RENDER_EXE)' --out kerr_4k.ppm --w 3840 --h 2160 --fov 55 --mass 10 --rcam 50 --theta 80 --spin 0.9 --spp-total 64 --samples 4 --checkpoint-every 4 --tile 256"
else
	./$(RENDER_EXE) \
	  --out kerr_4k.ppm --w 3840 --h 2160 --fov 55 \
	  --mass 10 --rcam 50 --theta 80 --spin 0.9 \
	  --spp-total 64 --samples 4 --checkpoint-every 4 --tile 256
endif
