# GEOKER bit‑exact parity roadmap

This plan drives the CUDA GeoKerr (GPU) implementation to bit‑exact parity with the reference (CPU) along a single, reproducible configuration. It breaks the work into strict determinism controls, math library alignment, algorithmic alignment, and verification hardening.

Success criteria:
- For the target platform (specified below), every validated field sample matches bit‑for‑bit: ULP distance = 0 across the entire grid and all cases.
- Validator passes with quantization disabled (quant_lsb = 0) and ULP check enabled, with atol = rtol = 0 in “bit‑exact mode”.

Scope and assumptions:
- Target a single deterministic hardware/software stack first; broaden later.
- We accept any performance loss during the bring‑up; we can re‑optimize under constraints once parity is achieved.

Target baseline (locked):
- GPU arch: one specific SM (e.g., sm_86 on the current machine; adjust as needed).
- CUDA: the version installed in this repo’s environment (record exact version).
- Compiler(s): GCC/Clang exact version; NVCC exact version.
- OS/driver: current WSL/driver versions recorded.

Deliverables:
- Deterministic build configurations (Makefile updates) for CPU and GPU.
- Shared deterministic math layer used by both CPU and GPU paths.
- Extended validator configuration: bit‑exact mode and detailed mismatch tracing.
- Documentation and CI check that enforces bit‑exact on the target stack.

Bit‑exact profile defaults (recommended for initial bring‑up):
- Grid size S: start small (e.g., 256–512) to speed iteration.
- Grid spacing delta: 1e-4 (tunable). 1e-3 is acceptable but coarser; smaller deltas help isolate early‑step drift.
- Quantization: 0 (disabled).
- Tolerances: atol=0, rtol=0, max_ulps=0.

---

## Phase 0 — Baseline capture and reproducibility

1) Inventory and freeze environment
- Record: GPU model + SM, NVIDIA driver, CUDA toolkit, NVCC, GCC/Clang.
- Lock compiler flags via Makefile to a new "deterministic" config.

2) Determinism audit (current code)
- Ensure no parallel reductions or data‑dependent iteration counts that differ between CPU/GPU.
- Ensure branch thresholds/epsilons identical and centralized.
- Ensure no non‑deterministic atomics or unordered writes are used in validated paths.

3) Validation harness hardening
- Re‑enable ULP check in `gkval_core.c` and add a "bit‑exact" profile:
  - tolerances: atol = 0, rtol = 0, max_ulps = 0
  - quantization: quant_lsb = 0.
- Add option to dump first mismatch with: index, field, a_bits, b_bits, ulp, and the expression lineage if available.

Exit criteria: repeatable runs (same hashes) on the target machine with today’s code/flags, even if not yet bit‑exact.

---

## Phase 1 — Strict FP semantics alignment (toolchain & flags)

Goal: eliminate compiler/runtime FP variance between CPU and GPU.

A) Host compiler flags (CPU reference build)
- Disable fast‑math and contraction; standard rounding and precision:
  - `-O2 -fno-fast-math -fno-associative-math -fno-unsafe-math-optimizations -ffp-contract=off`
  - `-fexcess-precision=standard -frounding-math`
  - On x86‑64, force SSE (avoid x87 extended precision): `-mfpmath=sse -msse2`
- Do not use `-ffloat-store` unless a last‑resort (hurts perf and not usually needed with SSE2).

B) Device compiler flags (CUDA build)
- Enforce IEEE‑like semantics and disable fused transformations unless used explicitly:
  - `--fmad=false` (or use explicit `fma` on both sides; see Phase 2)
  - `--prec-div=true --prec-sqrt=true`
  - `-ftz=false` (preserve subnormals)
- Do not use `--use_fast_math`.
- Target a single SM: `-arch=sm_XX -code=sm_XX` where XX is the chosen baseline.

C) Runtime FP state
- CPU: ensure FE_TONEAREST and that FTZ/DAZ are disabled (already implemented in `gkval_validate_fp_env`).
- GPU: rely on compile options above (no fast‑math, no FTZ), and avoid intrinsics with unspecified rounding.

Exit criteria: micro‑kernels using only +, −, ×, ÷, sqrt produce identical bit patterns CPU vs GPU for controlled inputs.

---

## Phase 2 — Deterministic math layer (shared CPU/GPU)

Goal: remove libm/libdevice drift; unify the exact operations and rounding points used in the algorithm.

1) Create `include/geokerr_math.h` + `src/geokerr_math.cu`
- Provide a minimal set of primitives as inline `__host__ __device__` functions for double:
  - add/sub/mul with explicit rounding on device: `__dadd_rn`, `__dsub_rn`, `__dmul_rn` (CUDA intrinsics).
  - `fma_rn(a,b,c)`: on CPU call `fma(a,b,c)`, on GPU call `__fma_rn(a,b,c)`.
  - `div_rn(a,b)`: device relies on precise division; host uses `/`.
  - `sqrt_rn(x)`: device uses `sqrt(x)` with `--prec-sqrt=true`; host uses `sqrt`.
- Replace compound expressions in hot paths with these primitives to lock operation order and rounding.

2) Transcendentals inventory
- Catalog all uses of `pow`, `exp`, `log`, `sin`, `cos`, `atan2`, etc.
- Prefer algorithmic rewrites to avoid transcendentals where possible (e.g., fixed powers via multiplies; `exp(log(x)*n)` replaced by repeated mul for small integer n).
- If any transcendental remains in the validated path, consider a shared, deterministic implementation:
  - Option A (best): integrate a correctly‑rounded library usable on CPU, and port a numerically identical path to CUDA (limited surface, e.g., `log`, `exp`).
  - Option B (pragmatic): implement the specific approximations we need in header‑only form and use identically on CPU/GPU.

3) Constants hygiene
- Express critical constants as hexadecimal floating literals to fix exact binary values (e.g., `0x1.921fb54442d18p+1` for π).
- Centralize all epsilons/tolerances in one header; avoid duplicated literals across files.

4) Operation ordering
- Parenthesize expressions to pin evaluation order.
- Replace fragile algebraic rewrites with explicit sequencing via the math wrappers.

Exit criteria: unit tests for the math layer show CPU/GPU identity at the bit level across a representative input corpus.

---

## Phase 3 — Algorithmic alignment and control flow stability

Goal: ensure both implementations traverse identical branches and iteration counts.

1) Carlson elliptic integrals path
- Ensure both CPU and GPU call the same source for Carlson RF/RC/RD/RJ with identical iteration limits and termination criteria.
- Use the shared math primitives for all steps to lock rounding and sequence.

2) Iterative solvers / steppers
- Fix maximum iteration counts and stopping conditions identically.
- Avoid early‑exit conditions that depend on sub‑ULP differences; compare against centralized epsilons.

3) Summations and reductions
- Replace naive accumulation with Kahan or Neumaier compensated summation on both sides if accumulation length > O(10^3).
- Keep summation order fixed (no parallel tree reductions in validated code paths).

Exit criteria: trace logs of branch decisions (optional debug mode) show identical sequences CPU vs GPU for representative cases.

---

## Phase 4 — Validator: bit‑exact mode and diagnostics

1) Re‑enable and enforce ULP checks
- Set `max_ulps = 0` for bit‑exact runs; set `atol = rtol = 0` and `quant_lsb = 0`.
- Add a CLI preset in `gkval_cli` (e.g., `--profile bitexact`) that configures the manifest accordingly.

2) Mismatch forensics
- On first mismatch, dump:
  - sample index, field, `a_bits`, `b_bits`, `a_hex`, `b_hex`, ULP distance.
  - optional: a compact provenance record (which high‑level function produced the value) via tagged checkpoints in code.

3) Reproducibility gate
- Add a CI job on the target runner that builds CPU and GPU in deterministic mode and runs the validator. Failing parity blocks merges to main.

---

## Phase 5 — Performance recapture (post‑parity)

- Gradually re‑enable optimizations that don’t change results (prove with validator):
  - Allow explicit `fma_rn` where mathematically equivalent; keep same on CPU using `fma`.
  - Evaluate `-Xptxas` optimizations that preserve rounding (validate each change).
  - Consider mixed precision only for non‑validated paths.

---

## Work breakdown and timeline (first pass)

Week 1
- Implement deterministic build targets (Makefile updates for flags above).
- Re‑enable ULP check and add bit‑exact profile in validator.
- Add math layer scaffolding and micro‑tests (CPU vs GPU) for +, −, ×, ÷, fma, sqrt.

Week 2
- Port Carlson integrals to use math layer; replace literals with hex‑floats.
- Instrument branch traces (debug) and align iteration thresholds.
- Achieve bit‑exact on a single end‑to‑end case with reduced grid.

Week 3
- Scale to full test suite; address any residual mismatches with targeted diffs.
- Add CI gate and documentation. Start performance recapture experiments.

---

## Makefile notes (to implement)

CPU (reference):
- CFLAGS (add): `-O2 -fno-fast-math -fno-associative-math -fno-unsafe-math-optimizations -ffp-contract=off -fexcess-precision=standard -frounding-math -mfpmath=sse -msse2`

CUDA:
- NVCC flags (device): `--fmad=false --prec-div=true --prec-sqrt=true -ftz=false -arch=sm_XX -code=sm_XX`
- NVCC host pass‐through: add CPU flags via `-Xcompiler` as appropriate.

Validator bit‑exact preset:
- `atol=0`, `rtol=0`, `max_ulps=0`, `quant_lsb=0`, grid unchanged.

---

## Risks and mitigations

- Libm vs libdevice differences: mitigate by using our shared math layer and avoiding generic libm/libdevice in validated paths.
- Subnormal/rounding behavior: enforce via flags (`-ftz=false`, SSE DAZ/FTZ off) and wrappers.
- Compiler re‑association: disabled via flags and by explicit parentheses and `ffp-contract=off`.
- Hardware variance: start with one SM; later add per‑SM golden baselines and/or per‑SM builds.
- Effort creep for transcendentals: minimize by avoiding them or isolating to a small, shared implementation.

---

## Definition of Done (bit‑exact)

- Validator run in bit‑exact mode passes across all selected cases and fields with ULP=0 and zero tolerances.
- Repeatable across reruns on the locked target machine and SM.
- CI job enforces the above on every change.
