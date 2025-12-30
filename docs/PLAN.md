# NullPath — Research-Grade Accuracy Implementation Plan

Purpose
- Make NullPath a scientifically rigorous, reproducible null‑geodesic tracer and renderer that matches or exceeds current gold standards (e.g., GeoKerr for geodesics; Gyoto/RAPTOR class for ray tracing) on accuracy, with clearly quantified error bars and validation artifacts.

Scope
- Primary focus: scientific accuracy ("research" mode). Non‑research/fast modes are out of scope for this plan.
- Spacetimes: Schwarzschild (now), Kerr (next). Only vacuum metrics; matter effects and polarization are out of scope for this phase.

Outcomes
- A verified geodesic solver with tight, documented accuracy targets; reproducible comparison datasets vs. literature and established codes; images and tables that can be cited.

--------------------------------------------------------------------------------

## 1) Definitions, Units, and Conventions

- Units
  - Geometric units (G = c = 1). Report both M and Rs=2M units in logs and artifacts.
- Impact parameter
  - b ≡ L/E at infinity. Always additionally report b/M.
- Critical scales (Schwarzschild)
  - Photon sphere: r_ph = 3M = 1.5 Rs.
  - Critical impact: b_crit = 3√3 M = (3√3/2) Rs.
- r_min
  - Turning radius (if it exists): the largest positive real root of the radial turning-point equation.
  - For captured rays (no turning point), r_turning is undefined; report NaN and a separate closest_radius_reached for diagnostics.
- Deflection and lensing angles
  - Store total Δφ along the path; define lensing bend α = Δφ − π for escaped rays.
- Observer tetrads
  - Schwarzschild: static observer.
  - Kerr: ZAMO by default (static only where timelike, i.e., g_tt < 0). Document the choice in metadata.

--------------------------------------------------------------------------------

## 2) Accuracy Targets (Research Mode)

Schwarzschild
- r_min(b): relative error ≤ 1e−12 for b ≥ 1.05 b_crit; ≤ 1e−6 down to 1.01 b_crit.
- b_crit: |b_num − 3√3 M| / (3√3 M) ≤ 1e−10.
- Deflection α(b):
  - Weak field (b ≥ 20 M): |α_num − α_PN| ≤ 1e−3.
  - Near-critical fit α(b) ≈ −c1 ln(b/b_crit − 1) + c2: c1 within 1% of the known analytic value; c2 within 2%.
- Invariants drift along accepted steps: |H|max ≤ 1e−12, |ΔE/E|, |ΔL/E| ≤ 1e−12.

Kerr
- Match reference tool (GeoKerr) over a modest grid of (a/M, b/M, θ_obs):
  - Δφ agreement ≤ 1e−6 rad (equatorial); ≤ 1e−5 off‑equatorial.
  - Identical fate (escape/capture/spherical) classification.
- Invariants with Carter constant Q: |H|max ≤ 1e−12, |ΔE/E|, |ΔLz/E|, |ΔQ/Q| ≤ 1e−10.

Renderer (when relevant)
- Equatorial thin disk intersections reproducible to ≤ 1e−9 r.
- Golden images PSNR ≥ 60 dB relative to research baseline (identical seeds/settings).

--------------------------------------------------------------------------------

## 3) Architecture Choices (Research Mode)

- Precision: FP64 in all geodesic math and event detection; FP32 allowed only for shading and final buffer packing.
- Integrator
  - Schwarzschild: Hamiltonian RK with analytic turning radius and potential‑aware events; optional exact quadrature for φ(r) to audit α(b).
  - Kerr: Prefer separated form in Mino time (λ̃) using Carter potentials R(r), Θ(θ) for research mode; keep Hamiltonian form available for cross‑checks.
- Coordinates
  - Use BL for far field; near the horizon or when g_tt → 0 (Kerr), prefer horizon‑regular forms (Eddington–Finkelstein or Kerr–Schild) or terminate with clear classification before ill‑conditioning.
- Adaptive error control
  - Embedded high‑order method (RKF78 / DOP853) or step‑doubling with strict rejection; invariant‑gated acceptance (H, E, Lz, Q).
- Events and classification
  - Escape: asymptotic check (monotone r beyond r_esc and positive radial potential), not just r threshold.
  - Turning points: explicit detection from potentials; analytic turning radius in Schwarzschild, robust root‑finding in Kerr.
  - Disk intersections: root‑find θ − π/2 between brackets with interpolation safe-guards; sort multi‑hits; self‑occlusion optional.

--------------------------------------------------------------------------------

## 4) Workstreams and Deliverables

WS‑0 — Research Mode Baseline and Governance (1–2 days)
- Add a single “research” mode switch that forces FP64, disables fast‑math, and sets strict integrator tolerances and invariant gates.
- CLI: `--mode research` is the only documented mode in this plan.
- Metadata JSON sidecar for every product: commit, compile flags, GPU, CC, driver, integrator kind/order, tolerances, scene params, RNG seed, timing.
- Deliverable: example run produces image/table + JSON; deterministic across runs on the same GPU/driver.

WS‑1 — Schwarzschild Accuracy Pack (4–6 days)
1. Exact turning radius and classification
   - Implement closed‑form r_min(b) from cubic r^3 − b^2 r + b^2 Rs = 0 for b > b_crit.
   - Output arrays: turning_radius (NaN if none), closest_radius_reached (for captured/diagnostics).
2. Escape and deflection
   - Replace “r>threshold” with potential‑aware escape: require r≥r_esc and monotone increase with R(r)>0 for k consecutive steps.
   - Validate α(b) with PN series at large b and strong‑deflection fit near b_crit (Bozza): report residuals.
3. Initial condition independence
   - Demonstrate α(b) invariance for r0∈{500,1000,2000} Rs (Δα ≤ 1e−6). Prefer turning‑point init for audit path.
4. Convergence & invariants
   - Adaptive stepping + invariant gates (H, E, L) with reject/halve on violation; plot error vs. step.
5. Tests & artifacts
   - Unit tests: b_crit, r_min(b) grid, α(b) vs PN/strong‑deflection; convergence slopes.
   - Artifacts: CSV tables, plots (saved via scripts), and a “Schwarzschild validation” README with figures.

WS‑2 — Kerr Accuracy Pack (8–12 days)
1. Separated equations in Mino time
   - Implement R(r), Θ(θ) potentials, Carter constant Q; integrate dφ/dλ̃, dt/dλ̃ alongside dr/dλ̃, dθ/dλ̃.
   - Robust root‑finding at turning points; event functions for spherical photon orbits (R=0, dR/dr=0).
2. Observer and emitter frames
   - Use ZAMO tetrad by default. Provide static only where g_tt < 0 (assert otherwise).
   - Disk inner radius r_in = r_isco(a) (analytic expression); expose override.
3. Validation grid
   - Build a modest set of cases (e.g., a/M ∈ {0,0.5,0.9,0.99}, b/M spanning near‑critical to weak‑field; θ_obs ∈ {π/2, π/3}).
   - Cross‑validate Δφ, turning radii, and fates vs GeoKerr.
4. Invariants and convergence
   - Track H, E, Lz, Q drift; require research‑mode thresholds. Show p‑order slopes.
5. Tests & artifacts
   - GPU tests gated by env; reference CSV from GeoKerr runs (stored under `geokerr_port/validation/test_data/`).
   - Plots and summary tables committed.

WS‑3 — Research‑Grade Disk Rendering (Optional, 4–6 days)
1. Disk geometry and intersections
   - Multi‑hit traversal; sort by affine parameter; self‑occlusion enabled.
   - Robust equatorial crossing root‑find with angle guarding.
2. Radiometry
   - I_obs = g^3 I_em with I_em(r) = r^−p; record p in metadata. Optionally add banded outputs (grayscale per ν band).
3. Golden images
   - Fixed seeds, metadata, PSNR ≥ 60 dB across reruns. Provide reference PNG/EXR assets.

WS‑4 — Reproducibility & Tooling (2–3 days)
1. Sidecar JSON metadata for all outputs (tables, images).
2. Experiment manifests (YAML/JSON): end‑to‑end reproducible configurations to regenerate figures.
3. Scripts under `scripts/` for: validation runs, plot generation, and GeoKerr harness.

WS‑5 — Repository Restructure and Documentation (2–3 days)
1. Layout: `src/`, `include/`, `tests/`, `scripts/`, `bin/` (as per AGENTS.md).
2. Developer docs: theory notes (b_crit, turning radius, Carter separation, ZAMO, Iν invariance), validation methodology, usage.
3. Build presets: `make research` (strict flags; FP64; no fast‑math), test targets; README badges for “research mode ready”.

--------------------------------------------------------------------------------

## 5) Validation Plan and Datasets

Analytic (Schwarzschild)
- b_crit exact value; r_min(b) closed form; α(b) PN series at large b; strong‑deflection expansion near b_crit.
- Datasets: CSV of (b/M, r_min/Rs, α, flags) for several grids; plot residuals.

Cross‑tool (Kerr & Schwarzschild)
- Use `geokerr_port/reference/*.f` to generate reference tables (M=1 units): Δφ, r_turning, and fates for curated cases.
- Store reference outputs under `geokerr_port/validation/test_data/` (already present); version them with SHA and parameter stamps.
- Comparison scripts compute Δ metrics and assert acceptance thresholds.

Convergence & Invariants
- For a set of rays: run at multiple tolerances/step sizes; fit slope p; assert expected order within ±0.2.
- Track invariant drifts; assert maxima within thresholds; save histograms.

Rendering
- Golden images (PNG/EXR) with sidecar JSON; PSNR/SSIM versus baseline.

Artifacts
- Store plots (PDF/PNG) and CSV/JSON summaries in `tests/artifacts/` with auto‑generated index.

--------------------------------------------------------------------------------

## 6) Numerical & GPU Considerations

Precision and math
- FP64 kernels for geodesics; require precise division/sqrt; prefer fused FMA.
- Avoid FTZ/DAZ for FP64 segments. Clamp only where mathematically justified (e.g., sin θ → sign‑preserving small value).

Step control
- Research mode: mandatory adaptive with strict tolerances (e.g., rel 1e−12, abs 1e−14) and invariant gates.
- Quantize accepted h only if it does not inflate error; otherwise leave continuous.

Event robustness
- Use bracketed root‑finders (bisection + secant) for events (turning points, θ=π/2 crossing) to guarantee convergence.

Parallelism and determinism
- One ray per thread; avoid race conditions; diagnostical atomics are okay but not used to alter physics decisions.
- Deterministic sampling: fixed seeds from pixel coords and pass index.

Memory and performance
- Keep pinned host buffers and tiling; overlap compute and transfers. Performance is secondary to accuracy but we should avoid pathological slowdowns.

--------------------------------------------------------------------------------

## 7) Build & Run Modes

- Research build flags (NVCC)
  - `-O3 -std=c++17 -Xcompiler -Wall -arch=sm_XX` with precise math; no `-use_fast_math`.
  - Ensure `--prec-div=true --prec-sqrt=true` if applicable; keep FMA enabled.
- Make targets
  - `make research` → builds app and tests with strict flags and FP64 for geodesics.
  - `make run-test` → full research validations on current GPU (time‑boxed by env variables).
  - `make quick-test` → small subset for smoke.

--------------------------------------------------------------------------------

## 8) File‑Level Changes (Outline Only; to be implemented under PRs)

Core geodesics
- `include/bh_common.h`: helpers for units, constants; invariant thresholds; research mode flags.
- `include/bh_schwarzschild.h`: FP64 variant; analytic turning radius helper; potential‑aware events; optional EF/Kerr–Schild note.
- `include/bh_kerr.h`: Mino‑time pathway (R, Θ, Carter Q); analytic derivatives already present; ZAMO observer default.
- `src/schwarzschild_precompute.cu` (refactor from single file): tables with turning_radius, closest_radius, deflection, bending, status.

Renderer
- `src/render.cu`: robust disk crossing root‑finding; multi‑hit; r_in=r_isco(a); invariant‑gated adaptive in research mode.

Tests & scripts
- `tests/`: schwarzschild_accuracy_tests.cu, kerr_grid_tests.cu, invariants_convergence.cu.
- `scripts/`: geokerr_harness.py, plot_schwarzschild.py, plot_kerr_grid.py, gen_manifests.py.

--------------------------------------------------------------------------------

## 9) Milestones, Timeline, and Success Criteria

M0 — Governance & Research Mode (Day 1–2)
- Done when: metadata sidecars, strict flags, deterministic runs; a minimal analytic check passes.

M1 — Schwarzschild Accuracy (Day 3–8)
- Done when: r_min(b) and α(b) datasets meet targets; convergence and invariants plots committed; b_crit passes to 1e−10.

M2 — Kerr Accuracy (Day 9–20)
- Done when: Kerr grid matches GeoKerr within targets; invariants with Q meet thresholds; spherical photon orbit detection validated.

M3 — Rendering Research Features (Day 21–26)
- Done when: disk multi‑hit with robust root‑finding; golden images with PSNR ≥ 60 dB; metadata present.

M4 — Reproducibility & Docs (Day 27–30)
- Done when: experiment manifests regenerate artifacts; docs include theory notes and validation methods; repository reorganized.

Stretch Goals
- Mino‑time + adaptive RKF78 implementation optimized for GPU; optional multi‑GPU tile distribution.

--------------------------------------------------------------------------------

## 10) Risks and Mitigations

- Near‑critical stiffness (b → b_crit)
  - Mitigate with analytic turning radii (Schwarzschild) and separated Mino‑time integration (Kerr).
- Horizon coordinate issues in BL
  - Use ZAMO observers; terminate before ill‑conditioning; consider EF/Kerr–Schild.
- Cross‑tool discrepancies (convention mismatches)
  - Normalize units; document definitions (b, r_min vs. closest radius); write comparison adapters.
- GPU nondeterminism
  - Avoid physics decisions based on atomics; log diagnostics only; fixed seeds.

--------------------------------------------------------------------------------

## 11) Review Checklist (Per PR)

- Physics
  - Equations match references; units consistent; events defined clearly; invariants bounded.
- Numerics
  - Adaptive control present; convergence evidence; tolerances justified.
- Validation
  - Tests updated; artifacts regenerated; comparisons vs. references included.
- Reproducibility
  - Metadata present; manifests updated; seeds fixed.
- Docs
  - Theory and method notes; CLI options; acceptance criteria.

--------------------------------------------------------------------------------

## 12) References (non‑exhaustive; for implementation guidance)

- S. Chandrasekhar, The Mathematical Theory of Black Holes, OUP (1983).
- V. Bozza, "Gravitational lensing in the strong field limit," Phys. Rev. D (2002) and follow‑ups.
- E. Teo, "Spherical photon orbits around a Kerr black hole," Gen. Rel. Grav. 52, 5 (2020).
- Bardeen, Press, Teukolsky (1972), rotating black holes and LNRF/ZAMO frames.
- Carter (1968), constant of motion and separability in Kerr.
- GeoKerr (Hamilton + collaborators): semi‑analytic Kerr geodesics; use as reference for Δφ and r_turning.
- Gyoto, RAPTOR: ray‑tracing implementations for reference behavior and features.

--------------------------------------------------------------------------------

## 13) Immediate Next Actions (No Code Changes Until Approved)

1) Confirm r_min convention (turning_radius NaN for captured; keep closest_radius_reached separately).
2) Approve WS‑0 and WS‑1 as first implementation batches.
3) Identify the initial Schwarzschild b grid for validation datasets (e.g., 200 logarithmically spaced b values from 1.01 b_crit to 100 M).
4) Decide on preferred analytic sources for PN and strong‑deflection coefficients (Bozza parameters) for inclusion in tests.

End of PLAN.md

