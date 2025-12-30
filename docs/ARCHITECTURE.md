# NullPath Architecture

This document describes the architecture, design decisions, and physics conventions used in the NullPath black hole ray tracer.

## Overview

NullPath is a CUDA-accelerated null geodesic integrator and black hole renderer. It traces photon paths through curved spacetime around Schwarzschild (non-rotating) and Kerr (rotating) black holes, producing physically accurate visualizations of gravitational lensing, photon spheres, and accretion disk emission.

```
+------------------+     +-------------------+     +------------------+
|   Application    |     |    Renderer       |     |   Precompute     |
| (nullpath demo)  |     |   (render.cu)     |     |   (lookup tables)|
+--------+---------+     +---------+---------+     +--------+---------+
         |                         |                        |
         v                         v                        v
+--------+---------+     +---------+---------+     +--------+---------+
|  bh_schwarzschild.h   |  bh_kerr.h        |     |  bh_schwarzschild.h
|  (Schwarzschild       |  (Kerr metric,    |     |  (geodesic        |
|   geodesics)          |   geodesics)      |     |   integration)    |
+--------+---------+     +---------+---------+     +--------+---------+
         |                         |                        |
         +-----------+-------------+------------------------+
                     |
                     v
         +-----------+-----------+
         |      bh_types.h       |
         |  (GeodesicState,      |
         |   SchwarzschildMetric,|
         |   GeodesicLookupTable)|
         +-----------+-----------+
                     |
                     v
         +-----------+-----------+
         |     bh_common.h       |
         |  (constants, helpers, |
         |   RayStatus, macros)  |
         +------------------------+
```

## Module Hierarchy

### Layer 0: Foundation (`bh_common.h`)

The base layer providing:
- **Physical constants**: Pi, epsilon values for numerical stability
- **Unit system constants**: Schwarzschild radius multiplier (Rs = 2M), photon sphere multiplier (r_ph = 1.5 Rs)
- **Integrator parameters**: Default step sizes, RK tolerances, step limits
- **Utility functions**: Clamping, safe division, safe sin^2 for pole handling
- **Ray classification**: `RayStatus` enum (Escaped, Captured, PhotonSphere, Inconclusive)
- **CUDA error handling**: `BH_CUDA_CHECK` macro

All other modules depend on this layer.

### Layer 1: Type Definitions (`bh_types.h`)

Core data structures used throughout:

- **`SchwarzschildMetric`**: Stores black hole mass M and Schwarzschild radius Rs. Provides covariant metric component accessors `g_tt(r)`, `g_rr(r)`, etc.

- **`GeodesicState`**: The 8-component phase space state (t, r, theta, phi, pt, pr, ptheta, pphi) representing a photon's position and covariant 4-momentum.

- **`GeodesicLookupTable`**: Host-side storage for precomputed geodesic data including impact parameters, deflection angles, bending angles, closest approaches, and ray status arrays.

### Layer 2: Metric-Specific Physics

#### Schwarzschild (`bh_schwarzschild.h`)

Implements geodesic integration in Schwarzschild spacetime:
- **Hamiltonian**: `hamiltonian_schwarzschild()` computes H = (1/2) g^{uv} p_u p_v (should be ~0 for null geodesics)
- **Derivatives**: `computeGeodesicDerivatives()` computes dx^u/dlambda and dp_u/dlambda using Hamilton's equations
- **Integrator**: `integrateGeodesic()` performs one RK4 step

#### Kerr (`bh_kerr.h`)

Implements geodesic integration in Kerr spacetime (Boyer-Lindquist coordinates):
- **Metric blocks**: `kerr_blocks_cov()` and `kerr_blocks_contra()` compute covariant and contravariant metric components
- **Analytic derivatives**: `kerr_blocks_contra_derivs()` provides closed-form partial derivatives of inverse metric w.r.t. r and theta
- **Hamiltonian**: `hamiltonian_kerr()` for constraint monitoring
- **Derivatives and integrators**: Both FP32 and FP64 implementations
- **Camera initialization**: `init_ray_from_camera_general()` (static observer) and `init_ray_from_camera_zamo()` (ZAMO frame)
- **Doppler factor**: `doppler_g_general()` computes redshift/blueshift for disk emission
- **Adaptive stepping**: `rk4_adaptive_step_d()` with error control via step-doubling

## Data Flow

### Precomputation Path (schwarzschild_blackhole.cu)

```
1. Generate impact parameter array b[] (logarithmic spacing around b_crit)
2. Copy b[] to GPU
3. Launch precomputeGeodesicsKernel:
   For each ray:
   a. Initialize GeodesicState at r_start (effective infinity, >=1000 and >~1.2*b_max)
   b. Set p_t = -E, p_phi = bE, solve for p_r from H=0
   c. Integrate inward, tracking:
      - Minimum radius (refined at turning point)
      - Total phi change (deflection)
      - Termination condition (escape/capture)
   d. Store deflection_angle, closest_approach, ray_status
4. Copy results to host
5. Compute bending_angle = deflection - pi for escaped rays
6. Build GeodesicLookupTable
```

### Rendering Path (render.cu)

```
1. Parse command-line arguments into RenderParams
2. For each tile (rows processed in batches):
   For each pixel (x, y):
   a. Compute ray direction in camera frame
   b. Transform to local tetrad (static or ZAMO)
   c. Initialize GeodesicState with covariant momenta
   d. Integrate geodesic (backward ray tracing):
      - Check for disk crossings (equatorial plane)
      - Accumulate disk emission with Doppler beaming
      - Track minimum radius
   e. If ray escapes: sample background
   f. Apply exposure, tone mapping, gamma
   g. Write pixel to output buffer
3. Accumulate across SPP passes
4. Write final PPM image
```

## Physics Conventions

### Unit System

NullPath uses **geometric units** where G = c = 1:
- Mass M is the fundamental scale
- Schwarzschild radius: Rs = 2M
- Photon sphere radius: r_ph = 3M = 1.5 Rs
- Critical impact parameter: b_crit = 3*sqrt(3)*M = (3*sqrt(3)/2)*Rs

All lengths are measured in units of M (or equivalently Rs/2). Time has the same units as length.

### Coordinate Systems

**Schwarzschild**: Standard Schwarzschild coordinates (t, r, theta, phi)
- Line element: ds^2 = -(1-Rs/r)dt^2 + (1-Rs/r)^{-1}dr^2 + r^2(dtheta^2 + sin^2(theta)dphi^2)

**Kerr**: Boyer-Lindquist coordinates (t, r, theta, phi)
- Sigma = r^2 + a^2 cos^2(theta)
- Delta = r^2 - 2Mr + a^2
- A = (r^2 + a^2)^2 - a^2 Delta sin^2(theta)
- The metric has off-diagonal g_{t,phi} term encoding frame dragging

### Sign Conventions

- **Metric signature**: (-,+,+,+) (timelike intervals have ds^2 < 0)
- **4-momentum**: Covariant components (p_t, p_r, p_theta, p_phi) are stored
- **Energy**: E = -p_t > 0 for future-directed photons
- **Angular momentum**: L = p_phi > 0 for prograde orbits
- **Impact parameter**: b = L/E = p_phi / (-p_t)
- **Ingoing rays**: p_r < 0 (decreasing radius)
- **Affine parameter**: lambda increases along the ray (backward tracing uses negative effective steps)

### Hamiltonian Formulation

The null geodesic equation is derived from the Hamiltonian:

```
H = (1/2) g^{uv} p_u p_v = 0
```

Hamilton's equations give:
```
dx^u / dlambda = dH/dp_u = g^{uv} p_v
dp_u / dlambda = -dH/dx^u = -(1/2) (d g^{ab}/dx^u) p_a p_b
```

For cyclic coordinates (t in stationary spacetimes, phi in axisymmetric spacetimes), the corresponding momenta are conserved: dp_t/dlambda = 0 and dp_phi/dlambda = 0.

## Integrator Design

### RK4 (Fixed Step)

Both Schwarzschild and Kerr modules provide a classic 4th-order Runge-Kutta integrator:

```
k1 = f(y_n)
k2 = f(y_n + h/2 * k1)
k3 = f(y_n + h/2 * k2)
k4 = f(y_n + h * k3)
y_{n+1} = y_n + h/6 * (k1 + 2*k2 + 2*k3 + k4)
```

Step size is scaled with radius to maintain accuracy near the black hole while allowing large steps in the far field:
```
scale = clamp(r / (0.5 * Rs), 1, 100)
h_effective = h_base * scale
```

### Adaptive RK4 (FP64 Kerr)

For high-precision work, `rk4_adaptive_step_d()` implements step-doubling error estimation:

1. Take one full step of size h: y_full
2. Take two half steps of size h/2: y_half
3. Estimate error as |y_full - y_half|
4. Accept step if error < tolerance; otherwise shrink h and retry
5. Adjust h for next step based on error ratio

Tolerances are controlled by `BH_RK_REL_TOL` (1e-10) and `BH_RK_ABS_TOL` (1e-12).

### Termination Conditions

Integration stops when:
- Ray escapes: r > r_start and p_r > 0 (outgoing at large radius)
- Ray captured: r <= Rs + epsilon (crosses event horizon)
- Theta bounds violated: theta <= 0 or theta >= pi (numerical instability at poles)
- Step limit reached: MAX_INTEGRATION_STEPS exceeded

## Camera Models

### Static Observer

Valid outside the ergosphere (where g_{tt} < 0). The observer's 4-velocity is proportional to the timelike Killing vector:

```
u^t = 1/sqrt(-g_{tt}), u^r = u^theta = u^phi = 0
```

An orthonormal tetrad is constructed with:
- e_(t) along u (time direction)
- e_(r) in the radial direction
- e_(theta) in the polar direction
- e_(phi) in the azimuthal direction

Ray directions in the camera frame are converted to covariant 4-momenta using this tetrad.

### ZAMO (Zero Angular Momentum Observer)

Valid everywhere outside the horizon. ZAMOs are dragged by the black hole's rotation but have zero angular momentum (L = 0). Their 4-velocity is:

```
u^t = 1/alpha, u^phi = omega/alpha
```

where alpha is the lapse and omega is the frame-dragging angular velocity. The ZAMO tetrad provides a locally non-rotating reference frame.

## Key Design Decisions

### Header-Only Implementation

All physics code is in header files with inline functions, enabling:
- Single compilation unit simplicity
- Automatic inlining by CUDA compiler
- Easy `__host__ __device__` dual compilation
- No link-time dependencies

### Separate FP32 and FP64 Paths

The Kerr module provides both single and double precision implementations:
- FP32 (`GeodesicState`, `integrateGeodesicKerr`): Faster, sufficient for visualization
- FP64 (`GeodesicStateD`, `integrateGeodesicKerrD`): Higher precision for validation and near-critical rays

The renderer uses FP64 by default (`BH_FP64_GEODESIC=1`) for accuracy.

### Analytic Metric Derivatives

Rather than using finite differences, `bh_kerr.h` provides closed-form expressions for all partial derivatives of the inverse metric. This:
- Eliminates finite-difference errors
- Reduces numerical noise in momentum evolution
- Enables more aggressive step sizes

### Backward Ray Tracing

Rays are traced backward from the camera into the scene (standard in ray tracing). This naturally handles:
- Gravitational lensing (multiple images map to same sky position)
- Disk emission (intersections detected during integration)
- Shadow computation (rays that never escape = black hole silhouette)

### Logarithmic Impact Parameter Sampling

The precomputation kernel uses logarithmic spacing for impact parameters, concentrating resolution near the critical value b_crit where deflection diverges. This provides better accuracy for strong-field lensing without excessive samples in the weak-field regime.

## File Reference

| File | Purpose |
|------|---------|
| `include/bh_common.h` | Constants, utilities, macros |
| `include/bh_types.h` | Core data structures |
| `include/bh_schwarzschild.h` | Schwarzschild geodesics |
| `include/bh_kerr.h` | Kerr geodesics, camera models |
| `schwarzschild_blackhole.cu` | Precomputation demo |
| `render.cu` | Full black hole renderer |
| `tests/bh_sanity_test.cu` | Geodesic validation tests |
| `tests/test_kerr_derivs.cu` | Kerr derivative verification |

## References

- Misner, Thorne, Wheeler (1973): Gravitation (Hamiltonian geodesics)
- Chandrasekhar (1983): The Mathematical Theory of Black Holes
- Bardeen, Press, Teukolsky (1972): Rotating Black Holes (LNRF, circular orbits)
- Luminet (1979): Image of a spherical black hole with thin accretion disk
- Dexter & Agol (2009): A Fast New Public Code for Computing Photon Orbits in a Kerr Spacetime
