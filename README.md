# NullPath — CUDA/C++ Null‑Geodesic Tracer and Renderer

NullPath is a CUDA/C++ prototype that integrates photon paths (null geodesics) in curved spacetime and renders black‑hole imagery. It includes an RK4 integrator, precomputed deflection/bending tables, a thin‑disk renderer with progressive accumulation, and sanity tests. Schwarzschild is the validated core; Kerr support is experimental (finite‑difference metric derivatives pending analytic forms).

## Quickstart
- Build renderer: `make render-fast`
- Render a still: `./bin/render --w 1280 --h 720 --samples 4 --out frame.ppm`
- Kerr preview: `./bin/render --spin 0.9 --w 1920 --h 1080 --samples 4 --spp-total 16 --checkpoint-every 4 --out kerr.ppm`
- Tests: `make quick-test` (honors `BH_NUM_RAYS`)

## Citations & References (with arXiv IDs where available)
- Bardeen, Press, Teukolsky (1972), “Rotating Black Holes: Locally Nonrotating Frames, Energy Extraction, and Scalar Synchrotron Radiation,” ApJ 178, 347. Classic LNRF and circular orbits in Kerr. (no arXiv; ADS bibcode 1972ApJ...178..347B)
- Chandrasekhar (1983), The Mathematical Theory of Black Holes, Oxford. Canonical text for Kerr/Schwarzschild geodesics. (no arXiv)
- Luminet (1979), “Image of a spherical black hole with thin accretion disk,” A&A 75, 228–235. Thin‑disk appearance; Doppler beaming. (no arXiv; ADS 1979A&A....75..228L)
- Misner, Thorne, Wheeler (1973), Gravitation, Freeman. Hamiltonian geodesics H = ½ g^{μν}p_μ p_ν; conserved quantities. (no arXiv)
- Lindquist (1966), “Relativistic Transport Theory,” Ann. Phys. 37, 487. Invariance of I_ν/ν³ ⇒ I_obs = g³ I_em. (no arXiv)
- Teo (2020), “Spherical photon orbits around a Kerr black hole,” Gen. Rel. Grav. 52, 5. arXiv:2007.04022 — practical formulas for Kerr photon orbits.
- Perlick & Tsupko (2021), “Calculating black hole shadows: Review,” Phys. Rep. 947, 1–39. arXiv:2105.07101 — comprehensive review of lensing/shadows in Kerr/Schwarzschild.
- Gralla & Lupsasca (2019), “Null geodesics of the Kerr exterior,” Phys. Rev. D 101, 044032. arXiv:1910.12881 — geodesic structure and photon rings; see also arXiv:1910.12873.
- Takahashi (2004), “Shapes and positions of black hole shadows...,” ApJ 611, 996. arXiv:astro-ph/0405099 — observational aspects tied to lensing/shadows.
- Kerr metric (Boyer–Lindquist): Σ=r²+a²cos²θ, Δ=r²−2Mr+a², A=(r²+a²)²−a²Δ sin²θ; g^{tt}=−A/(ΣΔ), g^{tφ}=−2aMr/(ΣΔ), g^{φφ}=(Δ−a² sin²θ)/(ΣΔ sin²θ), g^{rr}=Δ/Σ, g^{θθ}=1/Σ. We use these with analytic ∂g^{μν}/∂(r,θ) for the Hamiltonian integrator.
- Schwarzschild null turning point: for b≡L/E, r_min solves b² = r³/(r−2M); photon sphere r=3M, critical impact b_crit=3√3 M.

Notes
- Deflection storage: `deflection_angles` stores total Δφ; `bending_angles` = Δφ − π for escaped rays (conventional lensing bend).
- Units: geometric (G=c=1). Rs=2M, so r_ph=1.5 Rs and b_crit=(3√3/2)Rs.
- Camera frame: static observer outside the ergosphere; consider ZAMO for deep Kerr placements (see TODOs in AGENTS.md).
