/**
 * @file bh_types.h
 * @brief Core data structures for the NullPath geodesic tracer.
 *
 * This header defines the fundamental types used throughout NullPath:
 * - SchwarzschildMetric: Schwarzschild spacetime metric accessor
 * - GeodesicState: Phase space state for geodesic integration (FP32)
 * - GeodesicLookupTable: Precomputed geodesic data for fast lookups
 *
 * All structures are designed for use in both host and device code.
 * Header-only with inline functions.
 *
 * @note Coordinate system: Schwarzschild coordinates (t, r, theta, phi).
 * @note Metric signature: (-,+,+,+) with timelike ds^2 < 0.
 * @note Units: Geometric units (G = c = 1). Mass M sets the scale.
 *
 * @see bh_common.h for constants and utilities
 * @see bh_schwarzschild.h for geodesic integration
 * @see bh_kerr.h for Kerr spacetime structures (GeodesicStateD)
 */

#ifndef BH_TYPES_H
#define BH_TYPES_H

#include "bh_common.h"

/* ============================================================================
 * SCHWARZSCHILD METRIC
 * ============================================================================
 */

/**
 * @brief Schwarzschild metric accessor for a non-rotating black hole.
 *
 * Provides covariant metric components in Schwarzschild coordinates
 * (t, r, theta, phi). The line element is:
 *
 *   ds^2 = g_tt dt^2 + g_rr dr^2 + g_theta_theta dtheta^2 + g_phi_phi dphi^2
 *        = -(1 - Rs/r) dt^2 + (1 - Rs/r)^{-1} dr^2 + r^2 dtheta^2 + r^2 sin^2(theta) dphi^2
 *
 * @note Units: Geometric units where G = c = 1.
 * @note The metric is singular at r = Rs (coordinate singularity at event horizon)
 *       and r = 0 (true singularity). Functions do not protect against these.
 *
 * Example usage:
 * @code
 *   SchwarzschildMetric metric(10.0f);  // M = 10 in geometric units
 *   float Rs = metric.Rs;               // = 20.0
 *   float g_tt = metric.g_tt(30.0f);    // = -(1 - 20/30) = -1/3
 * @endcode
 */
struct SchwarzschildMetric {
    /**
     * @brief Black hole mass in geometric units.
     *
     * This is the mass parameter M appearing in the metric.
     * In physical units, M_physical = M * (G/c^2).
     */
    float mass;

    /**
     * @brief Schwarzschild radius Rs = 2M.
     *
     * The coordinate radius of the event horizon for a non-rotating black hole.
     * Precomputed for efficiency.
     */
    float Rs;

    /**
     * @brief Construct a Schwarzschild metric for a given mass.
     * @param m Black hole mass M in geometric units.
     *
     * Initializes mass = m and Rs = 2M.
     */
    __host__ __device__
    explicit SchwarzschildMetric(float m) : mass(m), Rs(BH_SCHWARZSCHILD_MULTIPLIER * m) {}

    /**
     * @brief Covariant g_tt component.
     *
     * g_tt = -(1 - Rs/r)
     *
     * @param r Radial coordinate (must be > Rs for physical region).
     * @return g_tt at radius r.
     *
     * @note Returns 0 at r = Rs (horizon), positive for r < Rs (interior).
     * @warning No protection against r = 0 or r <= Rs.
     */
    __host__ __device__ inline float g_tt(float r) const {
        return -(1.0f - Rs / r);
    }

    /**
     * @brief Covariant g_rr component.
     *
     * g_rr = 1 / (1 - Rs/r)
     *
     * @param r Radial coordinate (must be > Rs for physical region).
     * @return g_rr at radius r.
     *
     * @note Diverges as r -> Rs (coordinate singularity at horizon).
     * @note Protected with BH_EPS floor to avoid division by zero at horizon.
     */
    __host__ __device__ inline float g_rr(float r) const {
        return 1.0f / fmaxf(BH_EPS, 1.0f - Rs / r);
    }

    /**
     * @brief Covariant g_theta_theta component.
     *
     * g_theta_theta = r^2
     *
     * @param r Radial coordinate.
     * @return g_theta_theta at radius r.
     */
    __host__ __device__ inline float g_theta_theta(float r) const {
        return r * r;
    }

    /**
     * @brief Covariant g_phi_phi component.
     *
     * g_phi_phi = r^2 sin^2(theta)
     *
     * @param r Radial coordinate.
     * @param theta Polar angle in radians (0 = north pole, pi = south pole).
     * @return g_phi_phi at (r, theta).
     *
     * @note Returns 0 at theta = 0 or pi (coordinate singularity at poles).
     */
    __host__ __device__ inline float g_phi_phi(float r, float theta) const {
        float s = sinf(theta);
        return r * r * s * s;
    }
};

/* ============================================================================
 * GEODESIC STATE (FP32)
 * ============================================================================
 */

/**
 * @brief Phase space state for null geodesic integration (single precision).
 *
 * Represents the 8-dimensional phase space state of a photon:
 * - Position: (t, r, theta, phi) in Schwarzschild or Boyer-Lindquist coordinates
 * - Momentum: (p_t, p_r, p_theta, p_phi) as covariant 4-momentum components
 *
 * The Hamiltonian constraint for null geodesics is:
 *   H = (1/2) g^{mu nu} p_mu p_nu = 0
 *
 * Conserved quantities (for stationary, axisymmetric spacetimes):
 * - Energy at infinity: E = -p_t
 * - Angular momentum: L = p_phi
 * - Impact parameter: b = L/E = -p_phi/p_t
 *
 * @note Coordinates: (t, r, theta, phi) with theta in [0, pi] and phi in [0, 2*pi).
 * @note Momenta are covariant (lower indices): p_mu = g_{mu nu} dx^nu/d(lambda).
 * @note For ingoing rays (toward black hole), p_r < 0.
 *
 * @see GeodesicStateD in bh_kerr.h for double precision version.
 */
struct GeodesicState {
    /* --------------------------------------------------------------------
     * Position coordinates
     * -------------------------------------------------------------------- */

    /**
     * @brief Coordinate time t.
     *
     * In Schwarzschild/Boyer-Lindquist coordinates, t is the time measured
     * by a static observer at infinity. Not directly observable for photons
     * but tracked for consistency checking.
     */
    float t;

    /**
     * @brief Radial coordinate r.
     *
     * The areal radius in Schwarzschild coordinates. The event horizon is
     * at r = Rs = 2M. The photon sphere is at r = 3M = 1.5*Rs.
     */
    float r;

    /**
     * @brief Polar angle theta in radians.
     *
     * Range: [0, pi] where theta = 0 is the north pole (along the spin axis
     * for Kerr) and theta = pi is the south pole. The equatorial plane is
     * at theta = pi/2.
     */
    float theta;

    /**
     * @brief Azimuthal angle phi in radians.
     *
     * Range: [0, 2*pi) wrapping around the symmetry axis. For axisymmetric
     * spacetimes, the metric is independent of phi.
     */
    float phi;

    /* --------------------------------------------------------------------
     * Covariant 4-momentum components
     * -------------------------------------------------------------------- */

    /**
     * @brief Covariant time component p_t = g_{t mu} p^mu.
     *
     * For stationary spacetimes, p_t is conserved along geodesics.
     * Energy at infinity: E = -p_t > 0 for future-directed photons.
     */
    float pt;

    /**
     * @brief Covariant radial component p_r = g_{r mu} p^mu.
     *
     * Sign convention:
     * - p_r < 0: ingoing (decreasing r)
     * - p_r > 0: outgoing (increasing r)
     *
     * At a turning point (closest approach), p_r = 0.
     */
    float pr;

    /**
     * @brief Covariant polar component p_theta = g_{theta mu} p^mu.
     *
     * For equatorial orbits (theta = pi/2 constant), p_theta = 0.
     * Otherwise oscillates as the ray moves between theta bounds.
     */
    float ptheta;

    /**
     * @brief Covariant azimuthal component p_phi = g_{phi mu} p^mu.
     *
     * For axisymmetric spacetimes, p_phi is conserved along geodesics.
     * This is the angular momentum L about the symmetry axis.
     * - L > 0: prograde (orbiting in direction of black hole spin)
     * - L < 0: retrograde
     */
    float pphi;

    /**
     * @brief Default constructor initializes all components to zero.
     *
     * A zero state is not physically meaningful; use this for
     * allocation and then set proper initial conditions.
     */
    __host__ __device__ GeodesicState()
        : t(0), r(0), theta(0), phi(0), pt(0), pr(0), ptheta(0), pphi(0) {}
};

/* ============================================================================
 * GEODESIC LOOKUP TABLE
 * ============================================================================
 */

/**
 * @brief Precomputed geodesic data for fast impact parameter lookups.
 *
 * Stores results from integrating many geodesics with different impact
 * parameters, enabling fast interpolation of deflection angles, bending
 * angles, and ray classification without re-integrating.
 *
 * This structure is host-only (not used on GPU). The arrays are allocated
 * on the host and filled by copying results from device computations.
 *
 * @note Memory management: Caller is responsible for allocating and
 *       deallocating all pointer members. Use new[]/delete[].
 *
 * Example usage:
 * @code
 *   GeodesicLookupTable* table = precomputeGeodesics(mass, num_rays);
 *   // Use table for lookups...
 *   delete[] table->impact_parameters;
 *   delete[] table->deflection_angles;
 *   // ... delete other arrays ...
 *   delete table;
 * @endcode
 */
struct GeodesicLookupTable {
    /**
     * @brief Array of impact parameters b = L/E.
     *
     * Size: num_entries elements.
     * Typically logarithmically spaced to concentrate resolution near b_crit.
     * Units: geometric units (same as r, M).
     */
    float* impact_parameters;

    /**
     * @brief Array of total deflection angles.
     *
     * Size: num_entries elements.
     * The total change in phi from start to end of integration.
     * For escaped rays, this is the total angle swept.
     * Units: radians.
     *
     * @note For a ray that comes from infinity, turns, and escapes,
     *       deflection > pi indicates strong lensing.
     */
    float* deflection_angles;

    /**
     * @brief Array of bending angles (conventional lensing definition).
     *
     * Size: num_entries elements.
     * bending_angle = deflection_angle - pi for escaped rays.
     * This is the angle by which the ray's direction changes from
     * the incident to the outgoing asymptote.
     *
     * Set to NaN for captured rays (no outgoing asymptote).
     * Units: radians.
     */
    float* bending_angles;

    /**
     * @brief Array of closest approach radii.
     *
     * Size: num_entries elements.
     * The minimum r reached during integration. For escaped rays with
     * b > b_crit, this is the turning point radius.
     * For captured rays, this is typically close to Rs.
     * Units: geometric units.
     *
     * @note For escaped rays, an analytic formula may provide higher
     *       accuracy than the numerically tracked minimum.
     */
    float* closest_approaches;

    /**
     * @brief Array of ray status codes.
     *
     * Size: num_entries elements.
     * Integer codes matching RayStatus enum values:
     * - 0: Escaped
     * - 1: Captured
     * - 2: PhotonSphere
     * - 3: Inconclusive
     *
     * @note Uses int rather than RayStatus enum for CUDA compatibility.
     */
    int* ray_status;

    /**
     * @brief Number of entries in each array.
     */
    int num_entries;

    /**
     * @brief Black hole mass used for this table.
     *
     * Stored for reference when interpreting impact parameters
     * and closest approaches in physical units.
     */
    float mass;
};

#endif // BH_TYPES_H
