#!/usr/bin/env python3
"""
Generate comprehensive test cases for geokerr validation
Creates systematic parameter sweeps covering critical regions
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import json

def estimate_kerr_critical_b(a):
    """
    Estimate critical impact parameter for Kerr black hole
    Based on approximate formula from literature
    """
    if a == 0:
        return 3.0 * np.sqrt(3) / 2  # Schwarzschild: 3√3/2 ≈ 2.598
    else:
        # Approximate formula for Kerr (Bardeen et al. 1972)
        # This is a rough estimate - exact value requires solving quartic
        Rs = 2.0  # M = 1 in geometric units
        return Rs * (2.598 - 0.5 * a + 0.1 * a**2)

def generate_elliptic_integral_tests(num_tests=5000):
    """
    Generate test cases for Carlson elliptic integrals RF, RD, RJ, RC
    """
    print(f"Generating {num_tests} elliptic integral test cases...")
    
    tests = []
    
    # RF(x,y,z) tests - requires x,y,z ≥ 0 with at most one zero
    for i in range(num_tests // 4):
        if i < 100:
            # Edge cases: one parameter near zero
            x = np.random.uniform(1e-12, 1e-8) if i % 3 == 0 else np.random.uniform(0.1, 10)
            y = np.random.uniform(1e-12, 1e-8) if i % 3 == 1 else np.random.uniform(0.1, 10)  
            z = np.random.uniform(1e-12, 1e-8) if i % 3 == 2 else np.random.uniform(0.1, 10)
        else:
            # General cases
            x = np.random.lognormal(0, 2)  # Log-normal for wide dynamic range
            y = np.random.lognormal(0, 2)
            z = np.random.lognormal(0, 2)
        
        tests.append({
            'type': 'RF',
            'x': x, 'y': y, 'z': z, 'p': 0.0,
            'expected': 0.0  # To be filled by FORTRAN
        })
    
    # RD(x,y,z) tests - requires x,y ≥ 0, z > 0, with at most one of x,y zero
    for i in range(num_tests // 4):
        x = np.random.lognormal(0, 2) if i % 2 == 0 else np.random.uniform(1e-12, 1e-8)
        y = np.random.lognormal(0, 2) if x > 1e-6 else np.random.lognormal(0, 2)
        z = np.random.lognormal(0, 2) + 1e-6  # Ensure z > 0
        
        tests.append({
            'type': 'RD',
            'x': x, 'y': y, 'z': z, 'p': 0.0,
            'expected': 0.0
        })
    
    # RJ(x,y,z,p) tests - more complex constraints
    for i in range(num_tests // 4):
        x = np.random.lognormal(0, 1) + 1e-6
        y = np.random.lognormal(0, 1) + 1e-6
        z = np.random.lognormal(0, 1) + 1e-6
        p = np.random.lognormal(0, 1) + 1e-6
        
        tests.append({
            'type': 'RJ', 
            'x': x, 'y': y, 'z': z, 'p': p,
            'expected': 0.0
        })
    
    # RC(x,y) tests - requires x ≥ 0, y ≠ 0
    for i in range(num_tests // 4):
        x = np.random.lognormal(0, 2)
        y = np.random.lognormal(0, 2) * (1 if np.random.rand() > 0.5 else -1)  # Can be negative
        
        tests.append({
            'type': 'RC',
            'x': x, 'y': y, 'z': 0.0, 'p': 0.0,
            'expected': 0.0
        })
    
    return tests

def generate_geodesic_tests(num_tests=2000):
    """
    Generate test cases for semi-analytic geodesic computation
    """
    print(f"Generating {num_tests} geodesic test cases...")
    
    tests = []
    
    # Kerr spin parameter sweep: a/M ∈ [0, 0.999]
    a_values = np.concatenate([
        [0.0],                          # Schwarzschild (critical test)
        np.linspace(0.01, 0.1, 10),    # Slow rotation  
        np.linspace(0.1, 0.9, 40),     # Moderate rotation
        np.linspace(0.9, 0.999, 30)    # Near-extremal (critical)
    ])
    
    for i, a in enumerate(a_values[:num_tests // len(a_values)]):
        M = 1.0  # Geometric units
        Rs = 2.0 * M
        
        # For each spin, test multiple impact parameters
        b_crit = estimate_kerr_critical_b(a)
        
        # Focus on critical region ±10% around b_crit
        b_values = np.concatenate([
            np.linspace(b_crit - 0.1 * b_crit, b_crit + 0.1 * b_crit, 20),
            np.logspace(np.log10(b_crit + 0.1 * b_crit), 2, 10)  # Large impact parameters
        ])
        
        for b in b_values[:num_tests // len(a_values) // len(b_values) + 1]:
            # Starting conditions: from "infinity"
            r0 = 1000.0 * Rs  # Start from large radius
            theta0 = np.pi / 2  # Equatorial plane (most important)
            phi0 = 0.0
            t0 = 0.0
            
            # Initial momentum for photon (null geodesic)  
            E = 1.0          # Energy at infinity
            L = b * E        # Angular momentum
            pt0 = -E         # p_t = -E (conserved)
            pphi0 = L        # p_φ = L (conserved)
            ptheta0 = 0.0    # Equatorial motion
            
            # Radial momentum from null condition (will be computed in FORTRAN)
            pr0 = 0.0  # Placeholder
            
            tests.append({
                'type': 'GEODESIC',
                'a': a, 'M': M, 'E': E, 'L': L, 'b': b,
                'r0': r0, 'theta0': theta0, 'phi0': phi0, 't0': t0,
                'pt0': pt0, 'pr0': pr0, 'ptheta0': ptheta0, 'pphi0': pphi0,
                'r_final': 0.0, 'theta_final': 0.0, 'phi_final': 0.0, 't_final': 0.0,
                'lambda_total': 0.0, 'ray_status': 0  # To be filled by FORTRAN
            })
    
    return tests

def save_test_cases(elliptic_tests, geodesic_tests, output_dir):
    """
    Save test cases in multiple formats for easy processing
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)
    
    # Save as JSON for Python processing
    all_tests = {
        'elliptic_tests': elliptic_tests,
        'geodesic_tests': geodesic_tests,
        'metadata': {
            'num_elliptic': len(elliptic_tests),
            'num_geodesic': len(geodesic_tests),
            'generated_by': 'generate_test_cases.py'
        }
    }
    
    with open(output_dir / 'test_cases.json', 'w') as f:
        json.dump(all_tests, f, indent=2)
    
    # Save elliptic integrals as FORTRAN-readable format
    with open(output_dir / 'elliptic_inputs.dat', 'w') as f:
        f.write(f"{len(elliptic_tests)}\n")
        for test in elliptic_tests:
            f.write(f"{test['type']:<2} {test['x']:20.12e} {test['y']:20.12e} "
                   f"{test['z']:20.12e} {test['p']:20.12e}\n")
    
    # Save geodesic tests as FORTRAN-readable format  
    with open(output_dir / 'geodesic_inputs.dat', 'w') as f:
        f.write(f"{len(geodesic_tests)}\n")
        for test in geodesic_tests:
            f.write(f"{test['a']:15.8e} {test['M']:15.8e} {test['b']:15.8e} "
                   f"{test['r0']:15.8e} {test['theta0']:15.8e} {test['phi0']:15.8e}\n")
    
    print(f"Test cases saved to {output_dir}")
    print(f"  - {len(elliptic_tests)} elliptic integral tests")
    print(f"  - {len(geodesic_tests)} geodesic tests")
    
    return output_dir / 'test_cases.json'

def plot_parameter_coverage(geodesic_tests, output_dir):
    """
    Create visualization of parameter space coverage
    """
    output_dir = Path(output_dir)
    
    # Extract parameters for plotting
    a_vals = [test['a'] for test in geodesic_tests]
    b_vals = [test['b'] for test in geodesic_tests]
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # Spin parameter distribution
    ax1.hist(a_vals, bins=50, alpha=0.7, edgecolor='black')
    ax1.axvline(0, color='red', linestyle='--', label='Schwarzschild')
    ax1.axvline(0.998, color='red', linestyle='--', label='Near-extremal')
    ax1.set_xlabel('Kerr spin parameter a/M')
    ax1.set_ylabel('Number of tests')
    ax1.set_title('Spin Parameter Coverage')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Impact parameter vs spin
    ax2.scatter(a_vals, b_vals, alpha=0.5, s=1)
    ax2.set_xlabel('Kerr spin parameter a/M')
    ax2.set_ylabel('Impact parameter b/Rs')
    ax2.set_title('Parameter Space Coverage')
    ax2.set_yscale('log')
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_dir / 'parameter_coverage.png', dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"Parameter coverage plot saved to {output_dir}/parameter_coverage.png")

def main():
    """
    Main function to generate all test cases
    """
    print("GEOKERR VALIDATION TEST CASE GENERATOR")
    print("=====================================")
    
    # Generate test cases
    elliptic_tests = generate_elliptic_integral_tests(5000)
    geodesic_tests = generate_geodesic_tests(2000)
    
    # Save test cases
    output_dir = Path(__file__).parent / 'test_data'
    test_file = save_test_cases(elliptic_tests, geodesic_tests, output_dir)
    
    # Create visualization
    plot_parameter_coverage(geodesic_tests, output_dir)
    
    print(f"\\nTest generation complete!")
    print(f"Next steps:")
    print(f"  1. Run FORTRAN geokerr on these inputs to generate reference data")
    print(f"  2. Implement CUDA Carlson elliptic integrals") 
    print(f"  3. Implement CUDA semi-analytic geodesic solver")
    print(f"  4. Validate CUDA against FORTRAN reference")
    
    return test_file

if __name__ == '__main__':
    main()