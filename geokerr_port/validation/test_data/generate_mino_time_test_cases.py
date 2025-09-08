#!/usr/bin/env python3
"""
Generate comprehensive test cases for Mino time integration validation
Focuses on GEOPHITIME subroutine parameters with extensive edge case coverage

Based on analysis of geokerr_original.f GEOPHITIME subroutine:
- Input: u0, uf, mu0, muf, a, l, l2, q2, tpm, tpr, su, sm
- Key: Mino time parameterization (lambda) is the fundamental output
- Critical: Need cases covering all orbit classifications (NCASE 1-6)
"""

import numpy as np
import json
import sys
from itertools import product

def generate_mino_time_test_cases():
    """
    Generate comprehensive test cases for Mino time integration validation.
    
    Key Parameters from GEOPHITIME:
    - u0, uf: Initial and final inverse radius (1/r)  
    - mu0, muf: Initial and final cos(theta)
    - a: Black hole spin [0, 1)
    - l: Angular momentum 
    - q2: Carter's constant
    - tpm, tpr: Turning point counts
    - su, sm: Initial velocities in Mino time (du/dlambda, dmu/dlambda)
    
    Edge cases critical for numerical stability:
    - Near-horizon (u large)
    - Polar regions (mu → ±1) 
    - High spin (a → 1)
    - Critical orbits (separatrix)
    - All NCASE classifications (1-6)
    """
    
    test_cases = []
    case_id = 1
    
    print("Generating Mino time integration test cases...")
    
    # =================================================================
    # 1. BASIC PARAMETER RANGES (Systematic Coverage)
    # =================================================================
    
    # Black hole spin: Include critical values
    spin_values = [
        0.0,           # Schwarzschild
        0.1, 0.3, 0.5, # Moderate spin
        0.7, 0.9,      # High spin  
        0.95, 0.99,    # Near-extremal (numerical challenges)
        0.999, 0.9999  # Extreme spin (stability test)
    ]
    
    # Radial coordinates (inverse radius u = M/r)
    u_values = [
        # Far field
        1e-6, 1e-5, 1e-4, 1e-3,
        # Intermediate  
        0.01, 0.05, 0.1, 0.2,
        # Near black hole
        0.3, 0.4, 0.45, 0.49,
        # Very close to horizon (numerical challenge)
        0.495, 0.499, 0.4999
    ]
    
    # Polar angles (mu = cos(theta))
    mu_values = sorted(list(set([0.0, 0.1, -0.1, 0.3, -0.3, 0.5, -0.5, 
                                0.7, -0.7, 0.9, -0.9, 0.95, -0.95, 
                                0.99, -0.99, 0.995, -0.995, 0.999, -0.999,
                                0.9999, -0.9999, 0.99999, -0.99999])))
    
    # Angular momentum (dimensionless)
    l_values = [
        # Prograde
        0.0, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0,
        # High angular momentum
        10.0, 20.0, 50.0,
        # Retrograde  
        -0.5, -1.0, -2.0, -5.0
    ]
    
    # Carter's constant Q^2 (dimensionless)
    q2_values = [
        # Equatorial (Q^2 = 0)
        0.0,
        # Small Q^2
        0.01, 0.1, 0.5, 1.0,
        # Medium Q^2  
        2.0, 5.0, 10.0,
        # Large Q^2 (highly inclined orbits)
        20.0, 50.0, 100.0, 200.0
    ]
    
    # =================================================================
    # 2. SYSTEMATIC PARAMETER COMBINATIONS
    # =================================================================
    
    print(f"Generating systematic combinations...")
    
    # Sample subset for systematic coverage (prevent combinatorial explosion)
    spin_sample = [0.0, 0.5, 0.9, 0.99]
    u_sample = [1e-4, 0.1, 0.3, 0.49] 
    mu_sample = [0.0, 0.5, 0.9, 0.99]
    l_sample = [0.0, 1.0, 3.0, -1.0]
    q2_sample = [0.0, 1.0, 10.0]
    
    for a, u0, uf, mu0, muf, l, q2 in product(
        spin_sample, u_sample, u_sample, mu_sample, mu_sample, l_sample, q2_sample
    ):
        if u0 != uf or mu0 != muf:  # Ensure non-trivial geodesic
            
            # Calculate derived parameters
            l2 = l * l
            
            # Initial velocities (signs determine direction)
            # su = ±1 for radial motion, sm for polar motion
            for su in [1, -1]:
                for sm in [0, 1, -1]:  # 0 for equatorial orbits
                    
                    # Estimate turning points (simplified)
                    tpr = 0 if abs(uf - u0) < 0.1 else 1
                    tpm = 0 if abs(muf - mu0) < 0.1 else (1 if abs(muf - mu0) < 0.5 else 2)
                    
                    test_case = {
                        'case_id': case_id,
                        'type': 'systematic',
                        'u0': u0,
                        'uf': uf, 
                        'mu0': mu0,
                        'muf': muf,
                        'a': a,
                        'l': l,
                        'l2': l2,
                        'q2': q2,
                        'tpm': tpm,
                        'tpr': tpr,
                        'su': su,
                        'sm': sm,
                        'description': f'Systematic: a={a:.3f}, u=({u0:.4f},{uf:.4f}), mu=({mu0:.3f},{muf:.3f}), l={l:.1f}, q2={q2:.1f}'
                    }
                    
                    test_cases.append(test_case)
                    case_id += 1
                    
                    if case_id > 5000:  # Prevent excessive cases
                        break
                if case_id > 5000:
                    break
            if case_id > 5000:
                break
    
    print(f"Generated {len(test_cases)} systematic cases")
    
    # =================================================================
    # 3. CRITICAL EDGE CASES (High Priority for Stability)
    # =================================================================
    
    edge_cases = [
        
        # ========================
        # HORIZON PROXIMITY TESTS
        # ========================
        {
            'type': 'horizon_approach',
            'a': 0.0, 'u0': 0.001, 'uf': 0.499, 'mu0': 0.0, 'muf': 0.0,
            'l': 2.0, 'q2': 0.0, 'su': 1, 'sm': 0, 'tpr': 1, 'tpm': 0,
            'description': 'Schwarzschild horizon approach (equatorial)'
        },
        {
            'type': 'horizon_approach', 
            'a': 0.9, 'u0': 0.001, 'uf': 0.499, 'mu0': 0.0, 'muf': 0.0,
            'l': 2.0, 'q2': 0.0, 'su': 1, 'sm': 0, 'tpr': 1, 'tpm': 0,
            'description': 'Kerr horizon approach (equatorial, high spin)'
        },
        {
            'type': 'horizon_approach',
            'a': 0.99, 'u0': 0.01, 'uf': 0.4999, 'mu0': 0.5, 'muf': 0.7,  
            'l': 3.0, 'q2': 2.0, 'su': 1, 'sm': 1, 'tpr': 1, 'tpm': 1,
            'description': 'Near-extremal Kerr horizon (inclined orbit)'
        },
        
        # ========================  
        # POLAR REGION TESTS
        # ========================
        {
            'type': 'polar_region',
            'a': 0.5, 'u0': 0.1, 'uf': 0.2, 'mu0': 0.999, 'muf': 0.99,
            'l': 1.0, 'q2': 10.0, 'su': 1, 'sm': -1, 'tpr': 0, 'tpm': 1,
            'description': 'Near north pole (high Carter constant)'  
        },
        {
            'type': 'polar_region',
            'a': 0.8, 'u0': 0.05, 'uf': 0.3, 'mu0': -0.9999, 'muf': 0.9999,
            'l': 4.0, 'q2': 50.0, 'su': 1, 'sm': 1, 'tpr': 1, 'tpm': 2,
            'description': 'Polar crossing (very high inclination)'
        },
        
        # ========================
        # CRITICAL ORBIT TYPES  
        # ========================
        {
            'type': 'circular_orbit',
            'a': 0.0, 'u0': 0.167, 'uf': 0.167, 'mu0': 0.0, 'muf': 0.0,
            'l': 3.464, 'q2': 0.0, 'su': 0, 'sm': 0, 'tpr': 0, 'tpm': 0,
            'description': 'Schwarzschild circular orbit (r=6M)'
        },
        {
            'type': 'plunging_orbit',
            'a': 0.0, 'u0': 0.1, 'uf': 0.45, 'mu0': 0.1, 'muf': 0.2,
            'l': 1.0, 'q2': 0.5, 'su': 1, 'sm': 1, 'tpr': 1, 'tpm': 1, 
            'description': 'Plunging orbit (insufficient angular momentum)'
        },
        {
            'type': 'bound_orbit',
            'a': 0.7, 'u0': 0.05, 'uf': 0.15, 'mu0': 0.3, 'muf': 0.7,
            'l': 5.0, 'q2': 8.0, 'su': 1, 'sm': 1, 'tpr': 1, 'tpm': 2,
            'description': 'Bound orbit with turning points'
        },
        
        # ========================
        # RETROGRADE ORBITS
        # ======================== 
        {
            'type': 'retrograde',
            'a': 0.9, 'u0': 0.01, 'uf': 0.2, 'mu0': 0.0, 'muf': 0.5,
            'l': -2.0, 'q2': 1.0, 'su': 1, 'sm': 1, 'tpr': 0, 'tpm': 1,
            'description': 'Retrograde orbit (high spin)'
        },
        
        # ========================
        # EXTREME PARAMETERS
        # ========================
        {
            'type': 'extreme_spin',
            'a': 0.9999, 'u0': 0.001, 'uf': 0.1, 'mu0': 0.1, 'muf': 0.5,
            'l': 2.0, 'q2': 1.0, 'su': 1, 'sm': 1, 'tpr': 0, 'tpm': 1,
            'description': 'Extreme spin (numerical stability test)'
        },
        {
            'type': 'high_energy',
            'a': 0.5, 'u0': 1e-5, 'uf': 0.01, 'mu0': 0.0, 'muf': 0.3,
            'l': 50.0, 'q2': 100.0, 'su': 1, 'sm': 1, 'tpr': 0, 'tpm': 1,
            'description': 'High energy orbit (large conserved quantities)'
        }
    ]
    
    # Add derived parameters to edge cases
    for i, case in enumerate(edge_cases):
        case['case_id'] = case_id + i
        case['l2'] = case['l'] ** 2
        if 'description' not in case:
            case['description'] = f"Edge case {case['type']}"
    
    test_cases.extend(edge_cases)
    print(f"Added {len(edge_cases)} critical edge cases")
    
    # =================================================================
    # 4. ORBIT CLASSIFICATION TESTS (NCASE 1-6 Coverage)
    # =================================================================
    
    # Based on Dexter & Agol Table 1, different NCASE values correspond to
    # different root structures of the radial potential
    orbit_classification_cases = [
        {
            'type': 'ncase_1', 'a': 0.0, 'u0': 0.001, 'uf': 0.1,
            'mu0': 0.0, 'muf': 0.0, 'l': 10.0, 'q2': 0.0,
            'su': 1, 'sm': 0, 'tpr': 0, 'tpm': 0,
            'description': 'NCASE=1 (escape orbit, no turning points)'
        },
        {
            'type': 'ncase_2', 'a': 0.5, 'u0': 0.1, 'uf': 0.3, 
            'mu0': 0.2, 'muf': 0.6, 'l': 3.0, 'q2': 2.0,
            'su': 1, 'sm': 1, 'tpr': 1, 'tpm': 1,
            'description': 'NCASE=2 (bound radial motion)'
        },
        {
            'type': 'ncase_3', 'a': 0.8, 'u0': 0.05, 'uf': 0.4,
            'mu0': 0.1, 'muf': 0.9, 'l': 2.0, 'q2': 10.0,
            'su': 1, 'sm': 1, 'tpr': 1, 'tpm': 2,
            'description': 'NCASE=3 (complex root structure)'
        }
    ]
    
    for i, case in enumerate(orbit_classification_cases):
        case['case_id'] = case_id + len(edge_cases) + i
        case['l2'] = case['l'] ** 2
    
    test_cases.extend(orbit_classification_cases) 
    print(f"Added {len(orbit_classification_cases)} orbit classification cases")
    
    # =================================================================
    # 5. NUMERICAL PRECISION TESTS
    # =================================================================
    
    precision_cases = [
        {
            'type': 'tiny_difference',
            'a': 0.5, 'u0': 0.1, 'uf': 0.1 + 1e-10, 
            'mu0': 0.3, 'muf': 0.3 + 1e-10,
            'l': 2.0, 'q2': 1.0, 'su': 1, 'sm': 1, 'tpr': 0, 'tpm': 0,
            'description': 'Tiny coordinate differences (precision test)'
        },
        {
            'type': 'large_numbers',
            'a': 0.1, 'u0': 1e-8, 'uf': 1e-6,
            'mu0': 1e-6, 'muf': 1e-4,
            'l': 1000.0, 'q2': 10000.0, 'su': 1, 'sm': 1, 'tpr': 0, 'tpm': 0,
            'description': 'Large parameter values (overflow test)'
        }
    ]
    
    for i, case in enumerate(precision_cases):
        case['case_id'] = case_id + len(edge_cases) + len(orbit_classification_cases) + i
        case['l2'] = case['l'] ** 2
    
    test_cases.extend(precision_cases)
    print(f"Added {len(precision_cases)} precision test cases")
    
    print(f"\nTotal test cases generated: {len(test_cases)}")
    
    return test_cases

def write_test_cases(test_cases, output_format='both'):
    """Write test cases in FORTRAN and JSON formats"""
    
    if output_format in ['fortran', 'both']:
        print("Writing FORTRAN input file...")
        with open('mino_time_inputs.dat', 'w') as f:
            f.write(f"# Mino time integration test cases\n")
            f.write(f"# Format: u0 uf mu0 muf a l l2 q2 tpm tpr su sm\n")
            f.write(f"{len(test_cases)}\n")
            
            for case in test_cases:
                f.write(f"{case['u0']:.12e} {case['uf']:.12e} ")
                f.write(f"{case['mu0']:.12e} {case['muf']:.12e} ")  
                f.write(f"{case['a']:.12e} {case['l']:.12e} {case['l2']:.12e} {case['q2']:.12e} ")
                f.write(f"{case['tpm']} {case['tpr']} {case['su']} {case['sm']}\n")
        
        print(f"Written {len(test_cases)} cases to mino_time_inputs.dat")
    
    if output_format in ['json', 'both']:
        print("Writing JSON metadata file...")
        with open('mino_time_test_cases.json', 'w') as f:
            json.dump({
                'description': 'Mino time integration test cases for GEOKERR validation',
                'total_cases': len(test_cases),
                'format_version': '1.0',
                'test_cases': test_cases
            }, f, indent=2)
        
        print(f"Written test case metadata to mino_time_test_cases.json")

def analyze_test_coverage(test_cases):
    """Analyze parameter space coverage"""
    
    print("\n" + "="*60)
    print("TEST CASE COVERAGE ANALYSIS")
    print("="*60)
    
    # Parameter ranges
    params = ['a', 'u0', 'uf', 'mu0', 'muf', 'l', 'q2']
    for param in params:
        values = [case[param] for case in test_cases]
        print(f"{param:>4}: [{min(values):>8.4f}, {max(values):>8.4f}] "
              f"(n={len(set(values)):>3})")
    
    # Test types
    types = {}
    for case in test_cases:
        t = case.get('type', 'systematic')
        types[t] = types.get(t, 0) + 1
    
    print(f"\nTest types:")
    for t, count in sorted(types.items()):
        print(f"  {t:>20}: {count:>4} cases")
    
    # Critical coverage
    critical_checks = {
        'Near horizon (u > 0.4)': len([c for c in test_cases if max(c['u0'], c['uf']) > 0.4]),
        'High spin (a > 0.9)': len([c for c in test_cases if c['a'] > 0.9]),
        'Polar regions (|mu| > 0.9)': len([c for c in test_cases if max(abs(c['mu0']), abs(c['muf'])) > 0.9]),
        'Retrograde (l < 0)': len([c for c in test_cases if c['l'] < 0]),
        'High Carter constant (q2 > 10)': len([c for c in test_cases if c['q2'] > 10])
    }
    
    print(f"\nCritical region coverage:")
    for check, count in critical_checks.items():
        print(f"  {check:>30}: {count:>4} cases")
    
    print(f"\nRecommended validation targets:")
    print(f"  - Elliptic integral accuracy: >99% pass at 1e-12 tolerance") 
    print(f"  - Mino time lambda computation: exact agreement with FORTRAN")
    print(f"  - All orbit types (NCASE 1-6): complete coverage")
    print(f"  - Edge cases: robust handling without NaN/overflow")

if __name__ == "__main__":
    print("GEOKERR Mino Time Integration Test Case Generator")
    print("="*55)
    
    # Generate comprehensive test cases
    test_cases = generate_mino_time_test_cases()
    
    # Write output files
    write_test_cases(test_cases, 'both')
    
    # Coverage analysis
    analyze_test_coverage(test_cases)
    
    print(f"\nNext steps:")
    print(f"1. Create FORTRAN batch tester for GEOPHITIME subroutine")
    print(f"2. Generate reference data: ./build/mino_time_batch")
    print(f"3. Implement CUDA Mino time integration")
    print(f"4. Validate: exact numerical agreement required")
    print(f"\nFocus: Mino time parameter λ is the key output to validate!")