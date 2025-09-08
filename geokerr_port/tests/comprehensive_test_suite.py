#!/usr/bin/env python3
"""
Comprehensive Test Suite for GEOKERR Semi-Analytic Components

This suite generates systematic test cases covering:
1. Standard photon orbits across parameter space
2. Edge cases: near-extremal spins, critical impact parameters
3. Boundary conditions: horizon, turning points, caustics
4. Numerical challenges: small/large values, degeneracies
5. Physical validation: conservation laws, geodesic equation
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import json
import sys

class GeodeticTestGenerator:
    """Generate systematic test cases for geodesic components"""
    
    def __init__(self, output_dir="test_cases"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True, parents=True)
        
        # Physical constants and limits
        self.a_max = 0.9999  # Nearly extremal
        self.r_horizon_factor = 1.001  # Just outside horizon
        
    def generate_spin_parameter_sweep(self, n_points=50):
        """Generate comprehensive Kerr spin parameter test cases"""
        print(f"Generating {n_points} spin parameter test cases...")
        
        # Logarithmic distribution focusing on critical regions
        a_values = np.concatenate([
            [0.0],  # Schwarzschild (critical test)
            np.logspace(-4, -2, 10),  # Very small spins
            np.linspace(0.01, 0.1, 10),  # Small spins
            np.linspace(0.1, 0.8, 20),   # Moderate spins  
            np.linspace(0.8, 0.98, 15),  # High spins
            np.linspace(0.98, 0.9999, 10)  # Near-extremal (critical)
        ])[:n_points]
        
        return a_values
    
    def generate_impact_parameter_sweep(self, a):
        """Generate impact parameters covering all orbit types"""
        
        # Estimate critical impact parameters for given spin
        r_plus = 1.0 + np.sqrt(1.0 - a*a)  # Outer horizon
        
        # Approximate critical impact parameters (Bardeen et al. 1972)
        if a == 0:
            b_crit = 3.0 * np.sqrt(3) / 2  # Schwarzschild photon sphere
            b_unstable = b_crit
        else:
            # Approximate formulas for Kerr critical orbits
            b_crit = 2.6 - 0.4*a + 0.1*a*a  # Rough approximation
            b_unstable = b_crit * (1 + 0.1*a)
        
        # Systematic impact parameter coverage
        b_values = np.concatenate([
            # Direct capture (small b)
            np.linspace(0.1, 0.8*b_crit, 10),
            
            # Critical region (near photon sphere)
            np.linspace(0.9*b_crit, 1.1*b_crit, 20),
            
            # Moderate scattering
            np.linspace(1.2*b_crit, 2.0*b_crit, 15),
            
            # Weak field (large b)  
            np.logspace(np.log10(2.0*b_crit), 2, 10)
        ])
        
        return b_values, b_crit, b_unstable
        
    def generate_observer_geometries(self, n_configs=20):
        """Generate diverse observer positions and orientations"""
        
        configs = []
        
        # Standard configurations
        for i in range(n_configs):
            config = {
                'r_obs': np.random.uniform(10, 1000),  # Observer distance
                'theta_obs': np.random.uniform(0.1, np.pi - 0.1),  # Inclination
                'phi_obs': np.random.uniform(0, 2*np.pi),  # Azimuth
                
                # Image plane parameters
                'alpha_range': np.random.uniform(5, 20),  # Impact parameter range
                'beta_range': np.random.uniform(5, 20),
                'image_size': np.random.choice([64, 128, 256])
            }
            configs.append(config)
            
        # Add special cases
        special_cases = [
            # Face-on (θ = π/2)
            {'r_obs': 100, 'theta_obs': np.pi/2, 'phi_obs': 0,
             'alpha_range': 10, 'beta_range': 10, 'image_size': 128},
             
            # Edge-on (θ ~ 0)
            {'r_obs': 100, 'theta_obs': 0.1, 'phi_obs': 0,
             'alpha_range': 10, 'beta_range': 10, 'image_size': 128},
             
            # Near black hole
            {'r_obs': 5, 'theta_obs': np.pi/3, 'phi_obs': 0,
             'alpha_range': 5, 'beta_range': 5, 'image_size': 64},
             
            # Distant observer  
            {'r_obs': 10000, 'theta_obs': np.pi/4, 'phi_obs': 0,
             'alpha_range': 20, 'beta_range': 20, 'image_size': 256}
        ]
        
        return configs + special_cases
    
    def generate_edge_cases(self):
        """Generate challenging edge cases for robust validation"""
        
        edge_cases = []
        
        # Near-extremal spin cases
        for a in [0.999, 0.9999, 0.99999]:
            for b in np.logspace(-2, 1, 10):  # Wide impact parameter range
                edge_cases.append({
                    'type': 'near_extremal',
                    'a': a,
                    'b': b,
                    'r0': 1000,
                    'theta0': np.pi/2,
                    'challenge': 'Near-extremal spin'
                })
        
        # Near-horizon initial conditions
        for a in [0.0, 0.5, 0.9]:
            r_horizon = 1.0 + np.sqrt(1.0 - a*a)
            for r_factor in [1.001, 1.01, 1.1]:  # Very close to horizon
                edge_cases.append({
                    'type': 'near_horizon',
                    'a': a,
                    'b': 4.0,
                    'r0': r_horizon * r_factor,
                    'theta0': np.pi/2,
                    'challenge': f'r0 = {r_factor:.3f} * r_horizon'
                })
        
        # Critical impact parameters (photon sphere region)
        for a in [0.0, 0.3, 0.6, 0.9]:
            b_crit = 3.0 * np.sqrt(3) / 2 if a == 0 else 2.6 - 0.4*a
            for b_factor in [0.99, 1.0, 1.01]:  # Near critical
                edge_cases.append({
                    'type': 'critical_impact',
                    'a': a, 
                    'b': b_crit * b_factor,
                    'r0': 100,
                    'theta0': np.pi/2,
                    'challenge': f'b = {b_factor:.3f} * b_critical'
                })
        
        # Polar orbits (θ₀ near 0 or π)
        for theta_factor in [0.001, 0.01, 0.1]:
            for sign in [1, -1]:
                theta0 = theta_factor if sign > 0 else np.pi - theta_factor
                edge_cases.append({
                    'type': 'polar_orbit',
                    'a': 0.5,
                    'b': 4.0,
                    'r0': 50,
                    'theta0': theta0,
                    'challenge': f'θ₀ = {theta0:.4f} (near pole)'
                })
        
        # Large impact parameters (weak field)
        for b in [50, 100, 500, 1000]:
            edge_cases.append({
                'type': 'weak_field',
                'a': 0.9,
                'b': b,
                'r0': 1000,
                'theta0': np.pi/2,
                'challenge': f'Large impact parameter b = {b}'
            })
            
        # Small impact parameters (strong field)
        for b in [0.1, 0.5, 1.0]:
            edge_cases.append({
                'type': 'strong_field',
                'a': 0.8,
                'b': b,
                'r0': 100,
                'theta0': np.pi/2,
                'challenge': f'Small impact parameter b = {b}'
            })
        
        return edge_cases
    
    def generate_conservation_tests(self, n_tests=100):
        """Generate tests specifically for conservation law validation"""
        
        tests = []
        
        for i in range(n_tests):
            # Random parameters
            a = np.random.uniform(0, 0.95)
            b = np.random.uniform(1, 20)
            r0 = np.random.uniform(10, 1000)
            theta0 = np.random.uniform(0.1, np.pi - 0.1)
            
            # Compute initial conserved quantities
            E = 1.0  # Photons
            L = b * E
            cos_theta0 = np.cos(theta0)
            sin_theta0 = np.sin(theta0)
            
            # Carter constant for equatorial motion (Q = 0)
            # For general motion: Q = β² + α²cos²θ₀ - a²cos²θ₀
            Q = 0.0  # Start with equatorial
            
            tests.append({
                'a': a, 'b': b, 'r0': r0, 'theta0': theta0,
                'E': E, 'L': L, 'Q': Q,
                'test_type': 'conservation',
                'expected_conserved': [E, L, Q]
            })
            
        return tests
    
    def generate_boundary_tests(self):
        """Generate tests for boundary conditions and special cases"""
        
        boundary_tests = []
        
        # Event horizon approach
        for a in [0.0, 0.5, 0.99]:
            r_horizon = 1.0 + np.sqrt(1.0 - a*a)
            boundary_tests.append({
                'type': 'horizon_approach',
                'a': a,
                'r_target': r_horizon * 1.0001,
                'expect_behavior': 'coordinate_time_divergence',
                'validation': 'Check t → ∞ as r → r₊'
            })
        
        # Turning points
        for orbit_type in ['bound', 'plunging', 'scattering']:
            boundary_tests.append({
                'type': 'turning_point',
                'orbit_type': orbit_type,
                'validation': 'Check dr/dλ = 0 at turning points'
            })
            
        # Infinity approach  
        boundary_tests.append({
            'type': 'infinity_approach',
            'r_start': 1000,
            'expect_behavior': 'flat_spacetime_limit',
            'validation': 'Check Minkowski limit at large r'
        })
        
        return boundary_tests
    
    def save_test_suite(self):
        """Generate and save complete test suite"""
        
        print("COMPREHENSIVE GEOKERR TEST SUITE GENERATOR")
        print("=" * 50)
        
        # Generate all test categories
        spin_sweep = self.generate_spin_parameter_sweep(50)
        observer_configs = self.generate_observer_geometries(20)  
        edge_cases = self.generate_edge_cases()
        conservation_tests = self.generate_conservation_tests(100)
        boundary_tests = self.generate_boundary_tests()
        
        # Systematic geodesic test matrix
        systematic_tests = []
        for a in spin_sweep[:10]:  # Sample of spins
            b_values, b_crit, b_unstable = self.generate_impact_parameter_sweep(a)
            for b in b_values[:5]:  # Sample of impact parameters
                systematic_tests.append({
                    'a': a, 'b': b, 'r0': 100, 'theta0': np.pi/2,
                    'b_critical': b_crit, 'b_unstable': b_unstable,
                    'E': 1.0, 'L': b, 'Q': 0.0
                })
        
        # Compile master test suite
        test_suite = {
            'metadata': {
                'generator': 'comprehensive_test_suite.py',
                'total_tests': len(systematic_tests) + len(edge_cases) + 
                             len(conservation_tests) + len(boundary_tests),
                'categories': {
                    'systematic': len(systematic_tests),
                    'edge_cases': len(edge_cases),
                    'conservation': len(conservation_tests),
                    'boundary': len(boundary_tests)
                }
            },
            'systematic_tests': systematic_tests,
            'edge_cases': edge_cases,
            'conservation_tests': conservation_tests,
            'boundary_tests': boundary_tests,
            'observer_configurations': observer_configs,
            'spin_parameter_sweep': spin_sweep.tolist()
        }
        
        # Convert numpy types to native Python types for JSON serialization
        def convert_numpy_types(obj):
            if isinstance(obj, np.ndarray):
                return obj.tolist()
            elif isinstance(obj, (np.int64, np.int32, np.int_)):
                return int(obj)
            elif isinstance(obj, (np.float64, np.float32, np.float_)):
                return float(obj)
            elif isinstance(obj, dict):
                return {key: convert_numpy_types(value) for key, value in obj.items()}
            elif isinstance(obj, list):
                return [convert_numpy_types(item) for item in obj]
            else:
                return obj
        
        test_suite = convert_numpy_types(test_suite)
        
        # Save test suite
        output_file = self.output_dir / 'master_test_suite.json'
        with open(output_file, 'w') as f:
            json.dump(test_suite, f, indent=2)
        
        # Generate FORTRAN input files for reference validation
        self.generate_fortran_inputs(systematic_tests, edge_cases)
        
        # Generate validation plots
        self.generate_test_coverage_plots(test_suite)
        
        print(f"\nTest suite generated:")
        print(f"  Total tests: {test_suite['metadata']['total_tests']}")
        print(f"  Systematic: {len(systematic_tests)}")
        print(f"  Edge cases: {len(edge_cases)}")
        print(f"  Conservation: {len(conservation_tests)}")
        print(f"  Boundary: {len(boundary_tests)}")
        print(f"  Saved to: {output_file}")
        
        return test_suite
    
    def generate_fortran_inputs(self, systematic_tests, edge_cases):
        """Generate FORTRAN-compatible input files for validation"""
        
        all_tests = systematic_tests + edge_cases
        
        # Save as FORTRAN input format
        fortran_file = self.output_dir / 'validation_inputs.dat'
        with open(fortran_file, 'w') as f:
            f.write(f"{len(all_tests)}\n")
            for test in all_tests:
                a = test['a']
                M = 1.0  # Geometric units
                b = test['b'] 
                r0 = test['r0']
                theta0 = test['theta0']
                phi0 = 0.0  # Standard
                f.write(f"{a:15.8e} {M:15.8e} {b:15.8e} {r0:15.8e} {theta0:15.8e} {phi0:15.8e}\n")
        
        print(f"FORTRAN input file: {fortran_file}")
        
    def generate_test_coverage_plots(self, test_suite):
        """Generate visualization of test parameter coverage"""
        
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        
        # Extract test data
        systematic = test_suite['systematic_tests']
        edges = test_suite['edge_cases']
        
        a_sys = [t['a'] for t in systematic]
        b_sys = [t['b'] for t in systematic]
        a_edge = [t['a'] for t in edges if 'a' in t]
        b_edge = [t['b'] for t in edges if 'b' in t]
        
        # Spin parameter distribution
        axes[0,0].hist(test_suite['spin_parameter_sweep'], bins=20, alpha=0.7, 
                      label='Spin sweep', edgecolor='black')
        axes[0,0].axvline(0, color='red', linestyle='--', label='Schwarzschild')
        axes[0,0].axvline(0.998, color='red', linestyle='--', label='Near-extremal')
        axes[0,0].set_xlabel('Kerr spin parameter a/M')
        axes[0,0].set_ylabel('Number of tests')
        axes[0,0].set_title('Spin Parameter Coverage')
        axes[0,0].legend()
        axes[0,0].grid(True, alpha=0.3)
        
        # Impact parameter distribution
        axes[0,1].scatter(a_sys, b_sys, alpha=0.6, s=20, label='Systematic')
        axes[0,1].scatter(a_edge, b_edge, alpha=0.8, s=30, color='red', label='Edge cases')
        axes[0,1].set_xlabel('Kerr spin parameter a/M')
        axes[0,1].set_ylabel('Impact parameter b/M')
        axes[0,1].set_title('Parameter Space Coverage')
        axes[0,1].set_yscale('log')
        axes[0,1].legend()
        axes[0,1].grid(True, alpha=0.3)
        
        # Edge case types
        edge_types = {}
        for edge in edges:
            etype = edge.get('type', 'unknown')
            edge_types[etype] = edge_types.get(etype, 0) + 1
        
        axes[1,0].bar(edge_types.keys(), edge_types.values())
        axes[1,0].set_xlabel('Edge Case Type')
        axes[1,0].set_ylabel('Count')
        axes[1,0].set_title('Edge Case Distribution')
        axes[1,0].tick_params(axis='x', rotation=45)
        
        # Test category summary
        categories = test_suite['metadata']['categories']
        axes[1,1].pie(categories.values(), labels=categories.keys(), autopct='%1.1f%%')
        axes[1,1].set_title('Test Suite Composition')
        
        plt.tight_layout()
        plt.savefig(self.output_dir / 'test_coverage.png', dpi=150, bbox_inches='tight')
        plt.close()
        
        print(f"Coverage plot: {self.output_dir / 'test_coverage.png'}")

def main():
    """Generate comprehensive test suite for GEOKERR port validation"""
    
    if len(sys.argv) > 1:
        output_dir = sys.argv[1]
    else:
        output_dir = "comprehensive_tests"
        
    generator = GeodeticTestGenerator(output_dir)
    test_suite = generator.save_test_suite()
    
    print(f"\n✅ COMPREHENSIVE TEST SUITE READY")
    print(f"Use this for validating CUDA port against FORTRAN reference")
    
    return test_suite

if __name__ == '__main__':
    main()