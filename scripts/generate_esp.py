#!/usr/bin/env python3
"""
Generate ESP cube file using Psi4

Usage:
    python generate_esp.py molecule.xyz

Requirements:
    pip install psi4

This script:
1. Reads XYZ geometry
2. Runs DFT calculation (B3LYP/6-31G*)
3. Generates ESP cube file
"""

import sys
import os

def generate_esp_psi4(xyz_file):
    """Generate ESP cube file using Psi4"""
    
    try:
        import psi4
    except ImportError:
        print("Error: Psi4 not installed")
        print("Install with: conda install -c conda-forge psi4")
        sys.exit(1)
    
    # Read XYZ file
    with open(xyz_file, 'r') as f:
        lines = f.readlines()
    
    # Parse XYZ format
    n_atoms = int(lines[0].strip())
    comment = lines[1].strip()
    
    # Build geometry string
    geom_lines = []
    for i in range(2, 2 + n_atoms):
        geom_lines.append(lines[i].strip())
    
    geom_str = '\n'.join(geom_lines)
    
    # Setup molecule
    mol = psi4.geometry(f"""
    0 1
    {geom_str}
    """)
    
    # Set options
    psi4.set_options({
        'basis': '6-31G*',
        'scf_type': 'df',
        'cubeprop_tasks': ['ESP'],
        'cubic_grid_spacing': [0.3, 0.3, 0.3],
        'cubic_grid_overage': [4.0, 4.0, 4.0]
    })
    
    # Run calculation
    print(f"Running DFT calculation on {xyz_file}...")
    print("This may take a few minutes...")
    
    energy, wfn = psi4.energy('B3LYP', return_wfn=True)
    
    print(f"SCF Energy: {energy:.6f} Hartree")
    
    # Generate cube file
    psi4.cubeprop(wfn)
    
    # Find and rename cube file
    base_name = os.path.splitext(xyz_file)[0]
    cube_file = f"{base_name}_esp.cube"
    
    # Psi4 generates ESP.cube by default
    if os.path.exists('ESP.cube'):
        os.rename('ESP.cube', cube_file)
        print(f"\nESP cube file created: {cube_file}")
        print(f"\nNow run:")
        print(f"  ./charge_optimizer {xyz_file} {cube_file}")
    else:
        print("\nError: ESP.cube not generated")
        sys.exit(1)

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage: python generate_esp.py molecule.xyz")
        sys.exit(1)
    
    xyz_file = sys.argv[1]
    
    if not os.path.exists(xyz_file):
        print(f"Error: File not found: {xyz_file}")
        sys.exit(1)
    
    generate_esp_psi4(xyz_file)
