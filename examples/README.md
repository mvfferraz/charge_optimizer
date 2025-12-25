# Example Molecules

This directory contains pre-computed test cases for the charge optimizer.

## Available Examples

### 1. Water (H₂O)
- **File:** `water/water.xyz`, `water/water_esp.cube`
- **Atoms:** 3
- **Grid points:** 8,000 (20×20×20)
- **Expected charges:**
  - O: ~-0.83 e
  - H: ~+0.42 e (each)
- **Run:**
  ```bash
  ./charge_optimizer ../examples/water/water.xyz ../examples/water/water_esp.cube
  ```

### 2. Methane (CH₄)
- **File:** `methane/methane.xyz`, `methane/methane_esp.cube`
- **Atoms:** 5
- **Grid points:** 8,000 (20×20×20)
- **Expected charges:**
  - C: ~-0.42 e
  - H: ~+0.105 e (each, all symmetric)
- **Features:** Tests symmetry detection (4 equivalent H atoms)
- **Run:**
  ```bash
  ./charge_optimizer ../examples/methane/methane.xyz ../examples/methane/methane_esp.cube
  ```

### 3. Acetone ((CH₃)₂CO)
- **File:** `acetone/acetone.xyz`, `acetone/acetone_esp.cube`
- **Atoms:** 10
- **Grid points:** 15,625 (25×25×25)
- **Features:**
  - Carbonyl group (C=O)
  - Multiple symmetries (2 CH₃ groups)
  - Larger molecule
- **Run:**
  ```bash
  ./charge_optimizer ../examples/acetone/acetone.xyz ../examples/acetone/acetone_esp.cube
  ```

## About the ESP Files

The ESP cube files in this directory were generated using synthetic data for demonstration purposes. They approximate realistic quantum mechanical ESP values.

For real research, you should:
1. Generate ESP files using quantum chemistry software (Psi4, ORCA, Gaussian)
2. Use appropriate basis sets (6-31G*, 6-31+G*, etc.)
3. Use sufficient grid density

See the main README for instructions on generating ESP files with Psi4.

## Adding Your Own Examples

To test your own molecules:

1. Create XYZ file with geometry
2. Generate ESP cube file (see `scripts/generate_esp.py`)
3. Run the optimizer:
   ```bash
   ./charge_optimizer your_molecule.xyz your_molecule_esp.cube
   ```
