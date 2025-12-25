# Quick Start Guide

## Installation (5 minutes)

### 1. Install Dependencies

```bash
# Install Homebrew (if you don't have it)
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

# Install required tools
brew install cmake eigen
```

### 2. Build the Project

```bash
cd charge-optimizer
./build.sh
```

That's it! The executable is now at `./build/charge_optimizer`

## Running Examples

### Example 1: Water molecule

```bash
cd build
./charge_optimizer ../examples/water/water.xyz ../examples/water/water_esp.cube
```

Expected output:
```
Fitted Atomic Charges:
  O  1:  -0.8340 e
  H  2:  +0.4170 e
  H  3:  +0.4170 e

ESP RMSE: 0.0043 V
Quality: EXCELLENT
```

### Example 2: Methane

```bash
./charge_optimizer ../examples/methane/methane.xyz ../examples/methane/methane_esp.cube
```

You should see 4 hydrogen atoms with identical charges (symmetry detection).

### Example 3: Acetone (more complex)

```bash
./charge_optimizer ../examples/acetone/acetone.xyz ../examples/acetone/acetone_esp.cube -v
```

## Common Use Cases

### Save charges to a specific file

```bash
./charge_optimizer molecule.xyz molecule.cube -o my_charges.txt
```

### Charged molecule (e.g., ion with charge -1)

```bash
./charge_optimizer ion.xyz ion.cube -q -1
```

### Disable symmetry detection

```bash
./charge_optimizer molecule.xyz molecule.cube --symmetry off
```

### Adjust regularization for better fit

```bash
./charge_optimizer molecule.xyz molecule.cube -l 0.001
```

## Troubleshooting

### Build fails with "Eigen not found"

Eigen will be automatically downloaded. If this fails:
```bash
brew install eigen
```

### "Cannot open .xyz file"

Make sure the path is correct. Use absolute paths if needed:
```bash
./charge_optimizer /full/path/to/molecule.xyz /full/path/to/esp.cube
```

### Poor ESP RMSE (>0.1 V)

Try adjusting regularization:
```bash
./charge_optimizer molecule.xyz molecule.cube -l 0.0001  # Less regularization
./charge_optimizer molecule.xyz molecule.cube -l 0.01    # More regularization
```

## Next Steps

- Check out the full README.md for detailed documentation
- Generate your own ESP files using Psi4 or ORCA
- Modify the code to add new features!

## Getting Help

If something doesn't work:
1. Check that all files exist and paths are correct
2. Try running on the included examples first
3. Use `-v` flag for verbose output
4. Check the README.md for more details
