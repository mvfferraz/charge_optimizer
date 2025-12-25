# Mathematical Details

## Problem Formulation

### Objective

Given:
- Molecular geometry: atom positions **r**₁, **r**₂, ..., **r**ₙ
- Quantum mechanical electrostatic potential at grid points: V_QM(**r**ᵢ) for i = 1...m

Find: atomic partial charges q₁, q₂, ..., qₙ that minimize:

```
L(q) = Σᵢ [V_QM(**r**ᵢ) - V_charges(**r**ᵢ; q)]² + λ Σⱼ qⱼ²
```

where:
- V_charges(**r**; q) = Σⱼ qⱼ / |**r** - **r**ⱼ| (Coulomb potential)
- λ is regularization parameter (prevents overfitting)

### Constraints

1. **Charge conservation**: Σⱼ qⱼ = Q_total
2. **Symmetry**: qᵢ = qⱼ for chemically equivalent atoms i, j

## Quadratic Programming Formulation

The problem can be written as a QP:

**Minimize:** ½ **q**ᵀ H **q** + **f**ᵀ **q**

**Subject to:** A_eq **q** = **b**_eq

### Matrix Construction

1. **Build matrix A** (m × n):
   ```
   A_ij = 1 / |**r**ᵢ - **r**ⱼ|
   ```
   where i indexes grid points, j indexes atoms

2. **Get target vector:**
   ```
   **v** = [V_QM(**r**₁), V_QM(**r**₂), ..., V_QM(**r**ₘ)]ᵀ
   ```

3. **Expand objective:**
   ```
   L = (A**q** - **v**)ᵀ(A**q** - **v**) + λ **q**ᵀ**q**
     = **q**ᵀ(AᵀA)**q** - 2**v**ᵀA**q** + **v**ᵀ**v** + λ **q**ᵀ**q**
     = **q**ᵀ(AᵀA + λI)**q** - 2**v**ᵀA**q** + const
   ```

4. **Therefore:**
   ```
   H = 2(AᵀA + λI)
   **f** = -2Aᵀ**v**
   ```

### Constraint Matrices

**Charge conservation:**
```
A_eq = [1, 1, 1, ..., 1]    (row vector of ones)
b_eq = [Q_total]
```

**Symmetry between atoms i and j:**
```
A_eq row = [0, ..., 1, ..., -1, ..., 0]
           (1 at position i, -1 at position j)
b_eq = 0
```

## Solution Method: KKT System

For equality-constrained QP, we solve the KKT system:

```
[ H    Aᵀ ] [ **q** ]   [ -**f** ]
[ A    0  ] [ **λ** ] = [ **b**  ]
```

where:
- **q** are the primal variables (charges we want)
- **λ** are Lagrange multipliers
- A, **b** are constraint matrix and vector

### Algorithm

```python
def solve_equality_qp(H, f, A, b):
    n = H.shape[0]  # number of variables
    m = A.shape[0]  # number of constraints
    
    # Build KKT matrix
    KKT = zeros(n + m, n + m)
    KKT[0:n, 0:n] = H
    KKT[0:n, n:n+m] = A.T
    KKT[n:n+m, 0:n] = A
    
    # Build RHS
    rhs = zeros(n + m)
    rhs[0:n] = -f
    rhs[n:n+m] = b
    
    # Solve linear system
    solution = solve(KKT, rhs)
    
    # Extract charges (first n components)
    q = solution[0:n]
    
    return q
```

## Regularization

### L2 Regularization (Ridge)

The term λ Σⱼ qⱼ² adds to the diagonal of H:

```
H_regularized = H + 2λI
```

**Effect:**
- Prevents overfitting to noisy ESP data
- Penalizes large charges
- Makes problem better conditioned

**Typical values:** λ ∈ [0.0001, 0.01]

### Choosing λ

- **Small λ** (0.0001): Tight fit to ESP, risk of overfitting
- **Large λ** (0.01): Smooth charges, may underfit
- **Default** (0.0005): Good balance for most molecules

## Symmetry Detection

### Algorithm

For each atom i:
1. Find atoms j with same element
2. Compute distance vectors to all other atoms
3. Sort distances
4. If sorted distances match (within tolerance), atoms are equivalent

### Example: Water (H₂O)

```
O at (0, 0, 0.12)
H1 at (0, 0.76, -0.47)
H2 at (0, -0.76, -0.47)

Distances from H1: [0.96 to O, 1.52 to H2]
Distances from H2: [0.96 to O, 1.52 to H1]

→ H1 and H2 are equivalent
→ Constraint: q_H1 = q_H2
```

## Validation Metrics

### ESP RMSE

```
RMSE = sqrt( (1/m) Σᵢ [V_QM(**r**ᵢ) - V_fit(**r**ᵢ)]² )
```

**Interpretation:**
- < 0.01 V: Excellent
- < 0.05 V: Good
- < 0.10 V: Acceptable
- \> 0.10 V: Poor (adjust parameters)

### Dipole Moment

```
**μ** = Σⱼ qⱼ **r**ⱼ    (in e·Å)

|**μ**| (in Debye) = |**μ**| / 0.2082
```

Compare with experimental or high-level QM values.

## Computational Complexity

- **Matrix A construction:** O(mn) where m = grid points, n = atoms
- **H = AᵀA:** O(mn²)
- **KKT solution:** O((n+k)³) where k = number of constraints
- **Typical:** n ~ 10-100 atoms, m ~ 10,000-100,000 grid points
- **Runtime:** Seconds to minutes for small molecules

## Numerical Stability

### Issues and Solutions

1. **Ill-conditioned H:**
   - Add regularization (λI term)
   - Use LDLT instead of Cholesky if needed

2. **Near-singular constraints:**
   - Check for linearly dependent constraints
   - Use full-pivot LU for KKT system

3. **Grid points too close to atoms:**
   - Clip distances: r = max(r, r_min) with r_min ~ 0.01 Å

## Extensions

### Multi-Conformer Fitting

Minimize ESP error across multiple conformations:

```
L = Σₖ wₖ Σᵢ [V_QM,k(**r**ᵢ,ₖ) - V_charges(**r**ᵢ,ₖ; **q**)]² + λ ||**q**||²
```

Ensures charges are transferable between conformations.

### Hyperbolic Restraints

Instead of L2 regularization, use RESP hyperbolic restraints:

```
L_restraint = Σⱼ aⱼ [(qⱼ² + b²)^(1/2) - b]
```

Allows larger charges on buried atoms while restraining surface atoms.

## References

1. Bayly et al., "A well-behaved electrostatic potential based method using charge restraints for deriving atomic charges" (1993)
2. Cornell et al., "A second generation force field for the simulation of proteins, nucleic acids, and organic molecules" (1995)
3. Nocedal & Wright, "Numerical Optimization" (2006)
