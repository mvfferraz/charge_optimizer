#include "active_set.hpp"
#include <iostream>
#include <cmath>

namespace chargeopt {

Eigen::VectorXd ActiveSetSolver::solve_unconstrained(const Eigen::MatrixXd& H, 
                                                     const Eigen::VectorXd& f) {
    // Solve H * x = -f using Cholesky decomposition
    Eigen::LLT<Eigen::MatrixXd> llt(H);
    if (llt.info() != Eigen::Success) {
        std::cerr << "Warning: Cholesky decomposition failed, using LDLT instead" << std::endl;
        Eigen::LDLT<Eigen::MatrixXd> ldlt(H);
        return ldlt.solve(-f);
    }
    return llt.solve(-f);
}

Eigen::VectorXd ActiveSetSolver::solve_equality_constrained(const Eigen::MatrixXd& H,
                                                            const Eigen::VectorXd& f,
                                                            const Constraints& constraints) {
    // Solve equality-constrained QP using KKT system:
    // [H   A^T] [x]   [-f]
    // [A    0 ] [λ] = [b ]
    //
    // Where A*x = b are the equality constraints
    
    const Eigen::MatrixXd& A = constraints.A_eq();
    const Eigen::VectorXd& b = constraints.b_eq();
    
    const int n = H.rows();  // Number of variables
    const int m = A.rows();  // Number of constraints
    
    if (m == 0) {
        // No constraints - solve unconstrained
        return solve_unconstrained(H, f);
    }
    
    // Build KKT system
    Eigen::MatrixXd KKT(n + m, n + m);
    KKT.setZero();
    
    KKT.topLeftCorner(n, n) = H;
    KKT.topRightCorner(n, m) = A.transpose();
    KKT.bottomLeftCorner(m, n) = A;
    
    Eigen::VectorXd rhs(n + m);
    rhs.head(n) = -f;
    rhs.tail(m) = b;
    
    // Solve KKT system
    Eigen::VectorXd solution = KKT.fullPivLu().solve(rhs);
    
    // Extract primal variables (charges)
    return solution.head(n);
}

QPSolution ActiveSetSolver::solve(const Eigen::MatrixXd& H,
                                 const Eigen::VectorXd& f,
                                 const Constraints& constraints) {
    
    QPSolution result;
    
    if (verbose_) {
        std::cout << "Active-Set QP Solver" << std::endl;
        std::cout << "  Variables: " << H.rows() << std::endl;
        std::cout << "  Constraints: " << constraints.num_constraints() << std::endl;
    }
    
    // For equality-constrained QP, we can solve directly using KKT conditions
    result.charges = solve_equality_constrained(H, f, constraints);
    
    // Check convergence
    result.converged = constraints.is_satisfied(result.charges, tol_);
    result.iterations = 1;
    
    // Compute objective value: 0.5 * x^T * H * x + f^T * x
    result.objective_value = 0.5 * result.charges.dot(H * result.charges) + f.dot(result.charges);
    
    if (verbose_) {
        std::cout << "  Converged: " << (result.converged ? "Yes" : "No") << std::endl;
        std::cout << "  Objective: " << result.objective_value << std::endl;
        
        // Check constraint satisfaction
        if (constraints.num_constraints() > 0) {
            const Eigen::VectorXd& b = constraints.b_eq();
            Eigen::VectorXd residual = constraints.A_eq() * result.charges - b;
            std::cout << "  Constraint residual: " << residual.norm() << std::endl;
        }
    }
    
    return result;
}

} // namespace chargeopt
