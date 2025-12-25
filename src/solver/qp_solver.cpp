#include "qp_solver.hpp"
#include "active_set.hpp"
#include <iostream>
#include <cmath>

namespace chargeopt {

void QPSolver::build_esp_matrices(const Molecule& mol,
                                  const ESPGrid& grid,
                                  Eigen::MatrixXd& H,
                                  Eigen::VectorXd& f) {
    const int n_atoms = mol.num_atoms();
    const int n_points = grid.num_points();
    
    // Build matrix A: ESP contribution matrix
    // A(i,j) = 1/r_ij where r_ij is distance from atom j to grid point i
    Eigen::MatrixXd A(n_points, n_atoms);
    
    for (int i = 0; i < n_points; ++i) {
        const auto& grid_point = grid.point(i);
        
        for (int j = 0; j < n_atoms; ++j) {
            const auto& atom = mol.atom(j);
            double r = (grid_point.position - atom.position).norm();
            
            // Coulomb potential: V = q/r (in atomic units, k=1)
            // Avoid division by zero
            if (r < 1e-10) r = 1e-10;
            
            A(i, j) = 1.0 / r;
        }
    }
    
    // Get target ESP values
    Eigen::VectorXd V_target = grid.potentials();
    
    // QP formulation: min 0.5 * q^T * H * q + f^T * q
    // where the objective is: ||A*q - V_target||^2 + lambda * ||q||^2
    //
    // Expanding: (A*q - V)^T * (A*q - V) + lambda * q^T * q
    //          = q^T * A^T * A * q - 2 * V^T * A * q + V^T * V + lambda * q^T * q
    //          = q^T * (A^T * A + lambda * I) * q - 2 * V^T * A * q + const
    //
    // Therefore:
    // H = 2 * (A^T * A + lambda * I)
    // f = -2 * A^T * V_target
    
    H = 2.0 * (A.transpose() * A);
    f = -2.0 * A.transpose() * V_target;
}

QPSolution QPSolver::solve(const Eigen::MatrixXd& H,
                          const Eigen::VectorXd& f,
                          const Constraints& constraints) {
    
    // Add regularization to H
    Eigen::MatrixXd H_reg = H + 2.0 * config_.regularization * Eigen::MatrixXd::Identity(H.rows(), H.cols());
    
    // Use active-set method for constrained QP
    ActiveSetSolver solver(config_.tolerance, config_.max_iterations, config_.verbose);
    return solver.solve(H_reg, f, constraints);
}

double QPSolver::compute_esp(const Eigen::Vector3d& grid_point,
                            const Molecule& mol,
                            const Eigen::VectorXd& charges) {
    double esp = 0.0;
    
    for (size_t i = 0; i < mol.num_atoms(); ++i) {
        const auto& atom = mol.atom(i);
        double r = (grid_point - atom.position).norm();
        
        if (r < 1e-10) r = 1e-10;
        
        esp += charges(i) / r;
    }
    
    return esp;
}

} // namespace chargeopt
