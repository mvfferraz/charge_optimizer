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
            
            // Coulomb potential: V = q/r (in atomic units)
            // Avoid division by zero
            if (r < 1e-10) r = 1e-10;
            
            A(i, j) = 1.0 / r;
        }
    }
    
    // Get target ESP values
    Eigen::VectorXd V_target = grid.potentials();
    
    // NORMALIZE A for better conditioning
    Eigen::VectorXd scale(n_atoms);
    for (int j = 0; j < n_atoms; ++j) {
        scale(j) = A.col(j).norm();
        if (scale(j) > 1e-10) {
            A.col(j) /= scale(j);
        }
    }
    
    // QP formulation with normalized A
    H = 2.0 * (A.transpose() * A);
    f = -2.0 * A.transpose() * V_target;
    
    // Scale f back to account for normalization
    for (int j = 0; j < n_atoms; ++j) {
        if (scale(j) > 1e-10) {
            f(j) /= scale(j);
        }
    }
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
