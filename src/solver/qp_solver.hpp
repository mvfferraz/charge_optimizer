#pragma once

#include "../core/molecule.hpp"
#include "../core/esp_grid.hpp"
#include "constraints.hpp"
#include <Eigen/Dense>

namespace chargeopt {

struct QPSolution {
    Eigen::VectorXd charges;
    double objective_value;
    bool converged;
    int iterations;
    
    QPSolution() : objective_value(0), converged(false), iterations(0) {}
};

class QPSolver {
public:
    struct Config {
        double tolerance = 1e-6;
        double regularization = 0.0005;  // Lambda for L2 regularization
        int max_iterations = 1000;
        bool verbose = false;
        
        Config() {}
    };
    
    QPSolver(const Config& config = Config()) : config_(config) {}
    
    // Solve: min 0.5 * x^T * H * x + f^T * x
    //        subject to: A_eq * x = b_eq
    QPSolution solve(const Eigen::MatrixXd& H, 
                     const Eigen::VectorXd& f,
                     const Constraints& constraints);
    
    // Build QP problem from molecule and ESP grid
    static void build_esp_matrices(const Molecule& mol,
                                   const ESPGrid& grid,
                                   Eigen::MatrixXd& H,
                                   Eigen::VectorXd& f);

private:
    Config config_;
    
    // Compute ESP at grid point from point charges
    static double compute_esp(const Eigen::Vector3d& grid_point,
                            const Molecule& mol,
                            const Eigen::VectorXd& charges);
};

} // namespace chargeopt
