#pragma once

#include "constraints.hpp"
#include "qp_solver.hpp"
#include <Eigen/Dense>

namespace chargeopt {

class ActiveSetSolver {
public:
    ActiveSetSolver(double tolerance = 1e-6, int max_iter = 1000, bool verbose = false)
        : tol_(tolerance), max_iter_(max_iter), verbose_(verbose) {}
    
    QPSolution solve(const Eigen::MatrixXd& H,
                    const Eigen::VectorXd& f,
                    const Constraints& constraints);

private:
    double tol_;
    int max_iter_;
    bool verbose_;
    
    // Solve unconstrained QP: min 0.5 * x^T * H * x + f^T * x
    Eigen::VectorXd solve_unconstrained(const Eigen::MatrixXd& H, const Eigen::VectorXd& f);
    
    // Solve equality-constrained QP using null-space method
    Eigen::VectorXd solve_equality_constrained(const Eigen::MatrixXd& H,
                                               const Eigen::VectorXd& f,
                                               const Constraints& constraints);
};

} // namespace chargeopt
