#pragma once

#include <Eigen/Dense>
#include <vector>
#include <set>

namespace chargeopt {

class Constraints {
public:
    Constraints() {}
    
    // Add equality constraint: sum of charges = total_charge
    void add_charge_constraint(int num_atoms, double total_charge) {
        Eigen::VectorXd a = Eigen::VectorXd::Ones(num_atoms);
        A_eq_.conservativeResize(A_eq_.rows() + 1, num_atoms);
        A_eq_.row(A_eq_.rows() - 1) = a;
        
        b_eq_.conservativeResize(b_eq_.size() + 1);
        b_eq_(b_eq_.size() - 1) = total_charge;
    }
    
    // Add symmetry constraint: charge[i] = charge[j]
    void add_symmetry_constraint(int i, int j, int num_atoms) {
        Eigen::VectorXd a = Eigen::VectorXd::Zero(num_atoms);
        a(i) = 1.0;
        a(j) = -1.0;
        
        A_eq_.conservativeResize(A_eq_.rows() + 1, num_atoms);
        A_eq_.row(A_eq_.rows() - 1) = a;
        
        b_eq_.conservativeResize(b_eq_.size() + 1);
        b_eq_(b_eq_.size() - 1) = 0.0;
    }
    
    const Eigen::MatrixXd& A_eq() const { return A_eq_; }
    const Eigen::VectorXd& b_eq() const { return b_eq_; }
    
    int num_constraints() const { return A_eq_.rows(); }
    
    // Check if constraints are satisfied
    bool is_satisfied(const Eigen::VectorXd& x, double tol = 1e-6) const {
        if (A_eq_.rows() == 0) return true;
        Eigen::VectorXd residual = A_eq_ * x - b_eq_;
        return residual.norm() < tol;
    }

private:
    Eigen::MatrixXd A_eq_;  // Equality constraint matrix (m x n)
    Eigen::VectorXd b_eq_;  // Equality constraint vector (m)
};

} // namespace chargeopt
