#pragma once

#include "atom.hpp"
#include <vector>
#include <string>
#include <Eigen/Dense>

namespace chargeopt {

class Molecule {
public:
    Molecule() : total_charge_(0.0) {}
    
    void add_atom(const Atom& atom) {
        atoms_.push_back(atom);
        atoms_.back().index = atoms_.size() - 1;
    }
    
    size_t num_atoms() const { return atoms_.size(); }
    
    const Atom& atom(size_t i) const { return atoms_[i]; }
    Atom& atom(size_t i) { return atoms_[i]; }
    
    const std::vector<Atom>& atoms() const { return atoms_; }
    std::vector<Atom>& atoms() { return atoms_; }
    
    void set_total_charge(double charge) { total_charge_ = charge; }
    double total_charge() const { return total_charge_; }
    
    // Get positions as Nx3 matrix
    Eigen::MatrixXd positions() const {
        Eigen::MatrixXd pos(atoms_.size(), 3);
        for (size_t i = 0; i < atoms_.size(); ++i) {
            pos.row(i) = atoms_[i].position;
        }
        return pos;
    }
    
    // Get current charges as vector
    Eigen::VectorXd charges() const {
        Eigen::VectorXd q(atoms_.size());
        for (size_t i = 0; i < atoms_.size(); ++i) {
            q(i) = atoms_[i].charge;
        }
        return q;
    }
    
    // Set charges from vector
    void set_charges(const Eigen::VectorXd& charges) {
        for (size_t i = 0; i < atoms_.size(); ++i) {
            atoms_[i].charge = charges(i);
        }
    }
    
    // Compute molecular center of mass
    Eigen::Vector3d center_of_mass() const {
        Eigen::Vector3d com(0, 0, 0);
        double total_mass = 0.0;
        
        for (const auto& atom : atoms_) {
            double mass = atom.atomic_number();  // Simplified
            com += mass * atom.position;
            total_mass += mass;
        }
        
        return com / total_mass;
    }
    
    // Compute dipole moment from current charges (Debye)
    double dipole_moment() const {
        Eigen::Vector3d dipole(0, 0, 0);
        
        for (const auto& atom : atoms_) {
            dipole += atom.charge * atom.position;
        }
        
        // Convert from e*Angstrom to Debye (1 D = 0.2081943 e*Angstrom)
        return dipole.norm() / 0.2081943;
    }

private:
    std::vector<Atom> atoms_;
    double total_charge_;
};

} // namespace chargeopt
