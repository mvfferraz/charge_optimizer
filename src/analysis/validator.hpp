#pragma once

#include "../core/molecule.hpp"
#include "../core/esp_grid.hpp"
#include <Eigen/Dense>
#include <cmath>
#include <iostream>

namespace chargeopt {

class Validator {
public:
    struct ValidationResults {
        double esp_rmse;
        double esp_max_error;
        double dipole_moment;
        double total_charge;
        
        std::string quality() const {
            if (esp_rmse < 0.01) return "EXCELLENT";
            if (esp_rmse < 0.05) return "GOOD";
            if (esp_rmse < 0.10) return "ACCEPTABLE";
            return "POOR";
        }
    };
    
    static ValidationResults validate(const Molecule& mol, const ESPGrid& grid) {
        ValidationResults results;
        
        // Compute ESP RMSE (all in atomic units)
        double sum_sq_error = 0.0;
        double max_error = 0.0;
        
        for (size_t i = 0; i < grid.num_points(); ++i) {
            const auto& point = grid.point(i);
            
            // Compute ESP from fitted charges (atomic units: q/r with r in Bohr)
            double esp_fitted = compute_esp_at_point(mol, point.position);
            double error = std::abs(esp_fitted - point.potential);
            
            sum_sq_error += error * error;
            if (error > max_error) max_error = error;
        }
        
        // RMSE in atomic units (1 a.u. ≈ 27.2 eV ≈ 27.2 V)
        results.esp_rmse = std::sqrt(sum_sq_error / grid.num_points());
        results.esp_max_error = max_error;
        
        // Compute dipole moment
        // In atomic units: μ (a.u.) = Σ q_i * r_i (with r in Bohr)
        // Convert to Debye: 1 a.u. = 2.5417464 Debye
        Eigen::Vector3d dipole(0, 0, 0);
        for (size_t i = 0; i < mol.num_atoms(); ++i) {
            const auto& atom = mol.atom(i);
            dipole += atom.charge * atom.position;  // Both in a.u.
        }
        
        double dipole_au = dipole.norm();
        results.dipole_moment = dipole_au * 2.5417464;  // Convert to Debye
        
        // Total charge
        results.total_charge = 0.0;
        for (size_t i = 0; i < mol.num_atoms(); ++i) {
            results.total_charge += mol.atom(i).charge;
        }
        
        return results;
    }
    
    static void print_results(const ValidationResults& results, bool verbose = false) {
        std::cout << "\n=== Validation Results ===" << std::endl;
        std::cout << "  ESP RMSE:       " << results.esp_rmse << " a.u." << std::endl;
        std::cout << "  ESP max error:  " << results.esp_max_error << " a.u." << std::endl;
        std::cout << "  Dipole moment:  " << results.dipole_moment << " D" << std::endl;
        std::cout << "  Total charge:   " << results.total_charge << " e" << std::endl;
        std::cout << "  Quality:        " << results.quality() << std::endl;
        
        if (verbose) {
            std::cout << "\nInterpretation:" << std::endl;
            std::cout << "  RMSE < 0.01 a.u.  : Excellent fit" << std::endl;
            std::cout << "  RMSE < 0.05 a.u.  : Good fit" << std::endl;
            std::cout << "  RMSE < 0.10 a.u.  : Acceptable" << std::endl;
            std::cout << "  RMSE > 0.10 a.u.  : Poor fit" << std::endl;
            std::cout << "\n  Note: 1 a.u. ≈ 27.2 eV" << std::endl;
        }
    }

private:
    static double compute_esp_at_point(const Molecule& mol, const Eigen::Vector3d& point) {
        double esp = 0.0;
        
        // Coulomb potential in atomic units: V = q/r (with r in Bohr)
        for (size_t i = 0; i < mol.num_atoms(); ++i) {
            const auto& atom = mol.atom(i);
            double r = (point - atom.position).norm();  // Distance in Bohr
            
            if (r < 1e-10) r = 1e-10;
            
            esp += atom.charge / r;  // Both in a.u.
        }
        
        return esp;
    }
};

} // namespace chargeopt
