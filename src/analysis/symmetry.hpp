#pragma once

#include "../core/molecule.hpp"
#include <vector>
#include <set>
#include <map>

namespace chargeopt {

class SymmetryDetector {
public:
    // Find groups of symmetrically equivalent atoms
    // Returns: vector of sets, where each set contains indices of equivalent atoms
    static std::vector<std::set<int>> detect_equivalent_atoms(const Molecule& mol, 
                                                               double tolerance = 0.1) {
        std::vector<std::set<int>> groups;
        std::vector<bool> assigned(mol.num_atoms(), false);
        
        for (size_t i = 0; i < mol.num_atoms(); ++i) {
            if (assigned[i]) continue;
            
            std::set<int> group;
            group.insert(i);
            assigned[i] = true;
            
            const auto& atom_i = mol.atom(i);
            
            // Find atoms with same element and similar local environment
            for (size_t j = i + 1; j < mol.num_atoms(); ++j) {
                if (assigned[j]) continue;
                
                const auto& atom_j = mol.atom(j);
                
                // Must be same element
                if (atom_i.element != atom_j.element) continue;
                
                // Check if local environments are similar
                if (is_equivalent_environment(mol, i, j, tolerance)) {
                    group.insert(j);
                    assigned[j] = true;
                }
            }
            
            if (group.size() > 1) {
                groups.push_back(group);
            }
        }
        
        return groups;
    }
    
private:
    // Simple heuristic: atoms are equivalent if they have similar distances to other atoms
    static bool is_equivalent_environment(const Molecule& mol, int i, int j, double tol) {
        const auto& atom_i = mol.atom(i);
        const auto& atom_j = mol.atom(j);
        
        // Compute distances to all other atoms
        std::vector<double> dist_i, dist_j;
        
        for (size_t k = 0; k < mol.num_atoms(); ++k) {
            if (k == static_cast<size_t>(i) || k == static_cast<size_t>(j)) continue;
            
            double d_i = (atom_i.position - mol.atom(k).position).norm();
            double d_j = (atom_j.position - mol.atom(k).position).norm();
            
            dist_i.push_back(d_i);
            dist_j.push_back(d_j);
        }
        
        // Sort distances
        std::sort(dist_i.begin(), dist_i.end());
        std::sort(dist_j.begin(), dist_j.end());
        
        // Compare sorted distances
        if (dist_i.size() != dist_j.size()) return false;
        
        for (size_t k = 0; k < dist_i.size(); ++k) {
            if (std::abs(dist_i[k] - dist_j[k]) > tol) {
                return false;
            }
        }
        
        return true;
    }
};

} // namespace chargeopt
