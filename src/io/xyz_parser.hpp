#pragma once

#include "../core/molecule.hpp"
#include <string>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace chargeopt {

class XYZParser {
public:
    static Molecule parse(const std::string& filename) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open file: " + filename);
        }
        
        Molecule mol;
        std::string line;
        int num_atoms = 0;
        int line_num = 0;
        
        // Read number of atoms (first line)
        if (std::getline(file, line)) {
            std::istringstream iss(line);
            if (!(iss >> num_atoms) || num_atoms <= 0) {
                throw std::runtime_error("Invalid number of atoms in XYZ file");
            }
            line_num++;
        } else {
            throw std::runtime_error("Empty XYZ file");
        }
        
        // Skip comment line (second line)
        if (std::getline(file, line)) {
            line_num++;
        }
        
        // Read atoms
        int atoms_read = 0;
        while (std::getline(file, line) && atoms_read < num_atoms) {
            line_num++;
            std::istringstream iss(line);
            
            std::string element;
            double x, y, z;
            
            if (iss >> element >> x >> y >> z) {
                Eigen::Vector3d pos(x, y, z);
                mol.add_atom(Atom(element, pos, atoms_read));
                atoms_read++;
            } else if (!line.empty()) {
                throw std::runtime_error("Invalid atom line " + std::to_string(line_num) + " in XYZ file");
            }
        }
        
        if (atoms_read != num_atoms) {
            throw std::runtime_error("Expected " + std::to_string(num_atoms) + 
                                   " atoms but read " + std::to_string(atoms_read));
        }
        
        // Auto-detect total charge (neutral by default)
        mol.set_total_charge(0.0);
        
        return mol;
    }
};

} // namespace chargeopt
