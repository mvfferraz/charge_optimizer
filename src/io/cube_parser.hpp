#pragma once

#include "../core/esp_grid.hpp"
#include <string>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <iostream>
#include <cmath>
#include <algorithm>

namespace chargeopt {

class CubeParser {
public:
    static ESPGrid parse(const std::string& filename) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open CUBE file: " + filename);
        }
        
        ESPGrid grid;
        std::string line;
        
        // Line 1-2: Comments
        std::getline(file, line);
        std::getline(file, line);
        
        // Line 3: num_atoms, origin
        std::getline(file, line);
        std::istringstream iss1(line);
        int num_atoms;
        double origin_x, origin_y, origin_z;
        iss1 >> num_atoms >> origin_x >> origin_y >> origin_z;
        
        Eigen::Vector3d origin(origin_x, origin_y, origin_z);
        
        // Lines 4-6: Grid vectors
        int nx, ny, nz;
        Eigen::Vector3d vx, vy, vz;
        
        std::getline(file, line);
        std::istringstream iss2(line);
        iss2 >> nx >> vx(0) >> vx(1) >> vx(2);
        
        std::getline(file, line);
        std::istringstream iss3(line);
        iss3 >> ny >> vy(0) >> vy(1) >> vy(2);
        
        std::getline(file, line);
        std::istringstream iss4(line);
        iss4 >> nz >> vz(0) >> vz(1) >> vz(2);
        
        // Convert from Bohr to Angstrom if needed
        bool is_bohr = (std::abs(vx(0)) > 0.1 || std::abs(vy(1)) > 0.1 || std::abs(vz(2)) > 0.1);
        const double bohr_to_angstrom = 0.529177249;
        
        if (is_bohr) {
            origin *= bohr_to_angstrom;
            vx *= bohr_to_angstrom;
            vy *= bohr_to_angstrom;
            vz *= bohr_to_angstrom;
        }
        
        std::cout << "  Grid dimensions: " << nx << " x " << ny << " x " << nz << std::endl;
        
        // Read and store atom positions (for filtering later)
        std::vector<Eigen::Vector3d> atom_positions;
        for (int i = 0; i < std::abs(num_atoms); ++i) {
            std::getline(file, line);
            std::istringstream iss_atom(line);
            int atomic_num;
            double charge, ax, ay, az;
            iss_atom >> atomic_num >> charge >> ax >> ay >> az;
            
            Eigen::Vector3d atom_pos(ax, ay, az);
            if (is_bohr) {
                atom_pos *= bohr_to_angstrom;
            }
            atom_positions.push_back(atom_pos);
        }
        
        // Read volumetric data
        std::vector<double> values;
        values.reserve(nx * ny * nz);
        
        double val;
        while (file >> val) {
            values.push_back(val);
        }
        
        std::cout << "  Values read: " << values.size() << " (expected: " << (nx*ny*nz) << ")" << std::endl;
        
        if (values.empty()) {
            throw std::runtime_error("No ESP values read from CUBE file!");
        }
        
        // Build grid, filtering points too close to nuclei
        size_t idx = 0;
        int filtered_count = 0;
        int extreme_value_count = 0;
        
        for (int i = 0; i < nx && idx < values.size(); ++i) {
            for (int j = 0; j < ny && idx < values.size(); ++j) {
                for (int k = 0; k < nz && idx < values.size(); ++k) {
                    Eigen::Vector3d pos = origin + i * vx + j * vy + k * vz;
                    double esp_val = values[idx];
                    
                    // Check if point is too close to any nucleus
                    bool too_close = false;
                    for (const auto& atom_pos : atom_positions) {
                        double dist = (pos - atom_pos).norm();
                        if (dist < 0.8) {  // Within 0.8 Angstrom of nucleus
                            too_close = true;
                            filtered_count++;
                            break;
                        }
                    }
                    
                    // Also filter extreme ESP values (divergent near nuclei)
                    if (std::abs(esp_val) > 1.0) {
                        extreme_value_count++;
                        too_close = true;
                    }
                    
                    // Only add point if it's reasonable
                    if (!too_close) {
                        grid.add_point(pos, esp_val);
                    }
                    
                    idx++;
                }
            }
        }
        
        std::cout << "  Grid points accepted: " << grid.num_points() << std::endl;
        std::cout << "  Filtered (too close to nuclei): " << filtered_count << std::endl;
        std::cout << "  Filtered (extreme ESP values): " << extreme_value_count << std::endl;
        
        if (grid.num_points() == 0) {
            throw std::runtime_error("No valid ESP points after filtering!");
        }
        
        // Report final ESP range
        double min_val = grid.min_potential();
        double max_val = grid.max_potential();
        std::cout << "  Final ESP range: [" << min_val << ", " << max_val << "] a.u." << std::endl;
        
        return grid;
    }
};

} // namespace chargeopt
