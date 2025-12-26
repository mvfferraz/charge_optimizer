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
        
        // Line 3: num_atoms, origin (in Bohr)
        std::getline(file, line);
        std::istringstream iss1(line);
        int num_atoms;
        double origin_x, origin_y, origin_z;
        iss1 >> num_atoms >> origin_x >> origin_y >> origin_z;
        
        Eigen::Vector3d origin(origin_x, origin_y, origin_z);
        
        // Lines 4-6: Grid vectors (in Bohr)
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
        
        std::cout << "  Grid dimensions: " << nx << " x " << ny << " x " << nz << std::endl;
        std::cout << "  Grid spacing: " << vx.norm() << " Bohr (keeping atomic units)" << std::endl;
        
        // Read and store atom positions (KEEP IN BOHR!)
        std::vector<Eigen::Vector3d> atom_positions;
        std::vector<int> atomic_numbers;
        
        for (int i = 0; i < std::abs(num_atoms); ++i) {
            std::getline(file, line);
            std::istringstream iss_atom(line);
            int atomic_num;
            double charge, ax, ay, az;
            iss_atom >> atomic_num >> charge >> ax >> ay >> az;
            
            // CRITICAL: Keep positions in Bohr (atomic units)
            Eigen::Vector3d atom_pos(ax, ay, az);
            atom_positions.push_back(atom_pos);
            atomic_numbers.push_back(atomic_num);
        }
        
        std::cout << "  Atom positions stored in Bohr (atomic units)" << std::endl;
        
        // Read volumetric data (ESP in atomic units)
        std::vector<double> values;
        values.reserve(nx * ny * nz);
        
        double val;
        while (file >> val) {
            values.push_back(val);
        }
        
        std::cout << "  ESP values read: " << values.size() << " (expected: " << (nx*ny*nz) << ")" << std::endl;
        
        if (values.empty()) {
            throw std::runtime_error("No ESP values read from CUBE file!");
        }
        
        // AUTO-DETECT SIGN CONVENTION
        bool should_flip_sign = false;
        
        // Check if molecule has electronegative atoms
        bool has_electroneg = false;
        for (int Z : atomic_numbers) {
            if (Z >= 6) {  // C, N, O, F, etc.
                has_electroneg = true;
                break;
            }
        }
        
        if (has_electroneg) {
            // Sample ESP in the "shell" around the molecule
            // Range: 2.0-5.0 Bohr from nearest atom
            double sum_esp = 0.0;
            int count = 0;
            
            size_t idx = 0;
            for (int i = 0; i < nx && idx < values.size(); ++i) {
                for (int j = 0; j < ny && idx < values.size(); ++j) {
                    for (int k = 0; k < nz && idx < values.size(); ++k) {
                        Eigen::Vector3d pos = origin + i * vx + j * vy + k * vz;
                        
                        // Find distance to nearest atom (in Bohr)
                        double min_dist = 1000.0;
                        for (const auto& atom_pos : atom_positions) {
                            double dist = (pos - atom_pos).norm();
                            if (dist < min_dist) min_dist = dist;
                        }
                        
                        // Sample in the "shell" region (2-5 Bohr)
                        if (min_dist >= 2.0 && min_dist <= 5.0 && std::abs(values[idx]) < 5.0) {
                            sum_esp += values[idx];
                            count++;
                        }
                        idx++;
                    }
                }
            }
            
            if (count > 100) {
                double avg_esp = sum_esp / count;
                std::cout << "  Sign detection: sampled " << count << " points" << std::endl;
                std::cout << "  Average ESP in molecular shell: " << avg_esp << " a.u." << std::endl;
                
                // For molecules with electronegative atoms, avg ESP should be negative
                if (avg_esp > 0.001) {
                    should_flip_sign = true;
                    std::cout << "  ⚠️  INVERTED SIGN DETECTED - flipping ESP signs!" << std::endl;
                } else {
                    std::cout << "  ✓ Standard ESP sign convention" << std::endl;
                }
            }
        }
        
        // Build grid, filtering extreme points
        size_t idx = 0;
        int filtered_close = 0;
        int filtered_extreme = 0;
        
        for (int i = 0; i < nx && idx < values.size(); ++i) {
            for (int j = 0; j < ny && idx < values.size(); ++j) {
                for (int k = 0; k < nz && idx < values.size(); ++k) {
                    Eigen::Vector3d pos = origin + i * vx + j * vy + k * vz;
                    double esp_val = values[idx];
                    
                    // Find minimum distance to any nucleus (in Bohr)
                    bool too_close = false;
                    double min_dist = 1000.0;
                    
                    for (const auto& atom_pos : atom_positions) {
                        double dist = (pos - atom_pos).norm();
                        if (dist < min_dist) min_dist = dist;
                        
                        // Filter points VERY close to nuclei (< 1 Bohr ≈ 0.53 Å)
                        if (dist < 1.5) {
                            too_close = true;
                            filtered_close++;
                            break;
                        }
                    }
                    
                    // Filter extreme ESP values with distance-dependent threshold
                    double esp_limit = 20.0;  // Default for points far from nuclei
                    if (min_dist < 2.0) {
                        esp_limit = 50.0;  // More lenient for closer points
                    }
                    
                    if (std::abs(esp_val) > esp_limit) {
                        too_close = true;
                        filtered_extreme++;
                    }
                    
                    // Add point to grid if reasonable
                    if (!too_close) {
                        double final_esp = should_flip_sign ? -esp_val : esp_val;
                        
                        // CRITICAL: Store position in BOHR (atomic units)
                        // Store ESP in a.u. (Hartree/e)
                        grid.add_point(pos, final_esp);
                    }
                    
                    idx++;
                }
            }
        }
        
        std::cout << "  Grid points accepted: " << grid.num_points() << std::endl;
        std::cout << "  Filtered (too close to nuclei): " << filtered_close << std::endl;
        std::cout << "  Filtered (extreme ESP values): " << filtered_extreme << std::endl;
        
        if (grid.num_points() == 0) {
            throw std::runtime_error("No valid ESP points after filtering!");
        }
        
        // Report final ESP range
        double min_val = grid.min_potential();
        double max_val = grid.max_potential();
        std::cout << "  Final ESP range: [" << min_val << ", " << max_val << "] a.u." << std::endl;
        std::cout << "  ✓ All data in atomic units (Bohr, Hartree/e)" << std::endl;
        
        return grid;
    }
};

} // namespace chargeopt
