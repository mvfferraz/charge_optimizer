#include "core/molecule.hpp"
#include "core/esp_grid.hpp"
#include "io/xyz_parser.hpp"
#include "io/cube_parser.hpp"
#include "solver/qp_solver.hpp"
#include "solver/constraints.hpp"
#include "analysis/validator.hpp"
#include "analysis/symmetry.hpp"

#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>

using namespace chargeopt;

void print_usage(const char* prog_name) {
    std::cout << "\nCharge Optimizer - Atomic Partial Charge Fitting via QP\n" << std::endl;
    std::cout << "Usage: " << prog_name << " <geometry.xyz> <esp.cube> [options]\n" << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << "  -o, --output <file>    Output file for charges (default: charges.txt)" << std::endl;
    std::cout << "  -q, --total-charge <n> Total molecular charge (default: 0)" << std::endl;
    std::cout << "  -t, --tolerance <val>  Convergence tolerance (default: 1e-6)" << std::endl;
    std::cout << "  -l, --lambda <val>     Regularization parameter (default: 0.0005)" << std::endl;
    std::cout << "  -s, --symmetry <on|off> Auto-detect symmetry (default: on)" << std::endl;
    std::cout << "  -v, --verbose          Verbose output" << std::endl;
    std::cout << "  -h, --help             Show this help message" << std::endl;
    std::cout << "\nExamples:" << std::endl;
    std::cout << "  " << prog_name << " water.xyz water_esp.cube" << std::endl;
    std::cout << "  " << prog_name << " molecule.xyz molecule.cube -q -1 -o my_charges.txt" << std::endl;
    std::cout << "  " << prog_name << " complex.xyz complex.cube -l 0.001 -v\n" << std::endl;
}

int main(int argc, char** argv) {
    // Parse command-line arguments
    if (argc < 3) {
        print_usage(argv[0]);
        return 1;
    }
    
    std::string xyz_file = argv[1];
    std::string cube_file = argv[2];
    std::string output_file = "charges.txt";
    double total_charge = 0.0;
    double tolerance = 1e-6;
    double lambda = 0.0005;
    bool use_symmetry = true;
    bool verbose = false;
    
    // Parse options
    for (int i = 3; i < argc; ++i) {
        std::string arg = argv[i];
        
        if (arg == "-h" || arg == "--help") {
            print_usage(argv[0]);
            return 0;
        }
        else if ((arg == "-o" || arg == "--output") && i + 1 < argc) {
            output_file = argv[++i];
        }
        else if ((arg == "-q" || arg == "--total-charge") && i + 1 < argc) {
            total_charge = std::stod(argv[++i]);
        }
        else if ((arg == "-t" || arg == "--tolerance") && i + 1 < argc) {
            tolerance = std::stod(argv[++i]);
        }
        else if ((arg == "-l" || arg == "--lambda") && i + 1 < argc) {
            lambda = std::stod(argv[++i]);
        }
        else if ((arg == "-s" || arg == "--symmetry") && i + 1 < argc) {
            std::string val = argv[++i];
            use_symmetry = (val == "on" || val == "true" || val == "1");
        }
        else if (arg == "-v" || arg == "--verbose") {
            verbose = true;
        }
        else {
            std::cerr << "Unknown option: " << arg << std::endl;
            print_usage(argv[0]);
            return 1;
        }
    }
    
    try {
        // Banner
        std::cout << "\n╔════════════════════════════════════════════╗" << std::endl;
        std::cout << "║  Charge Optimizer v1.0                     ║" << std::endl;
        std::cout << "║  QP-based Atomic Charge Fitting            ║" << std::endl;
        std::cout << "╚════════════════════════════════════════════╝\n" << std::endl;
        
        // Load molecule
        std::cout << "Loading molecule from: " << xyz_file << std::endl;
        Molecule mol = XYZParser::parse(xyz_file);
        mol.set_total_charge(total_charge);
        std::cout << "  Atoms: " << mol.num_atoms() << std::endl;
        std::cout << "  Total charge: " << total_charge << " e\n" << std::endl;
        
        // Load ESP grid
        std::cout << "Loading ESP grid from: " << cube_file << std::endl;
        ESPGrid grid = CubeParser::parse(cube_file);
        std::cout << "  Grid points: " << grid.num_points() << std::endl;
        
        // DEBUG: Check ESP range immediately after loading
        std::cout << "\n  DEBUG MAIN: Verifying ESP range after loading..." << std::endl;
        double debug_min = grid.min_potential();
        double debug_max = grid.max_potential();
        std::cout << "  DEBUG MAIN: min_potential() = " << std::scientific << debug_min << std::endl;
        std::cout << "  DEBUG MAIN: max_potential() = " << std::scientific << debug_max << std::endl;
        
        // DEBUG: Manually check first/last points
        if (grid.num_points() > 0) {
            std::cout << "  DEBUG MAIN: First point ESP = " << grid.point(0).potential << std::endl;
            std::cout << "  DEBUG MAIN: Last point ESP = " << grid.point(grid.num_points()-1).potential << std::endl;
            
            // Manually scan for actual min/max
            double manual_min = grid.point(0).potential;
            double manual_max = grid.point(0).potential;
            size_t max_idx = 0;
            
            for (size_t i = 0; i < grid.num_points(); i++) {
                double val = grid.point(i).potential;
                if (val < manual_min) manual_min = val;
                if (val > manual_max) {
                    manual_max = val;
                    max_idx = i;
                }
            }
            
            std::cout << "  DEBUG MAIN: Manual scan - min = " << manual_min << ", max = " << manual_max << std::endl;
            std::cout << "  DEBUG MAIN: Max value found at grid point index " << max_idx << std::endl;
            std::cout << "  DEBUG MAIN: That point's position = (" 
                      << grid.point(max_idx).position(0) << ", "
                      << grid.point(max_idx).position(1) << ", "
                      << grid.point(max_idx).position(2) << ")" << std::endl;
        }
        std::cout << std::defaultfloat << std::endl;
        
        std::cout << "  ESP range: [" << grid.min_potential() << ", " 
                  << grid.max_potential() << "] V\n" << std::endl;
        
        // Build QP matrices
        std::cout << "Building QP problem..." << std::endl;
        Eigen::MatrixXd H;
        Eigen::VectorXd f;
        QPSolver::build_esp_matrices(mol, grid, H, f);
        
        // Setup constraints
        Constraints constraints;
        
        // Total charge constraint
        constraints.add_charge_constraint(mol.num_atoms(), total_charge);
        std::cout << "  Added total charge constraint\n" << std::endl;
        
        // Symmetry constraints
        if (use_symmetry) {
            auto equiv_groups = SymmetryDetector::detect_equivalent_atoms(mol);
            
            if (!equiv_groups.empty()) {
                std::cout << "Detected symmetry:" << std::endl;
                for (const auto& group : equiv_groups) {
                    std::cout << "  Equivalent atoms: ";
                    for (int idx : group) {
                        std::cout << mol.atom(idx).element << (idx + 1) << " ";
                    }
                    std::cout << std::endl;
                    
                    // Add constraints: all atoms in group have same charge
                    auto it = group.begin();
                    int first = *it;
                    ++it;
                    for (; it != group.end(); ++it) {
                        constraints.add_symmetry_constraint(first, *it, mol.num_atoms());
                    }
                }
                std::cout << std::endl;
            }
        }
        
        // Solve QP
        std::cout << "Solving QP..." << std::endl;
        QPSolver::Config config;
        config.tolerance = tolerance;
        config.regularization = lambda;
        config.verbose = verbose;
        
        QPSolver solver(config);
        QPSolution solution = solver.solve(H, f, constraints);
        
        if (!solution.converged) {
            std::cerr << "\nWarning: Optimization did not fully converge!" << std::endl;
        }
        
        std::cout << "  Converged: " << (solution.converged ? "Yes" : "No") << std::endl;
        std::cout << "  Iterations: " << solution.iterations << std::endl;
        std::cout << "  Objective value: " << std::scientific << solution.objective_value << std::defaultfloat << "\n" << std::endl;
        
        // Update molecule with fitted charges
        mol.set_charges(solution.charges);
        
        // Print charges
        std::cout << "=== Fitted Atomic Charges ===" << std::endl;
        std::cout << std::fixed << std::setprecision(4);
        double charge_sum = 0.0;
        for (size_t i = 0; i < mol.num_atoms(); ++i) {
            const auto& atom = mol.atom(i);
            std::cout << "  " << std::setw(3) << atom.element << std::setw(2) << (i + 1) 
                      << ":  " << std::setw(8) << std::showpos << atom.charge << std::noshowpos << " e" << std::endl;
            charge_sum += atom.charge;
        }
        std::cout << "  Sum:  " << std::showpos << charge_sum << std::noshowpos << " e\n" << std::endl;
        
        // Validate
        auto validation = Validator::validate(mol, grid);
        Validator::print_results(validation, verbose);
        
        // Write output
        std::cout << "\nWriting charges to: " << output_file << std::endl;
        std::ofstream out(output_file);
        if (!out.is_open()) {
            throw std::runtime_error("Cannot open output file: " + output_file);
        }
        
        out << "# Atomic partial charges fitted using QP optimization" << std::endl;
        out << "# Molecule: " << xyz_file << std::endl;
        out << "# Total charge: " << total_charge << std::endl;
        out << "# ESP RMSE: " << validation.esp_rmse << " V" << std::endl;
        out << "# Dipole moment: " << validation.dipole_moment << " D" << std::endl;
        out << "#" << std::endl;
        out << "# Atom  Element  Charge(e)" << std::endl;
        
        out << std::fixed << std::setprecision(6);
        for (size_t i = 0; i < mol.num_atoms(); ++i) {
            const auto& atom = mol.atom(i);
            out << std::setw(5) << (i + 1) << "  "
                << std::setw(7) << std::left << atom.element << std::right << "  "
                << std::setw(12) << atom.charge << std::endl;
        }
        
        out.close();
        
        std::cout << "\n✓ Optimization complete!\n" << std::endl;
        
        return 0;
        
    } catch (const std::exception& e) {
        std::cerr << "\nError: " << e.what() << std::endl;
        return 1;
    }
}
