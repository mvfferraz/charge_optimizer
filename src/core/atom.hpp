#pragma once

#include <string>
#include <Eigen/Dense>

namespace chargeopt {

struct Atom {
    std::string element;
    Eigen::Vector3d position;  // Angstroms
    double charge;             // Partial charge (what we're solving for)
    int index;                 // 0-based index
    
    Atom() : element(""), position(0, 0, 0), charge(0.0), index(-1) {}
    
    Atom(const std::string& elem, const Eigen::Vector3d& pos, int idx = -1)
        : element(elem), position(pos), charge(0.0), index(idx) {}
    
    // Atomic number from element symbol
    int atomic_number() const {
        if (element == "H") return 1;
        if (element == "C") return 6;
        if (element == "N") return 7;
        if (element == "O") return 8;
        if (element == "F") return 9;
        if (element == "P") return 15;
        if (element == "S") return 16;
        if (element == "Cl") return 17;
        return 0;
    }
    
    // Van der Waals radius (Angstroms) - for grid generation
    double vdw_radius() const {
        if (element == "H") return 1.20;
        if (element == "C") return 1.70;
        if (element == "N") return 1.55;
        if (element == "O") return 1.52;
        if (element == "F") return 1.47;
        if (element == "S") return 1.80;
        if (element == "P") return 1.80;
        if (element == "Cl") return 1.75;
        return 1.70;
    }
};

} // namespace chargeopt
