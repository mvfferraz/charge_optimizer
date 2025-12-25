#pragma once

#include <vector>
#include <Eigen/Dense>

namespace chargeopt {

struct GridPoint {
    Eigen::Vector3d position;  // Angstroms
    double potential;          // Electrostatic potential (V or a.u.)
    
    GridPoint() : position(0, 0, 0), potential(0.0) {}
    GridPoint(const Eigen::Vector3d& pos, double pot) 
        : position(pos), potential(pot) {}
};

class ESPGrid {
public:
    ESPGrid() {}
    
    void add_point(const GridPoint& point) {
        points_.push_back(point);
    }
    
    void add_point(const Eigen::Vector3d& pos, double potential) {
        points_.emplace_back(pos, potential);
    }
    
    size_t num_points() const { return points_.size(); }
    
    const GridPoint& point(size_t i) const { return points_[i]; }
    const std::vector<GridPoint>& points() const { return points_; }
    
    // Get all positions as Nx3 matrix
    Eigen::MatrixXd positions() const {
        Eigen::MatrixXd pos(points_.size(), 3);
        for (size_t i = 0; i < points_.size(); ++i) {
            pos.row(i) = points_[i].position;
        }
        return pos;
    }
    
    // Get all potentials as vector
    Eigen::VectorXd potentials() const {
        Eigen::VectorXd pot(points_.size());
        for (size_t i = 0; i < points_.size(); ++i) {
            pot(i) = points_[i].potential;
        }
        return pot;
    }
    
    // Grid statistics
    double min_potential() const {
        if (points_.empty()) return 0.0;
        double min_val = points_[0].potential;
        for (const auto& p : points_) {
            if (p.potential < min_val) min_val = p.potential;
        }
        return min_val;
    }
    
    double max_potential() const {
        if (points_.empty()) return 0.0;
        double max_val = points_[0].potential;
        for (const auto& p : points_) {
            if (p.potential > max_val) max_val = p.potential;
        }
        return max_val;
    }

private:
    std::vector<GridPoint> points_;
};

} // namespace chargeopt
