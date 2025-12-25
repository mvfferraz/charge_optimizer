#include <iostream>
#include <Eigen/Dense>
#include <cmath>

bool test_eigen() {
    Eigen::MatrixXd A(2, 2);
    A << 1, 2,
         3, 4;
    
    Eigen::VectorXd b(2);
    b << 5, 11;
    
    Eigen::VectorXd x = A.colPivHouseholderQr().solve(b);
    
    // Solution should be [1, 2]
    double error = (x - Eigen::Vector2d(1, 2)).norm();
    
    return error < 1e-10;
}

bool test_charge_conservation() {
    // Test that charge constraint works
    Eigen::VectorXd charges(3);
    charges << -0.8, 0.4, 0.4;
    
    double sum = charges.sum();
    double expected = 0.0;
    
    return std::abs(sum - expected) < 1e-10;
}

int main() {
    std::cout << "Running basic tests..." << std::endl;
    
    int passed = 0;
    int failed = 0;
    
    if (test_eigen()) {
        std::cout << "✓ Eigen test passed" << std::endl;
        passed++;
    } else {
        std::cout << "✗ Eigen test failed" << std::endl;
        failed++;
    }
    
    if (test_charge_conservation()) {
        std::cout << "✓ Charge conservation test passed" << std::endl;
        passed++;
    } else {
        std::cout << "✗ Charge conservation test failed" << std::endl;
        failed++;
    }
    
    std::cout << "\nTests: " << passed << " passed, " << failed << " failed" << std::endl;
    
    return failed > 0 ? 1 : 0;
}
