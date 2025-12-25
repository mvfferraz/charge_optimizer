#!/bin/bash

echo "╔════════════════════════════════════════════╗"
echo "║  Charge Optimizer - Build Script          ║"
echo "╚════════════════════════════════════════════╝"
echo ""

if ! command -v brew &> /dev/null; then
    echo "Homebrew not found!"
    echo "Install it from: https://brew.sh"
    echo ""
    echo "Run this command:"
    echo '/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"'
    exit 1
fi

echo "Homebrew found"

if ! command -v cmake &> /dev/null; then
    echo "Installing CMake..."
    brew install cmake
fi

echo "CMake found"

if brew list eigen &> /dev/null; then
    echo "✓ Eigen3 found (via Homebrew)"
else
    echo "⚠ Eigen3 not found, will be downloaded automatically during build"
fi

echo ""
echo "Creating build directory..."
rm -rf build
mkdir -p build
cd build

echo ""
echo "Running CMake..."
cmake -DCMAKE_BUILD_TYPE=Release ..

if [ $? -ne 0 ]; then
    echo ""
    echo "CMake configuration failed!"
    exit 1
fi

echo ""
echo "Building project..."
make -j$(sysctl -n hw.ncpu)

if [ $? -ne 0 ]; then
    echo ""
    echo "Build failed!"
    exit 1
fi

echo ""
echo "Running tests..."
ctest --output-on-failure

echo ""
echo "╔════════════════════════════════════════════╗"
echo "║  Build Complete!                          ║"
echo "╚════════════════════════════════════════════╝"
echo ""
echo "Executable created: ./build/charge_optimizer"
echo ""
echo "Try it out:"
echo "  ./build/charge_optimizer ../examples/water/water.xyz ../examples/water/water_esp.cube"
echo ""
