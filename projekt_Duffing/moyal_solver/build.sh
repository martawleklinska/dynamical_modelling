#!/bin/bash

# Build script for Structured Moyal Solver

echo "=== Building Structured Moyal Solver ==="

# Check if build directory exists
if [ ! -d "build" ]; then
    mkdir build
    echo "Created build directory"
fi

cd build

# Configure with CMake
echo "Configuring with CMake..."
cmake .. -DCMAKE_BUILD_TYPE=Release

# Build
echo "Building..."
make -j$(nproc)

if [ $? -eq 0 ]; then
    echo "Build successful!"
    echo "Executable: ./basic_example"
else
    echo "Build failed!"
    exit 1
fi
