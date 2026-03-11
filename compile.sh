#!/bin/bash

set -e

SRC_DIR=src
BUILD_DIR=build
EXE=gen

mkdir -p ${BUILD_DIR}

echo "Compiling..."

g++ -std=c++17 -O2 -Wall -Wextra ${SRC_DIR}/*.cpp -o ${BUILD_DIR}/${EXE}

echo "Done."