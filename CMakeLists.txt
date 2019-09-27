cmake_minimum_required(VERSION 3.14)
project(SubgraphCount)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(SubgraphCount.out main.cpp include/Graph.h include/JointSort.h Graph.cpp include/ErrorCode.h include/GraphIO.h GraphIO.cpp include/Digraph.h include/GetAllCounts.h include/Triadic.h include/TriangleEstimators.h include/EstimatorUtil.h include/BaselineEstimators.h include/Utils.h)

add_executable(sanitize.out sanitize.cpp sanitize.h)