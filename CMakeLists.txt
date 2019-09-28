cmake_minimum_required(VERSION 3.11)
project(SubgraphCount)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(SubgraphCount.out main.cpp Graph.cpp GraphIO.cpp include/Graph.h include/JointSort.h include/ErrorCode.h include/GraphIO.h  include/Digraph.h include/GetAllCounts.h include/Triadic.h include/TriangleEstimators.h include/EstimatorUtil.h include/BaselineEstimators.h include/ConfigReader.h)

add_executable(sanitize.out sanitize.cpp sanitize.h)

add_executable(ExactCount.out ExactCount.cpp Graph.cpp GraphIO.cpp include/Graph.h include/JointSort.h include/ErrorCode.h include/GraphIO.h  include/Digraph.h include/GetAllCounts.h include/Triadic.h include/TriangleEstimators.h include/EstimatorUtil.h include/BaselineEstimators.h)