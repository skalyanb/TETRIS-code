cmake_minimum_required(VERSION 3.11)
project(SubgraphCount)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(SubgraphCount.out main.cpp Graph.cpp GraphIO.cpp include/Graph.h include/JointSort.h include/ErrorCode.h
        include/GraphIO.h  include/Digraph.h include/GetAllCounts.h include/Triadic.h include/TriangleEstimators.h
        include/EstimatorUtil.h  include/util/ConfigReader.h include/EstimateEdgeCount.h
        include/FileIO.h include/EstimatorUtilStruct.h include/EstimatorUtilStats.h include/baseline/VertexMCMC.h
        include/baseline/TripleMCMC.h include/util/RandomWalkUtils.h include/baseline/SubgraohRandomWalk_SRW.h include/TETRIS.h)

add_executable(sanitize.out sanitize.cpp sanitize.h)

add_executable(ExactCount.out ExactCount.cpp Graph.cpp GraphIO.cpp include/Graph.h include/JointSort.h include/ErrorCode.h include/GraphIO.h  include/Digraph.h include/GetAllCounts.h include/Triadic.h include/TriangleEstimators.h include/EstimatorUtil.h include/BaselineEstimators.h)

add_executable(GraphProperties.out GraphProperties.cpp Graph.cpp GraphIO.cpp include/Graph.h include/GraphIO.h)

add_executable(FormatConverter.out GraphFormatConverter.cpp Graph.cpp GraphIO.cpp include/Graph.h include/GraphIO.h)