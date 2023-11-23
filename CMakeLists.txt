project ("NSFiniteDifference")

cmake_minimum_required(VERSION 3.10)

# Add C++ standard (e.g., C++11)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Define path to third-party libraries
set(EIGEN3_INCLUDE_DIR C:\\Users\\luis2\\eigen-3.4.0\\eigen-3.4.0)

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /O2")

# Include directories
include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/common) 
include_directories(${CMAKE_SOURCE_DIR}/src) 


# Add your source files
add_executable(
	SOLVER 
	"src/main.cpp"
	"common/Grid.cpp" 
	"common/Connectivity.cpp" 
	 "common/LaplacianOperator.cpp" "common/RHS.cpp" "common/Field.cpp")

