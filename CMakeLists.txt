cmake_minimum_required(VERSION 3.0)
project(MyProject)

# Add C++11 support
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add your source files
set(SOURCES main.cpp)

# Create an executable target
add_executable(my_app ${SOURCES})
