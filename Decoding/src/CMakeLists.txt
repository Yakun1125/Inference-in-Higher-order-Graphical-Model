cmake_minimum_required(VERSION 3.10)
project(LDPC_Decoder)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

###########################
# Configure your own gurobi path
###########################
set(GUROBI_HOME "")

set(SOURCES
    main.cpp
    Generator.cpp
    Optimizer.cpp
)

set(HEADERS
    Generator.h
    Optimizer.h
)

add_executable(${PROJECT_NAME} ${SOURCES} ${HEADERS})
target_include_directories(${PROJECT_NAME} PRIVATE 
    include
    ${GUROBI_HOME}/include
)
target_link_libraries(${PROJECT_NAME} PRIVATE 
    ${GUROBI_HOME}/lib/libgurobi_c++.a
    ${GUROBI_HOME}/lib/libgurobi110.so
)