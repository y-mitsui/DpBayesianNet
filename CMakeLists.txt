cmake_minimum_required (VERSION 2.6)

project (DpBayesianNet)

SET(CMAKE_C_FLAGS "-Wall -g")

ADD_LIBRARY(dpBayesianNet SHARED dpBayesianNet.c laa.c util.c)

add_executable (example example.c dpBayesianNet.c laa.c util.c)
target_link_libraries(example m)
