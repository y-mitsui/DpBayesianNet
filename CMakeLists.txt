cmake_minimum_required (VERSION 2.6)

project (DpBayesianNet)

add_executable (example example.c dpBayesianNet.c laa.c)
target_link_libraries(example m)
