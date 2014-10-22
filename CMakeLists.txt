cmake_minimum_required (VERSION 2.6)

project (DpBayesianNet)

SET(CMAKE_C_FLAGS "-Wall -g")

ADD_LIBRARY(dpBayesianNet STATIC dpBayesianNet.c laa.c util.c)


INSTALL(TARGETS dpBayesianNet ARCHIVE DESTINATION lib)
INSTALL(FILES dpBayesianNetwork.h DESTINATION include)

add_executable (example example.c)
target_link_libraries(example dpBayesianNet m)
