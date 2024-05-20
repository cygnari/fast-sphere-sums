add_library(FastSphereSums
    general_utils.cpp
    initial_conditions.cpp
    initialize_tree.cpp
    interp_utils.cpp
    io_utils.cpp
    mpi_utils.cpp)

target_link_libraries(FastSphereSums PUBLIC MPI::MPI_C LAPACK::LAPACK)
target_include_directories(FastSphereSums PUBLIC
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  $<INSTALL_INTERFACE:include>)

install(FILES
    # all FastSphereSums header/include files go here
    ${CMAKE_BINARY_DIR}/fast-sphere-sums-config.h
    initial_conditions.hpp
    initialize_tree.hpp
    interp_utils.hpp
    io_utils.hpp
    mpi_utils.hpp
    structs.hpp
    DESTINATION ${PROJECT_SOURCE_DIR}/include)

#install(TARGETS FastSphereSums
#        EXPORT FastSphereSumsTargets
#        DESTINATION ${CMAKE_INSTALL_LIBDIR})
