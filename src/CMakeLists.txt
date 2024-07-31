add_library(FastSphereSums
    direct_sum_funcs.cpp
    ./fast_sum_interactions/inverse_laplacian_interactions.cpp
    ./fast_sum_interactions/inverse_biharmonic_interactions.cpp
    ./fast_sum_interactions/laplacian_interactions.cpp
    ./fast_sum_interactions/sal_interactions.cpp
    ./fast_sum_interactions/sal_interactions_lat_deriv.cpp
    ./fast_sum_interactions/sal_interactions_lon_deriv.cpp
    fast_sum_funcs.cpp
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
    direct_sum_funcs.hpp
    ./fast_sum_interactions/inverse_laplacian_interactions.hpp
    ./fast_sum_interactions/inverse_biharmonic_interactions.hpp
    ./fast_sum_interactions/laplacian_interactions.hpp
    ./fast_sum_interactions/sal_interactions.hpp
    ./fast_sum_interactions/sal_interactions_lat_deriv.hpp
    ./fast_sum_interactions/sal_interactions_lon_deriv.hpp
    fast_sum_funcs.hpp
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
