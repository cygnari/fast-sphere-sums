add_libary(FastSphereSums
    mpi_utils.cpp)

target_link_libraries(FastSphereSums PUBLIC MPI::MPI_C LAPACK::LAPACK)
target_include_directories(FastSphereSums PUBLIC
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  $<INSTALL_INTERFACE:include>)

install(FILES
    # all FastSphereSums header/include files go here
    ${CMAKE_BINARY_DIR}/fast-sphere-sums-config.h
    mpi_utils.hpp
    structs.hpp
    DESTINATION ${PROJECT_SOURCE_DIR}/include)

#install(TARGETS FastSphereSums
#        EXPORT FastSphereSumsTargets
#        DESTINATION ${CMAKE_INSTALL_LIBDIR})
