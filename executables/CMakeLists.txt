add_executable(inverse_laplacian inverse_laplacian.cpp)
target_link_libraries(inverse_laplacian PUBLIC FastSphereSums)
target_include_directories(inverse_laplacian PUBLIC
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
  $<INSTALL_INTERFACE:include>)

add_executable(inverse_biharmonic inverse_biharmonic.cpp)
target_link_libraries(inverse_biharmonic PUBLIC FastSphereSums)
target_include_directories(inverse_biharmonic PUBLIC
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
  $<INSTALL_INTERFACE:include>)

add_executable(icos_grid_gen icos_grid_gen.cpp)
target_link_libraries(icos_grid_gen PUBLIC FastSphereSums)
target_include_directories(icos_grid_gen PUBLIC
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
  $<INSTALL_INTERFACE:include>)

#install(TARGETS fast_sum
#        EXPORT FastSphereSumsTargets
#        DESTINATION ${CMAKE_INSTALL_BINDIR})
