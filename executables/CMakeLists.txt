add_executable(fish_delta fish_delta.cpp)
target_link_libraries(fish_delta PUBLIC FastSphereSums)
target_include_directories(fish_delta PUBLIC
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
  $<INSTALL_INTERFACE:include>)

add_executable(gradient gradient.cpp)
target_link_libraries(gradient PUBLIC FastSphereSums)
target_include_directories(fish_delta PUBLIC
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
  $<INSTALL_INTERFACE:include>)

add_executable(inverse_laplacian inverse_laplacian.cpp)
target_link_libraries(inverse_laplacian PUBLIC FastSphereSums)
target_include_directories(inverse_laplacian PUBLIC
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
  $<INSTALL_INTERFACE:include>)

add_executable(inverse_laplacian_reg inverse_laplacian_reg.cpp)
target_link_libraries(inverse_laplacian_reg PUBLIC FastSphereSums)
target_include_directories(inverse_laplacian_reg PUBLIC
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
  $<INSTALL_INTERFACE:include>)

add_executable(laplacian laplacian.cpp)
target_link_libraries(laplacian PUBLIC FastSphereSums)
target_include_directories(laplacian PUBLIC
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
  $<INSTALL_INTERFACE:include>)

add_executable(inverse_biharmonic inverse_biharmonic.cpp)
target_link_libraries(inverse_biharmonic PUBLIC FastSphereSums)
target_include_directories(inverse_biharmonic PUBLIC
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
  $<INSTALL_INTERFACE:include>)

add_executable(sal sal.cpp)
target_link_libraries(sal PUBLIC FastSphereSums)
target_include_directories(sal PUBLIC
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
  $<INSTALL_INTERFACE:include>)

add_executable(sal_derivs sal_derivs.cpp)
target_link_libraries(sal_derivs PUBLIC FastSphereSums)
target_include_directories(sal_derivs PUBLIC
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
  $<INSTALL_INTERFACE:include>)

add_executable(icos_grid_gen icos_grid_gen.cpp)
target_link_libraries(icos_grid_gen PUBLIC FastSphereSums)
target_include_directories(icos_grid_gen PUBLIC
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
  $<INSTALL_INTERFACE:include>)

add_executable(icos_grid_gen_centroid icos_grid_gen_centroid.cpp)
target_link_libraries(icos_grid_gen_centroid PUBLIC FastSphereSums)
target_include_directories(icos_grid_gen PUBLIC
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
  $<INSTALL_INTERFACE:include>)

#install(TARGETS fast_sum
#        EXPORT FastSphereSumsTargets
#        DESTINATION ${CMAKE_INSTALL_BINDIR})
