add_executable(fast_sum fast_sum.cpp)
target_link_libraries(fast_sum PUBLIC FastSphereSums)
target_include_directories(fast_sum PUBLIC
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
  $<INSTALL_INTERFACE:include>)

add_executable(icos_grid_gen icos_grid_gen.cpp)
target_link_libraries(icos_grid_gen PUBLIC FastSphereSums)
target_include_directories(icos_grid_gen PUBLIC
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
  $<INSTALL_INTERFACE:include>)

add_executable(direct_sum direct_sum.cpp)
target_link_libraries(direct_sum PUBLIC FastSphereSums)
target_include_directories(direct_sum PUBLIC
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
  $<INSTALL_INTERFACE:include>)

#install(TARGETS fast_sum
#        EXPORT FastSphereSumsTargets
#        DESTINATION ${CMAKE_INSTALL_BINDIR})
