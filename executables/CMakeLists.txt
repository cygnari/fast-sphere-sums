add_executable(fast_sum fast_sum.cpp)
target_link_libraries(fast_sum PUBLIC FastSphereSums)
target_include_directories(fast_sum PUBLIC
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
  $<INSTALL_INTERFACE:include>)

#install(TARGETS fast_sum
#        EXPORT FastSphereSumsTargets
#        DESTINATION ${CMAKE_INSTALL_BINDIR})
