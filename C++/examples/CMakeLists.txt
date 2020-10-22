# SE-Sync command-line driver
add_executable(SE-Sync main.cpp)
target_link_libraries(SE-Sync SESync)
if(${ENABLE_PROFILING})
target_compile_definitions(SE-Sync PRIVATE GPERFTOOLS)
target_link_libraries(SE-Sync profiler)
endif()


message(STATUS "Building main SE-Sync command-line executable in directory ${EXECUTABLE_OUTPUT_PATH}\n")


# SE-Sync visualizer
if(${ENABLE_VISUALIZATION})
  add_executable(SE-SyncViz mainviz.cpp)
  target_link_libraries(SE-SyncViz SESync SESyncViz stdc++fs)
endif()
