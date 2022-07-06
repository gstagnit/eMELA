# build examples
file(GLOB examplecodes *.cc )
foreach(examplesource ${examplecodes})
  GET_FILENAME_COMPONENT(filename ${examplesource} NAME_WE)
  add_executable(${filename} ${filename}.cc)
  set_target_properties(${filename} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/example")
  target_link_libraries(${filename} eMELA)
endforeach()
