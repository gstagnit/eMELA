# build examples
file(GLOB examplecodes *.cc )
if (WITH_LHAPDF)
  foreach(examplecode ${examplecodes})
    get_filename_component(filename ${examplecode} NAME_WE)
    add_executable(${filename} ${examplecode})
    target_link_libraries(${filename} eMELA)
  endforeach()
endif()
