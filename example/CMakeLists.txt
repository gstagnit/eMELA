# build examples
file(GLOB examplecodes *.cc )
if (WITH_LHAPDF)
  foreach(examplecode ${examplecodes})
    get_filename_component(examplename ${examplecode} NAME_WE)
    add_executable(${examplename} ${examplecode})
    target_link_libraries(${examplename} eMELA)
  endforeach()
endif()
