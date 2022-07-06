# define minimum version of cmake
cmake_minimum_required (VERSION 3.5)

# define project name and its language
project(eMELA CXX Fortran)

# define c++ standard
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}  -fPIC -Wunused")
set(CMAKE_SHARED_LINKER_FLAGS -w)

if (NOT DEFINED CMAKE_MACOSX_RPATH)
   set(CMAKE_MACOSX_RPATH 0)
endif()

# set version
set(eMELA_VERSION 1.0)

# export version to file
configure_file(
  "${PROJECT_SOURCE_DIR}/inc/eMELA/version.h.in"
  "${PROJECT_SOURCE_DIR}/inc/eMELA/version.h"
)

# Configuration script
set(prefix ${CMAKE_INSTALL_PREFIX})
set(exec_prefix "${prefix}")
set(includedir "${prefix}/include")
set(libdir "${prefix}/lib")
configure_file(
  "${PROJECT_SOURCE_DIR}/bin/eMELA-config.in"
  "${PROJECT_SOURCE_DIR}/bin/eMELA-config"
)

# include directory
include_directories(${PROJECT_SOURCE_DIR}/inc)

# generate list of source files
file(GLOB_RECURSE  source_files src/*f src/*cc)

# build test codes
add_subdirectory(example)

# LHAPDF optional
option(WITH_LHAPDF "Use LHAPDF" OFF)

if(WITH_LHAPDF)
  message("-- LHAPDF enabled")
  add_definitions(-DUSE_LHAPDF)
  find_program(LHAPDF_CONFIG lhapdf-config REQUIRED)
  if (LHAPDF_CONFIG)
     exec_program(${LHAPDF_CONFIG}
     ARGS --cflags
     OUTPUT_VARIABLE LHAPDF_CXX_FLAGS
     )
     set(LHAPDF_CXX_FLAGS ${LHAPDF_CXX_FLAGS} CACHE STRING INTERNAL)
     exec_program(${LHAPDF_CONFIG}
     ARGS --libs
     OUTPUT_VARIABLE LHAPDF_LIBRARIES
     )
     set(LHAPDF_LIBRARIES ${LHAPDF_LIBRARIES} CACHE STRING INTERNAL)
  endif(LHAPDF_CONFIG)
endif(WITH_LHAPDF)

if(WITH_LHAPDF)
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${LHAPDF_CXX_FLAGS} -std=c++11 -Wall -Wextra -Wpedantic -fPIC")
else()
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -Wextra -Wpedantic -fPIC")
endif()

# build target library
# add_subdirectory(src)
# add_subdirectory(run)

option(SHARED "Build shared-libray instead of static-libray" ON)
if(SHARED)
  add_library(eMELA SHARED ${source_files})
else(SHARED)
  add_library(eMELA STATIC ${source_files})
endif(SHARED)

target_link_libraries(eMELA ${LHAPDF_LIBRARIES} ${MELA_LIBRARIES})
install(DIRECTORY ${PROJECT_SOURCE_DIR}/inc/eMELA DESTINATION include)
install(TARGETS eMELA DESTINATION lib)

install(FILES ${PROJECT_SOURCE_DIR}/bin/eMELA-config
	   DESTINATION bin
	   PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

# install(DIRECTORY ${PROJECT_SOURCE_DIR}/grids/ DESTINATION share/eMELA)
