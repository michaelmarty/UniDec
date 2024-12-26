cmake_minimum_required(VERSION 3.22.1)
include(GenerateExportHeader)

project(isodeclib C CXX)

set(CMAKE_C_STANDARD 99)
set(CMAKE_SHARED_LIBRARY_PREFIX "")

get_filename_component(INSTALL_DIR ../ ABSOLUTE)
message(STATUS "INSTALL_DIR: ${INSTALL_DIR}")


include_directories("C:/Python/UniDec3/unidec/IsoDec/src_cmake")

add_library(isogenmass SHARED isogenmass.c isogenmass.h)
set_target_properties(isogenmass PROPERTIES LINKER_LANGUAGE C)

add_executable(isogen isogen.cpp)
target_link_libraries(isogen PUBLIC isogenmass)
set_target_properties(isogen PROPERTIES LINKER_LANGUAGE C)

install(TARGETS isogenmass DESTINATION ${INSTALL_DIR})


add_library(isodeclib SHARED isodeclib.c isodeclib.h phase_model_4.h phase_model_8.h)
set_target_properties(isodeclib PROPERTIES LINKER_LANGUAGE C)
target_link_libraries(isodeclib PUBLIC isogenmass)

include(find_avx.cmake)
CHECK_FOR_AVX()


if(HAVE_AVX2_EXTENSIONS)
    message(STATUS "Using AVX2")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${AVX_FLAGS}")
else()
    message(STATUS "Not Using AVX2")
    find_package(OpenMP REQUIRED)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3")
endif()

#generate_export_header(isodeclib
#        BASE_NAME isodeclib
#        EXPORT_MACRO_NAME ISODEC_EXPORTS
#        EXPORT_FILE_NAME isodec_exports.h
#        DEFINE_NO_DEPRECATED
#)

install(TARGETS isodeclib DESTINATION ${INSTALL_DIR})

include_directories("${INSTALL_DIR}")
add_executable(isodec_test isodec_test.cpp)
target_link_libraries(isodec_test PUBLIC isodeclib)
set_target_properties(isodec_test PROPERTIES LINKER_LANGUAGE C)






