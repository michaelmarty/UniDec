# ======================== Setup ========================
cmake_minimum_required(VERSION 3.22.1)
include(GenerateExportHeader)

project(isodeclib C CXX)

set(CMAKE_C_STANDARD 99)
set(CMAKE_SHARED_LIBRARY_PREFIX "")

get_filename_component(INSTALL_DIR ../ ABSOLUTE)
message(STATUS "INSTALL_DIR: ${INSTALL_DIR}")

# Set SRC_DIR to the path of the IsoDec source code
get_filename_component(SRC_DIR ./ ABSOLUTE)
message(STATUS "SRC_DIR: ${SRC_DIR}")

# FFTW directory is src_dir/fftw3
set(FFTW_DIR "${SRC_DIR}/fftw")
message(STATUS "FFTW_DIR: ${FFTW_DIR}")

include_directories("${SRC_DIR}")
include_directories("${FFTW_DIR}")

# ======================== IsoGen Mass Builds ========================

add_library(isogenmass SHARED isogenmass.c isogenmass.h)
set_target_properties(isogenmass PROPERTIES LINKER_LANGUAGE C)

add_executable(isogen isogen.cpp)
target_link_libraries(isogen PUBLIC isogenmass)
set_target_properties(isogen PROPERTIES LINKER_LANGUAGE C)

install(TARGETS isogenmass DESTINATION ${INSTALL_DIR})


# ======================== IsoFFT Builds ========================

add_library(isofft SHARED isofft.c isofft.h ${FFTW_DIR}/fftw3.h)
set_target_properties(isofft PROPERTIES LINKER_LANGUAGE C)
target_include_directories(isofft PUBLIC ${FFTW_DIR})
add_compile_definitions(isofft "ISOFFTlib_EXPORTS")

if(MSVC OR CMAKE_C_COMPILER_ID STREQUAL "GNU")
    add_library(fftw3 STATIC IMPORTED)
    set_target_properties(fftw3 PROPERTIES IMPORTED_LOCATION ${FFTW_DIR}/fftw3.lib)
else()
    message(STATUS "Not MSVC")
    message(STATUS "CMAKE_C_COMPILER_ID: ${CMAKE_C_COMPILER_ID}")
    #    set(FFTW_LIB "/usr/lib/x86_64-linux-gnu/libfftw3.so")
    #target_link_libraries(isofft PUBLIC ${FFTW_LIB})
    add_library(fftw3 SHARED IMPORTED)
    set_target_properties(fftw3 PROPERTIES IMPORTED_LOCATION /usr/lib/x86_64-linux-gnu/libfftw3.so)
    target_link_libraries(isofft PUBLIC m)
endif()

#add_library(fftw3 STATIC IMPORTED)
#set_target_properties(fftw3 PROPERTIES IMPORTED_LOCATION ${FFTW_DIR}/fftw3.lib)
target_link_libraries(isofft PUBLIC fftw3)

# isofft_test
add_executable(isofft_test isofft_test.c)
target_link_libraries(isofft_test PUBLIC isofft)
set_target_properties(isofft_test PROPERTIES LINKER_LANGUAGE C)

#install(TARGETS isofft_test DESTINATION ${INSTALL_DIR})
install(TARGETS isofft DESTINATION ${INSTALL_DIR})


# ======================== IsoDec Builds ========================

add_library(isodeclib SHARED isodeclib.c isodeclib.h phase_model_4.h phase_model_8.h)
set_target_properties(isodeclib PROPERTIES LINKER_LANGUAGE C)
#target_link_libraries(isodeclib PUBLIC isogenmass)
target_link_libraries(isodeclib PUBLIC isofft)

# Compiler flags
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

install(TARGETS isodeclib DESTINATION ${INSTALL_DIR})

# isodec_test
add_executable(isodec_test isodec_test.cpp)
target_link_libraries(isodec_test PUBLIC isodeclib)
set_target_properties(isodec_test PROPERTIES LINKER_LANGUAGE C)


#generate_export_header(isodeclib
#        BASE_NAME isodeclib
#        EXPORT_MACRO_NAME ISODEC_EXPORTS
#        EXPORT_FILE_NAME isodec_exports.h
#        DEFINE_NO_DEPRECATED
#)










