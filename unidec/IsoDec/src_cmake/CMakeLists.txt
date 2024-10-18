cmake_minimum_required(VERSION 3.22.1)
include(GenerateExportHeader)

project(isodeclib C CXX)

set(CMAKE_C_STANDARD 99)
set(CMAKE_SHARED_LIBRARY_PREFIX "")

get_filename_component(INSTALL_DIR ../ ABSOLUTE)
message(STATUS "INSTALL_DIR: ${INSTALL_DIR}")

add_library(isogenmass SHARED isogenmass.c isogenmass.h)
set_target_properties(isogenmass PROPERTIES LINKER_LANGUAGE C)

add_executable(isogen isogen.cpp)
target_link_libraries(isogen PUBLIC isogenmass)
set_target_properties(isogen PROPERTIES LINKER_LANGUAGE C)

install(TARGETS isogenmass DESTINATION ${INSTALL_DIR})


add_library(isodeclib SHARED isodeclib.c isodeclib.h phase_model_4.h phase_model_8.h)
set_target_properties(isodeclib PROPERTIES LINKER_LANGUAGE C)
target_link_libraries(isodeclib PUBLIC isogenmass)

if(MSVC)
    add_definitions(/arch:AVX2)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /arch:AVX2")
#    set(MKL_DIR "C:/Program Files (x86)/Intel/oneAPI/mkl/latest/lib/cmake/mkl")
    #    find_package(MKL CONFIG REQUIRED)
#    message(STATUS "${MKL_IMPORTED_TARGETS}")
#    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -openmp:experimental")
#    target_compile_options(isodeclib PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>)
#    target_include_directories(isodeclib PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)
#    target_link_libraries(isodeclib PUBLIC $<LINK_ONLY:MKL::MKL>)
#    target_link_libraries(isodeclib PUBLIC OpenMP::OpenMP_C)
else()
    find_package(OpenMP REQUIRED)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mavx")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3")
endif()
#set_target_properties(isodeclib PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${LIBRARY_OUTPUT_PATH})
#set_target_properties(isodeclib PROPERTIES LIBRARY_OUTPUT_DIRECTORY_DEBUG ${LIBRARY_OUTPUT_PATH} )
#set_target_properties(isodeclib PROPERTIES LIBRARY_OUTPUT_DIRECTORY_RELEASE ${LIBRARY_OUTPUT_PATH} )
#
#message(STATUS "LIBRARY_OUTPUT_PATH: ${LIBRARY_OUTPUT_PATH}")
#message(STATUS "LIBRARY_OUTPUT_DIRECTORY_RELEASE: ${LIBRARY_OUTPUT_DIRECTORY_RELEASE}")
#message(STATUS "CMAKE_LIBRARY_OUTPUT_DIRECTORY: ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")
##
#generate_export_header(isodeclib
#        BASE_NAME isodeclib
#        EXPORT_MACRO_NAME ISODEC_EXPORTS
#        EXPORT_FILE_NAME isodec_exports.h
#        DEFINE_NO_DEPRECATED
#)


install(TARGETS isodeclib DESTINATION ${INSTALL_DIR})

#find_library(ISODECLIB isodeclib PATHS ${INSTALL_DIR} ${CMAKE_CURRENT_BINARY_DIR} ${LIBRARY_OUTPUT_PATH} ${LIBRARY_OUTPUT_DIRECTORY_RELEASE})
#message(STATUS "ISODECLIB: ${ISODECLIB}")
#set(LIBRARY_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR})
#message(STATUS "CMAKE_CURRENT_BINARY_DIR: ${CMAKE_CURRENT_BINARY_DIR}")
include_directories("${INSTALL_DIR}")
add_executable(isodec_test isodec_test.cpp)
#target_link_libraries(isodec_test PUBLIC ${ISODECLIB})
target_link_libraries(isodec_test PUBLIC isodeclib)
set_target_properties(isodec_test PROPERTIES LINKER_LANGUAGE C)






