cmake_minimum_required(VERSION 3.22.1)
project(isofft C)

set(CMAKE_C_STANDARD 11)

get_filename_component(INSTALL_DIR ./ ABSOLUTE)
message(STATUS "INSTALL_DIR: ${INSTALL_DIR}")

add_library(isofft SHARED isofft.c)
set_target_properties(isofft PROPERTIES LINKER_LANGUAGE C)
add_compile_definitions(isofft "ISOFFTlib_EXPORTS")

#add_library(libisofft STATIC isofft.c)
#set_target_properties(libisofft PROPERTIES LINKER_LANGUAGE C)
#add_compile_definitions(libisofft "ISOFFTLIB_STATIC_DEFINE")

if(MSVC)
    add_library(fftw3 STATIC IMPORTED)
    set_target_properties(fftw3 PROPERTIES IMPORTED_LOCATION ${INSTALL_DIR}/fftw3.lib)
#    set(MKL_DIR "C:/Program Files (x86)/Intel/oneAPI/mkl/latest/lib/cmake/mkl")
#    message(STATUS "MKL_DIR: ${MKL_DIR}")
##    set(MKL_LINK "static")
#    find_package(MKL CONFIG REQUIRED)
#    message(STATUS "${MKL_IMPORTED_TARGETS}")
#    target_compile_options(isofft PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>)
#    target_include_directories(isofft PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)
#    target_include_directories(isofft PUBLIC "C:/Program Files (x86)/Intel/oneAPI/mkl/latest/include/fftw")
#    target_link_libraries(isofft PUBLIC $<LINK_ONLY:MKL::MKL>)

else()
    message(STATUS "Not MSVC")
    #    set(FFTW_LIB "/usr/lib/x86_64-linux-gnu/libfftw3.so")
    #target_link_libraries(isofft PUBLIC ${FFTW_LIB})
    add_library(fftw3 SHARED IMPORTED)
    set_target_properties(fftw3 PROPERTIES IMPORTED_LOCATION /usr/lib/x86_64-linux-gnu/libfftw3.so)
    target_link_libraries(isofft PUBLIC m)
endif()

target_link_libraries(isofft PUBLIC fftw3)


add_executable(isofft_test isofft_test.c)
target_link_libraries(isofft_test PUBLIC isofft)
set_target_properties(isofft_test PROPERTIES LINKER_LANGUAGE C)

install(TARGETS isofft_test DESTINATION ${INSTALL_DIR})
install(TARGETS isofft DESTINATION ${INSTALL_DIR})



