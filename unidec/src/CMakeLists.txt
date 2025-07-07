cmake_minimum_required(VERSION 3.31)
project(UniDec C)

set(CMAKE_C_STANDARD 11)

if(MSVC)
    add_definitions(-D_CRT_SECURE_NO_WARNINGS)
endif()

#include_directories
# Absolute path
#set(UNIDEC_DIR "${CMAKE_CURRENT_SOURCE_DIR}../")
cmake_path(GET CMAKE_CURRENT_SOURCE_DIR PARENT_PATH UNIDEC_DIR)

set(FFTW_DIR "${UNIDEC_DIR}/IsoDec/src_cmake/fftw")
set(INSTALL_DIR "${UNIDEC_DIR}/bin")

# Set the HDF5 directory
set(HDF5_DIR "C:/Program Files/HDF_Group/HDF5/1.14.6/")

message(STATUS "UniDec Dir: ${UNIDEC_DIR}")
message(STATUS "FFTW_DIR: ${FFTW_DIR}")
message(STATUS "HDF5_DIR: ${HDF5_DIR}")
message(STATUS "Install Dir: ${INSTALL_DIR}")

# Add FFTW3 and HDF5 directories to the include path and link directories
include_directories("${FFTW_DIR}")
link_directories("${FFTW_DIR}")

include_directories("${HDF5_DIR}/include")
link_directories("${HDF5_DIR}/lib")

#add_library(unideclib INTERFACE)
#target_include_directories(unideclib INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
add_executable(unidec
        "${FFTW_DIR}/fftw3.h"
        UniDec.h UD_conv.h UD_dataproc.h UD_H5_IO.h
        UD_IO.h
        UD_sg.h UD_peak_width.h UD_score.h
        UD_match.h UniDec.c UD_match.h
        UniDec_Main.h UniDecLC_Main.h UniDecIM.h UniDecIM_Main.h
        udtools.c
        udtools.h
        udstruct.c
        udstruct.h)

#target_link_libraries(unidec PRIVATE unideclib)
#target_include_directories(unideclib INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
#target_compile_definitions(unidec PRIVATE UNIDEC_BUILD_DLL)

target_link_libraries(unidec fftw3)
target_link_libraries(unidec hdf5 hdf5_hl)


find_package(OpenMP)
if (OpenMP_C_FOUND)
    target_link_libraries(unidec OpenMP::OpenMP_C)
else()
    message(STATUS "OpenMP not found, proceeding without it.")
endif()

# Optimization
target_compile_options(unidec PRIVATE
        $<$<AND:$<C_COMPILER_ID:MSVC>,$<CONFIG:Release>>:/O3>
        $<$<AND:$<NOT:$<C_COMPILER_ID:MSVC>>,$<CONFIG:Release>>:-O3>
)

# Install the executable
install(TARGETS unidec
        RUNTIME DESTINATION ${INSTALL_DIR}
        LIBRARY DESTINATION ${INSTALL_DIR}
        ARCHIVE DESTINATION ${INSTALL_DIR})

