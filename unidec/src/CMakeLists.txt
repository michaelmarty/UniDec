cmake_minimum_required(VERSION 3.25)
project(UniDec C)

set(CMAKE_C_STANDARD 11)

# Removes annoying warnings to define this
if(MSVC)
    add_definitions(-D_CRT_SECURE_NO_WARNINGS)
endif()

# Set the UniDec directory
#set(UNIDEC_DIR "${CMAKE_CURRENT_SOURCE_DIR}../")
cmake_path(GET CMAKE_CURRENT_SOURCE_DIR PARENT_PATH UNIDEC_DIR)
message(STATUS "UniDec Dir: ${UNIDEC_DIR}")

# Set the installation directory
set(INSTALL_DIR "${UNIDEC_DIR}/bin")
message(STATUS "Install Dir: ${INSTALL_DIR}")

# Add HDF5 directories to the include path and link directories
if(MSVC)
    set(HDF5_DIR "C:/Program Files/HDF_Group/HDF5/1.14.6/")
    message(STATUS "HDF5_DIR: ${HDF5_DIR}")
    include_directories("${HDF5_DIR}/include")
    link_directories("${HDF5_DIR}/lib")
else()
    set(HDF5_INCLUDE_DIR "/usr/include/hdf5/serial")
    message(STATUS "HDF5_DIR: ${HDF5_INCLUDE_DIR}")
    include_directories("${HDF5_INCLUDE_DIR}")
#    set(HDF5_LINK_DIR "/usr/lib/x86_64-linux-gnu/hdf5/serial")
#    link_directories("${HDF5_LINK_DIR}")
    link_directories("/usr/lib/x86_64-linux-gnu")
    link_directories("/usr/lib/x86_64-linux-gnu/hdf5/serial/lib")


endif()

add_library(unideclib SHARED udmain.c udmain.h
        udtools.c
        udtools.h
        udstruct.c
        udstruct.h
        udio.c
        udio.h
        udcore.c
        udcore.h
        h5io.c
        h5io.h
        UD_score.c
        UD_score.h
        UD_analysis.c
        UD_analysis.h
        UD_dataproc.c
        UD_dataproc.h
        UD_conv.h
        UD_conv.c
        udinterface.cpp
        udinterface.hpp
)
#target_include_directories(unideclib INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
add_executable(unidec
        UniDec.c
        UniDec.h
        udtools.c
        udtools.h
        udstruct.c
        udstruct.h
        udio.c
        udio.h
        udmain.c
        udmain.h
        udcore.c
        udcore.h
        h5io.c
        h5io.h
        UD_dataproc.c
        UD_dataproc.h
        UD_analysis.c
        UD_analysis.h
        UD_score.c
        UD_score.h
        UD_conv.h
        UD_conv.c
        UniDecIM.c
        UniDecIM.h
        UniDecIM_Main.c
        UniDecIM_Main.h
        UD_charge.h
        UD_charge.c
        UniDecCD_Main.c
        UniDecCD_Main.h
        MetaUniDec_Main.c
        MetaUniDec_Main.h
)

#target_link_libraries(unidec PRIVATE unideclib)
#target_include_directories(unideclib INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_definitions(unideclib PRIVATE UNIDEC_BUILD_DLL)
set_target_properties(unideclib PROPERTIES LINKER_LANGUAGE C)

# Set linking for HDF5 libraries
target_link_libraries(unidec PUBLIC hdf5 hdf5_hl)
target_link_libraries(unideclib PUBLIC hdf5 hdf5_hl)
# IF not msvc add this
if(NOT MSVC)
    target_link_libraries(unidec PUBLIC m)
    target_link_libraries(unideclib PUBLIC m)
    # If linux
    if (UNIX AND NOT APPLE)
        set_target_properties(unidec PROPERTIES OUTPUT_NAME "unideclinux")
    elseif (APPLE)
        set_target_properties(unidec PROPERTIES OUTPUT_NAME "unidecmac")
    endif()
endif()

#target_link_libraries(unidec PUBLIC unideclib)

# Set linking for FFTW libraries
# Find MKL (Math Kernel Library)
find_package(MKL)
if (MKL_FOUND)
    target_link_libraries(unidec PUBLIC ${MKL_LIBRARIES})
    target_link_libraries(unidec PUBLIC mkl_rt)
    target_link_libraries(unideclib PUBLIC mkl_rt)
    include_directories(${MKL_INCLUDE_DIRS})
    set(FFTW_DIR "C:/Program Files (x86)/Intel/oneAPI/mkl/2025.1/include/fftw")
    include_directories("${FFTW_DIR}")
    message(STATUS "MKL Include: ${MKL_INCLUDE_DIRS}")
else()
    message("MKL not found. Looking for FFTW elsewhere.")
    set(FFTW_DIR "${UNIDEC_DIR}/IsoDec/src_cmake/fftw")

    message(STATUS "FFTW_DIR: ${FFTW_DIR}")
    include_directories("${FFTW_DIR}")
    link_directories("${FFTW_DIR}")
    target_link_libraries(unidec PUBLIC fftw3f)
    target_link_libraries(unideclib PUBLIC fftw3f)
endif()


find_package(OpenMP)
if (OpenMP_C_FOUND)
    target_link_libraries(unidec PUBLIC OpenMP::OpenMP_C)
    target_link_libraries(unideclib PUBLIC OpenMP::OpenMP_C)
else()
    message(STATUS "OpenMP not found, proceeding without it.")
endif()

# Optimization
target_compile_options(unidec PRIVATE
        $<$<AND:$<C_COMPILER_ID:MSVC>,$<CONFIG:Release>>:/O3>
        $<$<AND:$<NOT:$<C_COMPILER_ID:MSVC>>,$<CONFIG:Release>>:-O3>
)
# Optimization
target_compile_options(unideclib PRIVATE
        $<$<AND:$<C_COMPILER_ID:MSVC>,$<CONFIG:Release>>:/O3>
        $<$<AND:$<NOT:$<C_COMPILER_ID:MSVC>>,$<CONFIG:Release>>:-O3>
)


# Install the executable
install(TARGETS unidec unideclib
        RUNTIME DESTINATION ${INSTALL_DIR}
        LIBRARY DESTINATION ${INSTALL_DIR}
        ARCHIVE DESTINATION ${INSTALL_DIR})

