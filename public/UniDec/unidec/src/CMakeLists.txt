cmake_minimum_required(VERSION 3.25)
project(UniDec C)

set(CMAKE_C_STANDARD 11)

option(UNIDEC_REQUIRE_INTELLLVM "Fail configuration unless IntelLLVM is the C compiler" OFF)
if(UNIDEC_REQUIRE_INTELLLVM AND NOT CMAKE_C_COMPILER_ID STREQUAL "IntelLLVM")
    message(FATAL_ERROR
            "UNIDEC_REQUIRE_INTELLLVM is enabled, but CMake selected ${CMAKE_C_COMPILER_ID} "
            "(${CMAKE_C_COMPILER}).")
endif()

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

# Discover HDF5 through CMake so Homebrew, Linux distributions, and custom
# prefixes all work without hard-coded include or library directories.
if(MSVC AND NOT DEFINED HDF5_ROOT AND EXISTS "C:/Program Files/HDF_Group/HDF5/2.1.1")
    set(HDF5_ROOT "C:/Program Files/HDF_Group/HDF5/2.1.1" CACHE PATH "HDF5 installation directory")
endif()
find_package(HDF5 REQUIRED COMPONENTS C HL)
message(STATUS "HDF5 version: ${HDF5_VERSION}")
message(STATUS "HDF5 include directories: ${HDF5_INCLUDE_DIRS}")

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
        UCCD_Main.c
        UCCD_Main.h
        MetaUniDec_Main.c
        MetaUniDec_Main.h
)

target_include_directories(unidec PRIVATE ${HDF5_INCLUDE_DIRS})
target_include_directories(unideclib PRIVATE ${HDF5_INCLUDE_DIRS})

# The Python package launches the executable from unidec/bin. Keep that copy
# synchronized with normal IDE/CMake builds; the install target is not run by
# CLion's regular Build action.
add_custom_command(TARGET unidec POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E make_directory "${INSTALL_DIR}"
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
                "$<TARGET_FILE:unidec>"
                "${INSTALL_DIR}/$<TARGET_FILE_NAME:unidec>"
        COMMENT "Deploying unidec executable to ${INSTALL_DIR}"
)

# The ctypes wrapper loads the shared library from unidec/bin. Deploy it next
# to the platform executable for normal IDE builds and CI builds alike.
add_custom_command(TARGET unideclib POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E make_directory "${INSTALL_DIR}"
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
                "$<TARGET_FILE:unideclib>"
                "${INSTALL_DIR}/$<TARGET_FILE_NAME:unideclib>"
        COMMENT "Deploying UniDec shared library to ${INSTALL_DIR}"
)

#target_link_libraries(unidec PRIVATE unideclib)
#target_include_directories(unideclib INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_definitions(unideclib PRIVATE UNIDEC_BUILD_DLL)
set_target_properties(unideclib PROPERTIES LINKER_LANGUAGE C)

# Set linking for HDF5 libraries
target_link_libraries(unidec PUBLIC ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES})
target_link_libraries(unideclib PUBLIC ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES})
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
find_package(MKL QUIET)
if (MKL_FOUND)
    target_compile_definitions(unidec PRIVATE UNIDEC_USE_MKL)
    target_include_directories(unidec PRIVATE "${MKL_ROOT}/include")
    target_link_libraries(unidec PUBLIC ${MKL_LIBRARIES})
    if (WIN32)
        target_link_libraries(unidec PUBLIC "${MKL_ROOT}/lib/mkl_rt.lib")
        target_link_libraries(unideclib PUBLIC "${MKL_ROOT}/lib/mkl_rt.lib")
    else()
        target_link_libraries(unidec PUBLIC mkl_rt)
        target_link_libraries(unideclib PUBLIC mkl_rt)
    endif()
    include_directories(${MKL_INCLUDE_DIRS})
    set(FFTW_DIR "${MKL_ROOT}/include/fftw")
    include_directories("${FFTW_DIR}")
    message(STATUS "MKL Include: ${MKL_INCLUDE_DIRS}")
else()
    message("MKL not found. Looking for FFTW elsewhere.")
    find_package(PkgConfig QUIET)
    if(PkgConfig_FOUND)
        pkg_check_modules(FFTW3F QUIET IMPORTED_TARGET fftw3f)
    endif()

    if(TARGET PkgConfig::FFTW3F)
        target_link_libraries(unidec PUBLIC PkgConfig::FFTW3F)
        target_link_libraries(unideclib PUBLIC PkgConfig::FFTW3F)
    else()
        find_path(FFTW3F_INCLUDE_DIR NAMES fftw3.h HINTS "${FFTW_ROOT}" ENV FFTW_ROOT PATH_SUFFIXES include)
        find_library(FFTW3F_LIBRARY NAMES fftw3f HINTS "${FFTW_ROOT}" ENV FFTW_ROOT PATH_SUFFIXES lib lib64)
        if(NOT FFTW3F_INCLUDE_DIR OR NOT FFTW3F_LIBRARY)
            message(FATAL_ERROR "FFTW single-precision library not found. Install fftw3/fftw or set FFTW_ROOT.")
        endif()
        target_include_directories(unidec PRIVATE "${FFTW3F_INCLUDE_DIR}")
        target_include_directories(unideclib PRIVATE "${FFTW3F_INCLUDE_DIR}")
        target_link_libraries(unidec PUBLIC "${FFTW3F_LIBRARY}")
        target_link_libraries(unideclib PUBLIC "${FFTW3F_LIBRARY}")
    endif()

    find_library(FFTW3F_THREADS_LIBRARY NAMES fftw3f_threads HINTS "${FFTW_ROOT}" ENV FFTW_ROOT PATH_SUFFIXES lib lib64)
    if (FFTW3F_THREADS_LIBRARY)
        target_compile_definitions(unidec PRIVATE UNIDEC_USE_FFTW_THREADS)
        target_link_libraries(unidec PUBLIC ${FFTW3F_THREADS_LIBRARY})
    endif()
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

