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

#set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,-rpath='$ORIGIN'")
set(CMAKE_INSTALL_RPATH ${INSTALL_DIR})
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)



# ======================== IsoGen Builds ========================

add_library(isogen SHARED
                isogendep.c isogendep.h
                ${FFTW_DIR}/fftw3.h
                isogenpep.c isogenpep.h
                isogenrna.c isogenrna.h
        isogen_test.cpp
)
set_target_properties(isogen PROPERTIES LINKER_LANGUAGE C)
target_include_directories(isogen PUBLIC ${FFTW_DIR})
target_compile_definitions(isogen PUBLIC ISOGEN_BUILD_DLL)


if(MSVC)
    # x64 keeps using the bundled import library exactly as before, so nothing has to be installed for
    # the existing Windows build. Other Windows architectures cannot: fftw/fftw3.lib is x86-64 only
    # (every one of its 429 COFF members is AMD64), so arm64 has to get FFTW from somewhere else.
    # vcpkg's fftw3 port supports all triplets, so "vcpkg install fftw3:arm64-windows" plus its
    # toolchain file is enough.
    if(CMAKE_SYSTEM_PROCESSOR MATCHES "^(AMD64|amd64|x86_64)$")
        set(FFTW3_LIBRARY "${FFTW_DIR}/fftw3.lib")
        message(STATUS "FFTW3_LIBRARY (bundled x64): ${FFTW3_LIBRARY}")
    else()
        find_library(FFTW3_LIBRARY NAMES fftw3 libfftw3)
        if(NOT FFTW3_LIBRARY)
            message(FATAL_ERROR
                "FFTW3 not found for ${CMAKE_SYSTEM_PROCESSOR}, and the bundled fftw/fftw3.lib is "
                "x86-64 only. Install it for this architecture and point CMake at it, for example:\n"
                "  vcpkg install fftw3:arm64-windows\n"
                "  cmake -DCMAKE_TOOLCHAIN_FILE=<vcpkg>/scripts/buildsystems/vcpkg.cmake "
                "-DVCPKG_TARGET_TRIPLET=arm64-windows ...")
        endif()
        message(STATUS "FFTW3_LIBRARY (found): ${FFTW3_LIBRARY}")
    endif()
    add_library(fftw3 STATIC IMPORTED)
    set_target_properties(fftw3 PROPERTIES IMPORTED_LOCATION ${FFTW3_LIBRARY})
else()
    message(STATUS "Not MSVC")
    message(STATUS "CMAKE_C_COMPILER_ID: ${CMAKE_C_COMPILER_ID}")
    # Locate FFTW rather than hardcoding a path. The previous absolute
    # /usr/lib/x86_64-linux-gnu/libfftw3.so is Debian-and-x86-64 only, so it misses macOS (.dylib,
    # and a Homebrew prefix that differs between Apple Silicon and Intel) and arm64 Linux
    # (/usr/lib/aarch64-linux-gnu).
    find_library(FFTW3_LIBRARY NAMES fftw3 libfftw3
        HINTS ENV FFTW3_DIR
        PATH_SUFFIXES lib)
    find_path(FFTW3_INCLUDE_DIR NAMES fftw3.h
        HINTS ENV FFTW3_DIR
        PATH_SUFFIXES include)

    if(NOT FFTW3_LIBRARY)
        message(FATAL_ERROR
            "FFTW3 not found. Install it (apt: libfftw3-dev, brew: fftw, dnf: fftw-devel) "
            "or set FFTW3_DIR to its prefix.")
    endif()

    message(STATUS "FFTW3_LIBRARY: ${FFTW3_LIBRARY}")
    add_library(fftw3 SHARED IMPORTED)
    set_target_properties(fftw3 PROPERTIES IMPORTED_LOCATION ${FFTW3_LIBRARY})
    if(FFTW3_INCLUDE_DIR)
        target_include_directories(isogen PUBLIC ${FFTW3_INCLUDE_DIR})
    endif()
    target_link_libraries(isogen PUBLIC m)
endif()

target_link_libraries(isogen PUBLIC fftw3)
install(TARGETS isogen DESTINATION ${INSTALL_DIR})


# ======================== IsoDec Builds ========================

add_library(isodeclib SHARED isodeclib.c isodeclib.h phase_model_4.h phase_model_8.h isogenpep.h)
set_target_properties(isodeclib PROPERTIES LINKER_LANGUAGE C)

target_link_libraries(isodeclib PUBLIC isogen)

target_include_directories(isodeclib PUBLIC ${SRC_DIR})

target_compile_definitions(isodeclib PRIVATE ISODEC_BUILD_DLL)


# Compiler flags
include(find_avx.cmake)
CHECK_FOR_AVX()
if(HAVE_AVX2_EXTENSIONS)
    message(STATUS "Using AVX2")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${AVX_FLAGS}")
else()
    message(STATUS "Not Using AVX2")
    # Every arm64 build lands here, since there is no AVX2 to detect. OpenMP is optional rather than
    # required because Apple Clang does not ship it: an unrecognised "#pragma omp" is ignored, so the
    # library still builds and runs, with that one loop single-threaded.
    find_package(OpenMP)
    if(OpenMP_C_FOUND)
        message(STATUS "OpenMP found")
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    else()
        message(STATUS "OpenMP not found - building single-threaded")
    endif()
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3")
endif()

install(TARGETS isodeclib DESTINATION ${INSTALL_DIR})

add_executable(isodec_test isodec_test.cpp)
target_include_directories(isodec_test PUBLIC ${SRC_DIR})
target_link_libraries(isodec_test PUBLIC isogen isodeclib)

add_executable(isogen_test isogen_test.cpp)
target_include_directories(isogen_test PUBLIC ${SRC_DIR})
target_link_libraries(isogen_test PUBLIC isogen isogen)

install(TARGETS isodec_test DESTINATION ${INSTALL_DIR})









