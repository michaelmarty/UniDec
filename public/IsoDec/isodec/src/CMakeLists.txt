cmake_minimum_required(VERSION 3.22.1)

# Keep standalone installs beside the source directory instead of using the
# platform-wide default (for example, Program Files on Windows). Packaging
# tools can still supply their own staging prefix via CMAKE_INSTALL_PREFIX.
get_filename_component(
    ISODEC_DEFAULT_INSTALL_PREFIX "${CMAKE_CURRENT_SOURCE_DIR}/.." ABSOLUTE
)
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
    set(
        CMAKE_INSTALL_PREFIX "${ISODEC_DEFAULT_INSTALL_PREFIX}"
        CACHE PATH "Default IsoDec installation directory"
    )
endif()

project(isodec LANGUAGES C CXX)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_SHARED_LIBRARY_PREFIX "")

option(ISODEC_ENABLE_AVX2 "Compile IsoDec's x86 native core with AVX2" ON)
option(ISODEC_STATIC_GNU_RUNTIME "Statically link the GNU runtime on Linux" ON)
option(ISODEC_REQUIRE_INTELLLVM "Fail unless Intel LLVM is the active compiler" OFF)

if(ISODEC_REQUIRE_INTELLLVM AND NOT CMAKE_C_COMPILER_ID STREQUAL "IntelLLVM")
    message(FATAL_ERROR
        "ISODEC_REQUIRE_INTELLLVM is ON, but the C compiler is ${CMAKE_C_COMPILER_ID}"
    )
endif()

# IsoGen is maintained independently and pinned by the repository submodule.
# Its target supplies the isotope-distribution routines used by isodeclib.
set(ISOGEN_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../extern/IsoGen/src")
if(NOT EXISTS "${ISOGEN_SOURCE_DIR}/CMakeLists.txt")
    message(FATAL_ERROR
        "IsoGen submodule is missing. Run: git submodule update --init --recursive"
    )
endif()
add_subdirectory("${ISOGEN_SOURCE_DIR}" "${CMAKE_CURRENT_BINARY_DIR}/isogen")

add_library(isodeclib SHARED
    isodeclib.c
    isodeclib.h
    phase_model_4.h
    phase_model_8.h
)
if(APPLE)
    set(ISODEC_LIBRARY_ORIGIN "@loader_path")
else()
    set(ISODEC_LIBRARY_ORIGIN "$ORIGIN")
endif()
set_target_properties(isodeclib PROPERTIES
    LINKER_LANGUAGE C
    INSTALL_RPATH "${ISODEC_LIBRARY_ORIGIN}"
    BUILD_RPATH "${ISODEC_LIBRARY_ORIGIN}"
)
target_include_directories(isodeclib PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}"
    "${ISOGEN_SOURCE_DIR}"
)
target_compile_definitions(isodeclib PRIVATE ISODEC_BUILD_DLL)
target_link_libraries(isodeclib PRIVATE isogen)

if(MSVC)
    target_compile_definitions(isodeclib PRIVATE _CRT_SECURE_NO_WARNINGS)
elseif(UNIX)
    target_link_libraries(isodeclib PRIVATE m)
endif()

if(ISODEC_ENABLE_AVX2 AND CMAKE_SYSTEM_PROCESSOR MATCHES "^(x86_64|AMD64|amd64|i[3-6]86)$")
    if(MSVC)
        target_compile_options(isodeclib PRIVATE /arch:AVX2)
    elseif(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|IntelLLVM")
        target_compile_options(isodeclib PRIVATE -mavx2 -mfma)
    endif()
endif()

if(ISODEC_STATIC_GNU_RUNTIME AND UNIX AND NOT APPLE AND
   CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|IntelLLVM")
    target_link_options(isodeclib PRIVATE -static-libgcc)
endif()

install(
    TARGETS isodeclib
    RUNTIME DESTINATION . COMPONENT Runtime
    LIBRARY DESTINATION . COMPONENT Runtime
    ARCHIVE DESTINATION . COMPONENT Development
)

# A standalone install into the source package keeps all native artifacts next
# to the test executables. Explicit install prefixes (including wheel staging
# directories) retain the standard CMake layout expected by their caller.
install(CODE "
    if(CMAKE_INSTALL_PREFIX STREQUAL \"${ISODEC_DEFAULT_INSTALL_PREFIX}\")
        file(MAKE_DIRECTORY \"\${CMAKE_INSTALL_PREFIX}/bin\")
        file(GLOB ISODEC_NATIVE_ARTIFACTS LIST_DIRECTORIES FALSE
            \"\${CMAKE_INSTALL_PREFIX}/*.dll\"
            \"\${CMAKE_INSTALL_PREFIX}/*.lib\"
            \"\${CMAKE_INSTALL_PREFIX}/*.so\"
            \"\${CMAKE_INSTALL_PREFIX}/*.so.*\"
            \"\${CMAKE_INSTALL_PREFIX}/*.dylib\"
            \"\${CMAKE_INSTALL_PREFIX}/*.a\"
        )
        foreach(ISODEC_NATIVE_ARTIFACT IN LISTS ISODEC_NATIVE_ARTIFACTS)
            get_filename_component(ISODEC_NATIVE_NAME
                \"\${ISODEC_NATIVE_ARTIFACT}\" NAME
            )
            file(RENAME
                \"\${ISODEC_NATIVE_ARTIFACT}\"
                \"\${CMAKE_INSTALL_PREFIX}/bin/\${ISODEC_NATIVE_NAME}\"
            )
        endforeach()
    endif()
" ALL_COMPONENTS)

add_executable(isodec_test isodec_test.cpp)
target_include_directories(isodec_test PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(isodec_test PRIVATE isodeclib)
if(ISODEC_STATIC_GNU_RUNTIME AND UNIX AND NOT APPLE AND
   CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|IntelLLVM")
    target_link_options(isodec_test PRIVATE -static-libstdc++ -static-libgcc)
endif()
install(TARGETS isodec_test RUNTIME DESTINATION bin COMPONENT Development)
