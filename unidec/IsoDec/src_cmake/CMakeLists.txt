cmake_minimum_required(VERSION 3.29)
include(GenerateExportHeader)

project(isodeclib C)

set(CMAKE_C_STANDARD 11)

include_directories("${CMAKE_CURRENT_BINARY_DIR}")
add_library(isodeclib SHARED isodeclib.c isodeclib.h isodec_exports.h)
generate_export_header(isodeclib
        BASE_NAME isodeclib
        EXPORT_MACRO_NAME ISODEC_EXPORTS
        EXPORT_FILE_NAME isodec_exports.h
        STATIC_DEFINE SHARED_EXPORTS_BUILT_AS_STATIC)

install(TARGETS isodeclib DESTINATION ./)
