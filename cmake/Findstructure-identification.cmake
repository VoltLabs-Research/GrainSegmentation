if(TARGET structure-identification::structure-identification)
    set(structure-identification_FOUND TRUE)
    return()
endif()

find_package(coretoolkit REQUIRED)
find_package(spdlog REQUIRED)
find_package(nlohmann_json REQUIRED)

if(NOT DEFINED STRUCTURE_IDENTIFICATION_SOURCE_DIR OR NOT EXISTS "${STRUCTURE_IDENTIFICATION_SOURCE_DIR}/CMakeLists.txt")
    message(FATAL_ERROR "structure-identification not found. Set STRUCTURE_IDENTIFICATION_SOURCE_DIR to a valid source tree.")
endif()

if(NOT DEFINED STRUCTURE_IDENTIFICATION_INSTALL_PREFIX OR STRUCTURE_IDENTIFICATION_INSTALL_PREFIX STREQUAL "")
    set(STRUCTURE_IDENTIFICATION_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/_deps/install/structure-identification")
endif()

set(_structure_identification_library "${STRUCTURE_IDENTIFICATION_INSTALL_PREFIX}/lib/libstructure-identification.a")
set(_structure_identification_build_dir "${CMAKE_BINARY_DIR}/_deps/structure-identification-build")

if(NOT EXISTS "${_structure_identification_library}")
    execute_process(
        COMMAND cmake -S "${STRUCTURE_IDENTIFICATION_SOURCE_DIR}" -B "${_structure_identification_build_dir}" -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=${STRUCTURE_IDENTIFICATION_INSTALL_PREFIX} -DCMAKE_MODULE_PATH=${CMAKE_SOURCE_DIR}/cmake -DCMAKE_FIND_PACKAGE_PREFER_CONFIG=FALSE -DCORETOOLKIT_SOURCE_DIR=${CORETOOLKIT_SOURCE_DIR} -DCORETOOLKIT_INSTALL_PREFIX=${CORETOOLKIT_INSTALL_PREFIX}
        RESULT_VARIABLE _structure_identification_configure_result
    )
    if(NOT _structure_identification_configure_result EQUAL 0)
        message(FATAL_ERROR "Failed to configure StructureIdentification from ${STRUCTURE_IDENTIFICATION_SOURCE_DIR}")
    endif()

    execute_process(
        COMMAND cmake --build "${_structure_identification_build_dir}"
        RESULT_VARIABLE _structure_identification_build_result
    )
    if(NOT _structure_identification_build_result EQUAL 0)
        message(FATAL_ERROR "Failed to build StructureIdentification from ${STRUCTURE_IDENTIFICATION_SOURCE_DIR}")
    endif()

    execute_process(
        COMMAND cmake --install "${_structure_identification_build_dir}"
        RESULT_VARIABLE _structure_identification_install_result
    )
    if(NOT _structure_identification_install_result EQUAL 0)
        message(FATAL_ERROR "Failed to install StructureIdentification into ${STRUCTURE_IDENTIFICATION_INSTALL_PREFIX}")
    endif()
endif()

add_library(structure-identification STATIC IMPORTED GLOBAL)
set_target_properties(structure-identification PROPERTIES
    IMPORTED_LOCATION "${_structure_identification_library}"
    INTERFACE_INCLUDE_DIRECTORIES "${STRUCTURE_IDENTIFICATION_INSTALL_PREFIX}/include"
    INTERFACE_LINK_LIBRARIES "coretoolkit::coretoolkit;spdlog::spdlog;nlohmann_json::nlohmann_json"
)

add_library(structure-identification::structure-identification ALIAS structure-identification)

set(structure-identification_FOUND TRUE)
