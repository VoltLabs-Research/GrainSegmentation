if(TARGET polyhedral-template-matching::polyhedral-template-matching)
    set(polyhedral-template-matching_FOUND TRUE)
    return()
endif()

find_package(coretoolkit REQUIRED)
find_package(structure-identification REQUIRED)
find_package(Threads REQUIRED)
find_package(TBB REQUIRED)
find_package(Boost REQUIRED)
find_package(spdlog REQUIRED)
find_package(nlohmann_json REQUIRED)

if(NOT DEFINED POLYHEDRAL_TEMPLATE_MATCHING_SOURCE_DIR OR NOT EXISTS "${POLYHEDRAL_TEMPLATE_MATCHING_SOURCE_DIR}/CMakeLists.txt")
    message(FATAL_ERROR "polyhedral-template-matching not found. Set POLYHEDRAL_TEMPLATE_MATCHING_SOURCE_DIR to a valid source tree.")
endif()

if(NOT DEFINED POLYHEDRAL_TEMPLATE_MATCHING_INSTALL_PREFIX OR POLYHEDRAL_TEMPLATE_MATCHING_INSTALL_PREFIX STREQUAL "")
    set(POLYHEDRAL_TEMPLATE_MATCHING_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/_deps/install/polyhedral-template-matching")
endif()

set(_polyhedral_template_matching_library "${POLYHEDRAL_TEMPLATE_MATCHING_INSTALL_PREFIX}/lib/libpolyhedral-template-matching_lib.a")
set(_polyhedral_template_matching_build_dir "${CMAKE_BINARY_DIR}/_deps/polyhedral-template-matching-build")

if(NOT EXISTS "${_polyhedral_template_matching_library}")
    execute_process(
        COMMAND cmake -S "${POLYHEDRAL_TEMPLATE_MATCHING_SOURCE_DIR}" -B "${_polyhedral_template_matching_build_dir}" -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=${POLYHEDRAL_TEMPLATE_MATCHING_INSTALL_PREFIX} -DCMAKE_MODULE_PATH=${CMAKE_SOURCE_DIR}/cmake -DCMAKE_FIND_PACKAGE_PREFER_CONFIG=FALSE -DCORETOOLKIT_SOURCE_DIR=${CORETOOLKIT_SOURCE_DIR} -DCORETOOLKIT_INSTALL_PREFIX=${CORETOOLKIT_INSTALL_PREFIX} -DSTRUCTURE_IDENTIFICATION_SOURCE_DIR=${STRUCTURE_IDENTIFICATION_SOURCE_DIR} -DSTRUCTURE_IDENTIFICATION_INSTALL_PREFIX=${STRUCTURE_IDENTIFICATION_INSTALL_PREFIX}
        RESULT_VARIABLE _polyhedral_template_matching_configure_result
    )
    if(NOT _polyhedral_template_matching_configure_result EQUAL 0)
        message(FATAL_ERROR "Failed to configure PolyhedralTemplateMatching from ${POLYHEDRAL_TEMPLATE_MATCHING_SOURCE_DIR}")
    endif()

    execute_process(
        COMMAND cmake --build "${_polyhedral_template_matching_build_dir}"
        RESULT_VARIABLE _polyhedral_template_matching_build_result
    )
    if(NOT _polyhedral_template_matching_build_result EQUAL 0)
        message(FATAL_ERROR "Failed to build PolyhedralTemplateMatching from ${POLYHEDRAL_TEMPLATE_MATCHING_SOURCE_DIR}")
    endif()

    execute_process(
        COMMAND cmake --install "${_polyhedral_template_matching_build_dir}"
        RESULT_VARIABLE _polyhedral_template_matching_install_result
    )
    if(NOT _polyhedral_template_matching_install_result EQUAL 0)
        message(FATAL_ERROR "Failed to install PolyhedralTemplateMatching into ${POLYHEDRAL_TEMPLATE_MATCHING_INSTALL_PREFIX}")
    endif()
endif()

set(_polyhedral_template_matching_boost_target "")
if(TARGET Boost::headers)
    set(_polyhedral_template_matching_boost_target Boost::headers)
elseif(TARGET Boost::boost)
    set(_polyhedral_template_matching_boost_target Boost::boost)
elseif(TARGET boost::headers)
    set(_polyhedral_template_matching_boost_target boost::headers)
elseif(TARGET boost::boost)
    set(_polyhedral_template_matching_boost_target boost::boost)
else()
    message(FATAL_ERROR "PolyhedralTemplateMatching requires a Boost headers target")
endif()

add_library(polyhedral-template-matching STATIC IMPORTED GLOBAL)
set_target_properties(polyhedral-template-matching PROPERTIES
    IMPORTED_LOCATION "${_polyhedral_template_matching_library}"
    INTERFACE_INCLUDE_DIRECTORIES "${POLYHEDRAL_TEMPLATE_MATCHING_INSTALL_PREFIX}/include;${CORETOOLKIT_SOURCE_DIR}/dependencies/ptm"
    INTERFACE_LINK_LIBRARIES "${_polyhedral_template_matching_boost_target};Threads::Threads;TBB::tbb;structure-identification::structure-identification;coretoolkit::coretoolkit;spdlog::spdlog;nlohmann_json::nlohmann_json"
)

add_library(polyhedral-template-matching::polyhedral-template-matching ALIAS polyhedral-template-matching)

set(polyhedral-template-matching_FOUND TRUE)
