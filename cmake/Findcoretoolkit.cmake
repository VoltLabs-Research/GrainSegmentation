if(TARGET coretoolkit::coretoolkit)
    set(coretoolkit_FOUND TRUE)
    return()
endif()

find_package(Threads REQUIRED)
find_package(TBB REQUIRED)
find_package(spdlog REQUIRED)
find_package(nlohmann_json REQUIRED)

if(NOT DEFINED CORETOOLKIT_SOURCE_DIR OR NOT EXISTS "${CORETOOLKIT_SOURCE_DIR}/CMakeLists.txt")
    message(FATAL_ERROR "coretoolkit not found. Set CORETOOLKIT_SOURCE_DIR to a valid source tree.")
endif()

if(NOT DEFINED CORETOOLKIT_INSTALL_PREFIX OR CORETOOLKIT_INSTALL_PREFIX STREQUAL "")
    set(CORETOOLKIT_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/_deps/install/coretoolkit")
endif()

set(_coretoolkit_library "${CORETOOLKIT_INSTALL_PREFIX}/lib/libcoretoolkit.a")
set(_coretoolkit_build_dir "${CMAKE_BINARY_DIR}/_deps/coretoolkit-build")

if(NOT EXISTS "${_coretoolkit_library}")
    execute_process(
        COMMAND cmake -S "${CORETOOLKIT_SOURCE_DIR}" -B "${_coretoolkit_build_dir}" -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=${CORETOOLKIT_INSTALL_PREFIX} -DCMAKE_MODULE_PATH=${CMAKE_SOURCE_DIR}/cmake -DCMAKE_FIND_PACKAGE_PREFER_CONFIG=FALSE
        RESULT_VARIABLE _coretoolkit_configure_result
    )
    if(NOT _coretoolkit_configure_result EQUAL 0)
        message(FATAL_ERROR "Failed to configure CoreToolkit from ${CORETOOLKIT_SOURCE_DIR}")
    endif()

    execute_process(
        COMMAND cmake --build "${_coretoolkit_build_dir}"
        RESULT_VARIABLE _coretoolkit_build_result
    )
    if(NOT _coretoolkit_build_result EQUAL 0)
        message(FATAL_ERROR "Failed to build CoreToolkit from ${CORETOOLKIT_SOURCE_DIR}")
    endif()

    execute_process(
        COMMAND cmake --install "${_coretoolkit_build_dir}"
        RESULT_VARIABLE _coretoolkit_install_result
    )
    if(NOT _coretoolkit_install_result EQUAL 0)
        message(FATAL_ERROR "Failed to install CoreToolkit into ${CORETOOLKIT_INSTALL_PREFIX}")
    endif()
endif()

add_library(coretoolkit STATIC IMPORTED GLOBAL)
set_target_properties(coretoolkit PROPERTIES
    IMPORTED_LOCATION "${_coretoolkit_library}"
    INTERFACE_INCLUDE_DIRECTORIES "${CORETOOLKIT_INSTALL_PREFIX}/include"
    INTERFACE_LINK_LIBRARIES "Threads::Threads;TBB::tbb;spdlog::spdlog;nlohmann_json::nlohmann_json;${CORETOOLKIT_INSTALL_PREFIX}/lib/libptm.a;${CORETOOLKIT_INSTALL_PREFIX}/lib/libmwm_csp.a;${CORETOOLKIT_INSTALL_PREFIX}/lib/libgeogram.a"
)

add_library(coretoolkit::coretoolkit ALIAS coretoolkit)

set(coretoolkit_FOUND TRUE)
