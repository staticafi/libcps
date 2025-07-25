set(THIS_TARGET_NAME cps)

find_package(Eigen3 CONFIG REQUIRED)

file(GLOB_RECURSE "${THIS_TARGET_NAME}_HPP" "./include/${THIS_TARGET_NAME}/*.hpp")
file(GLOB_RECURSE "${THIS_TARGET_NAME}_CPP" "./src/*.cpp")

add_library(${THIS_TARGET_NAME}
    "${${THIS_TARGET_NAME}_HPP}"
    "${${THIS_TARGET_NAME}_CPP}"
    )

target_link_libraries(${THIS_TARGET_NAME} PRIVATE Eigen3::Eigen utility)
target_include_directories(${THIS_TARGET_NAME} PUBLIC include)

if(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    target_compile_options(${THIS_TARGET_NAME} PUBLIC "/wd4996")
else()
    target_compile_options(${THIS_TARGET_NAME} PUBLIC "-Wno-deprecated")
endif()

set_target_properties(${THIS_TARGET_NAME} PROPERTIES
    DEBUG_OUTPUT_NAME "${THIS_TARGET_NAME}_${CMAKE_SYSTEM_NAME}_Debug"
    RELEASE_OUTPUT_NAME "${THIS_TARGET_NAME}_${CMAKE_SYSTEM_NAME}_Release"
    RELWITHDEBINFO_OUTPUT_NAME "${THIS_TARGET_NAME}_${CMAKE_SYSTEM_NAME}_RelWithDebInfo"
    )

if(DEFINED LIBCPS_INSTALL_DIR)
    install(TARGETS ${THIS_TARGET_NAME} DESTINATION "${LIBCPS_INSTALL_DIR}")
endif(DEFINED LIBCPS_INSTALL_DIR)

