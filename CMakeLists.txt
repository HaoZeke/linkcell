cmake_minimum_required(VERSION 3.22)
project(linkcell VERSION 0.1.2 LANGUAGES C CXX)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
find_package(Threads REQUIRED)
find_program(CARGO_EXECUTABLE cargo REQUIRED)

option(LINKCELL_BUILD_EXAMPLES "Build C and C++ example programs" ON)
option(LINKCELL_INSTALL "Install headers, archive, pkg-config, and CMake package" ON)

if(WIN32)
  set(LINKCELL_ARCHIVE_NAME linkcell.lib)
else()
  set(LINKCELL_ARCHIVE_NAME liblinkcell.a)
endif()

set(LINKCELL_ARCHIVE "${CMAKE_CURRENT_BINARY_DIR}/${LINKCELL_ARCHIVE_NAME}")
set(LINKCELL_CARGO_DIR "${CMAKE_CURRENT_BINARY_DIR}/cargo-target")

set(LINKCELL_RUST_SOURCES
  "${CMAKE_CURRENT_SOURCE_DIR}/Cargo.toml"
  "${CMAKE_CURRENT_SOURCE_DIR}/Cargo.lock"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/lib.rs"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/cell.rs"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/error.rs"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/knearest.rs"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/capi.rs"
)

add_custom_command(
  OUTPUT "${LINKCELL_ARCHIVE}"
  COMMAND "${CARGO_EXECUTABLE}" build --release --features capi
          --manifest-path "${CMAKE_CURRENT_SOURCE_DIR}/Cargo.toml"
          --target-dir "${LINKCELL_CARGO_DIR}"
  COMMAND "${CMAKE_COMMAND}" -E copy
          "${LINKCELL_CARGO_DIR}/release/${LINKCELL_ARCHIVE_NAME}"
          "${LINKCELL_ARCHIVE}"
  DEPENDS ${LINKCELL_RUST_SOURCES}
  COMMENT "Building liblinkcell via cargo"
  VERBATIM
)
add_custom_target(linkcell_cargo ALL DEPENDS "${LINKCELL_ARCHIVE}")

add_library(linkcell INTERFACE)
add_library(linkcell::linkcell ALIAS linkcell)
add_dependencies(linkcell linkcell_cargo)
target_include_directories(linkcell INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
set(_linkcell_sys Threads::Threads)
if(UNIX AND NOT APPLE)
  list(APPEND _linkcell_sys ${CMAKE_DL_LIBS} m)
elseif(APPLE)
  list(APPEND _linkcell_sys m)
elseif(WIN32)
  list(APPEND _linkcell_sys ws2_32 userenv bcrypt ntdll)
endif()
target_link_libraries(linkcell INTERFACE
  $<BUILD_INTERFACE:${LINKCELL_ARCHIVE}>
  ${_linkcell_sys}
)

if(LINKCELL_BUILD_EXAMPLES)
  add_executable(two_points examples/two_points.c)
  target_link_libraries(two_points PRIVATE linkcell::linkcell)
  add_dependencies(two_points linkcell_cargo)
  add_executable(two_points_cpp examples/two_points.cpp)
  target_link_libraries(two_points_cpp PRIVATE linkcell::linkcell)
  add_dependencies(two_points_cpp linkcell_cargo)
  enable_testing()
  add_test(NAME two-points-c COMMAND two_points)
  add_test(NAME two-points-cpp COMMAND two_points_cpp)
endif()

if(LINKCELL_INSTALL)
  if(UNIX AND NOT APPLE)
    set(LINKCELL_PKG_LIBS "-llinkcell -ldl -lpthread -lm")
  elseif(APPLE)
    set(LINKCELL_PKG_LIBS "-llinkcell -lpthread -lm")
  else()
    set(LINKCELL_PKG_LIBS "-llinkcell")
  endif()

  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/linkcell.pc.in"
    "${CMAKE_CURRENT_BINARY_DIR}/linkcell.pc"
    @ONLY
  )
  configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/linkcellConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/linkcellConfig.cmake"
    INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/linkcell"
  )
  write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/linkcellConfigVersion.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMinorVersion
  )

  install(FILES "${LINKCELL_ARCHIVE}" DESTINATION "${CMAKE_INSTALL_LIBDIR}")
  install(FILES
    "${CMAKE_CURRENT_SOURCE_DIR}/include/linkcell.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/linkcell.hpp"
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
  )
  install(FILES "${CMAKE_CURRENT_BINARY_DIR}/linkcell.pc"
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/pkgconfig"
  )
  install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/linkcellConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/linkcellConfigVersion.cmake"
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/linkcell"
  )
endif()
