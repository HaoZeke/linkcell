cmake_minimum_required(VERSION 3.22)
project(linkcell-consumer LANGUAGES CXX)

# Installed-prefix contract: find_package(linkcell) yields linkcell::linkcell
# with INTERFACE_INCLUDE_DIRECTORIES and IMPORTED_LOCATION.
find_package(linkcell 0.3 REQUIRED CONFIG)
if(NOT TARGET linkcell::linkcell)
  message(FATAL_ERROR "find_package(linkcell) did not define linkcell::linkcell")
endif()

add_executable(find_lc "${LINKCELL_EXAMPLE}")
target_link_libraries(find_lc PRIVATE linkcell::linkcell)

enable_testing()
add_test(NAME installed-linkcell-consumer COMMAND find_lc)
