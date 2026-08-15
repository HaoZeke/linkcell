cmake_minimum_required(VERSION 3.22)
project(linkcell-consumer LANGUAGES CXX)

find_package(linkcell 0.1 REQUIRED)
add_executable(find_lc "${LINKCELL_EXAMPLE}")
target_link_libraries(find_lc PRIVATE linkcell::linkcell)
