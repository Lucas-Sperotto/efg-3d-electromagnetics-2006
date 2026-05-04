cmake_minimum_required(VERSION 3.16)

project(efg_3d_electromagnetics_2006
  VERSION 0.1.0
  LANGUAGES C
)

# This repository is intentionally C11-only in the first didactic version.
# No external libraries are used at this stage.
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

add_library(efg_core
  src/analytical.c
  src/assembly.c
  src/error_norm.c
  src/gauss.c
  src/gmres.c
  src/mls.c
  src/nodes.c
  src/vec3.c
  src/weight.c
)

target_include_directories(efg_core
  PUBLIC
    ${PROJECT_SOURCE_DIR}/include
)

add_executable(smoke_test
  apps/smoke_test.c
)

target_link_libraries(smoke_test
  PRIVATE
    efg_core
)

add_executable(test_weight
  tests/test_weight.c
)

target_link_libraries(test_weight
  PRIVATE
    efg_core
)

if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(efg_core PRIVATE -Wall -Wextra -Wpedantic)
  target_compile_options(smoke_test PRIVATE -Wall -Wextra -Wpedantic)
  target_compile_options(test_weight PRIVATE -Wall -Wextra -Wpedantic)
endif()

enable_testing()
add_test(NAME smoke_test COMMAND smoke_test)
add_test(NAME test_weight COMMAND test_weight)
