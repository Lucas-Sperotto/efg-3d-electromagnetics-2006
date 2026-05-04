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
  src/dense_solver.c
  src/dirichlet.c
  src/error_norm.c
  src/gauss.c
  src/global_stiffness.c
  src/gmres.c
  src/lagrange_system.c
  src/mat4.c
  src/mls.c
  src/nodes.c
  src/stiffness.c
  src/support.c
  src/vec3.c
  src/weight.c
)

target_include_directories(efg_core
  PUBLIC
    ${PROJECT_SOURCE_DIR}/include
)

if(UNIX)
  target_link_libraries(efg_core PUBLIC m)
endif()

add_executable(smoke_test
  apps/smoke_test.c
)

target_link_libraries(smoke_test
  PRIVATE
    efg_core
)

add_executable(reproduce_analytical_slice
  apps/reproduce_analytical_slice.c
)

target_link_libraries(reproduce_analytical_slice
  PRIVATE
    efg_core
)

add_executable(test_gauss_cube
  apps/test_gauss_cube.c
)

target_link_libraries(test_gauss_cube
  PRIVATE
    efg_core
)

add_executable(test_dirichlet_constraints
  tests/test_dirichlet_constraints.c
)

target_link_libraries(test_dirichlet_constraints
  PRIVATE
    efg_core
)

add_executable(test_dense_solver
  tests/test_dense_solver.c
)

target_link_libraries(test_dense_solver
  PRIVATE
    efg_core
)

add_executable(test_gauss
  tests/test_gauss.c
)

target_link_libraries(test_gauss
  PRIVATE
    efg_core
)

add_executable(test_global_stiffness
  tests/test_global_stiffness.c
)

target_link_libraries(test_global_stiffness
  PRIVATE
    efg_core
)

add_executable(test_lagrange_system
  tests/test_lagrange_system.c
)

target_link_libraries(test_lagrange_system
  PRIVATE
    efg_core
)

add_executable(test_mat4
  tests/test_mat4.c
)

target_link_libraries(test_mat4
  PRIVATE
    efg_core
)

add_executable(test_mls_partition
  tests/test_mls_partition.c
)

target_link_libraries(test_mls_partition
  PRIVATE
    efg_core
)

add_executable(test_mls_linear_reproduction
  tests/test_mls_linear_reproduction.c
)

target_link_libraries(test_mls_linear_reproduction
  PRIVATE
    efg_core
)

add_executable(test_mls_gradient_partition
  tests/test_mls_gradient_partition.c
)

target_link_libraries(test_mls_gradient_partition
  PRIVATE
    efg_core
)

add_executable(test_mls_gradient_linear_reproduction
  tests/test_mls_gradient_linear_reproduction.c
)

target_link_libraries(test_mls_gradient_linear_reproduction
  PRIVATE
    efg_core
)

add_executable(test_stiffness_cell
  tests/test_stiffness_cell.c
)

target_link_libraries(test_stiffness_cell
  PRIVATE
    efg_core
)

add_executable(test_error_norm
  tests/test_error_norm.c
)

target_link_libraries(test_error_norm
  PRIVATE
    efg_core
)

add_executable(test_support
  tests/test_support.c
)

target_link_libraries(test_support
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

add_executable(test_weight_derivatives
  tests/test_weight_derivatives.c
)

target_link_libraries(test_weight_derivatives
  PRIVATE
    efg_core
)

add_executable(test_weight_tensor3d
  tests/test_weight_tensor3d.c
)

target_link_libraries(test_weight_tensor3d
  PRIVATE
    efg_core
)

if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(efg_core PRIVATE -Wall -Wextra -Wpedantic)
  target_compile_options(smoke_test PRIVATE -Wall -Wextra -Wpedantic)
  target_compile_options(reproduce_analytical_slice PRIVATE -Wall -Wextra -Wpedantic)
  target_compile_options(test_gauss_cube PRIVATE -Wall -Wextra -Wpedantic)
  target_compile_options(test_dirichlet_constraints PRIVATE -Wall -Wextra -Wpedantic)
  target_compile_options(test_dense_solver PRIVATE -Wall -Wextra -Wpedantic)
  target_compile_options(test_gauss PRIVATE -Wall -Wextra -Wpedantic)
  target_compile_options(test_global_stiffness PRIVATE -Wall -Wextra -Wpedantic)
  target_compile_options(test_lagrange_system PRIVATE -Wall -Wextra -Wpedantic)
  target_compile_options(test_mat4 PRIVATE -Wall -Wextra -Wpedantic)
  target_compile_options(test_mls_partition PRIVATE -Wall -Wextra -Wpedantic)
  target_compile_options(test_mls_linear_reproduction PRIVATE -Wall -Wextra -Wpedantic)
  target_compile_options(test_mls_gradient_partition PRIVATE -Wall -Wextra -Wpedantic)
  target_compile_options(test_mls_gradient_linear_reproduction PRIVATE -Wall -Wextra -Wpedantic)
  target_compile_options(test_stiffness_cell PRIVATE -Wall -Wextra -Wpedantic)
  target_compile_options(test_error_norm PRIVATE -Wall -Wextra -Wpedantic)
  target_compile_options(test_support PRIVATE -Wall -Wextra -Wpedantic)
  target_compile_options(test_weight PRIVATE -Wall -Wextra -Wpedantic)
  target_compile_options(test_weight_derivatives PRIVATE -Wall -Wextra -Wpedantic)
  target_compile_options(test_weight_tensor3d PRIVATE -Wall -Wextra -Wpedantic)
endif()

enable_testing()
add_test(NAME smoke_test COMMAND smoke_test)
add_test(NAME test_dirichlet_constraints COMMAND test_dirichlet_constraints)
add_test(NAME test_dense_solver COMMAND test_dense_solver)
add_test(NAME test_error_norm COMMAND test_error_norm)
add_test(NAME test_gauss COMMAND test_gauss)
add_test(NAME test_global_stiffness COMMAND test_global_stiffness)
add_test(NAME test_lagrange_system COMMAND test_lagrange_system)
add_test(NAME test_mat4 COMMAND test_mat4)
add_test(NAME test_mls_partition COMMAND test_mls_partition)
add_test(NAME test_mls_linear_reproduction COMMAND test_mls_linear_reproduction)
add_test(NAME test_mls_gradient_partition COMMAND test_mls_gradient_partition)
add_test(NAME test_mls_gradient_linear_reproduction COMMAND test_mls_gradient_linear_reproduction)
add_test(NAME test_stiffness_cell COMMAND test_stiffness_cell)
add_test(NAME test_support COMMAND test_support)
add_test(NAME test_weight COMMAND test_weight)
add_test(NAME test_weight_derivatives COMMAND test_weight_derivatives)
add_test(NAME test_weight_tensor3d COMMAND test_weight_tensor3d)
