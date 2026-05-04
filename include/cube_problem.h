#ifndef CUBE_PROBLEM_H
#define CUBE_PROBLEM_H

#include "dirichlet.h"
#include "nodes.h"

/*
 * Geometry and boundary-condition data for the electrostatic cube problem.
 *
 * This module only prepares regular nodes and Dirichlet points on the cube
 * boundary. It does not assemble K, G, the augmented Lagrange system, solve
 * linear systems, or call GMRES.
 */

#define CUBE_PROBLEM_OK 0
#define CUBE_PROBLEM_INVALID_ARGUMENT -1
#define CUBE_PROBLEM_OUTPUT_CAPACITY_TOO_SMALL -2

int cube_regular_node_count(int nx, int ny, int nz);

int cube_generate_regular_nodes(double L,
                                int nx,
                                int ny,
                                int nz,
                                Node3D *nodes,
                                int max_nodes,
                                int *node_count);

int cube_dirichlet_point_count(int nx, int ny, int nz);

int cube_generate_dirichlet_points(double L,
                                   double V0,
                                   int nx,
                                   int ny,
                                   int nz,
                                   DirichletPoint *points,
                                   int max_points,
                                   int *point_count);

#endif
