#ifndef MAT4_H
#define MAT4_H

/*
 * Small 4x4 matrix utilities.
 *
 * These helpers are intended for the local MLS moment matrix A(x) with a
 * linear 3-D basis p^T = [1, x, y, z]. They do not implement MLS yet; they only
 * provide the tiny dense linear algebra block that MLS will need later.
 */

#define MAT4_OK 0
#define MAT4_INVALID_ARGUMENT 1
#define MAT4_SINGULAR 2
#define MAT4_NEAR_SINGULAR 3

void mat4_zero(double A[4][4]);

void mat4_identity(double A[4][4]);

void mat4_copy(const double src[4][4], double dst[4][4]);

int mat4_solve(const double A[4][4], const double b[4], double x[4]);

#endif
