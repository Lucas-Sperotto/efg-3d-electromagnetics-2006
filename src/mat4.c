#include "mat4.h"

#include <math.h>
#include <stddef.h>

#define MAT4_SIZE 4
#define MAT4_PIVOT_TOLERANCE 1.0e-12

void mat4_zero(double A[4][4])
{
    if (A == NULL) {
        return;
    }

    for (int row = 0; row < MAT4_SIZE; ++row) {
        for (int col = 0; col < MAT4_SIZE; ++col) {
            A[row][col] = 0.0;
        }
    }
}

void mat4_identity(double A[4][4])
{
    if (A == NULL) {
        return;
    }

    mat4_zero(A);
    for (int i = 0; i < MAT4_SIZE; ++i) {
        A[i][i] = 1.0;
    }
}

void mat4_copy(const double src[4][4], double dst[4][4])
{
    if (src == NULL || dst == NULL) {
        return;
    }

    for (int row = 0; row < MAT4_SIZE; ++row) {
        for (int col = 0; col < MAT4_SIZE; ++col) {
            dst[row][col] = src[row][col];
        }
    }
}

static void swap_rows(double augmented[4][5], int row_a, int row_b)
{
    if (row_a == row_b) {
        return;
    }

    for (int col = 0; col < 5; ++col) {
        const double temporary = augmented[row_a][col];
        augmented[row_a][col] = augmented[row_b][col];
        augmented[row_b][col] = temporary;
    }
}

/*
 * Solve A x = b by Gaussian elimination with partial pivoting.
 *
 * This routine will be used later to solve the local MLS moment matrix A(x)
 * for a linear 3-D basis. It intentionally works on a local augmented copy, so
 * the caller's A and b are not modified.
 */
int mat4_solve(const double A[4][4], const double b[4], double x[4])
{
    double augmented[4][5];

    if (A == NULL || b == NULL || x == NULL) {
        return MAT4_INVALID_ARGUMENT;
    }

    for (int row = 0; row < MAT4_SIZE; ++row) {
        for (int col = 0; col < MAT4_SIZE; ++col) {
            augmented[row][col] = A[row][col];
        }
        augmented[row][MAT4_SIZE] = b[row];
    }

    for (int pivot_col = 0; pivot_col < MAT4_SIZE; ++pivot_col) {
        int pivot_row = pivot_col;
        double pivot_abs = fabs(augmented[pivot_col][pivot_col]);

        for (int row = pivot_col + 1; row < MAT4_SIZE; ++row) {
            const double candidate_abs = fabs(augmented[row][pivot_col]);

            if (candidate_abs > pivot_abs) {
                pivot_abs = candidate_abs;
                pivot_row = row;
            }
        }

        if (pivot_abs == 0.0) {
            return MAT4_SINGULAR;
        }

        if (pivot_abs < MAT4_PIVOT_TOLERANCE) {
            return MAT4_NEAR_SINGULAR;
        }

        swap_rows(augmented, pivot_col, pivot_row);

        for (int row = pivot_col + 1; row < MAT4_SIZE; ++row) {
            const double factor = augmented[row][pivot_col] / augmented[pivot_col][pivot_col];

            augmented[row][pivot_col] = 0.0;
            for (int col = pivot_col + 1; col <= MAT4_SIZE; ++col) {
                augmented[row][col] -= factor * augmented[pivot_col][col];
            }
        }
    }

    for (int row = MAT4_SIZE - 1; row >= 0; --row) {
        double sum = augmented[row][MAT4_SIZE];
        const double diagonal = augmented[row][row];

        if (diagonal == 0.0) {
            return MAT4_SINGULAR;
        }

        if (fabs(diagonal) < MAT4_PIVOT_TOLERANCE) {
            return MAT4_NEAR_SINGULAR;
        }

        for (int col = row + 1; col < MAT4_SIZE; ++col) {
            sum -= augmented[row][col] * x[col];
        }

        x[row] = sum / diagonal;
    }

    return MAT4_OK;
}

/*
 * Estimate condition number via max/min absolute pivot from Gaussian
 * elimination with partial pivoting (no rhs needed). The ratio gives a
 * cheap surrogate for the true condition number without requiring A^{-1}.
 */
int mat4_cond_estimate(const double A[4][4], double *cond_estimate)
{
    double work[4][4];
    double pivot_max = 0.0;
    double pivot_min = -1.0; /* sentinel */

    if (A == NULL || cond_estimate == NULL) {
        return MAT4_INVALID_ARGUMENT;
    }

    for (int row = 0; row < MAT4_SIZE; ++row) {
        for (int col = 0; col < MAT4_SIZE; ++col) {
            work[row][col] = A[row][col];
        }
    }

    for (int pivot_col = 0; pivot_col < MAT4_SIZE; ++pivot_col) {
        int pivot_row = pivot_col;
        double pivot_abs = fabs(work[pivot_col][pivot_col]);

        for (int row = pivot_col + 1; row < MAT4_SIZE; ++row) {
            const double cand = fabs(work[row][pivot_col]);
            if (cand > pivot_abs) {
                pivot_abs = cand;
                pivot_row = row;
            }
        }

        if (pivot_abs == 0.0) {
            return MAT4_SINGULAR;
        }

        /* swap rows */
        if (pivot_row != pivot_col) {
            for (int col = pivot_col; col < MAT4_SIZE; ++col) {
                double tmp = work[pivot_col][col];
                work[pivot_col][col] = work[pivot_row][col];
                work[pivot_row][col] = tmp;
            }
        }

        if (pivot_abs > pivot_max) {
            pivot_max = pivot_abs;
        }
        if (pivot_min < 0.0 || pivot_abs < pivot_min) {
            pivot_min = pivot_abs;
        }

        for (int row = pivot_col + 1; row < MAT4_SIZE; ++row) {
            const double factor = work[row][pivot_col] / work[pivot_col][pivot_col];
            work[row][pivot_col] = 0.0;
            for (int col = pivot_col + 1; col < MAT4_SIZE; ++col) {
                work[row][col] -= factor * work[pivot_col][col];
            }
        }
    }

    *cond_estimate = (pivot_min > 0.0) ? (pivot_max / pivot_min) : 0.0;
    return MAT4_OK;
}
