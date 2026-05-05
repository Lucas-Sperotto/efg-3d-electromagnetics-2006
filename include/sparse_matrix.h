#ifndef SPARSE_MATRIX_H
#define SPARSE_MATRIX_H

#include <stddef.h>

/*
 * Didactic sparse matrix in COO and CSR formats.
 *
 * This module is the first step toward replacing the dense O(n^2) stiffness
 * storage described in docs/TODO.md section 5. The design intentionally stays
 * simple: COO is used for assembly (easy to accumulate repeated entries),
 * then converted to CSR for the matrix-vector product required by GMRES.
 *
 * COO format: list of (row, col, value) triples, possibly with repeated
 *             (row,col) pairs that represent accumulated contributions.
 *
 * CSR format: standard compressed-sparse-row with row_ptr, col_ind, values.
 *             After conversion, duplicate entries are summed.
 */

/* ---------------------------------------------------------------- error codes */

#define SPARSE_OK                  0
#define SPARSE_INVALID_ARGUMENT    1
#define SPARSE_ALLOCATION_FAILED   2
#define SPARSE_INDEX_OUT_OF_RANGE  3

/* ---------------------------------------------------------------- COO */

typedef struct SparseCOO {
    int    *row;      /* row indices    [capacity] */
    int    *col;      /* column indices [capacity] */
    double *val;      /* values         [capacity] */
    int     count;    /* number of entries currently stored */
    int     capacity; /* allocated capacity */
    int     nrows;    /* number of rows (informational) */
    int     ncols;    /* number of columns (informational) */
} SparseCOO;

int sparse_coo_create(SparseCOO *coo, int nrows, int ncols, int initial_capacity);

void sparse_coo_destroy(SparseCOO *coo);

/*
 * Append a single (row, col, value) entry. Duplicate entries are allowed;
 * they will be summed when the matrix is converted to CSR.
 */
int sparse_coo_add(SparseCOO *coo, int row, int col, double value);

/* ---------------------------------------------------------------- CSR */

typedef struct SparseCSR {
    int    *row_ptr;  /* row pointer array [nrows+1] */
    int    *col_ind;  /* column indices    [nnz] */
    double *values;   /* non-zero values   [nnz] */
    int     nrows;
    int     ncols;
    int     nnz;      /* number of non-zeros after compression */
} SparseCSR;

/*
 * Convert a COO matrix to CSR. Entries with the same (row,col) are summed.
 * The COO matrix is not modified. *csr must be uninitialized (its fields will
 * be overwritten).
 */
int sparse_coo_to_csr(const SparseCOO *coo, SparseCSR *csr);

void sparse_csr_destroy(SparseCSR *csr);

/*
 * Compute y = A * x where A is in CSR format.
 * y must point to an array of length csr->nrows.
 * x must point to an array of length csr->ncols.
 */
int sparse_csr_matvec(const SparseCSR *csr, const double *x, double *y);

#endif
