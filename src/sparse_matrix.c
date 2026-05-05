#include "sparse_matrix.h"

#include <stddef.h>
#include <stdlib.h>
#include <string.h>

/* ---------------------------------------------------------------- COO */

int sparse_coo_create(SparseCOO *coo, int nrows, int ncols, int initial_capacity)
{
    if (coo == NULL || nrows <= 0 || ncols <= 0 || initial_capacity <= 0) {
        return SPARSE_INVALID_ARGUMENT;
    }

    coo->row = malloc((size_t)initial_capacity * sizeof(coo->row[0]));
    coo->col = malloc((size_t)initial_capacity * sizeof(coo->col[0]));
    coo->val = malloc((size_t)initial_capacity * sizeof(coo->val[0]));

    if (coo->row == NULL || coo->col == NULL || coo->val == NULL) {
        free(coo->row);
        free(coo->col);
        free(coo->val);
        coo->row = coo->col = NULL;
        coo->val = NULL;
        return SPARSE_ALLOCATION_FAILED;
    }

    coo->count    = 0;
    coo->capacity = initial_capacity;
    coo->nrows    = nrows;
    coo->ncols    = ncols;
    return SPARSE_OK;
}

void sparse_coo_destroy(SparseCOO *coo)
{
    if (coo == NULL) {
        return;
    }
    free(coo->row);
    free(coo->col);
    free(coo->val);
    coo->row = coo->col = NULL;
    coo->val = NULL;
    coo->count = coo->capacity = 0;
}

int sparse_coo_add(SparseCOO *coo, int row, int col, double value)
{
    if (coo == NULL) {
        return SPARSE_INVALID_ARGUMENT;
    }
    if (row < 0 || row >= coo->nrows || col < 0 || col >= coo->ncols) {
        return SPARSE_INDEX_OUT_OF_RANGE;
    }

    /* Grow by doubling when capacity is reached. */
    if (coo->count == coo->capacity) {
        int new_cap = coo->capacity * 2;
        int *nr = realloc(coo->row, (size_t)new_cap * sizeof(coo->row[0]));
        int *nc = realloc(coo->col, (size_t)new_cap * sizeof(coo->col[0]));
        double *nv = realloc(coo->val, (size_t)new_cap * sizeof(coo->val[0]));

        if (nr == NULL || nc == NULL || nv == NULL) {
            free(nr);
            free(nc);
            free(nv);
            return SPARSE_ALLOCATION_FAILED;
        }

        coo->row      = nr;
        coo->col      = nc;
        coo->val      = nv;
        coo->capacity = new_cap;
    }

    coo->row[coo->count] = row;
    coo->col[coo->count] = col;
    coo->val[coo->count] = value;
    ++coo->count;
    return SPARSE_OK;
}

/* ---------------------------------------------------------------- COO → CSR */

/*
 * Sort key: row-major order. Used by qsort to order the COO entries before
 * compression into CSR.
 */
typedef struct SortEntry {
    int    row;
    int    col;
    double val;
} SortEntry;

static int cmp_sort_entry(const void *a, const void *b)
{
    const SortEntry *ea = (const SortEntry *)a;
    const SortEntry *eb = (const SortEntry *)b;

    if (ea->row != eb->row) {
        return ea->row - eb->row;
    }
    return ea->col - eb->col;
}

int sparse_coo_to_csr(const SparseCOO *coo, SparseCSR *csr)
{
    SortEntry *entries = NULL;
    int nnz_unique = 0;
    int prev_row = -1;
    int prev_col = -1;

    if (coo == NULL || csr == NULL) {
        return SPARSE_INVALID_ARGUMENT;
    }

    /* Initialise output so partial failures are safe to destroy. */
    csr->row_ptr = NULL;
    csr->col_ind = NULL;
    csr->values  = NULL;
    csr->nrows   = coo->nrows;
    csr->ncols   = coo->ncols;
    csr->nnz     = 0;

    if (coo->count == 0) {
        csr->row_ptr = calloc((size_t)(coo->nrows + 1), sizeof(csr->row_ptr[0]));
        return (csr->row_ptr != NULL) ? SPARSE_OK : SPARSE_ALLOCATION_FAILED;
    }

    /* Copy entries to a sortable array. */
    entries = malloc((size_t)coo->count * sizeof(entries[0]));
    if (entries == NULL) {
        return SPARSE_ALLOCATION_FAILED;
    }

    for (int i = 0; i < coo->count; ++i) {
        entries[i].row = coo->row[i];
        entries[i].col = coo->col[i];
        entries[i].val = coo->val[i];
    }

    qsort(entries, (size_t)coo->count, sizeof(entries[0]), cmp_sort_entry);

    /* Count unique (row,col) pairs to determine CSR nnz. */
    for (int i = 0; i < coo->count; ++i) {
        if (entries[i].row != prev_row || entries[i].col != prev_col) {
            ++nnz_unique;
            prev_row = entries[i].row;
            prev_col = entries[i].col;
        }
    }

    /* Allocate CSR arrays. */
    csr->row_ptr = calloc((size_t)(coo->nrows + 1), sizeof(csr->row_ptr[0]));
    csr->col_ind = malloc((size_t)nnz_unique * sizeof(csr->col_ind[0]));
    csr->values  = calloc((size_t)nnz_unique, sizeof(csr->values[0]));

    if (csr->row_ptr == NULL || csr->col_ind == NULL || csr->values == NULL) {
        sparse_csr_destroy(csr);
        free(entries);
        return SPARSE_ALLOCATION_FAILED;
    }

    /* Fill CSR: accumulate duplicates. */
    int pos = -1; /* current position in col_ind / values */
    prev_row = -1;
    prev_col = -1;

    for (int i = 0; i < coo->count; ++i) {
        int r = entries[i].row;
        int c = entries[i].col;

        if (r != prev_row || c != prev_col) {
            ++pos;
            csr->col_ind[pos] = c;
            csr->values[pos]  = entries[i].val;
            prev_col = c;
            /* Fill row_ptr for every new row encountered. */
            for (int rr = prev_row + 1; rr <= r; ++rr) {
                csr->row_ptr[rr] = pos;
            }
            prev_row = r;
        } else {
            csr->values[pos] += entries[i].val;
        }
    }

    /* Close remaining rows. */
    for (int rr = prev_row + 1; rr <= coo->nrows; ++rr) {
        csr->row_ptr[rr] = nnz_unique;
    }

    csr->nnz = nnz_unique;
    free(entries);
    return SPARSE_OK;
}

void sparse_csr_destroy(SparseCSR *csr)
{
    if (csr == NULL) {
        return;
    }
    free(csr->row_ptr);
    free(csr->col_ind);
    free(csr->values);
    csr->row_ptr = NULL;
    csr->col_ind = NULL;
    csr->values  = NULL;
    csr->nnz     = 0;
}

/* ---------------------------------------------------------------- matvec */

int sparse_csr_matvec(const SparseCSR *csr, const double *x, double *y)
{
    if (csr == NULL || x == NULL || y == NULL) {
        return SPARSE_INVALID_ARGUMENT;
    }

    for (int i = 0; i < csr->nrows; ++i) {
        double sum = 0.0;
        for (int j = csr->row_ptr[i]; j < csr->row_ptr[i + 1]; ++j) {
            sum += csr->values[j] * x[csr->col_ind[j]];
        }
        y[i] = sum;
    }

    return SPARSE_OK;
}
