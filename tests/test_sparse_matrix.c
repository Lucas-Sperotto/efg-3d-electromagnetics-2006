#include "sparse_matrix.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static void fail(const char *msg)
{
    fprintf(stderr, "FAIL: %s\n", msg);
    exit(1);
}

/* ------------------------------------------------------------------ tests */

/* Insert three entries in a 3×3 matrix, convert, check CSR structure. */
static void test_basic_coo_to_csr(void)
{
    SparseCOO coo;
    SparseCSR csr;

    sparse_coo_create(&coo, 3, 3, 4);
    sparse_coo_add(&coo, 0, 0, 1.0);
    sparse_coo_add(&coo, 1, 1, 2.0);
    sparse_coo_add(&coo, 2, 2, 3.0);

    sparse_coo_to_csr(&coo, &csr);

    if (csr.nnz != 3) {
        fail("basic: expected nnz=3");
    }
    if (csr.row_ptr[0] != 0 || csr.row_ptr[1] != 1 ||
        csr.row_ptr[2] != 2 || csr.row_ptr[3] != 3) {
        fail("basic: wrong row_ptr");
    }

    sparse_coo_destroy(&coo);
    sparse_csr_destroy(&csr);
    printf("PASS test_basic_coo_to_csr\n");
}

/* Duplicate entries must be summed. */
static void test_duplicate_accumulation(void)
{
    SparseCOO coo;
    SparseCSR csr;

    sparse_coo_create(&coo, 2, 2, 4);
    sparse_coo_add(&coo, 0, 0, 1.0);
    sparse_coo_add(&coo, 0, 0, 2.0); /* duplicate */
    sparse_coo_add(&coo, 1, 1, 5.0);

    sparse_coo_to_csr(&coo, &csr);

    if (csr.nnz != 2) {
        fail("dup: expected nnz=2");
    }
    if (fabs(csr.values[0] - 3.0) > 1e-14) {
        fail("dup: (0,0) should be 3.0");
    }
    if (fabs(csr.values[1] - 5.0) > 1e-14) {
        fail("dup: (1,1) should be 5.0");
    }

    sparse_coo_destroy(&coo);
    sparse_csr_destroy(&csr);
    printf("PASS test_duplicate_accumulation\n");
}

/*
 * Matvec: compute y = A * x for a known 3×3 matrix.
 *   A = [[1,2,0],[0,3,4],[5,0,6]]  x = [1,1,1]  expected y = [3,7,11]
 */
static void test_matvec(void)
{
    SparseCOO coo;
    SparseCSR csr;
    double x[3] = {1.0, 1.0, 1.0};
    double y[3];

    sparse_coo_create(&coo, 3, 3, 6);
    sparse_coo_add(&coo, 0, 0, 1.0);
    sparse_coo_add(&coo, 0, 1, 2.0);
    sparse_coo_add(&coo, 1, 1, 3.0);
    sparse_coo_add(&coo, 1, 2, 4.0);
    sparse_coo_add(&coo, 2, 0, 5.0);
    sparse_coo_add(&coo, 2, 2, 6.0);

    sparse_coo_to_csr(&coo, &csr);
    sparse_csr_matvec(&csr, x, y);

    if (fabs(y[0] - 3.0) > 1e-12 ||
        fabs(y[1] - 7.0) > 1e-12 ||
        fabs(y[2] - 11.0) > 1e-12) {
        fprintf(stderr, "FAIL matvec: y = [%.6f, %.6f, %.6f]\n", y[0], y[1], y[2]);
        exit(1);
    }

    sparse_coo_destroy(&coo);
    sparse_csr_destroy(&csr);
    printf("PASS test_matvec\n");
}

/* Growth beyond initial capacity must work. */
static void test_capacity_growth(void)
{
    SparseCOO coo;
    SparseCSR csr;

    sparse_coo_create(&coo, 10, 10, 2); /* tiny initial capacity */

    for (int i = 0; i < 10; ++i) {
        sparse_coo_add(&coo, i, i, (double)(i + 1));
    }

    if (coo.count != 10) {
        fail("growth: expected 10 entries");
    }

    sparse_coo_to_csr(&coo, &csr);

    if (csr.nnz != 10) {
        fail("growth: expected nnz=10");
    }

    sparse_coo_destroy(&coo);
    sparse_csr_destroy(&csr);
    printf("PASS test_capacity_growth\n");
}

/* Empty COO must produce valid CSR with nnz=0. */
static void test_empty_matrix(void)
{
    SparseCOO coo;
    SparseCSR csr;
    double x[2] = {1.0, 2.0};
    double y[2];

    sparse_coo_create(&coo, 2, 2, 4);
    sparse_coo_to_csr(&coo, &csr);

    if (csr.nnz != 0) {
        fail("empty: expected nnz=0");
    }

    sparse_csr_matvec(&csr, x, y);
    if (fabs(y[0]) > 1e-14 || fabs(y[1]) > 1e-14) {
        fail("empty: y should be zero vector");
    }

    sparse_coo_destroy(&coo);
    sparse_csr_destroy(&csr);
    printf("PASS test_empty_matrix\n");
}

/* ------------------------------------------------------------------ main */

int main(void)
{
    test_basic_coo_to_csr();
    test_duplicate_accumulation();
    test_matvec();
    test_capacity_growth();
    test_empty_matrix();
    return 0;
}
