#include <assert.h>
#include <math.h>
#include <stdio.h>

#define main sparse_matrix_calculator_main
#include "../sparse_matrix_calculator.c"
#undef main

static int row_element_count(const SparseMatrix *matrix)
{
    int count = 0;
    Row_Node *row = matrix->rowHead;
    while(row != NULL)
    {
        Sm_Node *element = row->rowlist;
        while(element != NULL)
        {
            count++;
            element = element->right;
        }
        row = row->next;
    }
    return count;
}

static int column_element_count(const SparseMatrix *matrix)
{
    int count = 0;
    Col_Node *col = matrix->colHead;
    while(col != NULL)
    {
        Sm_Node *element = col->collist;
        while(element != NULL)
        {
            count++;
            element = element->down;
        }
        col = col->next;
    }
    return count;
}

static matrix_entry value_at(const SparseMatrix *matrix, int row, int col)
{
    Row_Node *row_node = matrix->rowHead;
    while(row_node != NULL && row_node->row < row)
    {
        row_node = row_node->next;
    }
    if(row_node == NULL || row_node->row != row)
    {
        return 0;
    }

    Sm_Node *element = row_node->rowlist;
    while(element != NULL && element->col < col)
    {
        element = element->right;
    }
    return element != NULL && element->col == col ? element->data : 0;
}

static void assert_index_consistent(const SparseMatrix *matrix)
{
    assert(row_element_count(matrix) == column_element_count(matrix));
}

static void test_clear_and_resize(void)
{
    SparseMatrix matrix;
    initializeMatrixWithSize(&matrix, 3, 3);
    assert(insertElement(0, 0, 1, &matrix) == SUCCESS);
    assert(insertElement(0, 2, 4, &matrix) == SUCCESS);
    assert(insertElement(1, 1, 2, &matrix) == SUCCESS);
    assert(insertElement(2, 2, 3, &matrix) == SUCCESS);

    resizeMatrix(&matrix, 2, 1);
    assert(matrix.rowCount == 2 && matrix.colCount == 1);
    assert(value_at(&matrix, 0, 0) == 1);
    assert(value_at(&matrix, 0, 2) == 0);
    assert(value_at(&matrix, 1, 1) == 0);
    assert_index_consistent(&matrix);

    resizeMatrix(&matrix, 1, 1);
    assert(matrix.rowCount == 1 && matrix.colCount == 1);
    assert_index_consistent(&matrix);

    clearMatrix(&matrix);
    assert(matrix.rowHead == NULL && matrix.colHead == NULL);
}

static void test_fractional_arithmetic(void)
{
    SparseMatrix left, right, result;
    initializeMatrixWithSize(&left, 1, 1);
    initializeMatrixWithSize(&right, 1, 1);
    assert(insertElement(0, 0, 1.3f, &left) == SUCCESS);
    assert(insertElement(0, 0, 0.3f, &right) == SUCCESS);
    assert(addMatrix(&left, &right, &result) == SUCCESS);
    assert(fabsf(value_at(&result, 0, 0) - 1.6f) < 0.001f);
    clearMatrix(&result);

    clearMatrix(&right);
    assert(insertElement(0, 0, 2.5f, &right) == SUCCESS);
    assert(multiplyMatrix(&left, &right, &result) == SUCCESS);
    assert(fabsf(value_at(&result, 0, 0) - 3.25f) < 0.001f);

    clearMatrix(&left);
    clearMatrix(&right);
    clearMatrix(&result);
}

static void test_update_and_inverse(void)
{
    SparseMatrix matrix, inverse;
    initializeMatrixWithSize(&matrix, 1, 1);
    assert(insertElement(0, 0, 5.0f, &matrix) == SUCCESS);
    assert(insertElement(0, 0, 9.0f, &matrix) == SUCCESS);
    assert(value_at(&matrix, 0, 0) == 9.0f);

    assert(inverseOfMatrix(&matrix, &inverse) == SUCCESS);
    assert(fabsf(value_at(&inverse, 0, 0) - (1.0f / 9.0f)) < 0.001f);
    assert_index_consistent(&inverse);

    clearMatrix(&matrix);
    clearMatrix(&inverse);
}

static void test_invalid_coordinates(void)
{
    SparseMatrix matrix;
    initializeMatrixWithSize(&matrix, 2, 2);
    assert(insertElement(-1, 0, 1, &matrix) == FAILURE);
    assert(insertElement(0, -1, 1, &matrix) == FAILURE);
    assert(insertElement(2, 0, 1, &matrix) == FAILURE);
    assert(insertElement(0, 2, 1, &matrix) == FAILURE);
    resizeMatrix(&matrix, -1, 2);
    assert(matrix.rowCount == 2 && matrix.colCount == 2);
    clearMatrix(&matrix);
}

static void test_zero_values_and_dimensions(void)
{
    SparseMatrix matrix, empty, other, result;
    initializeMatrixWithSize(&matrix, 2, 2);
    assert(insertElement(0, 0, 4.0f, &matrix) == SUCCESS);
    assert(insertElement(0, 0, 0.0f, &matrix) == SUCCESS);
    assert(value_at(&matrix, 0, 0) == 4.0f);
    scalarMultiplyMatrix(&matrix, 0.0f);
    assert(row_element_count(&matrix) == 0);
    assert(column_element_count(&matrix) == 0);

    initializeMatrixWithSize(&empty, 2, 2);
    initializeMatrixWithSize(&other, 3, 3);
    assert(addMatrix(&empty, &other, &result) == FAILURE);
    assert(multiplyMatrix(&empty, &other, &result) == FAILURE);
    assert(result.rowHead == NULL && result.colHead == NULL);

    initializeMatrixWithSize(&other, 2, 3);
    assert(multiplyMatrix(&empty, &other, &result) == SUCCESS);
    assert(result.rowCount == 2 && result.colCount == 3);
    assert(result.rowHead == NULL && result.colHead == NULL);

    clearMatrix(&matrix);
    clearMatrix(&empty);
    clearMatrix(&other);
    clearMatrix(&result);
}

static void test_scaled_singularity_tolerance(void)
{
    SparseMatrix small, mixed, inverse;
    float determinant_value;

    initializeMatrixWithSize(&small, 3, 3);
    assert(insertElement(0, 0, 0.001f, &small) == SUCCESS);
    assert(insertElement(1, 1, 0.001f, &small) == SUCCESS);
    assert(insertElement(2, 2, 0.001f, &small) == SUCCESS);
    assert(determinantOfMatrix(&small, &determinant_value) == SUCCESS);
    assert(fabsf(determinant_value - 1e-9f) < 1e-11f);
    assert(inverseOfMatrix(&small, &inverse) == SUCCESS);
    assert(fabsf(value_at(&inverse, 0, 0) - 1000.0f) < 0.1f);
    clearMatrix(&inverse);

    initializeMatrixWithSize(&mixed, 3, 3);
    assert(insertElement(0, 0, 1e6f, &mixed) == SUCCESS);
    assert(insertElement(1, 1, 1e-6f, &mixed) == SUCCESS);
    assert(insertElement(2, 2, 1e-6f, &mixed) == SUCCESS);
    assert(inverseOfMatrix(&mixed, &inverse) == SUCCESS);
    assert(fabsf(value_at(&inverse, 0, 0) - 1e-6f) < 1e-8f);
    assert(fabsf(value_at(&inverse, 1, 1) - 1e6f) < 1.0f);

    clearMatrix(&small);
    clearMatrix(&mixed);
    clearMatrix(&inverse);
}

int main(void)
{
    test_clear_and_resize();
    test_fractional_arithmetic();
    test_update_and_inverse();
    test_invalid_coordinates();
    test_zero_values_and_dimensions();
    test_scaled_singularity_tolerance();
    puts("All sparse matrix regression tests passed.");
    return 0;
}
