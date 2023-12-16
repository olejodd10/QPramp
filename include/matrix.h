#ifndef MATRIX_H
#define MATRIX_H

#include <stddef.h>

/**
 * @brief Compute matrix-vector product.
 *
 * @param[in] m Height of matrix
 * @param[in] n Width of matrix
 * @param[in] mat Flattened m by n array representing matrix
 * @param[in] vec Array with length n representing vector
 * @param[out] res Array with length m for storing result mat*vec
 */
void matrix_vector_product(size_t m, size_t n, const double mat[m][n], const double vec[n], double res[m]);

/**
 * @brief Compute matrix product.
 *
 * @param[in] m Height of m1
 * @param[in] n Width of m1
 * @param[in] p Width of m2
 * @param[in] m1 Flattened m by n array representing matrix
 * @param[in] m2t Flattened p by n array representing the transpose of matrix m2
 * @param[out] res Flattened m by p array for storing result m1*m2
 */
void matrix_product(size_t m, size_t n, size_t p, const double m1[m][n], const double m2t[p][n], double res[m][p]);

/**
 * @brief Negate a matrix.
 *
 * @param[in] m Height of matrix
 * @param[in] n Width of matrix
 * @param[in] mat Flattened m by n array representing matrix
 * @param[out] res Flattened m by n array for storing result -mat
 * @note mat and res can be the same arrays
 */
void matrix_negate(size_t m, size_t n, const double mat[m][n], double res[m][n]);

/**
 * @brief Transpose a matrix.
 *
 * @param[in] m Height of matrix
 * @param[in] n Width of matrix
 * @param[in] mat Flattened m by n array representing matrix
 * @param[out] res Flattened n by m array for storing result mat^T
 */
void matrix_transpose(size_t m, size_t n, const double mat[m][n], double res[n][m]);

#endif
