#ifndef VECTOR_H
#define VECTOR_H

#include <stddef.h>

/**
 * @brief Compute inner product of two vectors.
 *
 * @param[in] n Vector length
 * @param[in] v1 Array with length n representing a vector
 * @param[in] v2 Array with length n representing a vector
 * @return The inner product of v1 and v2
 */
double vector_inner_product(size_t n, const double v1[n], const double v2[n]);

/**
 * @brief Negate a vector.
 *
 * @param[in] n Vector length
 * @param[in] vec Array with length n representing a vector
 * @param[out] res Array with length n for storing result -vec
 * @note vec and res can be the same arrays
 */
void vector_negate(size_t n, const double vec[n], double res[n]);

/**
 * @brief Scale a vector.
 *
 * @param[in] n Vector length
 * @param[in] vec Array with length n representing a vector
 * @param[in] c Factor to scale vec by
 * @param[out] res Array with length n for storing result c*vec
 * @note vec and res can be the same arrays
 */
void vector_scale(size_t n, const double vec[n], double c, double res[n]);

/**
 * @brief Add two vectors.
 *
 * @param[in] n Vector length
 * @param[in] v1 Array with length n representing a vector
 * @param[in] v2 Array with length n representing a vector
 * @param[out] res Array with length n for storing result v1+v2
 * @note v1, v2 and res can be the same arrays
 */
void vector_add(size_t n, const double v1[n], const double v2[n], double res[n]);

/**
 * @brief Add a vector and a scaled vector.
 *
 * @param[in] n Vector length
 * @param[in] v1 Array with length n representing a vector
 * @param[in] v2 Array with length n representing a vector
 * @param[in] c Factor to scale v2 by
 * @param[out] res Array with length n for storing result v1+c*v2
 * @note v1, v2 and res can be the same arrays
 */
void vector_add_scaled(size_t n, const double v1[n], const double v2[n], double c, double res[n]);

#endif
