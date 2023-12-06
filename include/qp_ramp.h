#ifndef QP_RAMP_H
#define QP_RAMP_H

#include <stddef.h>

#define QPRAMP_ERROR_INFEASIBLE -1
#define QPRAMP_ERROR_RANK_2_UPDATE -2

/**
 * @brief Allocate internal memory for the solver.
 *
 * @param[in] c The number of constraints
 * @param[in] p The number of optimization variables
 */
void qp_ramp_init(size_t c, size_t p);

/**
 * @brief Deallocate internal memory used by the solver.
 */
void qp_ramp_cleanup(void);

/**
 * @brief Enable checks for infeasibility that use specified limits.
 *
 * Checks are done by examining the value of the divisor q_0+v_i, 
 * with q_0 and v_i as defined in Hovd and Valmorbida (2023).
 * Enabling the error makes the solver return QPRAMP_ERROR_INFEASIBLE 
 * if the limits are breached.
 *
 * @param[in] min Minimum absolute value of divisor
 * @param[in] max Maximum absolute value of divisor
 */
void qp_ramp_enable_infeasibility_error(double min, double max);

/**
 * @brief Solve QP.
 *
 * @param[in] c The number of constraints
 * @param[in] n The number of system states
 * @param[in] p The number of optimization variables
 * @param[in] neg_g_invh_gt Flattened c by c array representing -G*H^{-1}*G^T
 * @param[in] neg_s Flattened c by n array representing -S
 * @param[in] neg_w Array with length c representing -w
 * @param[in] neg_g_invh Flattened c by p array representing -G*H^{-1}
 * @param[in] x Array with length n representing current system state
 * @param[out] z Array with length p representing optimization variable
 * @return 0 on success, 
 *         QPRAMP_ERROR_INFEASIBLE if infeasibility errors are enabled and detected,
 *         QPRAMP_ERROR_RANK_2_UPDATE if a rank 2 update fails
 */
int qp_ramp_solve(size_t c, size_t n, size_t p, const double neg_g_invh_gt[c][c], const double neg_s[c][n], const double neg_w[c], const double neg_g_invh[c][p], const double x[n], double z[p]);

/**
 * @brief Solve QP for MPC.
 *
 * @param[in] c The number of constraints
 * @param[in] n The number of system states
 * @param[in] m The number of system inputs
 * @param[in] p The number of optimization variables
 * @param[in] neg_g_invh_gt Flattened c by c array representing -G*H^{-1}*G^T
 * @param[in] neg_s Flattened c by n array representing -S
 * @param[in] neg_w Array with length c representing -w
 * @param[in] neg_invh_f Flattened m by n array representing -H^{-1}*F
 * @param[in] neg_g_invh Flattened c by m array representing -G*H^{-1}
 * @param[in] x Array with length n representing current system state
 * @param[out] u Array with length m representing optimal system input
 * @return 0 on success, 
 *         QPRAMP_ERROR_INFEASIBLE if infeasibility errors are enabled and detected,
 *         QPRAMP_ERROR_RANK_2_UPDATE if a rank 2 update fails
 */
int qp_ramp_solve_mpc(size_t c, size_t n, size_t m, size_t p, const double neg_g_invh_gt[c][c], const double neg_s[c][n], const double neg_w[c], const double neg_invh_f[m][n], const double neg_g_invh[c][m], const double x[n], double u[m]);

#endif
