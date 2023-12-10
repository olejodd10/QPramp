#ifndef TIMING_H
#define TIMING_H

/**
 * @brief Print the precision of the timer.
 *
 */
void timing_print_precision(void);

/**
 * @brief Reset the timer.
 *
 */
void timing_reset(void);

/**
 * @brief Get the number of nanoseconds since the timer was last reset.
 *
 * @return The number of nanoseconds since the timer was last reset
 */
long timing_elapsed(void);

#endif
