#include <stdio.h>
#include <stdint.h>
#include <stddef.h>
#include <string.h>
#include <stdlib.h>

#include "qp_ramp.h"
#include "vector.h"
#include "matrix.h"
#include "csv.h"
#include "lti.h"
#include "timing.h"

int main(int argc, char *argv[]) {
    char *input_dir = argv[1];
    char *output_dir = argv[2];
    size_t simulation_timesteps = atoi(argv[3]);
    printf("Input directory: \"%s\" \nOutput directory: \"%s\" \nSimulation timesteps: %ld \n", input_dir, output_dir, simulation_timesteps);

    char input_path[80];
    timing_reset();
    sprintf(input_path, "%s/x0.csv", input_dir);
	size_t initial_conditions = csv_parse_matrix_height(input_path);
    sprintf(input_path, "%s/a.csv", input_dir);
	size_t n_dim = csv_parse_matrix_width(input_path);
    sprintf(input_path, "%s/b.csv", input_dir);
	size_t m_dim = csv_parse_matrix_width(input_path);
    sprintf(input_path, "%s/g.csv", input_dir);
	size_t c_dim = csv_parse_matrix_height(input_path);
	size_t p_dim  = csv_parse_matrix_width(input_path);
    printf("Input dimension parsing time: %ld us\n", timing_elapsed()/1000);

    timing_reset();
    double *a = (double*)malloc(n_dim*n_dim*sizeof(double));
    double *b = (double*)malloc(n_dim*m_dim*sizeof(double));
    double *x0 = (double*)malloc(initial_conditions*n_dim*sizeof(double));
    double *invh = (double*)malloc(p_dim*p_dim*sizeof(double));
	double *w = (double*)malloc(c_dim*sizeof(double));
    double *g = (double*)malloc(c_dim*p_dim*sizeof(double));
    double *s = (double*)malloc(c_dim*n_dim*sizeof(double));
    double *f = (double*)malloc(p_dim*n_dim*sizeof(double));

    double *ft = (double*)malloc(n_dim*p_dim*sizeof(double));
	double *neg_w = (double*)malloc(c_dim*sizeof(double));
    double *neg_s = (double*)malloc(c_dim*n_dim*sizeof(double));
    double *neg_invh = (double*)malloc(p_dim*p_dim*sizeof(double)); // Only used once for initial setup
    double *neg_invh_f = (double*)malloc(m_dim*n_dim*sizeof(double)); // Actually p_dim x n_dim, but we only need the m_dim first rows
    double *neg_g_invh = (double*)malloc(c_dim*p_dim*sizeof(double));
    double *neg_g_invh_gt = (double*)malloc(c_dim*c_dim*sizeof(double));
        
    double *x = (double*)malloc((simulation_timesteps+1)*n_dim*sizeof(double));
    double *u = (double*)malloc(simulation_timesteps*m_dim*sizeof(double));
	double *t = (double*)malloc(simulation_timesteps*sizeof(double));

#ifdef REFERENCE_DIR
    double *x_ref = (double*)malloc((simulation_timesteps+1)*n_dim*sizeof(double));
    double *u_ref = (double*)malloc(simulation_timesteps*m_dim*sizeof(double));
#endif
    printf("Allocation time: %ld us\n", timing_elapsed()/1000);

    timing_reset();
    sprintf(input_path, "%s/a.csv", input_dir);
    if (csv_parse_matrix(input_path, n_dim, n_dim, (double(*)[])a)) { 
        printf("Error while parsing input matrix a.\n"); 
        return 1;
    }
    sprintf(input_path, "%s/b.csv", input_dir);
    if (csv_parse_matrix(input_path, n_dim, m_dim, (double(*)[])b)) { 
        printf("Error while parsing input matrix b.\n"); 
        return 1;
    }
    sprintf(input_path, "%s/x0.csv", input_dir);
    if (csv_parse_matrix(input_path, initial_conditions, n_dim, (double(*)[])x0)) { 
        printf("Error while parsing input matrix x0.\n"); 
        return 1;
    }
    sprintf(input_path, "%s/invh.csv", input_dir);
    if (csv_parse_matrix(input_path, p_dim, p_dim, (double(*)[])invh)) { 
        printf("Error while parsing input matrix invh.\n"); 
        return 1; 
    }
    sprintf(input_path, "%s/w.csv", input_dir);
    if (csv_parse_vector(input_path, c_dim, w)) { 
        printf("Error while parsing input vector w.\n"); 
        return 1; 
    }
    sprintf(input_path, "%s/g.csv", input_dir);
    if (csv_parse_matrix(input_path, c_dim, p_dim, (double(*)[])g)) { 
        printf("Error while parsing input matrix g.\n"); 
        return 1; 
    }
    sprintf(input_path, "%s/s.csv", input_dir);
    if (csv_parse_matrix(input_path, c_dim, n_dim, (double(*)[])s)) { 
        printf("Error while parsing input matrix s.\n"); 
        return 1; 
    }
    sprintf(input_path, "%s/f.csv", input_dir);
    if (csv_parse_matrix(input_path, p_dim, n_dim, (double(*)[])f)) { 
        printf("Error while parsing input matrix f.\n"); 
        return 1; 
    }
    printf("Input parsing time: %ld us\n", timing_elapsed()/1000);

    // Other initialization
    timing_reset();
    negate_vector(c_dim, w, neg_w);
    negate_matrix(c_dim, n_dim, (double(*)[])s, (double(*)[])neg_s);
    transpose(p_dim, n_dim, (double(*)[])f, (double(*)[])ft);
    negate_matrix(p_dim, p_dim, (double(*)[])invh, (double(*)[])neg_invh);
    matrix_product(m_dim, p_dim, n_dim, (double(*)[])neg_invh, (double(*)[])ft, (double(*)[])neg_invh_f);
    matrix_product(c_dim, p_dim, p_dim, (double(*)[])g, (double(*)[])neg_invh, (double(*)[])neg_g_invh); // Exploiting the fact that invh is symmetric
    matrix_product(c_dim, p_dim, c_dim, (double(*)[])g, (double(*)[])neg_g_invh, (double(*)[])neg_g_invh_gt);
    matrix_product(c_dim, p_dim, m_dim, (double(*)[])g, (double(*)[])neg_invh, (double(*)[])neg_g_invh); // Make sure memory layout is correct for later use
    // Free memory that's no longer needed
    free(invh);
	free(w);
    free(g);
    free(s);
    free(f);
    free(ft);
    free(neg_invh);
    qp_ramp_init(c_dim, p_dim);
    // qp_ramp_enable_infeasibility_error(1e-12, 1e12); // Remove comment to enable infeasibility errors
    printf("Initialization time: %ld us\n", timing_elapsed()/1000);

    // Simulation
    double total_time = 0.0;
    for (uint16_t i = 0; i < initial_conditions; ++i) {
        // Initial state
        memcpy(&x[0], &x0[i*n_dim], sizeof(double)*n_dim);

        double test_case_time = 0.0;
        for (uint16_t j = 0; j < simulation_timesteps; ++j) {
            timing_reset();
            int err = qp_ramp_solve_mpc(c_dim, n_dim, m_dim, p_dim, (double(*)[])neg_g_invh_gt, (double(*)[])neg_s, neg_w, (double(*)[])neg_invh_f, (double(*)[])neg_g_invh, &x[j*n_dim], &u[j*m_dim]);
            if (err) {
                printf("ERROR: %d\n", err);
            }
            simulate(n_dim, m_dim, (double(*)[])a, &x[j*n_dim], (double(*)[])b, &u[j*m_dim], &x[(j+1)*n_dim]); 
            t[j] = (double)timing_elapsed();
            test_case_time += t[j];
        }
        total_time += test_case_time;

        // Test case result printing
        // printf("Total simulation time for %d iterations of test case %d: %.0f us\n", simulation_timesteps, i, test_case_time/1000);
        // print_vector(n_dim, x[simulation_timesteps]);

        char save_path[80];
        // Save test case results
        sprintf(save_path, "%s/xout%d.csv", output_dir, i);
        if (csv_save_matrix(save_path, simulation_timesteps+1, n_dim, (double(*)[])x) < 0) {
            printf("Error while saving x.\n");
        }
        sprintf(save_path, "%s/uout%d.csv", output_dir, i);
        if (csv_save_matrix(save_path, simulation_timesteps, m_dim, (double(*)[])u) < 0) {
            printf("Error while saving u.\n");
        }
        sprintf(save_path, "%s/tout%d.csv", output_dir, i);
        if (csv_save_vector(save_path, simulation_timesteps, t) < 0) {
            printf("Error while saving t.\n");
        }
    }

    // Summary
    printf("Average time per test case running %ld timesteps: %.0f us\n", simulation_timesteps, total_time/1000/initial_conditions);

    // Free the rest of the memory
    free(a);
    free(b);
    free(x0);
	free(neg_w);
    free(neg_s);
    free(neg_invh_f);
    free(neg_g_invh);
    free(neg_g_invh_gt);
    free(x);
    free(u);
	free(t);
    qp_ramp_cleanup();

    return 0;
}
