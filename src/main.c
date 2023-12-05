#include <stdio.h>
#include <stddef.h>
#include <string.h>
#include <stdlib.h>

#include "qp_ramp.h"
#include "vector.h"
#include "matrix.h"
#include "csv.h"
#include "lti.h"
#include "timing.h"

// Convenience macros to cast pointers to variable length array type, avoiding compiler warnings
#define VLA_CAST(ptr) (double(*)[])ptr
#define CONST_VLA_CAST(ptr) (const double(*)[])ptr

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
    printf("Allocation time: %ld us\n", timing_elapsed()/1000);

    timing_reset();
    sprintf(input_path, "%s/a.csv", input_dir);
    if (csv_parse_matrix(input_path, n_dim, n_dim, VLA_CAST(a))) { 
        printf("Error while parsing input matrix a.\n"); 
        return 1;
    }
    sprintf(input_path, "%s/b.csv", input_dir);
    if (csv_parse_matrix(input_path, n_dim, m_dim, VLA_CAST(b))) { 
        printf("Error while parsing input matrix b.\n"); 
        return 1;
    }
    sprintf(input_path, "%s/x0.csv", input_dir);
    if (csv_parse_matrix(input_path, initial_conditions, n_dim, VLA_CAST(x0))) { 
        printf("Error while parsing input matrix x0.\n"); 
        return 1;
    }
    sprintf(input_path, "%s/invh.csv", input_dir);
    if (csv_parse_matrix(input_path, p_dim, p_dim, VLA_CAST(invh))) { 
        printf("Error while parsing input matrix invh.\n"); 
        return 1; 
    }
    sprintf(input_path, "%s/w.csv", input_dir);
    if (csv_parse_vector(input_path, c_dim, w)) { 
        printf("Error while parsing input vector w.\n"); 
        return 1; 
    }
    sprintf(input_path, "%s/g.csv", input_dir);
    if (csv_parse_matrix(input_path, c_dim, p_dim, VLA_CAST(g))) { 
        printf("Error while parsing input matrix g.\n"); 
        return 1; 
    }
    sprintf(input_path, "%s/s.csv", input_dir);
    if (csv_parse_matrix(input_path, c_dim, n_dim, VLA_CAST(s))) { 
        printf("Error while parsing input matrix s.\n"); 
        return 1; 
    }
    sprintf(input_path, "%s/f.csv", input_dir);
    if (csv_parse_matrix(input_path, p_dim, n_dim, VLA_CAST(f))) { 
        printf("Error while parsing input matrix f.\n"); 
        return 1; 
    }
    printf("Input parsing time: %ld us\n", timing_elapsed()/1000);

    // Other initialization
    timing_reset();
    vector_negate(c_dim, w, neg_w);
    matrix_negate(c_dim, n_dim, CONST_VLA_CAST(s), VLA_CAST(neg_s));
    matrix_transpose(p_dim, n_dim, CONST_VLA_CAST(f), VLA_CAST(ft));
    matrix_negate(p_dim, p_dim, CONST_VLA_CAST(invh), VLA_CAST(neg_invh));
    matrix_product(m_dim, p_dim, n_dim, CONST_VLA_CAST(neg_invh), CONST_VLA_CAST(ft), VLA_CAST(neg_invh_f));
    matrix_product(c_dim, p_dim, p_dim, CONST_VLA_CAST(g), CONST_VLA_CAST(neg_invh), VLA_CAST(neg_g_invh)); // Exploiting the fact that invh is symmetric
    matrix_product(c_dim, p_dim, c_dim, CONST_VLA_CAST(g), CONST_VLA_CAST(neg_g_invh), VLA_CAST(neg_g_invh_gt));
    matrix_product(c_dim, p_dim, m_dim, CONST_VLA_CAST(g), CONST_VLA_CAST(neg_invh), VLA_CAST(neg_g_invh)); // Make sure memory layout is correct for later use
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
    for (size_t i = 0; i < initial_conditions; ++i) {
        // Initial state
        memcpy(&x[0], &x0[i*n_dim], sizeof(double)*n_dim);

        double test_case_time = 0.0;
        for (size_t j = 0; j < simulation_timesteps; ++j) {
            timing_reset();
            int err = qp_ramp_solve_mpc(c_dim, n_dim, m_dim, p_dim, CONST_VLA_CAST(neg_g_invh_gt), CONST_VLA_CAST(neg_s), neg_w, CONST_VLA_CAST(neg_invh_f), CONST_VLA_CAST(neg_g_invh), &x[j*n_dim], &u[j*m_dim]);
            if (err) {
                printf("ERROR: %d\n", err);
            }
            lti_simulate(n_dim, m_dim, CONST_VLA_CAST(a), &x[j*n_dim], CONST_VLA_CAST(b), &u[j*m_dim], &x[(j+1)*n_dim]); 
            t[j] = (double)timing_elapsed();
            test_case_time += t[j];
        }
        total_time += test_case_time;

        // Save test case results
        char save_path[80];
        sprintf(save_path, "%s/xout%ld.csv", output_dir, i);
        if (csv_save_matrix(save_path, simulation_timesteps+1, n_dim, CONST_VLA_CAST(x)) < 0) {
            printf("Error while saving x.\n");
        }
        sprintf(save_path, "%s/uout%ld.csv", output_dir, i);
        if (csv_save_matrix(save_path, simulation_timesteps, m_dim, CONST_VLA_CAST(u)) < 0) {
            printf("Error while saving u.\n");
        }
        sprintf(save_path, "%s/tout%ld.csv", output_dir, i);
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
