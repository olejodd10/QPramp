#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>

#include "qp_ramp.h"
#include "iterable_set.h"
#include "vector.h"
#include "matrix.h"
#include "csv.h"
#include "lti.h"
#include "timing.h"
#include "qpOASES_wrapper.h"

#define SIMULATION_TIMESTEPS 100

#define INPUT_DIR "../examples/osqp/n20"
#define OUTPUT_DIR INPUT_DIR "/out"
#define REFERENCE_DIR INPUT_DIR "/reference"

#define A_PATH INPUT_DIR "/a.csv"
#define B_PATH INPUT_DIR "/b.csv"
#define X0_PATH INPUT_DIR "/x0.csv"

#define QPOASES_H_PATH INPUT_DIR "/osqp_p.csv"
#define QPOASES_G_PATH INPUT_DIR "/osqp_q.csv"
#define QPOASES_LBA_PATH INPUT_DIR "/osqp_l.csv"
#define QPOASES_A_PATH INPUT_DIR "/osqp_a.csv"
#define QPOASES_UBA_PATH INPUT_DIR "/osqp_u.csv"

#define X_DIR OUTPUT_DIR
#define U_DIR OUTPUT_DIR
#define T_DIR OUTPUT_DIR

#ifdef REFERENCE_DIR
#define X_REF_DIR REFERENCE_DIR
#define U_REF_DIR REFERENCE_DIR
#define EQ_EPS 1e-8
#endif

int main() {
    timing_print_precision();
    timing_reset();
	size_t initial_conditions = csv_parse_matrix_height(X0_PATH);
	size_t n_dim = csv_parse_matrix_width(A_PATH);
	size_t m_dim = csv_parse_matrix_width(B_PATH);
	size_t c_dim = csv_parse_matrix_height(QPOASES_A_PATH);
	size_t p_dim  = csv_parse_matrix_width(QPOASES_A_PATH);
    size_t horizon = (p_dim - n_dim)/(n_dim + m_dim);
    printf("Input dimension parsing time: %ld us\n", timing_elapsed()/1000);

    timing_reset();
    double *a = (double*)malloc(n_dim*n_dim*sizeof(double));
    double *b = (double*)malloc(n_dim*m_dim*sizeof(double));
    double *x0 = (double*)malloc(initial_conditions*n_dim*sizeof(double));

    real_t *qpoases_h = (real_t*)malloc(p_dim*p_dim*sizeof(real_t));
    real_t *qpoases_a = (real_t*)malloc(c_dim*p_dim*sizeof(real_t));
    real_t *qpoases_g = (real_t*)malloc(p_dim*sizeof(real_t));
    real_t *qpoases_lba = (real_t*)malloc(c_dim*sizeof(real_t));
    real_t *qpoases_uba = (real_t*)malloc(c_dim*sizeof(real_t));
    real_t *xOpt = (real_t*)malloc(p_dim*sizeof(real_t));
    real_t *yOpt = (real_t*)malloc((c_dim+p_dim)*sizeof(real_t));

    double *x = (double*)malloc((SIMULATION_TIMESTEPS+1)*n_dim*sizeof(double));
    double *u = (double*)malloc(SIMULATION_TIMESTEPS*m_dim*sizeof(double));
	double *t = (double*)malloc(SIMULATION_TIMESTEPS*sizeof(double));

#ifdef REFERENCE_DIR
    double *x_ref = (double*)malloc((SIMULATION_TIMESTEPS+1)*n_dim*sizeof(double));
    double *u_ref = (double*)malloc(SIMULATION_TIMESTEPS*m_dim*sizeof(double));
#endif
    printf("Allocation time: %ld us\n", timing_elapsed()/1000);

    timing_reset();
    if (csv_parse_matrix(A_PATH, n_dim, n_dim, a)) { 
        printf("Error while parsing input matrix a.\n"); 
        return 1;
    }
    if (csv_parse_matrix(B_PATH, n_dim, m_dim, b)) { 
        printf("Error while parsing input matrix b.\n"); 
        return 1;
    }
    if (csv_parse_matrix(X0_PATH, initial_conditions, n_dim, x0)) { 
        printf("Error while parsing input matrix x0.\n"); 
        return 1;
    }

    if (csv_parse_matrix(QPOASES_H_PATH, p_dim, p_dim, qpoases_h)) { 
		printf("Error while parsing input matrix qpoases_h.\n");
		return 1;
	}
    if (csv_parse_matrix(QPOASES_A_PATH, c_dim, p_dim, qpoases_a)) { 
		printf("Error while parsing input matrix qpoases_a.\n");
		return 1;
	}
    if (csv_parse_vector(QPOASES_G_PATH, p_dim, qpoases_g)) { 
		printf("Error while parsing input vector qpoases_g.\n");
		return 1;
	}
    if (csv_parse_vector(QPOASES_LBA_PATH, c_dim, qpoases_lba)) { 
		printf("Error while parsing input vector qpoases_lba.\n");
		return 1;
	}
    if (csv_parse_vector(QPOASES_UBA_PATH, c_dim, qpoases_uba)) { 
		printf("Error while parsing input vector qpoases_uba.\n");
		return 1;
	}
    printf("Input parsing time: %ld us\n", timing_elapsed()/1000);

	qpOASES_Options options;
	qpOASES_Options_init( &options,2 ); // 0 for default, 1 for reliable and 2 for MPC
	options.printLevel = PL_NONE;

	QProblem_setup(	p_dim,c_dim,HST_POSDEF); // TODO

    real_t obj;
    int nWSR;
    int status;
    int initialized = 0; // Whether hot start can be used or not

    // Simulation
    double total_time = 0.0;
    for (uint16_t i = 0; i < initial_conditions; ++i) {
        // Initial state
        memcpy(&x[0], &x0[i*n_dim], sizeof(double)*n_dim);
        // This is redundant when i == 0 if qpoases_lba and qpoasas_uba are given correctly
        for (uint16_t k = 0; k < n_dim; ++k) {
            qpoases_lba[k] = -x[k];
            qpoases_uba[k] = -x[k];
        }

        double test_case_time = 0.0;
        for (uint16_t j = 0; j < SIMULATION_TIMESTEPS; ++j) {
            timing_reset();

            nWSR = c_dim+1; // TODO
            if (initialized) {
                int ret = QProblem_hotstart(	qpoases_g,NULL,NULL,qpoases_lba,qpoases_uba,
                        (int_t* const)&nWSR,0,
                        xOpt,yOpt,&obj,(int_t* const)&status
                        );
                // if (status == 0) {
                //     printf("Case %d timestep %d er vellykket\n", i, j);
                // }
                // printf("Hot start ret val and status: %d %d\n", ret, status);
            } else {
                int ret = QProblem_init(	qpoases_h,qpoases_g,qpoases_a,NULL,NULL,qpoases_lba,qpoases_uba,
                        (int_t* const)&nWSR,0,&options,
                        xOpt,yOpt,&obj,(int_t* const)&status
                        );
                // if (status == 0) {
                //     printf("Case %d timestep %d er vellykket\n", i, j);
                // }
                // printf("Init ret val and status: %d %d\n", ret, status);
                initialized = 1;
            }
            
            // /* Print solution of first QP. */	
            // printf( "\nxOpt = [ %e, %e ];  yOpt = [ %e, %e, %e ];  objVal = %e\n\n", 
            // 		xOpt[0],xOpt[1],yOpt[0],yOpt[1],yOpt[2], obj );


            double *u_temp = (double*)(&xOpt[n_dim*(horizon+1)]);

            // Apply input
            simulate(n_dim, m_dim, a, &x[j*n_dim], b, u_temp, &x[(j+1)*n_dim]); 

            // Update constraints
            for (uint16_t k = 0; k < n_dim; ++k) {
                qpoases_lba[k] = -x[(j+1)*n_dim + k];
                qpoases_uba[k] = -x[(j+1)*n_dim + k];
            }

            t[j] = (double)timing_elapsed();
            test_case_time += t[j];

            memcpy(&u[j*m_dim], u_temp, m_dim*sizeof(double));
        }
        total_time += test_case_time;

        // Test case result printing
        // printf("Total simulation time for %d iterations of test case %d: %.0f us\n", SIMULATION_TIMESTEPS, i, test_case_time/1000);
        // print_vector(n_dim, x[SIMULATION_TIMESTEPS]);

        char path[80];
        #ifdef REFERENCE_DIR
        // Verify test case results
        sprintf(path, X_REF_DIR "/xout%d.csv", i);
        FILE* f = fopen(path, "r");
        if (csv_parse_matrix(path, SIMULATION_TIMESTEPS+1, n_dim, x_ref)) {
            printf("Error while parsing reference matrix x_ref.\n"); 
        }
        if (!matrix_eq(SIMULATION_TIMESTEPS+1, n_dim, x, x_ref, EQ_EPS)) {
            printf("WARNING: Verification failed for x in test case %d\n", i);
        }
        fclose(f);

        sprintf(path, U_REF_DIR "/uout%d.csv", i);
        f = fopen(path, "r");
        if (csv_parse_matrix(path, SIMULATION_TIMESTEPS, m_dim, u_ref)) {
            printf("Error while parsing reference matrix u_ref.\n"); 
        }
        if (!matrix_eq(SIMULATION_TIMESTEPS, m_dim, u, u_ref, EQ_EPS)) {
            printf("WARNING: Verification failed for u in test case %d\n", i);
        }
        fclose(f);
        #endif

        // Save test case results
        sprintf(path, X_DIR "/xout%d.csv", i);
        if (csv_save_matrix(path, SIMULATION_TIMESTEPS+1, n_dim, x) < 0) {
            printf("Error while saving x.\n");
        }
        sprintf(path, U_DIR "/uout%d.csv", i);
        if (csv_save_matrix(path, SIMULATION_TIMESTEPS, m_dim, u) < 0) {
            printf("Error while saving u.\n");
        }
        sprintf(path, T_DIR "/tout%d.csv", i);
        if (csv_save_vector(path, SIMULATION_TIMESTEPS, t) < 0) {
            printf("Error while saving t.\n");
        }
    }

	QProblem_cleanup();

    // Summary
    printf("Average time per test case running %d timesteps: %.0f us\n", SIMULATION_TIMESTEPS, total_time/1000/initial_conditions);

    return 0;
}
