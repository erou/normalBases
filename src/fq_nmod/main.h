/*
 * 	***    HEADERS    ***
 *
 * 	A small description of each algorithm is availaible in each file
 * 	`algo.c`. A more complete one should be availaible in the report. All
 * 	these algorithms are based on Gao's PhD thesis.
 */


// Standard headers

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>

// Flint headers

#include <fmpz.h>
#include <fq_nmod.h>
#include <fq_nmod_poly.h>
#include <fq_nmod_mat.h>
#include <profiler.h>

// Functions 
int is_normal(const fq_nmod_t, const fq_nmod_ctx_t);

void normal_random(fq_nmod_t, const fq_nmod_ctx_t);

void sigma_order(fq_nmod_poly_t, const fq_nmod_t, const fq_nmod_ctx_t);

void help(const char*);

void factor_refinement(fq_nmod_poly_factor_t, const fq_nmod_poly_factor_t, const fq_nmod_ctx_t);

void frobenius_composition(fq_nmod_t, const fq_nmod_poly_t, const fq_nmod_t, const fq_nmod_ctx_t);

slong highest_exp(const fq_nmod_poly_t, const fq_nmod_poly_t, const fq_nmod_ctx_t);

void luneburg(fq_nmod_t, const fq_nmod_ctx_t);

void lenstra(fq_nmod_t, const fq_nmod_ctx_t);

void naive(fq_nmod_t, const fq_nmod_ctx_t);
