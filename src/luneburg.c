#include "main.h"

/* 	LUNEBURG
 *  Compute a normal element using Lüneburg algorithm.
 *  This algorithm is deterministic.
 *
 *  We first compute the sigma order polynomials f_i of
 *  alpha^i (i = 0 .. d-1) where alpha is a generator of 
 *  the field F_{p^d}. Then we apply factor refinement to
 *  the f_i's, to obtain polynomials g_j (j = 1 .. r) pairwise
 *  coprimes, and such that for each i, we have
 *  f_i = \prod g_j^e_ij. Then for each j, we find i(j) such that
 *  e_i(j)j is maximal, and we compute h_j = f_i(j) / g_j^e_i(j)j
 *  and beta_j = h_j(f)(alpha^i(j)) (where f is the frobenius map).
 *  Then beta = \sum beta_j is a normal element.
 *
 *  References : Gao's PhD thesis.
 */

void luneburg(fq_t res, const fq_ctx_t field) {

	// We declare future indices and constants
	slong d, i, j, n, m, exp, expmax, ind;
	d = fq_ctx_degree(field);

	// We initialize some elements of the field
	fq_t X, a, tmp, rop;
	fq_init(X, field);
	fq_init(a, field);
	fq_init(tmp, field);
	fq_init(rop, field);

	// And we set X to the generator of field
	fq_gen(X, field);

	// We initialise some factor, this is just a 
	// way to represent the list of polynomials
	// f_i and g_j
	fq_poly_factor_t f, g;
	fq_poly_factor_init(f, field);	
	fq_poly_factor_init(g, field);	

	// We initialise the polynomials h and P (a temporary polynomial)
	fq_poly_t P, h;
	fq_poly_init(P, field);
	fq_poly_init(h, field);

    // We create a list of polynomials
	fq_poly_struct *F;
	F = (fq_poly_struct*)malloc(d*sizeof(fq_poly_struct));

	// We create f_0, and we insert it in f
	fq_poly_init(F, field);
	fq_one(a, field);
	sigma_order(F, a, field);
	fq_poly_factor_insert(f, F, 1, field);


	// And we do the same for f_i (i = 1 .. d-1)
	for (i = 1; i < d; i++) {
		fq_poly_init(F + i, field);
		fq_mul(a, a, X, field);
		sigma_order(F + i, a, field);
		fq_poly_factor_insert(f, F + i, 1, field);
	}

	// We create the g_j
	factor_refinement(g, f, field);
	n = g->num;

    fq_poly_factor_print_pretty(f, "Z", field);
    flint_printf("\n\n");
    fq_poly_factor_print_pretty(g, "Z", field);
    flint_printf("\n\n");

	/* In the loop, we compute the indice
	 * i(j) such that e_i(j)j is maximal and we
	 * add beta_j to rop
	 */
	for (j = 0; j < n; j++) {
		exp = 0;
		expmax = 0;
		ind = 0;

		// We find i(j)
		for (i = 0; i < d; i++) {
            exp = highest_exp(F + i, g->poly + j, field);
			if (expmax < exp) {
				expmax = exp;
				ind = i;
			}
		}

		// And we compute beta_j
		fq_poly_pow(P, g->poly + j, expmax, field);
		fq_poly_divides(h, F + ind, P, field);
		fq_pow_ui(tmp, X, ind, field);
		frobenius_composition(tmp, h, tmp, field);

		// Then we add it to rop
		fq_add(rop, rop, tmp, field);
	}

	fq_set(res, rop, field);

	// Then we clear all the variables
	fq_clear(X, field);
	fq_clear(a, field);
	fq_clear(tmp, field);
	fq_clear(rop, field);
	fq_poly_clear(P, field);
	fq_poly_clear(h, field);
	fq_poly_factor_clear(f, field);
	fq_poly_factor_clear(g, field);
    free(F);
}
