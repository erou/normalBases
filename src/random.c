#include "main.h"

/*
 * Find a random normal element. 
 *
 * We take u at random in F_{p^n}, we compute g(u) and
 * we test if g(u) is normal, where
 *
 * 	g(x) = f(x) / [ (x - r) * f'(r) ]
 * 	
 * and f is an irreducible polynomial over F_p, r is a
 * root of f in F_{p^n}, and f' is the derivative of f.
 *
 */

void find_normal_random(fq_t res, const fq_poly_t f, const fq_t root, const fq_ctx_t F_pn) {

	// We create f'
	fq_poly_t fprime;
	fq_poly_init(fprime, F_pn);
	fq_poly_derivative(fprime, f, F_pn);

	// We create a random state
	flint_rand_t state;
	flint_randinit(state);

	// We create and initialize some temporary variables and u
	fq_t u,tmp1,tmp2;
	fq_init(u, F_pn);
	fq_init(tmp1, F_pn);
	fq_init(tmp2, F_pn);

	// We take u at random
	fq_randtest(u, state, F_pn);

	// But not equal to root, not to divise by zero
	while (fq_equal(u, root, F_pn)) {
		fq_randtest(u, state, F_pn);
	}

	// tmp1 = g(u)
	fq_sub(tmp1, u, root, F_pn);
	fq_poly_evaluate_fq(tmp2, fprime, root, F_pn);
	fq_mul(tmp1,tmp1,tmp2,F_pn);
	fq_inv(tmp1, tmp1, F_pn);
	fq_poly_evaluate_fq(tmp2, f, u, F_pn);
	fq_mul(tmp1, tmp1, tmp2, F_pn);

	// while g(u) is not normal, we take an other u and we do the same
	while (!(is_normal(tmp1, F_pn))) {
		fq_randtest(u, state, F_pn);

		while (fq_equal(u, root, F_pn)) {
			fq_randtest(u, state, F_pn);
		}

		fq_sub(tmp1, u, root, F_pn);
		fq_poly_evaluate_fq(tmp2, fprime, root, F_pn);
		fq_mul(tmp1,tmp1,tmp2,F_pn);
		fq_inv(tmp1, tmp1, F_pn);
		fq_poly_evaluate_fq(tmp2, f, u, F_pn);
		fq_mul(tmp1, tmp1, tmp2, F_pn);
	}

	// Finally g(u) = tmp1 is normal, so we set res to tmp1
	fq_set(res, tmp1, F_pn);

	// and we clear our variables
	fq_poly_clear(fprime, F_pn);
	flint_randclear(state);
	fq_clear(u,F_pn);
	fq_clear(tmp1,F_pn);
	fq_clear(tmp2,F_pn);
}
