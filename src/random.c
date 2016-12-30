#include "main.h"

/*
 * Find a normal element.
 *
 * The algorithm is randomized.
 *
 */

void find_normal_random(fq_t res, fq_poly_t f, fq_t root, fq_ctx_t F_qn, fq_ctx_t F_q) {
	fq_poly_t fprime;
	fq_poly_init(fprime, F_qn);
	fq_poly_derivative(fprime, f, F_qn);

	flint_rand_t state;
	flint_randinit(state);

	fq_t u,tmp1,tmp2;
	fq_init(u, F_qn);
	fq_init(tmp1, F_qn);
	fq_init(tmp2, F_qn);
	fq_randtest(u, state, F_qn);

	while (fq_equal(u, root, F_qn)) {
		fq_randtest(u, state, F_qn);
	}


	fq_sub(tmp1, u, root, F_qn);
	fq_poly_evaluate_fq(tmp2, fprime, root, F_qn);
	fq_mul(tmp1,tmp1,tmp2,F_qn);
	fq_inv(tmp1, tmp1, F_qn);
	fq_poly_evaluate_fq(tmp2, f, u, F_qn);
	fq_mul(tmp1, tmp1, tmp2, F_qn);

	while (!(is_normal(tmp1, F_qn, F_q))) {
		fq_randtest(u, state, F_qn);

		while (fq_equal(u, root, F_qn)) {
			fq_randtest(u, state, F_qn);
		}

		fq_sub(tmp1, u, root, F_qn);
		fq_poly_evaluate_fq(tmp2, fprime, root, F_qn);
		fq_mul(tmp1,tmp1,tmp2,F_qn);
		fq_inv(tmp1, tmp1, F_qn);
		fq_poly_evaluate_fq(tmp2, f, u, F_qn);
		fq_mul(tmp1, tmp1, tmp2, F_qn);
	}

	fq_set(res, tmp1, F_qn);

	fq_poly_clear(fprime, F_qn);
	flint_randclear(state);
	fq_clear(u,F_qn);
	fq_clear(tmp1,F_qn);
	fq_clear(tmp2,F_qn);
}
