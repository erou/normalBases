#include "main.h"

/*
 * Find a normal element.
 *
 * The algorithm is randomized.
 *
 */

void find_normal_random(fq_t res, fq_poly_t f, fq_t root, fq_ctx_t ctx) {
	fq_poly_t fprime;
	fq_poly_init(fprime, ctx);
	fq_poly_derivative(fprime, f, ctx);

	flint_rand_t state;
	flint_randinit(state);

	fq_t u,tmp1,tmp2;
	fq_init(u, ctx);
	fq_init(tmp1, ctx);
	fq_init(tmp2, ctx);
	fq_randtest(u, state, ctx);

	while (fq_equal(u, root, ctx)) {
		fq_randtest(u, state, ctx);
	}


	fq_sub(tmp1, u, root, ctx);
	fq_poly_evaluate_fq(tmp2, fprime, root, ctx);
	fq_mul(tmp1,tmp1,tmp2,ctx);
	fq_inv(tmp1, tmp1, ctx);
	fq_poly_evaluate_fq(tmp2, f, u, ctx);
	fq_mul(tmp1, tmp1, tmp2, ctx);

	while (!(is_normal(tmp1, ctx))) {
		fq_randtest(u, state, ctx);

		while (fq_equal(u, root, ctx)) {
			fq_randtest(u, state, ctx);
		}

		fq_sub(tmp1, u, root, ctx);
		fq_poly_evaluate_fq(tmp2, fprime, root, ctx);
		fq_mul(tmp1,tmp1,tmp2,ctx);
		fq_inv(tmp1, tmp1, ctx);
		fq_poly_evaluate_fq(tmp2, f, u, ctx);
		fq_mul(tmp1, tmp1, tmp2, ctx);
	}

	fq_set(res, tmp1, ctx);

	fq_poly_clear(fprime, ctx);
	flint_randclear(state);
	fq_clear(u,ctx);
	fq_clear(tmp1,ctx);
	fq_clear(tmp2,ctx);
}
