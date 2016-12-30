#include "main.h"

/*
 * Test if the element `elem` is normal in F_{q^n}/F_q. 
 *
 * For that purpose, we will check if Q =  X^n - 1 and
 * P are relatively primes. P is a polynomial that
 * depends only on `elem`.
 *
 * P = \sum_{i=0}^{n-1} elem^{q^{i-1}}X^i
 *
 * The present implementation works ONLY for normal elements over
 * prime fields. It could be changed by creating a function taking
 * two contexts as input.
 */

int is_normal(fq_t elem, fq_ctx_t ctx) {
	fq_poly_t P,Q,gcd; 
	slong n,i;
	n = fq_ctx_degree(ctx);
	fq_poly_init2(P, n, ctx);
	fq_poly_init2(Q, n+1, ctx);

	fq_t coef;
	fq_init(coef, ctx);

	fq_one(coef, ctx);

	fq_poly_set_coeff(Q, 0, coef, ctx);
	fq_poly_set_coeff(Q, n, coef, ctx);

	fq_set(coef, elem, ctx);

	fq_poly_set_coeff(P, 0, elem, ctx);

	for (i = 1; i < n; i++) {
		fq_frobenius(coef, coef, 1, ctx);
		fq_poly_set_coeff(P, i, coef, ctx);
	}
	
	fq_poly_init(gcd, ctx);
	fq_poly_gcd(gcd, P, Q, ctx);

	fq_poly_clear(P, ctx);
	fq_poly_clear(Q, ctx);
	fq_clear(coef, ctx);

	return fq_poly_is_one(gcd, ctx);

	fq_poly_clear(gcd, ctx);
}
