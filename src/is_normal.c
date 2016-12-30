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
 */

int is_normal(fq_t elem, fq_ctx_t F_qn, fq_ctx_t F_q) {
	fq_poly_t P,Q,gcd; 
	slong n,i;
	n = fq_ctx_degree(F_qn) / fq_ctx_degree(F_q);
	fq_poly_init2(P, n, F_qn);
	fq_poly_init2(Q, n+1, F_qn);

	fq_t coef;
	fq_init(coef, F_qn);

	fq_one(coef, F_qn);

	fq_poly_set_coeff(Q, 0, coef, F_qn);
	fq_poly_set_coeff(Q, n, coef, F_qn);

	fq_set(coef, elem, F_qn);

	fq_poly_set_coeff(P, 0, elem, F_qn);

	for (i = 1; i < n; i++) {
		fq_frobenius(coef, coef, fq_ctx_degree(F_q), F_qn);
		fq_poly_set_coeff(P, i, coef, F_qn);
	}
	
	fq_poly_init(gcd, F_qn);
	fq_poly_gcd(gcd, P, Q, F_qn);

	fq_poly_clear(P, F_qn);
	fq_poly_clear(Q, F_qn);
	fq_clear(coef, F_qn);

	return fq_poly_is_one(gcd, F_qn);

	fq_poly_clear(gcd, F_qn);
}
