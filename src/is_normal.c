#include "main.h"

/*
 * Test if the element `elem` is normal in F_{p^n}/F_p. 
 *
 * For that purpose, we will check if Q =  X^n - 1 and
 * P are relatively primes. P is a polynomial that
 * depends only on `elem`.
 *
 * P = \sum_{i=0}^{n-1} elem^{p^{i-1}}X^i
 *
 */

int is_normal(const fq_t elem, const fq_ctx_t F_pn) {

	// We create the variables
	fq_poly_t P,Q,gcd; 
	slong n,i;
	n = fq_ctx_degree(F_pn);

	// We initialize the polynomiales, deg P = n - 1 and deg Q = n
	fq_poly_init2(P, n, F_pn);
	fq_poly_init2(Q, n+1, F_pn);

	// We create a variable coef belonging to F_{p^n}
	fq_t coef;
	fq_init(coef, F_pn);

	// coef = 1
	fq_one(coef, F_pn);

	// We set Q to X^n
	fq_poly_set_coeff(Q, n, coef, F_pn);

	// coef = -1
	fq_neg(coef, coef, F_pn);

	// We set Q to X^n - 1
	fq_poly_set_coeff(Q, 0, coef, F_pn);

	// coef = elem
	fq_set(coef, elem, F_pn);

	// We set P to elem
	fq_poly_set_coeff(P, 0, elem, F_pn);


	/* And using the Frobenius homomorphism, we set P to the polynomial
	 described at the beginning. */
	for (i = 1; i < n; i++) {
		fq_frobenius(coef, coef, 1, F_pn);
		fq_poly_set_coeff(P, i, coef, F_pn);
	}

	// We set gcd to gcd(P, Q)
	fq_poly_init(gcd, F_pn);
	fq_poly_gcd(gcd, P, Q, F_pn);

	// We clear the variables, except gcd that is needed
	fq_poly_clear(P, F_pn);
	fq_poly_clear(Q, F_pn);
	fq_clear(coef, F_pn);

	// We return true if gcd = 1
	return fq_poly_is_one(gcd, F_pn);

	fq_poly_clear(gcd, F_pn);
}
