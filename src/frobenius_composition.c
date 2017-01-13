#include "main.h"

/* 	FROBENIUS COMPOSITION
 *
 *  Compute h(f)(x) where f is the frobenius map,
 *  h = \sum h_i X^i is a polynomial, and x is 
 *  an element of the field.
 */

void frobenius_composition(fq_t res, const fq_poly_t h, const fq_t x, const fq_ctx_t field) {

	// We initialize some temporary variables and copy of x
	fq_t tmp, xcopy;
	fq_init(tmp, field);
	fq_init(xcopy, field);
	fq_set(xcopy, x, field);

	// We set res to 0
	fq_zero(res, field);

	// then we add h_0 * x to res
	fq_mul(tmp, h->coeffs, xcopy, field);
	fq_add(res, res, tmp, field);

	// and we add h_i(f^i(x)) to res (for i = 1 to i = degree(h))
	for (slong i = 1; i < h->length; i++) {
		fq_frobenius(xcopy, xcopy, 1, field);
		fq_mul(tmp, h->coeffs + i, xcopy, field);
		fq_add(res, res, tmp, field);
	}

	// And we clear the variables
	fq_clear(tmp, field);
	fq_clear(xcopy, field);
}
