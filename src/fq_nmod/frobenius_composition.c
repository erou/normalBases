#include "main.h"

/* 	FROBENIUS COMPOSITION
 *
 *  Compute h(σ)(x) where σ is the frobenius map,
 *  h = \sum h_i X^i is a polynomial, and x is 
 *  an element of the field.
 */

void frobenius_composition(fq_nmod_t res, const fq_nmod_poly_t h, const fq_nmod_t x, const fq_nmod_ctx_t field) {

	// We initialize some temporary variables and copy of x
	fq_nmod_t tmp, xcopy;
	fq_nmod_init(tmp, field);
	fq_nmod_init(xcopy, field);
	fq_nmod_set(xcopy, x, field);

	// We set res to 0
	fq_nmod_zero(res, field);

	// then we add h_0 * x to res
	fq_nmod_mul(tmp, h->coeffs, xcopy, field);
	fq_nmod_add(res, res, tmp, field);

	// and we add h_i(σ^i(x)) to res (for i = 1 to i = degree(h))
	for (slong i = 1; i < h->length; i++) {
		fq_nmod_frobenius(xcopy, xcopy, 1, field);
		fq_nmod_mul(tmp, h->coeffs + i, xcopy, field);
		fq_nmod_add(res, res, tmp, field);
	}

	// And we clear the variables
	fq_nmod_clear(tmp, field);
	fq_nmod_clear(xcopy, field);
}
