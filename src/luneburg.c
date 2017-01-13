#include "main.h"

void luneburg(fq_t res, const fq_ctx_t field) {
	fq_zero(res, field);
	slong d, i, j, n, exp, expmax, ind;
	d = fq_ctx_degree(field);
	fq_t X, a, tmp;
	fq_init(X, field);
	fq_init(a, field);
	fq_init(tmp, field);
	fq_gen(X, field);

	fq_poly_factor_t f, g, fcopy;
	fq_poly_factor_init(f, field);	
	fq_poly_factor_init(fcopy, field);	
	fq_poly_factor_init(g, field);	

	fq_poly_t P, h;
	fq_poly_init(P, field);
	fq_poly_init(h, field);

	fq_one(a, field);
	sigma_order(P, a, field);
	fq_poly_factor_insert(f, P, 1, field);

	for (i = 1; i < d; i++) {
		fq_mul(a, a, X, field);
		sigma_order(P, a, field);
		fq_poly_factor_insert(f, P, 1, field);
	}

	fq_poly_factor_set(fcopy, f, field);

	factor_refinement(g, f, field);
	n = g->num;

	for (j = 0; j < n; j++) {
		exp = 0;
		expmax = 0;
		ind = 0;
		for (i = 0; i < d; i++) {
			exp = fq_poly_remove(fcopy->poly + i, g->poly + j, field);
			if (expmax < exp) {
				expmax = exp;
				ind = i;
			}
		}
		fq_poly_pow(P, g->poly + j, expmax, field);
		fq_poly_divides(h, f->poly + ind, P, field);
		fq_pow_ui(tmp, X, ind, field);
		frobenius_composition(tmp, h, tmp, field);
		fq_add(res, res, tmp, field);
	}

	fq_clear(X, field);
	fq_clear(a, field);
	fq_clear(tmp, field);
	fq_poly_clear(P, field);
	fq_poly_clear(h, field);
	fq_poly_factor_clear(f, field);
	fq_poly_factor_clear(g, field);
	fq_poly_factor_clear(fcopy, field);
}
