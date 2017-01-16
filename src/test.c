#include "main.h"

void main() {
	fmpz_t p;
	fmpz_init(p);

	fmpz_set_ui(p,2);
	fq_ctx_t field;

	slong d = 2;
	
	fq_ctx_init(field, p, d, "X");

	fq_ctx_print(field);

	fq_t res, x;
	fq_init(res, field);
	fq_init(x, field);
	fq_one(x, field);

	fq_poly_t X, P, Q;
	fq_poly_init(X, field);
	fq_poly_one(X, field);

	fq_poly_init2(P, 3, field);
	fq_poly_init2(Q, 3, field);

	for (slong i = 0; i < 3; i++) {
		fq_poly_set_coeff(P, i, x, field);
	}

	for (slong i = 1; i < 3; i++) {
		fq_poly_set_coeff(Q, i, x, field);
	}

	fq_poly_factor_t f,g;
	fq_poly_factor_init(f, field);

	fq_poly_factor_insert(f, P, 1, field);
	fq_poly_factor_insert(f, Q, 1, field);
	fq_poly_factor_print_pretty(f, "X", field);

	fq_poly_factor_init(g, field);

	factor_refinement(g, f, field);
	fq_poly_factor_print_pretty(g, "X", field);

	flint_printf("\n");

	fq_clear(res, field);
	fmpz_clear(p);
	fq_ctx_clear(field);
}
