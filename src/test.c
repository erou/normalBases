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
	fq_gen(x, field);

	fq_poly_t X;
	fq_poly_init(X, field);
	fq_poly_one(X, field);

	frobenius_composition(x, X, x, field);
	fq_print_pretty(x, field);
	flint_printf("\n---------\n");

	luneburg(res, field);
	fq_print_pretty(res, field);
	flint_printf("\n");

	fq_clear(res, field);
	fmpz_clear(p);
	fq_ctx_clear(field);
}
