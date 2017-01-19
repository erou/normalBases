#include "main.h"

void main() {
	fmpz_t p;
	fmpz_init(p);

	fmpz_set_ui(p, 2);
	slong d = 4;
	fq_ctx_t field;

	fq_poly_t P;
	fq_poly_init(P, field);

	
	fq_ctx_init(field, p, d, "X");

	fq_ctx_print(field);

	fq_t res, x;
	fq_init(res, field);

	luneburg(res, field);
	fq_print_pretty(res, field);
	flint_printf("\n");
	sigma_order(P, res, field);
	fq_print_pretty(res, field);
	flint_printf("\n");
	fq_poly_print_pretty(P, "Y", field);
	flint_printf("\n");

	if (is_normal(res, field)) {
		flint_printf("is_normal répond oui\n");
	}

	fq_print_pretty(res, field);
	flint_printf("\n");

	luneburg(res, field);
	fq_print_pretty(res, field);
	flint_printf("\n");
	sigma_order(P, res, field);
	fq_print_pretty(res, field);
	flint_printf("\n");
	fq_poly_print_pretty(P, "Y", field);
	flint_printf("\n");

	if (is_normal(res, field)) {
		flint_printf("is_normal répond oui\n");
	}

	fq_clear(res, field);
	fmpz_clear(p);
	fq_poly_clear(P, field);
	fq_ctx_clear(field);
}
