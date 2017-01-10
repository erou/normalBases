#include "main.h"

void main() {
	fmpz_t p;
	fmpz_init(p);

	fmpz_set_ui(p,5);

	fq_ctx_t field;

	slong d = 7;
	
	fq_ctx_init(field, p, d, "X");

	fq_ctx_print(field);

	fq_poly_t P;
	fq_poly_init2(P, d+1, field);


	for (slong i = 0; i < d+1; i++) {
		fq_poly_set_coeff_fmpz(P, i, (field->modulus)->coeffs +i, field);
	}

	fq_poly_print_pretty(P, "X", field);

	fq_poly_clear(P, field);
	fmpz_clear(p);
	fq_ctx_clear(field);
}
