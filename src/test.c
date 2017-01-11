#include "main.h"

void main() {
	fmpz_t p, uno;
	fmpz_init(p);
	fmpz_init(uno);

	fmpz_set_ui(p,2);
	fmpz_set_ui(uno,1);

	fq_ctx_t field;

	slong d = 2;
	
	fq_ctx_init(field, p, d, "X");

	fq_ctx_print(field);

	fq_poly_factor_t fac;
	fq_poly_factor_init(fac, field);
	
	fq_poly_t P, Q, X;
	fq_poly_init(X, field);
	fq_poly_init(P, field);
	fq_poly_init(Q, field);
	fq_poly_gen(X, field);
	fq_poly_pow(Q, X, 2, field);
	fq_poly_add(P, X, Q, field);
	fq_poly_one(Q, field);
	fq_poly_add(Q, Q, X, field);

	fq_poly_struct *fact;
	fact = fac->poly;

	fq_poly_factor_insert(fac, P, 1, field);
	fq_poly_factor_insert(fac, Q, 1, field);
	fq_poly_factor_print_pretty(fac, "X", field);

	factor_refinement(fac, field);

	fq_poly_factor_print_pretty(fac, "X", field);

	fq_poly_clear(P, field);
	fmpz_clear(p);
	fq_ctx_clear(field);
}
