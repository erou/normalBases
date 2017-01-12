#include "main.h"

void main() {
	fmpz_t p, uno, dos;
	fmpz_init(p);
	fmpz_init(uno);
	fmpz_init(dos);

	fmpz_set_ui(p,3);
	fmpz_set_ui(uno,1);
	fmpz_set_ui(dos,2);

	fq_ctx_t field;

	slong d = 2;
	
	fq_ctx_init(field, p, d, "X");

	fq_ctx_print(field);

	fq_poly_factor_t fac;
	fq_poly_factor_init(fac, field);
	
	fq_poly_t P, Q;
	fq_poly_init2(P, 5, field);
	fq_poly_set_coeff_fmpz(P, 1, uno, field);
	fq_poly_set_coeff_fmpz(P, 2, dos, field);
	fq_poly_set_coeff_fmpz(P, 3, dos, field);
	fq_poly_set_coeff_fmpz(P, 4, uno, field);

	fq_poly_init2(Q, 3, field);
	fq_poly_set_coeff_fmpz(Q, 1, dos, field);
	fq_poly_set_coeff_fmpz(Q, 2, uno, field);

	fq_poly_struct *fact;
	fact = fac->poly;

	fq_poly_factor_insert(fac, P, 1, field);
	fq_poly_factor_insert(fac, Q, 1, field);
	fq_poly_factor_print_pretty(fac, "X", field);

	factor_refinement(fac, field);

	fq_poly_factor_print_pretty(fac, "X", field);

	fq_poly_struct *tab;
	tab = (fq_poly_struct*)malloc(5*sizeof(fq_poly_struct));

	fq_poly_clear(P, field);
	fmpz_clear(p);
	fq_ctx_clear(field);
}
