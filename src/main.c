#include "main.h"

void main() {
	fmpz_t p;
	fmpz_init(p);

	fmpz_set_ui(p,2);

	fq_ctx_t fp, fpn;

	slong d = 2;
	slong n = 4;
	
	fq_ctx_init(fp, p, d, "X");
	fq_ctx_init(fpn, p, n, "Y");

	fq_ctx_print(fp);
	fq_ctx_print(fpn);

	fq_t one, Y;

	fq_init(one, fpn);
	fq_init(Y, fpn);
	fq_one(one, fpn);
	fq_gen(Y, fpn);

	fmpz *t;
	t = Y->coeffs;

	slong i;
	for (i = 0; i < 10; i++) {
		fmpz_print(t+i);
		flint_printf("\n");
	}

	fq_poly_t P;
	fq_poly_init2(P, n+1, fpn);

	fq_poly_set_coeff(P, 0, one, fpn);
	fq_poly_set_coeff(P, 1, one, fpn);
	fq_poly_set_coeff(P, 4, one, fpn);

	fq_t res;
	fq_init(res, fpn);

	find_normal_random(res, P, Y, fpn, fp);
	fq_print_pretty(res, fpn);
	flint_printf("\n");

	fmpz_clear(p);
	fq_ctx_clear(fp);
	fq_ctx_clear(fpn);
	fq_poly_clear(P, fpn);
	fq_clear(res, fpn);
	fq_clear(one, fpn);
	fq_clear(Y, fpn);
}
