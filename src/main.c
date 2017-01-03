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
	fq_pow_ui(Y, Y, 3, fpn);
	fq_add(Y, Y, one, fpn);

	fq_poly_t P;
	fq_poly_init(P, fpn);

	sigma_order(P, Y, fpn);
	fq_poly_print_pretty(P, "X", fpn);
	flint_printf("\n");

	fmpz_clear(p);
	fq_ctx_clear(fp);
	fq_ctx_clear(fpn);
	fq_poly_clear(P, fpn);
	fq_clear(one, fpn);
	fq_clear(Y, fpn);
}
