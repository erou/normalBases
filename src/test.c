#include "main.h"

void main() {
	fmpz_t p;
	fmpz_init(p);

	fmpz_set_ui(p,3);

	fq_ctx_t fp, fpn;

	slong d = 1;
	slong n = 2;
	
	fq_ctx_init(fp, p, d, "X");
	fq_ctx_init(fpn, p, n, "Y");

	fq_ctx_print(fp);
	fq_ctx_print(fpn);

	fq_t one, Y;

	fq_init(one, fpn);
	fq_init(Y, fpn);
	fq_one(one, fpn);
	fq_gen(Y, fpn);

	fq_t alpha, beta, a;

	fq_init(alpha, fpn);
	fq_init(beta, fpn);
	fq_init(a, fpn);

	fq_add(a, Y, one, fpn);
	fq_add(a, a, one, fpn);

	fq_add(alpha, one, Y, fpn);
	fq_set(beta, Y, fpn);

	fq_t tmp;
	fq_init(tmp, fp);

	fmpz *calpha;
	fmpz *cbeta;
	fmpz *ca;
	calpha = alpha->coeffs;
	cbeta = beta->coeffs;
	ca = a->coeffs;

	fq_mat_t M;
	fq_mat_init(M, 2, 3, fp);

	slong i;
	for (i = 0; i < 2; i++) {
		fq_set_fmpz(tmp, calpha + i, fp);
		fq_mat_entry_set(M, i, 0, tmp, fp);
	}

	for (i = 0; i < 2; i++) {
		fq_set_fmpz(tmp, cbeta+ i, fp);
		fq_mat_entry_set(M, i, 1, tmp, fp);
	}

	for (i = 0; i < 2; i++) {
		fq_set_fmpz(tmp, ca + i, fp);
		fq_mat_entry_set(M, i, 2, tmp, fp);
	}

	fq_mat_print_pretty(M, fp);

	fq_mat_rref(M, fp);

	fq_mat_print_pretty(M, fp);

	fmpz_clear(p);
	fq_ctx_clear(fp);
	fq_ctx_clear(fpn);
	fq_clear(one, fpn);
	fq_clear(alpha, fpn);
	fq_clear(beta, fpn);
	fq_clear(a, fpn);
	fq_mat_clear(M, fpn);
}
