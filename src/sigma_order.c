#include "main.h"

/*
 * Compute the sigma order of x in F_q^n / F_q.
 *
 * Only available for prime fields extensions.
 *
 * References ::
 * - Gao's PhD thesis (3.2)
 */

void sigma_order(fq_poly_t Q, fq_t x, fq_ctx_t fpn) {
	slong k,i,n;
        k = 1;
	n = fq_ctx_degree(fpn);

	fq_mat_t M, Mcopy, Mtmp;
	fq_mat_init(M,n,1, fpn);
	fq_mat_init(Mtmp,n,1, fpn);

	fq_t tmp, xcopy;
	fq_init(tmp, fpn);
	fq_init(xcopy, fpn);
	fq_set(xcopy, x, fpn);

	fmpz *coef;
	coef = x->coeffs;

	for (i = 0; i < n; i++) {
		fq_set_fmpz(tmp, coef + i, fpn);
		fq_mat_entry_set(M, i, 0, tmp, fpn);
	}

	fq_mat_init_set(Mcopy, M, fpn);

	while (fq_mat_rref(Mcopy, fpn) == k) {

		fq_mat_init_set(Mcopy, M, fpn);
		fq_mat_init(M, n, k+1, fpn);

		fq_frobenius(xcopy, xcopy, 1, fpn);
		coef = xcopy->coeffs;

		for (i = 0; i < n; i++) {

			fq_set_fmpz(tmp, coef + i, fpn);
			fq_mat_entry_set(Mtmp, i, 0, tmp, fpn);
		}

		fq_mat_concat_horizontal(M, Mcopy, Mtmp, fpn);
		fq_mat_init_set(Mcopy, M, fpn);

		k++;
	}

	fq_poly_t P;
	fq_poly_init2(P, k, fpn);

	for (i = 0; i < k-1; i++) {
		fq_poly_set_coeff(P, i, fq_mat_entry(Mcopy,i,k-1), fpn);
	}


	fq_poly_neg(P, P, fpn);
	fq_one(tmp, fpn);
	fq_poly_set_coeff(P, k-1, tmp, fpn); 
	fq_poly_set(Q, P, fpn);

	fq_clear(tmp, fpn);
	fq_clear(xcopy, fpn);
	fq_poly_clear(P, fpn);
	fq_mat_clear(M, fpn);
	fq_mat_clear(Mtmp, fpn);
}
