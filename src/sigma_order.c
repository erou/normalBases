#include "main.h"

/*
 * Compute the sigma order of x in F_{p^n} / F_p 
 *
 * The sigma order of x is a polynomial X^k - \sum_{i=0}^{k-1} c_i X^i
 * where k is the least positive integer such that f^k(x) belongs
 * to the F_p linear span of {f^0(x), f^1(x), ... , f^{k-1}(x)}, f is the
 * frobenius homomorphism, and the c_i's are the coefficient such that
 * f^k(x) = \sum c_i f^i(x).
 *
 * References ::
 * - Gao's PhD thesis (3.2)
 */

void sigma_order(fq_poly_t res, const fq_t x, const fq_ctx_t field) {

	// We treat the degenerate case x = 0
	if (fq_is_zero(x, field)) {
		fq_poly_one(res, field);
	}

	else {
		// We define some variables
		slong k,i,n;
		k = 1;
		n = fq_ctx_degree(field);

		// We define matrices
		fq_mat_t M, Mcopy, Mtmp;
		fq_mat_init(M,n,1, field);
		fq_mat_init(Mtmp,n,1, field);

		// And we define elements of F_{p^n}
		fq_t tmp, xcopy;
		fq_init(tmp, field);
		fq_init(xcopy, field);

		//xcopy is a copy of x
		fq_set(xcopy, x, field);

		// coef is the list of the coefficients of x in
		// the basis E = {1, X, ..., X^{n-1}} of F_{p^n}
		fmpz *coef;
		coef = x->coeffs;

		// M is a column matrix containing the coefficient of x in E
		for (i = 0; i < n; i++) {
			fq_set_fmpz(tmp, coef + i, field);
			fq_mat_entry_set(M, i, 0, tmp, field);
		}

		// Mcopy is a copy of M used to compute the rank of M
		// without changing M
		fq_mat_init_set(Mcopy, M, field);

		/* In order to know the least k and the coefficient, 
		 * we will use a reduced row echelon form :
		 * if we have a matrix (C1 C2 C3 C4) with columns Ci that are
		 * coefficients of vi, and (C1 C2 C3) of rank 3,
		 * then the rank will when adding C4 will show if v4 is in the linear
		 * span of {v1, v2, v3}, and the coefficients are those of the last
		 * column after performing the reduced row echelon form */
		while (fq_mat_rref(Mcopy, field) == k) {

			// We copy M into Mcopy
			fq_mat_init_set(Mcopy, M, field);

			// And then we add a column to M
			fq_mat_init(M, n, k+1, field);

			// We compute f^k(x)
			fq_frobenius(xcopy, xcopy, 1, field);
			coef = xcopy->coeffs;

			// And set its coefficients in E to the column matrix Mtmp
			for (i = 0; i < n; i++) {

				fq_set_fmpz(tmp, coef + i, field);
				fq_mat_entry_set(Mtmp, i, 0, tmp, field);
			}

			// We concatenate Mcopy and Mtmp
			// In other words, we add a column to M with the
			// coefficients of f^k(x)
			fq_mat_concat_horizontal(M, Mcopy, Mtmp, field);

			// We copy the new M in Mcopy to compute its rank
			fq_mat_init_set(Mcopy, M, field);

			k++;
		}

		// We define P and allocate memory for enough coefficients
		fq_poly_t P;
		fq_poly_init2(P, k, field);

		// We set the coefficients of P to the c_i's
		for (i = 0; i < k-1; i++) {
			fq_poly_set_coeff(P, i, fq_mat_entry(Mcopy,i,k-1), field);
		}

		// We set P to X^k - \sum c_i X^i
		fq_poly_neg(P, P, field);
		fq_one(tmp, field);
		fq_poly_set_coeff(P, k-1, tmp, field); 

		// We set res to P (our result)
		fq_poly_set(res, P, field);

		// And we clear our variables
		fq_clear(tmp, field);
		fq_clear(xcopy, field);
		fq_poly_clear(P, field);
		fq_mat_clear(M, field);
		fq_mat_clear(Mtmp, field);
	}
}
