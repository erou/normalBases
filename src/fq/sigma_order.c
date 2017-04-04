#include "main.h"

/*
 * Compute the sigma order of x in F_{p^n} / F_p 
 *
 * The sigma order of x is a polynomial X^k - \sum_{i=0}^{k-1} c_i X^i
 * where k is the least positive integer such that σ^k(x) belongs
 * to the F_p linear span of {σ^0(x), σ^1(x), ... , σ^{k-1}(x)}, σ is the
 * frobenius homomorphism, and the c_i's are the coefficients such that
 * σ^k(x) = \sum c_i σ^i(x).
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
		slong k, i, d;
		k = 1;
		d = fq_ctx_degree(field);

		// We define matrices
		fq_mat_t M;
		fq_mat_init(M, d, d+1, field);

		// And we define elements of F_{p^d}
		fq_t tmp, xcopy;
		fq_init(tmp, field);
		fq_init(xcopy, field);

		//xcopy is a copy of x
		fq_set(xcopy, x, field);

		// M is a matrix containing the coefficient of x in the
		// polynomial basis E = 1, ... , X^(d-1)
		for (i = 0; i < xcopy->length; i++) {
			fq_set_fmpz(tmp, xcopy->coeffs + i, field);
			fq_mat_entry_set(M, i, 0, tmp, field);
		}

		/* In order to know the least k and the coefficients, 
		 * we will use a reduced row echelon form :
		 * if we have a matrix (C1 C2 C3 C4) with columns Ci that are
		 * coefficients of vi, and (C1 C2 C3) of rank 3,
		 * then the rank will when adding C4 will show if v4 is in the linear
		 * span of {v1, v2, v3}, and the coefficients are those of the last
		 * column after performing the reduced row echelon form */
		while (fq_mat_rank(M, field) == k) {

			// We compute σ^k(x)
			fq_frobenius(xcopy, xcopy, 1, field);

			// And set its coefficients in E to the column k of the matrix M
			for (i = 0; i < xcopy->length; i++) {

				fq_set_fmpz(tmp, xcopy->coeffs + i, field);
				fq_mat_entry_set(M, i, k, tmp, field);
			}

			k++;
		}

        // We compute the reduced row echelon form of M
        fq_mat_rref(M, field);

		// We define P and allocate memory for enough coefficients
		fq_poly_t P;
		fq_poly_init(P, field);

		// We set the coefficients of P to the c_i's
		for (i = 0; i < k-1; i++) {
			fq_poly_set_coeff(P, i, fq_mat_entry(M,i,k-1), field);
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
	}
}
