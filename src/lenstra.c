#include "main.h"

/* 	LENSTRA
 *
 * Lenstra's algorithm for computing a normal element.
 *
 * The algorithm follow 4 steps, the explanation of the
 * is in the references. 
 * 
 *  References : Gao's phD thesis.
 */

void lenstra(fq_t res, const fq_ctx_t field){

	// We set some indices and the degree of the extension
	slong d, i, j, k;
        d = fq_ctx_degree(field);

	// We set some temporary variables and constants of the fields
	fq_t theta, eta, X, tmp, tmp2;
	fq_init(theta, field);
	fq_init(eta, field);
	fq_init(tmp, field);
	fq_init(tmp2, field);
	fq_init(X, field);

	fq_gen(theta, field);
	fq_gen(X, field);
	fq_one(tmp, field);

	// We set some polynomials 
	fq_poly_t ordTheta, ordBeta, P, g;
	fq_poly_init(ordTheta, field);
	fq_poly_init(ordBeta, field);
	fq_poly_init(g, field);

	// P = X^n - 1
	fq_poly_init2(P, d+1, field);
	fq_poly_set_coeff(P, d, tmp, field); 
	fq_neg(tmp, tmp, field);
	fq_poly_set_coeff(P, 0, tmp, field); 

	/* We finally set an n x (n+1) matrix M, M will be
	 * the concatenation of the matrix of the linear map
	 * g(f), where f is the frobenius map, and of the 
	 * coordinates of theta in the basis 1, ..., X^(d-1)
	 */
	fq_mat_t M;
	fq_mat_init(M, d, d+1, field);

	sigma_order(ordTheta, theta, field);

	// if ordTheta = X^d - 1, then theta is normal
	// so the algorithm stops
	while (!fq_poly_equal(ordTheta, P, field)) {

		// g = (X^d - 1) / Ord_\theta
		fq_poly_divides(g, P, ordTheta, field);

		// We compute the matrix of g(f)
		for (i = 0; i < d; i++) {
			fq_pow_ui(tmp, X, i, field);
			frobenius_composition(tmp, g, tmp, field);
			for (j = 0; j < d; j++) {
				fq_set_fmpz(tmp2, tmp->coeffs + j, field);
				fq_mat_entry_set(M, j, i, tmp2, field);
			}
		}

		// And the coodinates of \theta 
		for (i = 0; i < d; i++) {
			fq_set_fmpz(tmp, theta->coeffs + i, field);
			fq_mat_entry_set(M, i, d, tmp, field);
		}

		// We solve MxB=T with unknown B, where T represents 
		// the coordinates of \theta, 
		k = fq_mat_rref(M, field);

		// We set tmp2 to the solution \beta
		// (with coordinates B) of MxB=T
		fq_zero(tmp2, field);

		for (i = 0; i < d; i++) {
			fq_pow_ui(tmp, X, i, field);
			fq_mul(tmp, fq_mat_entry(M, i, d), tmp, field);
			fq_add(tmp2, tmp2, tmp, field);
		}

		fq_t test;
		fq_init(test, field);
		frobenius_composition(test, g, tmp2, field);
		fq_print_pretty(test, field);
		flint_printf(" =? ");
		fq_print_pretty(theta, field);
		flint_printf("\n");

		// We compute Ord_\beta
		sigma_order(ordBeta, tmp2, field);

		/* If the degree of ordTheta is higher than the degree
		 * of ordBeta, we change beta for theta = tmp2 and we
		 * continue in the while loop. If not, we look for an 
		 * element \eta such that g(f)(\eta) = 0 and we set
		 * \beta to \theta + \eta.
		 */
		if (ordTheta->length < ordBeta->length) {
				fq_set(theta, tmp2, field);
			}

		else {
			for (i = 0; i < k; i++) {
				fq_pow_ui(tmp, X, i, field);
				fq_neg(tmp2, fq_mat_entry(M, i, d-1), field);
				fq_mul(tmp, tmp, tmp2, field);
				fq_add(eta, eta, tmp, field);
			}

			if (k < d) {
			fq_pow_ui(tmp, X, d-1, field);
				fq_add(eta, eta, tmp, field);
			}
			fq_add(theta, theta, eta, field);
		}
	}
	fq_set(res, theta, field);

	// We clear all the variables
	fq_clear(theta, field);
	fq_clear(eta, field);
	fq_clear(X, field);
	fq_clear(tmp, field);
	fq_clear(tmp2, field);
	fq_poly_clear(ordTheta, field);
	fq_poly_clear(ordBeta, field);
	fq_poly_clear(P, field);
	fq_poly_clear(g, field);
	fq_mat_clear(M, field);
}
