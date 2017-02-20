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

void lenstra(fq_t res, const fq_ctx_t field) {

	// We set some indices and the degree of the extension
	slong d, b, i, j, k;
    d = fq_ctx_degree(field);

	// We set some temporary variables and constants of the fields
	fq_t theta, zeta, X, tmp, tmp2;
	fq_init(theta, field);
	fq_init(zeta, field);
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

	// P = X^d - 1
	fq_poly_init(P, field);
	fq_poly_set_coeff(P, d, tmp, field); 
	fq_neg(tmp, tmp, field);
	fq_poly_set_coeff(P, 0, tmp, field); 

	/* We finally set an n × (n+1) matrix M, M will be
	 * the concatenation of the matrix of the linear map
	 * g(σ), where σ is the frobenius map, and of the 
	 * coordinates of theta in the basis 1, ..., X^(d-1)
	 */
	fq_mat_t M;
	fq_mat_init(M, d, d+1, field);

	sigma_order(ordTheta, theta, field);

	// if ordTheta = X^d - 1, then theta is normal
	// so the algorithm stops
	while (!fq_poly_equal(ordTheta, P, field)) {

		// g = (X^d - 1) / Ord_θ
		fq_poly_divides(g, P, ordTheta, field);

		// We compute the matrix of g(σ)
		for (i = 0; i < d; i++) {
			fq_pow_ui(tmp, X, i, field);
			frobenius_composition(tmp, g, tmp, field);
			for (j = 0; j < tmp->length; j++) {
				fq_set_fmpz(tmp2, tmp->coeffs + j, field);
				fq_mat_entry_set(M, j, i, tmp2, field);
			}
		}

		// And the coodinates of θ
		for (i = 0; i < theta->length; i++) {
			fq_set_fmpz(tmp, theta->coeffs + i, field);
			fq_mat_entry_set(M, i, d, tmp, field);
		}

		// We solve M×B=T with unknown B, where T represents 
		// the coordinates of θ, 
		k = fq_mat_rref(M, field);

		// We set tmp2 to the solution β
		// (with coordinates B) of M×B=T
		fq_zero(tmp2, field);

        for (i = 0; i < k; i++) {
            for (j = i; j < d; j++) {
                if (!fq_is_zero(fq_mat_entry(M, i, j), field)) {
			        fq_pow_ui(tmp, X, j, field);
			        fq_mul(tmp, fq_mat_entry(M, i, d), tmp, field);
			        fq_add(tmp2, tmp2, tmp, field);
                    break;
                }
            }
        }

		// We compute Ord_β
		sigma_order(ordBeta, tmp2, field);

		/* If the degree of ord_θ is higher than the degree
		 * of ord_β, we change θ for β = tmp2 and we
		 * continue in the while loop. If not, we look for an 
		 * element ζ such that g(σ)(ζ) = 0 and we set
		 * θ to θ + ζ.
		 */
		if (ordTheta->length < ordBeta->length) {
				fq_set(theta, tmp2, field);
			}

		else {
            b = 1;

            for (i = 0; i < k; i++) {
                if (fq_is_zero(fq_mat_entry(M, i, i), field)) {
				    fq_pow_ui(zeta, X, i, field);
                    b = 0;
                    break;
                }
            }


            if (b) {
	    		for (i = 0; i < k; i++) {
		    		fq_pow_ui(tmp, X, i, field);
			    	fq_neg(tmp2, fq_mat_entry(M, i, d-1), field);
				    fq_mul(tmp, tmp, tmp2, field);
    				fq_add(zeta, zeta, tmp, field);
	    		}
    
	    		if (k < d) {
		    	    fq_pow_ui(tmp, X, d-1, field);
			    	fq_add(zeta, zeta, tmp, field);
    			}

            }

	    	fq_add(theta, theta, zeta, field);
            fq_zero(zeta, field);
		}

        fq_mat_zero(M, field);
        sigma_order(ordTheta, theta, field);
	}

	fq_set(res, theta, field);

	// We clear all the variables
	fq_clear(theta, field);
	fq_clear(zeta, field);
	fq_clear(X, field);
	fq_clear(tmp, field);
	fq_clear(tmp2, field);
	fq_poly_clear(ordTheta, field);
	fq_poly_clear(ordBeta, field);
	fq_poly_clear(P, field);
	fq_poly_clear(g, field);
	fq_mat_clear(M, field);
}
