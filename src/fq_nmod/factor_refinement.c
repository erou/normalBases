#include "main.h"

/*
 * 	FACTOR REFINEMENT
 *
 * Given a factorisation M = \prod_i f_i^d_i, compute a new
 * factorisation M = \prod_j g_j^e_j where the g_j's are
 * pairwise coprime.
 *
 * This algorithm is not exactly the factor refinement, since
 * we do not care about the powers.
 *
 * Reference : Factor Refinement, Eric Bach, James Driscoll
 * and Jeffrey Shallit. Journal of Algorithms.
 */

void factor_refinement(fq_nmod_poly_factor_t res, const fq_nmod_poly_factor_t fac, const fq_nmod_ctx_t field) {

	// We set res to fac
	fq_nmod_poly_factor_set(res, fac, field);

	// We set some indices 
	slong i, j, k;

    // And a boolean value
    int b = 1;

	// We initialise the gcd and quotient polynomials
	fq_nmod_poly_t gcd, quo;
	fq_nmod_poly_init(gcd, field);
	fq_nmod_poly_init(quo, field);

	// We use an other factor to copy res
	fq_nmod_poly_factor_t copy, empty;
	fq_nmod_poly_factor_init(copy, field);
	fq_nmod_poly_factor_init(empty, field);

	/* We start the while loop : if the factors are
	 * some f_i, we check that f_i is coprime with
	 * any f_j for j>i. If they are not coprime, and
	 * their gcd is d, we delete f_i and f_j of the 
	 * factors, and we add f_i/d, f_j/d and d. We do
	 * that until f_i and f_j are coprime for any i,j.
	 */
	i = 0;
	while (i < res->num - 1) {
		for (j = i + 1; j < res->num; j++) {

			// We compute the gcd of f_i and f_j
			fq_nmod_poly_gcd(gcd, res->poly + i, res->poly + j, field);

			// If it is not 1, we add f_i/gcd, f_j/gcd and gcd to
			// the factors, unless f_i/gcd = 1
			if (!fq_nmod_poly_is_one(gcd, field)) {

                b = 0;

				if (!fq_nmod_poly_equal(gcd, res->poly + i, field)) {
					fq_nmod_poly_divides(quo, res->poly + i, gcd, field);
					fq_nmod_poly_factor_insert(res, quo, 1, field);
				}

				if (!fq_nmod_poly_equal(gcd, res->poly + j, field)) {
					fq_nmod_poly_divides(quo, res->poly + j, gcd, field);
					fq_nmod_poly_factor_insert(res, quo, 1, field);
				}

				// and we delete f_i and f_j
				// in fact : we copy all the list but f_i and
				// f_j : we do not know how to do better
				for (k = 0; k < res->num; k++) {
					if (k!=i && k!=j) {
						fq_nmod_poly_factor_insert(copy, res->poly + k, 1, field);
					}
				}

				fq_nmod_poly_factor_insert(copy, gcd, 1, field);

				// and we set res to copy, which is res
				// without f_i and f_j
				fq_nmod_poly_factor_set(res, copy, field);				
				fq_nmod_poly_factor_set(copy, empty, field);
				break;
			}
		}

		// if b = 1, it means we never entered in the previous if
        // so f_i is coprime with any f_j, therefore we do not touch
		// this f_i until the end
		if (b) {
            i++;
        }
        else {
            b = 1;
        }
	}

	// We clear the variables
	fq_nmod_poly_clear(gcd, field);
	fq_nmod_poly_clear(quo, field);
	fq_nmod_poly_factor_clear(copy, field);
	fq_nmod_poly_factor_clear(empty, field);
}
