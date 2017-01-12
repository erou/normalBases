#include "main.h"

void factor_refinement(fq_poly_factor_t fac, const fq_ctx_t field) {

	// we set some indices 
	slong i, j, k;

	// We initialise the gcd and quotient polynomials
	fq_poly_t gcd, quo;
	fq_poly_init(gcd, field);
	fq_poly_init(quo, field);

	// We use an other factor to copy fac
	fq_poly_factor_t copy;

	/* We start the while loop : if the factors are
	 * some f_i, we check that f_i is coprime with
	 * any f_j for j>i. If they are not coprime, and
	 * their gcd is d, we delete f_i and f_j of the 
	 * factors, and we add f_i/d, f_j/d and d. We do
	 * that until f_i and f_j are coprime for any i,j.
	 */
	i = 0;
	while (i < fac->num - 1) {
		for (j = i + 1; j < fac->num; j++) {

			// We compute the gcd of f_i and f_j
			fq_poly_gcd(gcd, fac->poly + i, fac->poly + j, field);
			// If it is not 1, we add f_i/gcd, f_j/gcd and gcd to
			// the factors, unless f_i/gcd = 1
			if (!fq_poly_is_one(gcd, field)) {

				if (!fq_poly_equal(gcd, fac->poly + i, field)) {
					fq_poly_divides(quo, fac->poly + i, gcd, field);
					fq_poly_factor_insert(fac, quo, 1, field);
				}

				if (!fq_poly_equal(gcd, fac->poly + j, field)) {
					fq_poly_divides(quo, fac->poly + j, gcd, field);
					fq_poly_factor_insert(fac, quo, 1, field);
				}

				fq_poly_factor_init(copy, field);
				fq_poly_factor_insert(copy, gcd, 1, field);
				// and we delete f_i and f_j
				// in fact : we copy all the list but f_i and
				// f_j : we do not know how to do better
				for (k = 0; k < fac->num; k++) {
					if (k!=i && k!=j) {
						fq_poly_factor_insert(copy, fac->poly + k, 1, field);
					}
				}

				// and we set fac to copy, which is fac
				// without f_i and f_j
				fq_poly_factor_set(fac, copy, field);				
				break;
			}
		}

		// if we went to j = n, this means that f_i
		// is coprime with any f_j, so we do not touch
		// this f_i until the end
		if (j == fac->num) {
			i++;
		}
	}
}
