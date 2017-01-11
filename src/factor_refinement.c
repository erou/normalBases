#include "main.h"

void factor_refinement(fq_poly_factor_t fac, const fq_ctx_t field) {

	// we set some indices and n = length of fac
	slong n, i, j, k;
	n = fac->num;

	fq_poly_struct *factors;
	factors = fac->poly;

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
	while (i < n) {
		flint_printf("i = %wd\n", i);
		for (j = i + 1; j < n; j++) {
			flint_printf("j = %wd\n", j);

			// We compute the gcd of f_i and f_j
			fq_poly_gcd(gcd, factors + i, factors + j, field);
			// If it is not 1, we add f_i/gcd, f_j/gcd and gcd to
			// the factors, unless f_i/gcd = 1
			if (!fq_poly_is_one(gcd, field)) {

				if (!fq_poly_equal(gcd, factors + i, field)) {
					fq_poly_divides(quo, factors + i, gcd, field);
					fq_poly_factor_insert(fac, quo, 1, field);
				}

				if (!fq_poly_equal(gcd, factors + j, field)) {
					fq_poly_divides(quo, factors + j, gcd, field);
					fq_poly_factor_insert(fac, quo, 1, field);
				}

				fq_poly_factor_init(copy, field);
				fq_poly_factor_insert(copy, gcd, 1, field);
				// and we delete f_i and f_j
				// in fact : we copy all the list but f_i and
				// f_j : we do not know how to do better
				for (k = 0; k < fac->num; k++) {
					if (k!=i && k!=j) {
						fq_poly_factor_insert(copy, factors + k, 1, field);
					}
				}

				// and we set fac to copy, which is fac
				// without f_i and f_j
				fq_poly_factor_set(fac, copy, field);				
				fq_poly_factor_print_pretty(fac, "Y", field);
				// we compute the new length of fac
				factors = fac->poly;
				n = fac->num;
				break;
			}
		}

		// if we went to j = n, this means that f_i
		// is coprime with any f_j, so we do not touch
		// this f_i until the end
		if (j == n) {
			i++;
		}
	}
}
