#include "main.h"

/*
 * Find a random normal element. 
 *
 * We take u at random in F_{p^n}, we compute g(u) and
 * we test if g(u) is normal, where
 *
 * 	g(x) = f(x) / [ (x - r) * f'(r) ]
 * 	
 * and f is an irreducible polynomial over F_p, r is a
 * root of f in F_{p^n}, and f' is the derivative of f.
 *
 * Currently we take f defining F_{p^n} and root
 * the image of X in F_p[X]/(f).
 *
 */

void normal_random(fq_t res, const fq_ctx_t field) {

	// We initialize root to X
	fq_t root;
	fq_init(root, field);
	fq_gen(root, field);

	// And we initialize f to the modulus defining the field
	slong d = fq_ctx_degree(field);
	fq_poly_t f;
	fq_poly_init2(f, d + 1, field);

	for (slong i = 0; i < d + 1; i++) {
		fq_poly_set_coeff_fmpz(f, i, (field->modulus)->coeffs + i, field);
	}

	// We create f'
	fq_poly_t fprime;
	fq_poly_init(fprime, field);
	fq_poly_derivative(fprime, f, field);

	// We create a random state
	flint_rand_t state;
	flint_randinit(state);

	// We create and initialize some temporary variables and u
	fq_t u,tmp1,tmp2;
	fq_init(u, field);
	fq_init(tmp1, field);
	fq_init(tmp2, field);

	// We take u at random in F_p
    fmpz_t r;
    fmpz_init(r);
    fmpz_randm(r, state, fq_ctx_prime(field));
    fq_set_fmpz(u, r, field);

	// tmp1 = g(u)
	fq_sub(tmp1, u, root, field);
	fq_poly_evaluate_fq(tmp2, fprime, root, field);
	fq_mul(tmp1,tmp1,tmp2,field);
	fq_inv(tmp1, tmp1, field);
	fq_poly_evaluate_fq(tmp2, f, u, field);
	fq_mul(tmp1, tmp1, tmp2, field);

	// while g(u) is not normal, we take an other u and we do the same
	while (!(is_normal(tmp1, field))) {

        fmpz_randm(r, state, fq_ctx_prime(field));
        fq_set_fmpz(u, r, field);

		fq_sub(tmp1, u, root, field);
		fq_poly_evaluate_fq(tmp2, fprime, root, field);
		fq_mul(tmp1,tmp1,tmp2,field);
		fq_inv(tmp1, tmp1, field);
		fq_poly_evaluate_fq(tmp2, f, u, field);
		fq_mul(tmp1, tmp1, tmp2, field);
	}

	// Finally g(u) = tmp1 is normal, so we set res to tmp1
	fq_set(res, tmp1, field);

	// and we clear our variables
	fq_poly_clear(fprime, field);
	flint_randclear(state);
	fq_clear(u,field);
	fq_clear(tmp1,field);
	fq_clear(tmp2,field);
    fmpz_clear(r);
}
