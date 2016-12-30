#include "main.h"

/*
 * Test if the element `elem` is normal in F_{q^n}/F_q. 
 *
 * For that purpose, we will check if Q =  X^n - 1 and
 * P are relatively primes. P is a polynomial that
 * depends only on `elem`.
 *
 * P = \sum_{i=0}^{n-1} elem^{q^{i-1}}X^i
 *
 * The present implementation works ONLY for normal elements over
 * prime fields. It could be changed by creating a function taking
 * two contexts as input.
 */

int is_normal_test(fq_t elem, fq_ctx_t F_qn, fq_ctx_t F_q) {
	fq_poly_t P,Q,gcd; 
	slong n,i;
	n = fq_ctx_degree(F_qn) / fq_ctx_degree(F_q);
	fq_poly_init2(P, n, F_qn);
	fq_poly_init2(Q, n+1, F_qn);

	fq_t coef;
	fq_init(coef, F_qn);

	fq_one(coef, F_qn);

	fq_poly_set_coeff(Q, 0, coef, F_qn);
	fq_poly_set_coeff(Q, n, coef, F_qn);

	fq_set(coef, elem, F_qn);

	fq_poly_set_coeff(P, 0, elem, F_qn);

	for (i = 1; i < n; i++) {
		fq_frobenius(coef, coef, fq_ctx_degree(F_q), F_qn);
		fq_poly_set_coeff(P, i, coef, F_qn);
	}
	
	fq_poly_init(gcd, F_qn);
	fq_poly_gcd(gcd, P, Q, F_qn);

	fq_poly_clear(P, F_qn);
	fq_poly_clear(Q, F_qn);
	fq_clear(coef, F_qn);

	return fq_poly_is_one(gcd, F_qn);

	fq_poly_clear(gcd, F_qn);
}

/*
 * Find a normal element.
 *
 * The algorithm is randomized.
 *
 */
void find_normal_random_test(fq_t res, fq_poly_t f, fq_t root, fq_ctx_t F_qn, fq_ctx_t F_q) {
	fq_poly_t fprime;
	fq_poly_init(fprime, F_qn);
	fq_poly_derivative(fprime, f, F_qn);

	flint_rand_t state;
	flint_randinit(state);

	fq_t u,tmp1,tmp2;
	fq_init(u, F_qn);
	fq_init(tmp1, F_qn);
	fq_init(tmp2, F_qn);
	fq_randtest(u, state, F_qn);

	while (fq_equal(u, root, F_qn)) {
		fq_randtest(u, state, F_qn);
	}


	fq_sub(tmp1, u, root, F_qn);
	fq_poly_evaluate_fq(tmp2, fprime, root, F_qn);
	fq_mul(tmp1,tmp1,tmp2,F_qn);
	fq_inv(tmp1, tmp1, F_qn);
	fq_poly_evaluate_fq(tmp2, f, u, F_qn);
	fq_mul(tmp1, tmp1, tmp2, F_qn);

	while (!(is_normal_test(tmp1, F_qn, F_q))) {
		fq_randtest(u, state, F_qn);

		while (fq_equal(u, root, F_qn)) {
			fq_randtest(u, state, F_qn);
		}

		fq_sub(tmp1, u, root, F_qn);
		fq_poly_evaluate_fq(tmp2, fprime, root, F_qn);
		fq_mul(tmp1,tmp1,tmp2,F_qn);
		fq_inv(tmp1, tmp1, F_qn);
		fq_poly_evaluate_fq(tmp2, f, u, F_qn);
		fq_mul(tmp1, tmp1, tmp2, F_qn);
	}

	fq_set(res, tmp1, F_qn);

	fq_poly_clear(fprime, F_qn);
	flint_randclear(state);
	fq_clear(u,F_qn);
	fq_clear(tmp1,F_qn);
	fq_clear(tmp2,F_qn);
}

void main() {
	fmpz_t p;
	fmpz_init(p);

	fmpz_set_ui(p,2);

	fq_ctx_t fp, fpn;

	slong d = 1;
	slong n = 4;
	
	fq_ctx_init(fp, p, d, "X");
	fq_ctx_init(fpn, p, n, "Y");

	fq_ctx_print(fp);
	fq_ctx_print(fpn);

	fq_t one, Y;

	fq_init(one, fpn);
	fq_init(Y, fpn);
	fq_one(one, fpn);
	fq_gen(Y, fpn);

	fmpz *t;
	t = Y->coeffs;

	slong i;
	for (i = 0; i < 10; i++) {
		fmpz_print(t+i);
		flint_printf("\n");
	}

	fq_poly_t P;
	fq_poly_init2(P, n+1, fpn);

	fq_poly_set_coeff(P, 0, one, fpn);
	fq_poly_set_coeff(P, 1, one, fpn);
	fq_poly_set_coeff(P, 4, one, fpn);

	fq_t res;
	fq_init(res, fpn);

	find_normal_random_test(res, P, Y, fpn, fp);
	fq_print_pretty(res, fpn);
	flint_printf("\n");

	fmpz_clear(p);
	fq_ctx_clear(fp);
	fq_ctx_clear(fpn);
	fq_poly_clear(P, fpn);
	fq_clear(res, fpn);
	fq_clear(one, fpn);
	fq_clear(Y, fpn);
}
