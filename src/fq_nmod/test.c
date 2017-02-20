#include "main.h"

void main() {
	fmpz_t p;
	fmpz_init(p);

	fmpz_set_ui(p, 2);
	slong d = 31;

	fq_nmod_ctx_t field;
	fq_nmod_ctx_init(field, p, d, "X");
	fq_nmod_ctx_print(field);
    fq_nmod_ctx_clear(field);

	fq_nmod_ctx_init(field, p, d, "X");
	fq_nmod_ctx_print(field);

	fq_nmod_poly_t P;
	fq_nmod_poly_init(P, field);

	fq_nmod_t res, X;
	fq_nmod_init(res, field);
	fq_nmod_init(X, field);
    fq_nmod_gen(X, field);

    ceil(1.1*d);

	fq_nmod_clear(res, field);
	fq_nmod_clear(X, field);
	fmpz_clear(p);
	fq_nmod_poly_clear(P, field);
	fq_nmod_ctx_clear(field);
}
