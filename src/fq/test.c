#include "main.h"

void main() {
	fmpz_t p;
	fmpz_init(p);

	fmpz_set_ui(p, 2);
	slong d = 31;

	fq_ctx_t field;
	fq_ctx_init(field, p, d, "X");
	fq_ctx_print(field);
    fq_ctx_clear(field);

	fq_ctx_init(field, p, d, "X");
	fq_ctx_print(field);

	fq_poly_t P;
	fq_poly_init(P, field);

	fq_t res, X;
	fq_init(res, field);
	fq_init(X, field);
    fq_gen(X, field);

    ceil(1.1*d);

	fq_clear(res, field);
	fq_clear(X, field);
	fmpz_clear(p);
	fq_poly_clear(P, field);
	fq_ctx_clear(field);
}
