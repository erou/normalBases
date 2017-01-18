#include "main.h"

void main() {
	fmpz_t p;
	fmpz_init(p);

	fmpz_set_ui(p,2);
	fq_ctx_t field;

	slong d = 2;
	
	fq_ctx_init(field, p, d, "X");

	fq_ctx_print(field);

	fq_t res, x;
	fq_init(res, field);

	lenstra(res, field);
	fq_print_pretty(res, field);
	flint_printf("\n");

	fq_clear(res, field);
	fmpz_clear(p);
	fq_ctx_clear(field);
}
