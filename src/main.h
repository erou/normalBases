#include <stdlib.h>
#include <stdio.h>
#include <fmpz.h>
#include <fq.h>
#include <fq_poly.h>

int is_normal(fq_t, fq_ctx_t, fq_ctx_t);
void find_normal_random(fq_t, fq_poly_t, fq_t, fq_ctx_t, fq_ctx_t);
