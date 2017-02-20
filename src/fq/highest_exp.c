#include "main.h"

/* 
 *      HIGHEST EXP
 *
 * Compute the highest possible exposant of g in f.
 *
 *   Note : this function was copied and adaptated from the
 *   flint function fq_poly_remove, which was not convenient
 *   for our use.
 */

slong highest_exp(const fq_poly_t f, const fq_poly_t g, const fq_ctx_t field) {

    fq_poly_t q, r, fcopy;
    slong i = 0;

    fq_poly_init(q, field);
    fq_poly_init(r, field);
    fq_poly_init(fcopy, field);
    fq_poly_set(fcopy, f, field);

    while (1)
    {
        if (f->length < g->length)
            break;
        fq_poly_divrem(q, r, fcopy, g, field);
        if (r->length == 0)
            fq_poly_swap(q, fcopy, field);
        else
            break;
        i++;
    }

    fq_poly_clear(q, field);
    fq_poly_clear(r, field);
    fq_poly_clear(fcopy, field);

    return i;
}
