#include "main.h"

/*   NAIVE
 *
 * Naive random algorithm to compute normal elements.
 *
 * We juste pick an element at random until it is a 
 * normal element.
 */

void naive(fq_nmod_t res, const fq_nmod_ctx_t field) {
	
	// We create a random state
	flint_rand_t state;
	flint_randinit(state);

	// We create and initialize a temporary variable
	fq_nmod_t tmp;
    fq_nmod_init(tmp, field);

    // We set it to a random element
    fq_nmod_randtest_not_zero(tmp, state, field);

    while (!is_normal(tmp, field)) {
        fq_nmod_randtest_not_zero(tmp, state, field);
    }

    // We set res to the normal element tmp
    fq_nmod_set(res, tmp, field);
    
    // We clear the variables
    fq_nmod_clear(tmp, field);
    flint_randclear(state);
}
