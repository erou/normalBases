/*
 * 	MAIN
 *
 * Code permitting to choose an algorithm and
 * to specify the arguments such as the characteristic, 
 * the degree of the extension.
 *
 */

#include "main.h"

// Help function
void help(const char* prog) {
	fprintf(stderr, "Usage : %s [-p characteristic] [-d degreeOfExtension] [-a (random|Luneburg|Lenstra)]\n", prog);
}

int main(int argc, char **argv) {

	// Initialisation of the variables
	fmpz_t p;
	fmpz_init(p);
	slong d = 0;

	/* Use of the arguments :
	 *
	 * 	-p : the characteristic of the field
	 * 	-s : the degree of the extension
	 * 	-a : the algorithm wanted
	 */
	
	int opt;
	char alg = 'd';

	while ((opt = getopt(argc, argv, "hp:d:a:")) != -1) {
		switch (opt) {
			case 'p':
				fmpz_set_ui(p, atol(optarg));
				break;
			case 'd':
				d = atol(optarg);
				break;
			case 'a':
				alg = optarg[1];
				break;
			default:
				help(argv[0]);
				return 1;
		}
	}

	// We check that characteristic and degree were given
	if (fmpz_cmp_ui(p,0) && (d != 0)) {

		// Creation of the field F_{p^d}
		fq_ctx_t field;
		fq_ctx_init(field, p, d, "X");

		// We print the situation
		printf("\t\t *****\nCurrently working with the field F_{p^d}\nrepresented as F_p[X]/(f(X))\nwhere :\n\n");
		fq_ctx_print(field);
		printf("\t\t *****\n");

		fq_t res;
		fq_init(res, field);
		fq_one(res, field);

		// We compute P = X^d - 1
		fq_poly_t P, ord;
		fq_poly_init2(P, d+1, field);
		fq_poly_set_coeff(P, d, res, field);
		fq_neg(res, res, field);
		fq_poly_set_coeff(P, 0, res, field);
		fq_poly_init(ord, field);

		fq_zero(res, field);	

		switch (alg) {
			case 'a':
			case 'A':
				normal_random(res, field);
				sigma_order(ord, res, field);
				if (is_normal(res, field)) {
					flint_printf("is_normal : OK\n");
				}
				
				fq_poly_print_pretty(ord, "Y", field);
				flint_printf("\n");
				fq_print_pretty(res, field);
				flint_printf("\n");
				break;
			case 'u':
			case 'U':
				luneburg(res, field);
				sigma_order(ord, res, field);
				if (is_normal(res, field)) {
					flint_printf("is_normal : OK\n");
				}
				
				fq_poly_print_pretty(ord, "Y", field);
				flint_printf("\n");
				fq_print_pretty(res, field);
				flint_printf("\n");
				break;
			case 'e':
			case 'E':
				lenstra(res, field);
				sigma_order(ord, res, field);
				if (is_normal(res, field)) {
					flint_printf("is_normal : OK\n");
				}
				
				fq_poly_print_pretty(ord, "Y", field);
				flint_printf("\n");
				fq_print_pretty(res, field);
				flint_printf("\n");
				break;
			default:
				printf("Please specify an algorithm among random, Luneburg, Lenstra.\n");
		}

		fq_clear(res, field);
		fq_poly_clear(P, field);
		fq_poly_clear(ord, field);
		fq_ctx_clear(field);
	}
	else {
		help(argv[0]);
	}

	fmpz_clear(p);
}
