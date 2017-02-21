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
	fprintf(stderr, "Usage : %s [-p characteristic] [-d degreeOfExtension] [-a (naive|random|Luneburg|Lenstra)]\n", prog);
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
				fmpz_set_str(p, optarg, 10);
				break;
			case 'd':
				d = atol(optarg);
				break;
			case 'a':
				alg = optarg[3];
				break;
			default:
				help(argv[0]);
				return 1;
		}
	}

	// We check that characteristic and degree were given
	if (fmpz_cmp_ui(p,0) && (d != 0)) {

		// Creation of the field F_{p^d}
		fq_nmod_ctx_t field;
		fq_nmod_ctx_init(field, p, d, "X");

		// We print the situation
		printf("\t\t *****\nCurrently working with the field F_{p^d}\nrepresented as F_p[X]/(f(X))\nwhere :\n\n");
		fq_nmod_ctx_print(field);
		printf("\t\t *****\n");

		fq_nmod_t res;
		fq_nmod_init(res, field);
		fq_nmod_one(res, field);

		switch (alg) {
			case 'd':
			case 'D':
				normal_random(res, field);
				fq_nmod_print_pretty(res, field);
				flint_printf("\n");
				break;
			case 'e':
			case 'E':
				luneburg(res, field);
				fq_nmod_print_pretty(res, field);
				flint_printf("\n");
				break;
			case 's':
			case 'S':
				lenstra(res, field);
				fq_nmod_print_pretty(res, field);
				flint_printf("\n");
				break;
			case 'v':
			case 'V':
				naive(res, field);
				fq_nmod_print_pretty(res, field);
				flint_printf("\n");
				break;
			default:
				printf("Please specify an algorithm among random, Luneburg, Lenstra.\n");
		}

		fq_nmod_clear(res, field);
		fq_nmod_ctx_clear(field);
	}
	else {
		help(argv[0]);
	}

	fmpz_clear(p);
}
