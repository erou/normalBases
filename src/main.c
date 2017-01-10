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
void help(char* prog) {
	fprintf(stderr, "Usage : %s [-p characteristic] [-d degreeOfExtension] [-a (random|Luneburg|Lenstra)]\n", prog);
}

int main(int argc, char **argv) {

	// Initialisation of the variables
	fmpz_t p;
	fmpz_init(p);
	slong d;

	/* Use of the arguments :
	 *
	 * 	-p : the characteristic of the field
	 * 	-s : the degree of the extension
	 * 	-a : the algorithm wanted
	 */
	
	int opt;
	char alg = 'r';

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

	if (optind > 0 && optind < argc) {
		// Creation of the field F_{p^d}
		fq_ctx_t field;
		fq_ctx_init(field, p, d, "X");

		printf("\t\t *****\nCurrently working with the field F_{p^d}\nrepresented as F_p[X]/(f(X))\nwhere :\n\n");
		fq_ctx_print(field);
		printf("\t\t *****\n");

		fmpz_clear(p);
		fq_ctx_clear(field);
	}
	else {
		help(argv[0]);
	}
}
