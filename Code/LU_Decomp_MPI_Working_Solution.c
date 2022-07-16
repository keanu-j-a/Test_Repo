#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>

void LUDecomposition(int rank, int size, int r, long double** arrayInFunction) {

	/* Define inital varaibles and arrays */
	int c = r, j = 0, i = 0;
	long double sum = 0, x = 0, y = 0, z = 0;

	/* Define the L and U arrays */
	long double** lowerTArray = (long double**)calloc(r, sizeof(long double*));
	long double** upperTArray = (long double**)calloc(r, sizeof(long double*));
	long double* rowk = (long double*)calloc(r, sizeof(long double*));

	/* Define the range of values in which processes will operate */

	long int end_val = abs(rank - (r - 1));
	/*

	end_val is a ranged ranking algorithim that distributes the matrix work according
	to the id of the reviving process. The range is as follows:

	size = 3
	P_0 --> end_val = 2
	P_1 --> end_val = 1
	P_2 --> end_val = 0

	*/

	////////////////////////////////////////////////////////////////////////////////////////////////////////////

				/* Further define the L and U matrix values */
				for (int n = 0; n < c; n++) {
					lowerTArray[n] = (long double*)calloc(r, sizeof(long double));
				}
				for (int n = 0; n < c; n++) {
					upperTArray[n] = (long double*)calloc(r, sizeof(long double));
				}

				/* Implement a 2D - coarse grain LU decomposition solving algorithim */
				for (int k = rank; k <= r - end_val; k++) {
					MPI_Bcast(&arrayInFunction[k][j], 1, MPI_LONG_DOUBLE, rank, MPI_COMM_WORLD);
					if (k <= r) {
						for (int i = k + 1; i < r; i++) {
							lowerTArray[i][k] = arrayInFunction[i][k] / arrayInFunction[k][k];
							printf("This is lower %lf \n", lowerTArray[i][k]);
						}
					}
					MPI_Bcast(&lowerTArray[i][k], 1, MPI_LONG_DOUBLE, rank, MPI_COMM_WORLD);
					for (int j = k + 1; j < r; j++) {
						for (int i = k + 1; i < r; i++) {
							arrayInFunction[i][j] = arrayInFunction[i][j] - lowerTArray[i][k] * arrayInFunction[k][j];
							printf("This is upper %lf \n", arrayInFunction[i][j]);
						}
					}
				}

				/* Assign 0s to the U matrix */
				for (int m = 0; m < r; m++) {

					/* n = row variable */
					for (int n = 1 + m; n < r; n++) {
						arrayInFunction[n][m] = 0;
					}
				}

				/* Assign 1s to the L matrix diag */
				for (int k = 0; k < r; k++) {
					lowerTArray[k][k] = 1;
				}

	////////////////////////////////////////////////////////////////////////////////////////////////////////////

				/* Print the contents of the L matrix */
				printf(" This is the L Matrix from rank %d \n ", rank);
				for (int i = 0; i < r; i++) {
					printf("[ ");
					for (int j = 0; j < r; j++) {
						printf(", %lf", lowerTArray[i][j]);
					}
					printf(" ]\n");
				}


				printf("\n");

				/* Print the contents of the U matrix */

				printf(" This is the U Matrix from rank %d \n ", rank);
				for (int i = 0; i < r; i++) {
					printf("[ ");
					for (int j = 0; j < r; j++) {
						printf(", %lf", arrayInFunction[i][j]);
					}
					printf(" ] \n");
				}

	////////////////////////////////////////////////////////////////////////////////////////////////////////////

}

int main(int argc, char* argv[]) {

	/* Initalize the MPI variables to be used in this code */
	MPI_Init(&argc, &argv);
	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	/* Define inital variables */
	int r = 3, c = r, u = 0;
	int* numberArray = (int*)calloc(r * c, sizeof(int));
	long double** arrayInFunction = (long double**)calloc(r, sizeof(int*));

				/* Produce randomly generated numbers */
				srand(100);
				for (int i = 0; i < r * c; i++) {
					numberArray[i] = rand() / 100;
				}

				/* Further define 2D memory heap */
				for (int i = 0; i < r; i++) {
					arrayInFunction[i] = (long double*)calloc(c, sizeof(long double*));
				}

				/* Assign random numbers to memory heap */
				for (int i = 0; i < r; i++) {
					for (int j = 0; j < c; j++) {
						arrayInFunction[i][j] = numberArray[u];
						u++;
					}
				}

				/* Print the contents of the origional matrix for validation */
				if (rank == 0) {
					printf("This is the origional matrix: \n");
					for (int i = 0; i < r; i++) {
						printf("[ ");
						for (int j = 0; j < r; j
							++) {
							printf("%lf ", arrayInFunction[i][j]);
						}
						printf(" ] \n");
					}
					printf("\n");
				}

	/* Call the function */
	LUDecomposition(rank, size, r, arrayInFunction);

	/* Finalize the program */
	MPI_Finalize();
	return 0;
}