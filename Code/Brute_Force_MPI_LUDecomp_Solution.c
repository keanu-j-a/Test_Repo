#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>

/*

Written By Keanu J. Ammons
Modified on: 07/8/2022
Idaho National Labratory (INL) Computational Nuclear Engineering

This code is a best attempt at implementing MPI into a n x n matrix LU decomposition solver. The solution for this algorithim
was sucessfully verified against a 3 x 3 matrix, although further testing is needed for much larger matrix values. The ultimate goal
is to apply this code against a 1000 x 1000 size matrix, validate the output, and observe speedup % after designating N processes
on the sawtooth supercomputing cluster. The code shown here will be compared to the determinate solving algorithim (found in the Github
repository).

*/

void LUDecomposition(int rank, int size, int r, long double** arrayInFunction) {

	/*
	

	// Github recognizes this change
	// and a few more! =)
	// why not added another line?

	The arguments of the LU decomposition function are as follows

	1.) 'rank'.......................The definition of each processes's individual designation (0 to size - 1)
	2.) 'size'.......................The overall number of processes in the system
	3.) 'r'..........................The row (and column) of the n x n square matrix
	4.) 'arrayInFunction'............The n x n (2D) heap memory allocation that will hold the numbers for the A matrix.

	*/

	/* Define inital varaibles and arrays */
	int c = r, j = 0, i = 0;
	long double sum = 0, x = 0, y = 0, z = 0;

	/* Define the L and U arrays */
	long double** lowerTArray = (long double**)calloc(r, sizeof(long double*));
	long double** upperTArray = (long double**)calloc(r, sizeof(long double*));
	long double* rowk = (long double*)calloc(r, sizeof(long double*));

	/* MPI Request */
	MPI_Request request;

	// Define the range of values in which processes will operate
	int upperLimit = r;
	int start_val = rank * ceil(upperLimit / size), end_val;

	// generate an algorithim that defines the range of
	// each process to handle for the fibb_sequence problem.

	if (rank == (size - 1)) {
		end_val = upperLimit - 1;
	}
	else {
		end_val = start_val + ceil(upperLimit / size);
	}


	// long int end_val = abs(rank - (r - 1));

	// end_val is a ranged ranking algorithim that distributes the matrix work according
	// to the id of the reviving process. The range is as follows:

	// size = 3
	// P_0 --> end_val = 2
	// P_1 --> end_val = 1
	// P_2 --> end_val = 0


	////////////////////////////////////////////////////////////////////////////////////////////////////////////

				/* Further define the L and U matrix values */
				for (int n = 0; n < c; n++) {
					lowerTArray[n] = (long double*)calloc(r, sizeof(long double));
				}
				for (int n = 0; n < c; n++) {
					upperTArray[n] = (long double*)calloc(r, sizeof(long double));
				}

				/* Implement a 2D - coarse grain LU decomposition solving algorithim */
				for (int k = 0; k < r; k++) {
					MPI_Bcast(&arrayInFunction[k][j], r, MPI_LONG_DOUBLE, rank, MPI_COMM_WORLD);
					if (k < r) {
						for (int i = k + 1; i < r; i++) {
							for (int p = start_val; p < end_val; p++) {
								lowerTArray[i][k] = arrayInFunction[i][k] / arrayInFunction[k][k];
								printf("This is lower %lf from rank %d \n", lowerTArray[i][k], rank);
							}
						}
					}
					MPI_Bcast(&lowerTArray[i][k], r, MPI_LONG_DOUBLE, rank, MPI_COMM_WORLD);
					for (int j = k + 1; j < r; j++) {
						for (int i = k + 1; i < r; i++) {
							for (int p = start_val; p < end_val; p++) {
								arrayInFunction[i][j] = arrayInFunction[i][j] - lowerTArray[i][k] * arrayInFunction[k][j];
								printf("This is upper %lf from rank %d \n", arrayInFunction[i][j], rank);
							}
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
	long double** arrayInFunction = (long double**)calloc(r, sizeof(int));

				/* Produce randomly generated numbers */
				srand(100);
				for (int i = 0; i < r * c; i++) {
					numberArray[i] = rand() / 100;
				}

				/* Further define 2D memory heap */
				for (int i = 0; i < r; i++) {
					arrayInFunction[i] = (long double*)calloc(c, sizeof(long double)); // other issue was here
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
						for (int j = 0; j < r; j++) {
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