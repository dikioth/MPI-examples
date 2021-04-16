/*********************************************************************************
*                              Author: Elvis Rodas                               *
*                File Name: matmul.c                                             *
*                     Creation Date: April 1, 2021 01:16 PM                      *
*                      Last Updated: April 2, 2021 06:02 PM                      *
*                               Source Language: c                               *
*               Repository: https://github.com/dikioth/PDP-A2.git                *
*                                                                                *
*                            --- Code Description ---                            *
*               MPI implementation of matrix-matrix multiplication of            *
*               size n x n. The number of processes p is assumed to              *
*               obey:                                                            * 
*                               n % sqrt(p) == 0                                 *
*********************************************************************************/

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matmul.h"

int main(int argc, char **argv)
{
	int rank, num_procs;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

	/* Matrix mult params */
	double *A, *B, *C, *mA, *mB, *mC;
	int split_size, mat_size;

	MPI_Datatype row_t, col_t, block_t;

	/* Time measurment vars*/
	double tstart, texec;

	/* Check arguments */
	if (!(argc == (NUM_ARGS + 1) || argc == (NUM_ARGS)) && rank == ROOT)
	{
		printf("Usage : ./matmul ifile ofile\n");
		printf("ifile : Input file\n");
		printf("ofile : Output file (optional)\n");
		return MPI_Abort(MPI_COMM_WORLD, 1);
	}

	char *ifile = argv[1];
	char *ofile = argc == (NUM_ARGS + 1) ? argv[2] : NULL;

	if (rank == ROOT)
	{
		mat_size = readfile(ifile, &A, &B);

		/* Check divisibility */
		if (mat_size % num_procs != 0)
		{
			perror("Matrix size must be divisible by num procs\n");
			return MPI_Abort(MPI_COMM_WORLD, 1);
		}

		/* Allocate mem for results */
		if (NULL == (C = (double *)malloc(mat_size * mat_size * sizeof(double))))
			perror("Error allocating memory for results\n");
	}

	/* Start timing */
	if (rank == ROOT)
		tstart = MPI_Wtime();

	/* Broadcast mat_size to all PEs */
	MPI_Bcast(&mat_size, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
	/* block size */
	split_size = mat_size / sqrt(num_procs);
	int block_size = mat_size * split_size;

	/* Allocate mem for matrix blocks */
	mA = (double *)malloc(block_size * sizeof(double));
	mB = (double *)malloc(block_size * sizeof(double));
	mC = (double *)calloc(split_size * split_size, sizeof(double));

	/* Define derived data types */
	MPI_Type_contiguous(block_size, MPI_DOUBLE, &row_t);

	MPI_Type_vector(mat_size, split_size, mat_size, MPI_DOUBLE, &col_t);
	MPI_Type_create_resized(col_t, 0, split_size * sizeof(double), &col_t);

	MPI_Type_vector(split_size, split_size, mat_size, MPI_DOUBLE, &block_t);
	MPI_Type_create_resized(block_t, 0, 2 * sizeof(double), &block_t);

	MPI_Type_commit(&col_t);
	MPI_Type_commit(&row_t);
	MPI_Type_commit(&block_t);

	/* Scatter submatrices to all PEs */
	int *scounts = (int *)malloc(num_procs * sizeof(int));
	int *displsA = (int *)malloc(num_procs * sizeof(int));
	int *displsB = (int *)malloc(num_procs * sizeof(int));
	for (int i = 0; i < num_procs; i++)
	{
		scounts[i] = 1;
		displsA[i] = i / split_size;
		displsB[i] = i % split_size;
	}

	MPI_Scatterv(A, scounts, displsA, row_t, mA, 1, row_t, ROOT, MPI_COMM_WORLD);
	MPI_Scatterv(B, scounts, displsB, col_t, mB, block_size, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);

	free(scounts);
	free(displsA);
	free(displsB);

	/*Matrix multiplication */
	for (int i = 0; i < split_size; i++)
		for (int k = 0; k < mat_size; k++)
			for (int j = 0; j < split_size; j++)
				mC[i * split_size + j] += mA[i * mat_size + k] * mB[split_size * k + j];

	/* Send back results to root */
	int *recvCounts = (int *)malloc(num_procs * sizeof(int));
	int *recvdispls = (int *)malloc(num_procs * sizeof(int));
	for (int i = 0; i < num_procs; i++)
	{
		recvCounts[i] = 1;
		recvdispls[i] = i + split_size * (i / split_size);
	}

	MPI_Gatherv(mC, split_size * split_size, MPI_DOUBLE, C, recvCounts, recvdispls, block_t, ROOT, MPI_COMM_WORLD);
	free(recvCounts);
	free(recvdispls);
	free(mA);
	// free(mB);
	free(mC);

	/* Write output file */
	if (rank == ROOT)
	{
		texec = MPI_Wtime() - tstart;
		printf("%f\n", texec);

		if (ofile)
			writefile(ofile, C, mat_size, mat_size);

		free(A);
		free(B);
		free(C);
	}

	MPI_Finalize();
}

int readfile(char *fname, double **x, double **y)
{
	FILE *file;
	int num_values;

	if (NULL == (file = fopen(fname, "r")) && (num_values = FAILURE))
		perror("Couldn't open input file");

	if (EOF == fscanf(file, "%d", &num_values) && (num_values = FAILURE))
		perror("Couldn't read element count from input file");

	if (NULL == (*x = malloc(num_values * num_values * sizeof(double))) && (num_values = FAILURE))
		perror("Couldn't allocate memory for input");

	if (NULL == (*y = malloc(num_values * num_values * sizeof(double))) && (num_values = FAILURE))
		perror("Couldn't allocate memory for input");

	for (int i = 0; i < num_values * num_values; i++)
	{
		if (EOF == fscanf(file, "%lf", &((*x)[i])) && (num_values = FAILURE))
			perror("Couldn't read elements from input file");
	}

	for (int i = 0; i < num_values * num_values; i++)
	{
		if (EOF == fscanf(file, "%lf", &((*y)[i])) && (num_values = FAILURE))
			perror("Couldn't read elements from input file");
	}

	if (EOF == fclose(file) && (num_values = FAILURE))
		perror("Warning: couldn't close input file");

	return num_values;
}

int writefile(char *fname, double *output, int rows, int cols)
{
	FILE *file;
	if (NULL == (file = fopen(fname, "w")))
	{
		perror("Couldn't open output file");
		return EXIT_FAILURE;
	}

	for (int i = 0; i < rows * cols; i++)
	{
		if (0 > fprintf(file, (i + 1) % cols == 0 ? "%.6f\n" : "%.6f ", output[i]))
			perror("Couldn't write to output file");
	}

	return EXIT_SUCCESS;
}
