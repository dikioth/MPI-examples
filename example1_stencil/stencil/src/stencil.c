#include <mpi.h>
#include <stdio.h>
#include "stencil.h"

int main(int argc, char **argv)
{
	MPI_Comm Comm1D;
	int num_values, split_size;
	double *input, *output, *split_buf;
	int rank, left_rank, source_rank, right_rank, num_proc;
	int TAGLEFT = 1, TAGRIGHT = 2;

	/* MPI init */
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	/* Create 1D cartesian */
	int ndims = 1;
	int dims[1] = {num_proc};
	int periods[2] = {1, 1};
	int reorder = 0;
	MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, reorder, &Comm1D);

	/* Check arguments */
	if (!(argc == 4 || argc == 3) && rank == 0)
	{
		printf("Usage : ./stencil ifile [ofile] nsteps\n");
		printf("ifile : Input file\n");
		printf("ofile : Output file (optional)\n");
		printf("nsteps: Number of stencil applications\n");
		return MPI_Abort(MPI_COMM_WORLD, 1);
	}

	/* Local store args */
	char *input_name = argv[1];
	char *output_name = argc == 4 ? argv[2] : NULL;
	int num_steps = argc == 4 ? atoi(argv[3]) : atoi(argv[2]);

	/* Read and distribute data */
	if (rank == ROOT)
	{
		/* Read input file and Check divisibility */
		if (0 > (num_values = read_input(input_name, &input)) || num_values % num_proc != 0)
			return MPI_Abort(MPI_COMM_WORLD, 1);

		// Allocate data for result
		if (NULL == (output = (double *)malloc(num_values * sizeof(double))))
		{
			perror("Couldn't allocate memory for output");
			return MPI_Abort(MPI_COMM_WORLD, 1);
		}
	}

	/* Send num_values to all processes */
	MPI_Bcast(&num_values, 1, MPI_INT, ROOT, Comm1D);
	split_size = num_values / num_proc;

	/* Stencil constants */
	double h = 2.0 * PI / num_values;
	const int STENCIL_WIDTH = 5;
	const int EXTENT = STENCIL_WIDTH / 2;
	const double STENCIL[] = {1.0 / (12 * h), -8.0 / (12 * h), 0.0, 8.0 / (12 * h), -1.0 / (12 * h)};

	/* Allocate mem for split buf */
	split_buf = (double *)malloc((split_size + 2 * EXTENT) * sizeof(double));

	/* Split input array into sub arrays */
	MPI_Scatter(input, split_size, MPI_DOUBLE, split_buf + EXTENT, split_size, MPI_DOUBLE, ROOT, Comm1D);

	/* Get rank from neighbours */
	MPI_Cart_shift(Comm1D, 0, -1, &source_rank, &left_rank);
	MPI_Cart_shift(Comm1D, 0, 1, &source_rank, &right_rank);

	MPI_Request sendreq[2], recvreq[2];

	double *tmpout = (double *)malloc((split_size + 2 * EXTENT) * sizeof(double));
	double leftbuf[2], rightbuf[2];

	// Stencil operation
	double start = MPI_Wtime();
	for (int s = 0; s < num_steps; s++)
	{
		/* Receive extensions from right & left neighbours */
		MPI_Isend(split_buf + split_size, EXTENT, MPI_DOUBLE, right_rank, TAGRIGHT, Comm1D, &sendreq[1]); // Send to right
		MPI_Irecv(split_buf, EXTENT, MPI_DOUBLE, left_rank, TAGRIGHT, Comm1D, &recvreq[1]);				  // Receive left

		/* Send extensions to left & right neighbours */
		MPI_Isend(split_buf + EXTENT, EXTENT, MPI_DOUBLE, left_rank, TAGLEFT, Comm1D, &sendreq[0]);				  // Send to left
		MPI_Irecv(split_buf + EXTENT + split_size, EXTENT, MPI_DOUBLE, right_rank, TAGLEFT, Comm1D, &recvreq[0]); // Receive right

		/* Center stencil operation (After sending edge elements) */
		for (int i = 2 * EXTENT; i < split_size; i++)
		{
			double result = 0;
			for (int j = 0; j < STENCIL_WIDTH; j++)
			{
				result += split_buf[i + j - EXTENT] * STENCIL[j];
			}
			tmpout[i] = result;
		}

		/* Left and right stencil operation (After receiving edge elements) */
		MPI_Waitall(2, recvreq, MPI_STATUS_IGNORE);
		for (int i = EXTENT; i < 2 * EXTENT; i++)
		{
			double result = 0;
			for (int j = 0; j < STENCIL_WIDTH; j++)
			{
				result += split_buf[i + j - EXTENT] * STENCIL[j];
			}

			tmpout[i] = result;
		}

		for (int i = split_size; i < split_size + EXTENT; i++)
		{
			double result = 0;
			for (int j = 0; j < STENCIL_WIDTH; j++)
			{
				result += split_buf[i + j - EXTENT] * STENCIL[j];
			}
			tmpout[i] = result;
		}

		MPI_Waitall(2, sendreq, MPI_STATUS_IGNORE);
		double *tmp = split_buf;
		split_buf = tmpout;
		tmpout = tmp;
	}
	double execution_time = MPI_Wtime() - start;
	double max_time;

	/* Send back to root  */
	MPI_Gather(split_buf + EXTENT, split_size, MPI_DOUBLE, output, split_size, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);

	/* Handle execution time */
	MPI_Reduce(&execution_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, ROOT, MPI_COMM_WORLD);

	if (rank == ROOT)
	{
		if (output_name)
		{
			if (0 != write_output(output_name, output, num_values))
			{
				return 2;
			}
		}
		printf("%f\n", execution_time);
		free(output);
	}
	MPI_Comm_free(&Comm1D);
	MPI_Finalize();
	return MPI_SUCCESS;
}

int read_input(const char *file_name, double **values)
{
	FILE *file;
	if (NULL == (file = fopen(file_name, "r")))
	{
		perror("Couldn't open input file");
		return -1;
	}
	int num_values;
	if (EOF == fscanf(file, "%d", &num_values))
	{
		perror("Couldn't read element count from input file");
		return -1;
	}
	if (NULL == (*values = malloc(num_values * sizeof(double))))
	{
		perror("Couldn't allocate memory for input");
		return -1;
	}
	for (int i = 0; i < num_values; i++)
	{
		if (EOF == fscanf(file, "%lf", &((*values)[i])))
		{
			perror("Couldn't read elements from input file");
			return -1;
		}
	}
	if (0 != fclose(file))
	{
		perror("Warning: couldn't close input file");
	}
	return num_values;
}

int write_output(char *file_name, const double *output, int num_values)
{
	FILE *file;
	if (NULL == (file = fopen(file_name, "w")))
	{
		perror("Couldn't open output file");
		return -1;
	}
	for (int i = 0; i < num_values; i++)
	{
		if (0 > fprintf(file, "%.4f ", output[i]))
		{
			perror("Couldn't write to output file");
		}
	}
	if (0 > fprintf(file, "\n"))
	{
		perror("Couldn't write to output file");
	}
	if (0 != fclose(file))
	{
		perror("Warning: couldn't close output file");
	}
	return EXIT_SUCCESS;
}
