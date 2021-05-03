#include "conjugate.h"

void createMat(int rows, int cols, mat_t *mat)
{
    (*mat).data = (double *)malloc(rows * cols * sizeof(double));
    (*mat).rows = rows;
    (*mat).cols = cols;
}

void stencil(mat_t in, mat_t out, MPI_Comm comm)
{
    reqn = 0;

    /* Left neighbour comm */
    if (left != MPI_PROC_NULL)
    {
        MPI_Irecv(leftbuff.data, in.rows, MPI_DOUBLE, left, 0, comm, &reqs[reqn++]);
        MPI_Send(in.data, 1, col_t, left, 0, comm);
    }

    /* Right neighbour comm */
    if (right != MPI_PROC_NULL)
    {
        MPI_Irecv(righbuff.data, in.rows, MPI_DOUBLE, right, 0, comm, &reqs[reqn++]);
        MPI_Send(in.data + (in.cols - 1), 1, col_t, right, 0, comm);
    }

    /* Top neighbour comm */
    if (top != MPI_PROC_NULL)
    {
        MPI_Irecv(topbuff.data, in.cols, MPI_DOUBLE, top, 0, comm, &reqs[reqn++]);
        MPI_Send(in.data, 1, row_t, top, 0, comm);
    }

    /* Down neighbour comm */
    if (down != MPI_PROC_NULL)
    {
        MPI_Irecv(downbuff.data, in.cols, MPI_DOUBLE, down, 0, comm, &reqs[reqn++]);
        MPI_Send(in.data + (in.rows - 1) * in.cols, 1, row_t, down, 0, comm);
    }

    MPI_Waitall(reqn, reqs, MPI_STATUSES_IGNORE);

    /* Perform stencil operations */
    for (int i = 0; i < in.rows; i++)
        for (int j = 0; j < in.cols; j++)
        {
            /* Center points */
            out.data[i * in.cols + j] = STENCIL[2] * in.data[i * in.cols + j];

            /* Top stencil */
            out.data[i * in.cols + j] += (i == 0)
                                             ? ((top != MPI_PROC_NULL) ? STENCIL[1] * topbuff.data[j] : 0)
                                             : STENCIL[1] * in.data[(i - 1) * in.cols + j];

            /* Down stencil */
            out.data[i * in.cols + j] += (i == in.rows - 1)
                                             ? ((down != MPI_PROC_NULL) ? STENCIL[3] * downbuff.data[j] : 0)
                                             : STENCIL[3] * in.data[(i + 1) * in.cols + j];

            /* left stencil */
            out.data[i * in.cols + j] += (j == 0)
                                             ? ((left != MPI_PROC_NULL) ? STENCIL[0] * leftbuff.data[i] : 0)
                                             : STENCIL[0] * in.data[i * in.cols + j - 1];

            /* right stencil */
            out.data[i * in.cols + j] += (j == in.cols - 1)
                                             ? ((right != MPI_PROC_NULL) ? STENCIL[4] * righbuff.data[i] : 0)
                                             : STENCIL[4] * in.data[i * in.cols + j + 1];
        }
}

void scalarprod(double scalar, mat_t in, mat_t out)
{
    if (in.cols != out.cols || in.rows != out.rows)
        perror("scalarprod: Matrices must be same size\n");

    for (int i = 0; i < in.rows; i++)
        for (int j = 0; j < in.cols; j++)
            out.data[i * in.cols + j] = scalar * in.data[i * in.cols + j];
}

void matsum(mat_t v1, mat_t v2, mat_t out)
{
    if (v1.rows != v2.rows || v1.cols != v2.cols)
        perror("matsum: Matrices must be same size\n");

    for (int i = 0; i < v1.rows; i++)
        for (int j = 0; j < v1.cols; j++)
            out.data[i * out.cols + j] = v1.data[i * v1.cols + j] + v2.data[i * v2.cols + j];
}

void matsub(mat_t v1, mat_t v2, mat_t out)
{
    if (v1.rows != v2.rows || v1.cols != v2.cols)
        perror("matsub: Matrices must be same size\n");

    for (int i = 0; i < v1.rows; i++)
        for (int j = 0; j < v1.cols; j++)
            out.data[i * v1.cols + j] = v1.data[i * v1.cols + j] - v2.data[i * v1.cols + j];
}

void dotprod(mat_t v1, mat_t v2, double *out)
{
    if (v1.cols != v2.cols || v1.rows != v2.rows)
        perror("Dotprod: Matrices must be same size\n");

    *out = 0.0;
    for (int i = 0; i < v1.rows; i++)
        for (int j = 0; j < v1.cols; j++)
            *out += v1.data[i * v1.cols + j] * v2.data[i * v2.cols + j];
}

void vectorb(double x, double y, double h, double *out)
{
    *out = 2 * h * h * (x * (1.0 - x) + y * (1.0 - y));
}

int main(int argc, char *argv[])
{
    int rank, numprocs;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

    /* Handle args */
    if (argc != NARGS + 1 && rank == ROOT)
    {
        printf("usage: ./conjugate n\n");
        printf("n: The number of domain splits\n");
    }

    /* local vars */
    sqrtp = sqrt(numprocs);
    size = (n = atoi(argv[1])) - 1;
    h = 1.0 / n;

    /* Create 2D topological world */
    dims[0] = sqrtp;
    dims[1] = sqrtp;
    MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, reorder, &Comm2D);

    /* Compute local matrices size and offsets */
    splitsize = round(size / (float)sqrtp);
    mcols = ((rank + 1) % sqrtp == 0)
                ? size - splitsize * (sqrtp - 1)
                : splitsize;
    mrows = (rank + sqrtp >= numprocs)
                ? size - splitsize * (sqrtp - 1)
                : splitsize;

    xoffset = 1 + splitsize * (rank - rank % sqrtp) / sqrtp;
    yoffset = 1 + splitsize * (rank % sqrtp);

    /* Allocate memory for local matrices */
    createMat(mrows, mcols, &mu);
    createMat(mrows, mcols, &mb);
    createMat(mrows, mcols, &mg);
    createMat(mrows, mcols, &md);
    createMat(mrows, mcols, &mq);
    createMat(mrows, mcols, &temp1);
    createMat(mrows, mcols, &temp2);

    /* Get neighbours */
    MPI_Cart_shift(Comm2D, XAXIS, -1, &source, &top);
    MPI_Cart_shift(Comm2D, XAXIS, 1, &source, &down);
    MPI_Cart_shift(Comm2D, YAXIS, -1, &source, &left);
    MPI_Cart_shift(Comm2D, YAXIS, 1, &source, &right);

    /* Alloc memory for neighbours recv messages */
    if (left != MPI_PROC_NULL)
        createMat(mu.rows, 1, &leftbuff);
    if (right != MPI_PROC_NULL)
        createMat(mu.rows, 1, &righbuff);
    if (top != MPI_PROC_NULL)
        createMat(1, mu.cols, &topbuff);
    if (down != MPI_PROC_NULL)
        createMat(1, mu.cols, &downbuff);

    /* Derived datatypes (columns and rows) */
    MPI_Type_contiguous(mu.cols, MPI_DOUBLE, &row_t);
    MPI_Type_vector(mu.rows, 1, mu.cols, MPI_DOUBLE, &col_t);
    MPI_Type_commit(&row_t);
    MPI_Type_commit(&col_t);

    /* Start iteration */
    tstart = MPI_Wtime();

    /*************************************************************
     * CONJUGATE GRADIENT ALGORITHM. 
     ************************************************************/

    /* Step 1: Initialize local vars (u=0, g=-b, d=b) */
    for (int i = 0; i < mrows; i++)
        for (int j = 0; j < mcols; j++)
        {
            vectorb((i + xoffset) * h, (j + yoffset) * h, h, &mb.data[i * mcols + j]);
            mg.data[i * mcols + j] = -1 * mb.data[i * mcols + j];
            md.data[i * mcols + j] = mb.data[i * mcols + j];
            mu.data[i * mcols + j] = 0;
        }

    /* Step 2: q0 = g^T*g */
    dotprod(mg, mg, &q0);
    MPI_Allreduce(&q0, &q0, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    /* Step 3: Iterate */
    for (int iter = 0; iter < MAX_ITER; iter++)
    {
        /* Step 4: q = A*d */
        stencil(md, mq, Comm2D);

        /* Step 5: tau = q0/(d^T*q) */
        dotprod(md, mq, &tau);
        MPI_Allreduce(&tau, &tau, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        tau = q0 / tau;

        /* Step 6: u = u + tau * d */
        scalarprod(tau, md, temp1);
        matsum(mu, temp1, mu);

        /* Step 7: g = g + tau * q */
        scalarprod(tau, mq, temp2);
        matsum(mg, temp2, mg);

        /* Step 8: q1 = g^T*g  */
        dotprod(mg, mg, &q1);
        MPI_Allreduce(&q1, &q1, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        /* Step 9: beta = q1/q0 */
        beta = q1 / q0;

        /* Step 10: d = beta * d - g */
        scalarprod(beta, md, md);
        matsub(md, mg, md);

        /* Step 11: q0 = q1 */
        q0 = q1;
    }

    tend = MPI_Wtime() - tstart;
    MPI_Reduce(&tend, &tmax, 1, MPI_DOUBLE, MPI_MAX, ROOT, MPI_COMM_WORLD);

    /* Print results */
    if (rank == ROOT)
        printf("Time(s): %f, error: %f\n", tmax, sqrt(q1));

    MPI_Comm_free(&Comm2D);
    return MPI_Finalize();
}
