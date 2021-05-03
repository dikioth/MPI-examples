#ifndef __CONJUGATE_H__
#define __CONJUGATE_H__

#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/*****************************************************************************
 * GLOBAL VARIABLES 
 *****************************************************************************/

/* Custom variable (matrix) */
typedef struct mat_t
{
    double *data;
    int rows;
    int cols;
} mat_t;


/* Program params */
const int NARGS = 1; // NUmber of arguments (expecting only n)

/*Algorithm params */
const double STENCIL[] = {-1.0, -1.0, 4.0, -1.0, -1.0}; // Five-points stencil
const int MAX_ITER = 200;                               // Max number of interations
double h;                                               // Mesh size
mat_t mu, mb, mg, md, mq, temp1, temp2;                 // Local matrices
double q0, q1, tau, beta;                               // Local scalars

/*MPI params */
const int ROOT = 0; // Root rank
const int XAXIS = 0;
const int YAXIS = 1;

/*MPI mesh params*/
int n, size, sqrtp, splitsize, mcols, mrows, xoffset, yoffset;

int reqn;                                    // Number of requests
MPI_Request reqs[4];                         // Request array
mat_t leftbuff, righbuff, topbuff, downbuff; // Recvbuff from neighbours
int source, left, right, top, down;

/* 2D Cartesian topology */
MPI_Comm Comm2D;
int ndims = 2;             // Number of dims
int dims[2];               // size along each axis
int periods[2] = {0, 0};   // Non periodic boundary
int reorder = 0;           // Rank reordering
MPI_Datatype row_t, col_t; // Custom derived datatypes

/* Performance tests*/
double tstart, tend, tmax;


/*****************************************************************************
 * FUNCTIONS 
 *****************************************************************************/

/**
 * @brief Create a mat_t struct.
 * 
 * @param rows The number of matrix rows
 * @param cols The number of matrix columns
 * @param mat  Pointer to a mat_t struct.
 */
void createMat(int rows, int cols, mat_t *mat);

/**
 * @brief Matrix-vector operation using CG method.
 * The multiplication is in the form: 
 *             q = A * d
 * where A is the following five-points stencil 
 * 
 *              (-1)
 *                |
 *        (-1)--( 4 )--(-1)
 *                |
 *              (-1)
 * 
 * @param in   The mat_t struct respresenting vector d
 * @param out  The mat_t struct where results will be written
 * @param comm The 2D MPI cartesian topology
 */
void stencil(mat_t in, mat_t out, MPI_Comm comm);

/**
 * @brief Scalar-matrix operation in the form
 *         out = scalar * in
 * 
 * @param scalar The double scalar
 * @param in  The mat_t struct matrix
 * @param out A mat_t struct where results will be written
 */
void scalarprod(double scalar, mat_t in, mat_t out);

/**
 * @brief Matrix addition.
 * 
 * @param v1 A mxn matrix
 * @param v2 A mxn matrix
 * @param out A mat_t struct where reults will be written.
 */
void matsum(mat_t v1, mat_t v2, mat_t out);

/**
 * @brief Matrix substraction
 * 
 * @param v1 A mxn matrix
 * @param v2 A mxn matrix
 * @param out A mat_t struct where reults will be written.
 */
void matsub(mat_t v1, mat_t v2, mat_t out);

/**
 * @brief Dot product.
 * 
 * @param v1 A mxn matrix
 * @param v2 A mxn matrix
 * @param out A mat_t struct where reults will be written.
 */
void dotprod(mat_t v1, mat_t v2, double *out);

/**
 * @brief The load vector of the Possion equation.
 * Returns a function value at (x,y)
 * @param x The x coordinate (axis 0)
 * @param y The y coordinate (axis 1)
 * @param h The mesh size
 * @param out A double pointer 
 */
void vectorb(double x, double y, double h, double *out);

#endif // __CONJUGATE_H__
