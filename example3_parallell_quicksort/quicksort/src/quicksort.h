#ifndef __QUICKSORT_H__
#define __QUICKSORT_H__

#define PI 3.14159265358979323846

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

/* Expected 3 args (seq, len , piv) */
#define NUM_ARGS 3

/* Sequence distribution */
#define SEQ_UNIFORM 0
#define SEQ_EXPONENTIAL 1
#define SEQ_NORMAL 2
#define LAMBDA 10

/* Pivot methods */
int PIVOT;
#define PIVOT_METHOD1 0
#define PIVOT_METHOD2 1
#define PIVOT_METHOD3 2

/* MPI params */
#define ROOT 0
#define MSG_TAG 1

/* 1D vector (stores vector and length) */
typedef struct vector_t
{
    double *data;
    int len;
} vector_t;

/**
 * @brief Check if vector is sorted.
 * 
 * @param arr The vector_t vector to be checked.
 * @return int EXIT_SUCCESS if sorted else EXIT_FAILURE
 */
int is_sorted(vector_t arr);

/**
 * @brief Create a 1D vector with sequence distribution
 * specified by 'seq' and length 'len'
 * 
 * @param seq The sequence distrubition. 
 * 0 (SEQ_UNIFORM), 1 (SEQ_EXPONENTIAL) or 2 (SEQ_NORMAL)
 * @param len The length of vector
 * @return vector_t The created 1D vector
 */
vector_t dalloc(int seq, int len);

/**
 * @brief (For debugging purporse) Print the 1D vector.
 * 
 * @param arr The 1D vector to be printed.
 */
void print_array(vector_t arr);

/**
 * @brief Get the median of 1D array. Assuming arr exits.
 * 
 * @param arr The 1D array
 * @param median Pointer where results will be assigned.
 */
void get_median(vector_t arr, double *median);

/**
 * @brief Get the pivot value specified by 'PIVOT' method.
 * 
 * @param arr The current PE's 1D array.
 * @param pivot Pointer where results will be assigned.
 * @param comm The MPI_Comm world/group.
 */
void get_pivot(vector_t arr, double *pivot, MPI_Comm comm);

/**
 * @brief Local sequential quick-sort
 * 
 * @param data The 1D array to be sorted
 * @param left The left start index
 * @param right The right start index
 */
void local_qsort(double *data, int left, int right);

/**
 * @brief Global parallel quick sort.
 * 
 * @param mdata The <vector_t> 1D vector. 
 * @param comm The MPI_Comm world/group.
 */
void global_qsort(vector_t *mdata, MPI_Comm comm);

/**
 * @brief Parition 1D array into elements less and greather
 * than the pivot value (at ipivot).
 * 
 * @param data The 1D array to be partioned
 * @param left The initial left index
 * @param right The initial right index
 * @param ipivot The index of the pivot
 * @return int The index where data splits less and greather partitions.
 */
int partition(double *data, int left, int right, int ipivot);

/**
 * @brief Merge two vectors into one
 * 
 * @param v1 The left 1D vector pointer
 * @param n1 The left 1D vector length
 * @param v2 The right 1D vector pointer
 * @param n2 The right 1D vector length
 * @return double* Pointer to the new merged vector
 */
double *merge(double *v1, int n1, double *v2, int n2);

#endif // __QUICKSORT_H__