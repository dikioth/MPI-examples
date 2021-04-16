/**
 * This program applies a five point stencil on a vector of function values
 * f(x) for 0<=x<2*PI in order to approximate the derivative of the function.
 * The program takes 3 arguments:
 * - the path to the input file.
 * - the path to the output file. Note that this file will be overwritten!
 * - an integer specifying how many times the stencil will be applied
 * The formats of the files are described in the documentation of the functions
 * read_input and read_output. The program also prints the number of seconds
 * needed for the stencil application to standard out.
 */

#ifndef _ASSIGNMENT1_STENCIL_H_
#define _ASSIGNMENT1_STENCIL_H_

#define PI 3.14159265358979323846
#define PRODUCE_OUTPUT_FILE
#define ROOT 0

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

/**
 * Read function data from an input file and store in an array. The input file
 * is supposed to contain an integer representing the number of function values,
 * followed by the function values. All values should be separated by white
 * spaces.
 * @param file_name Name of input file
 * @param values Pointer to array where the values are to be stored
 * @return 0 on success, -1 on error
 */
int read_input(const char *file_name, double **values);

/**
 * Write function data to a file, with 4 decimal places. The values are
 * separated by spaces.
 * @param file_name Name of output file
 * @param output Function values to write
 * @param num_values Number of values to print
 * @return 0 on success, -1 on error
 */
int write_output(char *file_name, const double *output, int num_values);

#endif /* _ASSIGNMENT1_STENCIL_H_ */
