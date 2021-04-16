#ifndef __MATMUL_H__
#define __MATMUL_H__

#define SUCCESS 0
#define FAILURE -1
#define ROOT 0
#define NUM_ARGS 2

struct person_t
{
    int age;
    double height;
    char name[10];
};

int try(void *call, int ret, char *msg);
// int main(int argc, char argv[]);

/**
 * @brief Read input file. 
 * File must contain 2*n^2+1 numbers separeted by whitespaces. 
 * Numbers are stored in the following sequence:
 * n A B
 * where n is the matrix size, A is the n^2 elements of first matrix in row-major order,
 * B is the n^2 elements of second matrix.
 * 
 * @param fname  Name of file to be read.
 * @param values Pointer to array to store data
 * @return int 0 ix sucess else -1
 */
int readfile(char *fname, double **x, double **y);

/**
 * @brief Write output file
 * All matrix elements are foating-point numbers 
 * and will  be stored with six (6) decimal places.
 * @param fname 
 * @param values 
 * @return int 
 */
int writefile(char *fname, double *output, int rows, int cols);
#endif // __MATMUL_H__