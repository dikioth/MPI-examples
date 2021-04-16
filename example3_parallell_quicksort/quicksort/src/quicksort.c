/**********************************************************************
 * Quick sort
 * Usage: ./a.out sequence length
 *
 **********************************************************************/
#include "quicksort.h"

int is_sorted(vector_t arr)
{
    for (int i = 0; i < arr.len - 1; i++)
        if (arr.data[i] > arr.data[i + 1])
        {
            printf("Wrong: data[%d] = %f, data[%d] = %f\n",
                   i, arr.data[i], i + 1, arr.data[i + 1]);
            return EXIT_FAILURE;
        }

    return EXIT_SUCCESS;
}

vector_t dalloc(int seq, int len)
{
    vector_t arr = {.data = (double *)malloc(len * sizeof(double)),
                    .len = len};

    /* Uniform random numbers */
    if (seq == SEQ_UNIFORM)
    {
        for (int i = 0; i < len; i++)
            arr.data[i] = drand48();
    }

    /* Exponential distribution */
    else if (seq == SEQ_EXPONENTIAL)
    {
        for (int i = 0; i < len; i++)
            arr.data[i] = -LAMBDA * log(1 - drand48());
    }

    /* Normal distribution*/
    else if (seq == SEQ_NORMAL)
    {
        for (int i = 0; i < len; i++)
        {
            double x = drand48();
            double y = drand48();
            arr.data[i] = sqrt(-2 * log(x)) * cos(2 * PI * y);
        }
    }

    return arr;
}

void print_array(vector_t arr)
{
    for (int i = 0; i < arr.len; i++)
        printf(i == 0 ? "[%.2f " : ((i == arr.len - 1) ? " %.2f]\n" : "%.2f "), arr.data[i]);
}

void global_qsort(vector_t *mdata, MPI_Comm comm)
{
    double pivot;
    int recv_size;
    MPI_Status recv_stat;
    MPI_Comm group_comm;
    int local_size, local_rank;

    MPI_Comm_rank(comm, &local_rank);
    MPI_Comm_size(comm, &local_size);

    /* Exit recursion when a group only has one (1) process */
    if (local_size <= 1)
        return;

    /* Step 2.1: Select pivot (withing group) and brocast to all PE's */
    get_pivot(*mdata, &pivot, comm);
    MPI_Bcast(&pivot, 1, MPI_DOUBLE, ROOT, comm);

    /* Step 2.2: Partition local array (already sorted at step 1)*/
    int ipivot = (*mdata).len;
    for (int i = 0; i < (*mdata).len; i++)
        if ((*mdata).data[i] > pivot)
        {
            ipivot = i;
            break;
        }

    vector_t mleft = {.data = (*mdata).data,
                      .len = ipivot};
    vector_t mright = {.data = (*mdata).data + ipivot,
                       .len = (*mdata).len - ipivot};

    /* Step 2.3: Exchange data pariwise (mgroup is 0 or 1) */
    int mgroup = local_rank / (local_size / 2);
    int dest = (local_rank + local_size / 2) % local_size;

    /* Send msg (Group 0 send right block, Group sends left block) */

    double *recv_buff;
    if (mgroup == 0)
    {

        MPI_Send((mgroup == 0) ? mright.data : mleft.data,
                 (mgroup == 0) ? mright.len : mleft.len,
                 MPI_DOUBLE, dest, MSG_TAG, comm);

        /* Recv msg (Probe, Get length msg, allocate memory and receive msg)*/
        MPI_Probe(dest, MSG_TAG, comm, &recv_stat);
        MPI_Get_count(&recv_stat, MPI_DOUBLE, &recv_size);
        recv_buff = (double *)malloc(recv_size * sizeof(double));
        MPI_Recv(recv_buff, recv_size, MPI_DOUBLE, dest, MSG_TAG, comm, &recv_stat);
    }
    else if (mgroup == 1)
    {

        /* Recv msg (Probe, Get length msg, allocate memory and receive msg)*/
        MPI_Probe(dest, MSG_TAG, comm, &recv_stat);
        MPI_Get_count(&recv_stat, MPI_DOUBLE, &recv_size);
        recv_buff = (double *)malloc(recv_size * sizeof(double));
        MPI_Recv(recv_buff, recv_size, MPI_DOUBLE, dest, MSG_TAG, comm, &recv_stat);

        MPI_Send((mgroup == 0) ? mright.data : mleft.data,
                 (mgroup == 0) ? mright.len : mleft.len,
                 MPI_DOUBLE, dest, MSG_TAG, comm);
    }

    /* 'Repoint' left and right pointers */
    mleft.data = (mgroup == 0) ? mleft.data : recv_buff;
    mleft.len = (mgroup == 0) ? mleft.len : recv_size;
    mright.data = (mgroup == 0) ? recv_buff : mright.data;
    mright.len = (mgroup == 0) ? recv_size : mright.len;

    /* Merge left and right data */
    double *temp = (*mdata).data;
    (*mdata).data = merge(mleft.data, mleft.len, mright.data, mright.len);
    (*mdata).len = mleft.len + mright.len;
    free(temp);

    /* Step 2.4: Divide processes two groups (evenly) and recursevly sort */
    int color = local_rank / (local_size / 2);
    MPI_Comm_split(comm, color, local_rank, &group_comm);
    global_qsort(mdata, group_comm);
}

void get_median(vector_t arr, double *median)
{
    *median = arr.len % 2 == 0
                  ? 0.5 * (arr.data[arr.len / 2 - 1] + arr.data[arr.len / 2])
                  : arr.data[arr.len / 2];
}

void get_pivot(vector_t arr, double *pivot, MPI_Comm comm)
{
    double median;
    int mrank, msize;
    MPI_Comm_rank(comm, &mrank);
    MPI_Comm_size(comm, &msize);

    /* Method 1: Median in one processor 
     * Brodcast the mediam of local array of ROOT rank */
    if (PIVOT == PIVOT_METHOD1)
    {
        if (mrank == ROOT)
            get_median(arr, pivot);
    }

    /* Mehod 2: Median of all medians 
     * Get local median. ROOT gathers all medians and gets median of medians. */
    else if (PIVOT == PIVOT_METHOD2)
    {
        vector_t medians;
        get_median(arr, &median);

        if (mrank == ROOT)
            medians = (vector_t){.data = malloc(msize * sizeof(double)),
                                 .len = msize};

        MPI_Gather(&median, 1, MPI_DOUBLE,
                   medians.data, 1, MPI_DOUBLE, ROOT, comm);

        if (mrank == ROOT)
        {
            local_qsort(medians.data, 0, medians.len - 1);
            get_median(medians, pivot);
        }
    }

    /* Method 3: Mean value of all medians 
     * Similar to method 2 but mean instaed of median */
    else if (PIVOT == PIVOT_METHOD3)
    {
        get_median(arr, &median);
        MPI_Reduce(&median, pivot, 1, MPI_DOUBLE, MPI_SUM, ROOT, comm);
        if (mrank == ROOT)
            *pivot /= msize;
    }
}

int partition(double *data, int left, int right, int ipivot)
{

    double pivot = data[ipivot];
    double temp = data[ipivot];
    int isplit = left;

    data[ipivot] = data[right];
    data[right] = temp;

    for (int i = left; i < right; i++)
        if (data[i] <= pivot)
        {
            temp = data[i];
            data[i] = data[isplit];
            data[isplit] = temp;
            isplit += 1;
        }
    temp = data[isplit];
    data[isplit] = data[right];
    data[right] = temp;

    return isplit;
}

void local_qsort(double *data, int left, int right)
{
    /* Index pivot, new index pivot */
    int ipivot, nipivot;

    if (right > left)
    {
        ipivot = left + (right - left) / 2;
        nipivot = partition(data, left, right, ipivot);
        local_qsort(data, left, nipivot - 1);
        local_qsort(data, nipivot + 1, right);
    }
}

double *merge(double *v1, int n1, double *v2, int n2)
{
    int i = 0, j = 0, k = 0;
    double *result = (double *)malloc((n1 + n2) * sizeof(double));

    while (i < n1 && j < n2)
        if (v1[i] < v2[j])
            result[k++] = v1[i++];
        else
            result[k++] = v2[j++];

    if (i == n1)
        while (j < n2)
            result[k++] = v2[j++];
    else
        while (i < n1)
            result[k++] = v1[i++];

    return result;
}

int main(int argc, char *argv[])
{
    vector_t data, mdata;
    int grank, gsize, mlen;
    double tstart, tend;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &grank);
    MPI_Comm_size(MPI_COMM_WORLD, &gsize);

    if (argc != NUM_ARGS + 1 && grank == ROOT)
    {
        printf("usage: ./quicksort sequence length pivot\n");
        printf("sequence: The random number sequence \n");
        printf("length  : The length of the sequence\n");
        printf("pivot   : The method of choosing the global pivot element. 0,1 or 2.\n");
        printf("\t   0: Median of a single PE in group\n");
        printf("\t   1: Median of medians of all PE's in group\n");
        printf("\t   2: Mean of medians of all PE's in group\n");
        return EXIT_FAILURE;
    };

    if (grank == ROOT)
    {
        int seq = atoi(argv[1]);
        int len = atoi(argv[2]);
        mlen = len / gsize;

        /* Create 1D data array */
        data = dalloc(seq, len);
    }

    /* Get chosen pivoting method */
    PIVOT = atoi(argv[3]);

    /* Brodcast sub-arrays size */
    MPI_Bcast(&mlen, 1, MPI_INT, ROOT, MPI_COMM_WORLD);

    /* Allocate memory for local arrays */
    mdata.len = mlen;

    if (NULL == (mdata.data = (double *)malloc(mlen * sizeof(double))))
    {
        printf("Could not allocate memory\n");
        return MPI_Abort(MPI_COMM_WORLD, 1);
    }

    /* START TIMER */
    if (grank == ROOT)
        tstart = MPI_Wtime();

    /* Split data to all PE's */
    MPI_Scatter(data.data, mdata.len, MPI_DOUBLE,
                mdata.data, mdata.len, MPI_DOUBLE,
                ROOT, MPI_COMM_WORLD);

    /***************************************
     * Step 1: Local sequential quick sort 
     ***************************************/
    local_qsort(mdata.data, 0, mdata.len - 1);

    /***************************************
     * Step 2: Global parallel quick srot 
     **************************************/
    global_qsort(&mdata, MPI_COMM_WORLD);

    /***************************************
     * Step 3: Global gather data (at ROOT)
     **************************************/
    int *recvCounts = (int *)malloc(gsize * sizeof(int));
    int *recvdispls = (int *)calloc(gsize, sizeof(int));

    MPI_Gather(&mdata.len, 1, MPI_INT,
               recvCounts, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
    for (int i = 1; i < gsize; i++)
        recvdispls[i] = recvdispls[i - 1] + recvCounts[i - 1];

    MPI_Gatherv(mdata.data, mdata.len, MPI_DOUBLE,
                data.data, recvCounts, recvdispls, MPI_DOUBLE,
                ROOT, MPI_COMM_WORLD);

    free(recvCounts);
    free(recvdispls);

    /* Check results */
    if (grank == ROOT)
    {
        if (is_sorted(data) == EXIT_SUCCESS)
        {
            tend = MPI_Wtime() - tstart;
            printf("%f\n", tend);
        }
        free(data.data);
    }

    MPI_Finalize();
    return EXIT_SUCCESS;
}