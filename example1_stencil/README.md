# 1D stencil

Introduction
============

Stencil application is a common operation in computer science. It is
performed many times on large number of arrays and the performance is
highly dependent on it. Researches have developed different methods for
improving the performance of the stencil application. One such method is
the use of parallelization techniques, such as Message Passing Interface
(`MPI`).

Theory
======

The stencil application can be seen as a weighted sum of vectors, matrix
or tensors. It is computed using the current node and its neighbours
(See Fig. below). The boundary is periodic, meaning that the
stencil operation will be applied to it’s neighbours as a torus.

![Stencil application on a vector, matrix and
tensor.](figs/stencil_123D.png "fig:")

One application of the stencil operation is approximating derivatives of
a function. The first derivative using five-points stencil is defined in
Eq. below . Notice that the center point `f(x)` is not involved.

![fprim.](figs/fprim.png)

where `h` is the space between points in the grid. The stencil elements
are

![](figs/stencils.png)

Method
======

The parallelization of the stencil operation is done using . Fig.
below  illustrates the parallelization of a 1D array. Let
`P` be the number of processes. The root process (rank 0 in this case)
reads the input file and distributes equally splitted `P` sub-arrays to
all the processes (including itself). Each process applies the stencil
operation on their local copy of the sub-arrays and the root process
collects back the result of each process and saves it in an output file.

![parallelization of stencil operation. $P$ is the number of processes
and $S$ is the number of stencil
applications.](figs/diagram_stencil.png ) 

The five-points stencil operation needs two (2) elements from the left
and right neighbours at the array edges. The method implemented is to
allocate memory for two extra elements at right and left edges. (See
Fig. below). Then the edge elements are sent and received
to and from its neighbours.

![Send and received extension elements using MPI
functions.](figs/rank_sendrecv.png )

The list of all important functions used in the project is shown in
Table below



|**MPI function**|   **Description/used to**|
|------------------|-----------------------------------------|
|`MPI_Cart_Create`|     Creates a new 1D left and right periodic cartesian topology |
|`MPI_Bcast`|   Broadcast the length of array to all processes|
|`MPI_Scatter`|  Split contiguous input array into subarrays to all processes|
|`MPI_Cart_Shift`| Get rank left and right neighbors.|
|`MPI_Irecv`|Non-blocking receive of edge elements from left and right neighbours.|
|`MPI_Isend`|Non-blocking send of edge elements to left and right neighbours|
|`MPI_Gather`|Gather all subarrays into a single array in root process|
|`MPI_Reduce`|Get the max value of execution time using the operation|



Result
======

The program was tested on the server with the following specifications

-   **CPU:** AMD Opteron (Bulldozer) 6282SE, 2.6 GHz, 16-core, dual
    socket

-   **Memory:** 128GB

-   **OS:** Scientific Linux 6.10

The results of running 1,000 stencil applications on a 4M input size
array with different number of processes are shown in Fig.
below. The strong and weak scalability tests are presented below

![Speedup. 4M input size. 1,000 stencil
applications.](figs/speedup.png "fig:") [fig:speedup]


### **Strong scale test**
| # Processes | # Stencils applications | Computation time | 
|  ----| -------| ---------|
|  1  |  10000  | 29.8895|
|  2  |  10000  | 13.1404|
|  4  |  10000  | 6.2794|
|  8  |  10000  | 3.2142|
|  16 |  10000  | 2.1151|
 

### **Weak scale test**
| # Processes | # Stencils applications | Computation time | 
|  ----| -------| ---------|
  1 |   10000  |  29.7540
  2  |  20000   | 26.5911
  4   | 40000   | 25.1541
  8   | 80000   | 24.3669
  16  | 160000  | 27.1398
  ---- -------- ---------

Discussion
==========

The parallelized method could be seen as three separate stencil
applications on the left, center and right side of the array. The center
stencil application is independent on the process neighbours and the
left and right stencil applications must communication with respectively
neighbour in order to proceed the stencil computation. It is thus
important that the program implements non-blocking communication for the
left and right stencil applications. Thus, A decreased performance would
have achieved if the blocking communication and were used instead of the
non-blocking communication and respectively.

Overall, the program is proven to be efficiently parallelized. It shows
a reasonable speedup behaviour as well as a strong and weak scalability.
However, there is more potential for further optimization which is
outside the scope of the assignment.
