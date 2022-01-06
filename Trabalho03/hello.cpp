/* To compile, 
$ mpic++ -O2 -Wall -Wextra -Wpedantic hello.cpp -o hello

to run, 
mpirun -np 2 ./hello

or to run more than one process for each core

mpirun -oversubscribe -np 5 ./hello
*/
#include <mpi.h>
#include <stdio.h>

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    int nprocs, rank;

    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    std::cout<< "Process " << rank << " of " << nprocs << " is running." << std::endl;
    
    MPI_Finalize();
    return 0;
}