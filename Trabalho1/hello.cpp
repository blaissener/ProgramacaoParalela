#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>

/*
Example to test the openmp enviroment and also timing and 
getting arguments from the command line
*/

enum Error {SUCCESS, BAD_ARGUMENT};

int main(char *argv[]) {
    
    //Testing if the number of arguments is one
    if(argc != 2){
        
        fprintf(stderr, "Wrong number of arguments. Need 1, recived %d", argc);
        return BAD_ARGUMENT;

    }



    int nprocs = omp_get_num_procs();

    printf("\n\n\nRunning on %d processor(s). \n\n\n", nprocs);
    
    omp_set_num_threads(4* nprocs);

    int nthreads, tid;
    


    #pragma omp parallel default(none) private(tid, nthreads) 
    {

        tid = omp_get_thread_num();
        printf("Testing from thread %d. \n", tid);

        if(tid ==0){

            nthreads = omp_get_num_threads();
            printf("The total current number of threads is %d.\n", nthreads);

        }

    }

    printf("Now there is only one thread.\n");

    return 0;
}