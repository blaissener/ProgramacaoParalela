/*
Bruno Sartorelli Laissener - 10391628
Trabalho 2 - Programação Paralela.

To Compile:
$ mpic++ -O2 -Wall -Wextra -Wpedantic parallelCodeMPI_rev1.cpp -o parallelCodeMPI_rev1

To run:
$ mpirun -oversubscribe -np 5 ./parallelCodeMPI_rev1 huge-2-009.net

or:
$ mpic++ -O2 parallelCodeMPI_rev1.cpp -o parallelCodeMPI_rev1 && mpirun -oversubscribe -np 5 ./parallelCodeMPI_rev1 ./huge-2-009.net
*/

#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <fstream>
#include <algorithm>
#include <iomanip>
#include <chrono>
#include <mpi.h>


#ifndef _SIMPLE_MATRIX_
#define _SIMPLE_MATRIX_

#define DEFMATRIX(type, name)                                                  \
  typedef struct {                                                             \
    type *_data;                                                               \
    size_t nrows, ncols;                                                       \
  } name;

DEFMATRIX(double, MatrixDouble)
DEFMATRIX(float, MatrixFloat)
DEFMATRIX(int, MatrixInt)

typedef MatrixDouble Matrix;

#define IDX(m, i, j) (m)._data[(size_t)(i) * (m).ncols + (size_t)(j)]

#define INIT_MATRIX(mat, nr, nc)                                               \
  {                                                                            \
    (mat).nrows = (nr);                                                        \
    (mat).ncols = (nc);                                                        \
    (mat)._data = malloc((nr) * (nc) * sizeof((mat)._data[0]));                \
  }

#define ZERO(m)                                                                \
  memset((m)._data, 0, (m).nrows *(m).ncols * sizeof((m)._data[0]))

#endif /* _SIMPLE_MATRIX_ */


/*
Functions used in this program
*/




std::vector<int> degree(std::vector<std::vector<int>> aL)
{
    /*
        Function that creates a vector with the degree of 
        all the entries of a adjacency list (aL). The degree
        of a node is the number of links it has. In our
        case, it's the size of each entry of the aL.

        Parameters:
            aL -> Adjacency List
        Return:
            k -> List with the degrees of all the nodes.

    */

    size_t len = aL.size();
    std::vector<int> k(len);

    int deg = 0;

    for (int i = 0; i < len; i++)
    {
        k[i] = (aL[i].size());
    }

    return k;
}

std::vector<std::vector<int>> gDegree(std::vector<std::vector<int>> aL, int maxDeg)
{
    /*
    Function that creates a vector with all nodes of a
    certain degree for all the entries of a adjacency list (aL).
    The degree of a node is the number of links it has. In our
    case, it's the size of each entry of the aL. The degree 
    of our list is the index of the entry.

    Parameters:
        aL -> Adjacency List
        maxDeg -> The max degree of our net.
    Return:
        graphDegree -> List of lists of nodes from each degree 
        (the index is the digree)

    */

    size_t len = aL.size();
    std::vector<std::vector<int>> graphDegree(maxDeg);

    int deg = 0;

    for (int i = 0; i < len; i++)
    {
        int index = aL[i].size();
        for (int j = 0; j < index; j++)
        {
            graphDegree[j].push_back(i);
        }
    }

    return graphDegree;
}


std::vector<double> faster(std::vector<std::vector<int>> graphDegree, std::vector<std::vector<int>> aL, std::vector<int> deg)
{
    /*
        Function that calculates the rich-club coeficient using the aL list and also the deg list.
        The calculation is based on setting a degree k, finding a node that has at least degree k+1
        and for all the links this node has, the ones that also have degree bigger than k+1.

        Note that this method is different from the one we used in the first work, because
        in that case, we used to find 2 nodes that had degree > k+1 and in the list
        of all links, we searched for a link between these two. This works fine for 
        small to medium files, but for large and huge, it could not even retur a  value
        in a resonable time. 

        This was noticed by the professor and now is fixed following the method above.
        
        As this probleme was fixed, we implemented a parallel code that runs every search 
        and comparison between the nodes, in a different thread. This speeds up the process
        and in the extreme case, we could gain 1 second from the sequential algorithm.

        Using the file huge-2-009.net, the sequential algorithm takes 5 seconds and the
        parallel one, 4 seconds. As it is that fast for the biggest file, no more changes
        were necessary.


    */

    std::vector<double> richCoef;

    for (int k = 0; k < graphDegree.size(); k++) //size of graphDegree - 1 (<) is the maximum degree
    {

        int n = 0;

        double sum = 0.0;

        double rK=0.0;

        /*
            Here is where the paralization occurs. We must set the shared variables and also the variables
            that must be kept in summation only in each thread, for that we use the "reduction" statement,
            that copies a new variable for every thread.
            
            We also tested the schedule(static, x) statement, but no resonable changes occurs for different 
            orders of x.

            Last, we tryed to change the method from parallel to sequential based on the degree k
            as the lowests degrees are bigger in size, we could perform better just using the parallel
            mode for k < threshold. Again, no substantial gains occured.
        */
        
        /*
            Note:
            The only function that it is usefull to use different threads is
            the one that calculates the rich club coeficient. Other funtions 
            takes less than one second in sequential mode.
        */

        
        for (int i = 0; i < deg.size(); i++) 
        {                           
            
            if (deg[i] > k) 
            {   
                n++;
                
                for (int j = 0; j < aL[i].size(); j++)
                {
                    if (deg[aL[i][j]] > k)
                        sum += 1;
                };
            }
        }
        
            
        if (n <= 1)
        {
            rK = 1.0;
        }
        else
        {
            double a = 1.0 / (n * (n - 1.0));
            rK = a * sum;
        }


        richCoef.push_back(rK);
    }

    return richCoef;
}



/*
Main Function
*/

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    int rank;

    //MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    int conf0=0;
    int conf1=0;

    if(rank==0)
    {
        std::fstream file(argv[1]);
        //std::vector<int> conf(2);
        file >> conf0 >> conf1;
        
    }
    MPI_Bcast( &conf0 , 1 , MPI_INT , 0 , MPI_COMM_WORLD);
    
    std::vector<int> deg;
    std::vector<std::vector<int>> aL(conf0);

    int degSize = -1;
    int alSize = -1;
    int maxDegSize=-1;
    int alCols = -1;

    if(rank == 0)
    {
        std::fstream file(argv[1]);
        std::vector<int> conf(2);
        file >> conf[0] >> conf[1];

        //std::vector<std::vector<int>> aL(conf[0]);

        for (int i = 0; i < conf[1]; i++)
        {
            //Creating the input graph list.
            int a, b;
            file >> a >> b;

            //Creating the aL.
            aL[a].push_back(b);
            aL[b].push_back(a);
        }

        deg = degree(aL);
        int maxDegree = *std::max_element(deg.begin(), deg.end());
        maxDegSize = gDegree(aL, maxDegree).size();
        degSize = deg.size();
        alSize = aL.size();
        alCols = aL[0].size();
        
    }
    //Every process that participates on the MPI_COMM_WORLD com., has to call the Bcast.
    //So all of them waits until the process "root" sends the message 
    MPI_Bcast( &degSize , 1 , MPI_INT , 0 , MPI_COMM_WORLD);
    MPI_Bcast( &alSize , 1 , MPI_INT , 0 , MPI_COMM_WORLD);
    MPI_Bcast( &maxDegSize , 1 , MPI_INT , 0 , MPI_COMM_WORLD);
    MPI_Bcast( &alCols , 1 , MPI_INT , 0 , MPI_COMM_WORLD);
    //std::cout<< "Process " << rank << " of " << nprocs << " is printing the deg size: "  << degSize << std::endl;
    //std::cout<< "Process " << rank << " of " << nprocs << " is printing the aL size: "  << alSize << std::endl;
    
    int hmProcs = 0;
    MPI_Comm_size( MPI_COMM_WORLD, &hmProcs);
    std::cout<< "Process " << rank << " of " << hmProcs << std::endl;
    
    //Scattering both the vector of degrees and also the aL matrix

    std::vector<int> degCount;
    std::vector<int> degOffset;
    int degOff = 0;

    std::vector<int> alCount;
    std::vector<int> alOffset;
    int alOff = 0;



    int degQ = degSize/hmProcs;
    int degR = degSize%hmProcs;
    
    int alQ = alSize/hmProcs;
    int alR = alSize%hmProcs;

    // if(rank == 0){
    //     std::cout<< q << " " << r << std::endl;
    // }
    
    
    for (int i=0; i< hmProcs; i++)
    {
        degCount.push_back(i<degR ? degQ+1 : degQ);
        
        degOffset.push_back(degOff);

        degOff += degCount[i];

        alCount.push_back(i<alR ? alQ+1 : alQ);
        
        alOffset.push_back(alOff);

        alOff += alCount[i];
        // if(rank == 0){
        //     std::cout<< "Count "<< count[i] << std::endl;
        // }
    }

    std::cout<< "al data def from process " << rank << std::endl;

    /* int *degData = (int *)malloc(degCount[rank]*sizeof(int));
    //int *alData = (int *)malloc(alCount[rank]*sizeof(int));
    int **alData = (int **)malloc(alCount[rank] * sizeof(int*));
    
    for(int i = 0; i < alCount[rank]; i++)
    {
        alData[i] = (int *)malloc(alCols * sizeof(int));
    } */

    Matrix alData;
    INIT_MATRIX(alData, int(alCount[rank]), int(alCols));

    std::cout<< "al data DEFINED from process " << rank << std::endl;

    MPI_Datatype line;
    MPI_Type_contiguous( alCols , MPI_INT , &line);
    MPI_Type_commit( &line);

    //MPI_Scatterv(&deg[0] , &degCount[0] , &degOffset[0] , MPI_INT , &degData[0] , degCount[rank] , MPI_INT , 0 , MPI_COMM_WORLD);
    std::cout<< "Starting scatterv from process " << rank << std::endl;
    //MPI_Scatterv(&aL[0][0] , &alCount[0] , &alOffset[0] , line , &IDX(alData, 0, 0) , alCount[rank] , line , 0 , MPI_COMM_WORLD);
    /* for(int i=0; i<sizeof(alData); i++){
        for (int j=0; j<sizeof(alData[i]);j++){
            std::cout<< "Al Data "<< alData[i][j] << " from process" << rank << std::endl;
        }
    } */
    std::cout<< "Scatterv DONE process " << rank << std::endl;
    //std::cout<< "Process " << rank << " of " << hmProcs << " is printing the deg count size: "  << degCount[rank] << std::endl;
    //std::cout<< "Process " << rank << " of " << hmProcs << " is printing the aL count size: "  << alCount[rank] << std::endl;
    
    // std::vector<double> richCoef;
    // for(int k = 0; k < maxDegSize; k++)
    // {
    //     int n = 0;

    //     double sum = 0.0;

    //     double rK=0.0;

    //     /*
    //         Here is where the paralization occurs. We must set the shared variables and also the variables
    //         that must be kept in summation only in each thread, for that we use the "reduction" statement,
    //         that copies a new variable for every thread.
            
    //         We also tested the schedule(static, x) statement, but no resonable changes occurs for different 
    //         orders of x.

    //         Last, we tryed to change the method from parallel to sequential based on the degree k
    //         as the lowests degrees are bigger in size, we could perform better just using the parallel
    //         mode for k < threshold. Again, no substantial gains occured.
    //     */
        
    //     /*
    //         Note:
    //         The only function that it is usefull to use different threads is
    //         the one that calculates the rich club coeficient. Other funtions 
    //         takes less than one second in sequential mode.
    //     */

        
    //     for (int i = 0; i < sizeof(degData); i++) 
    //     {                           
            
    //         if (degData[i] > k) 
    //         {   
    //             n++;
                
    //             for (int j = 0; j < sizeof(alData[i]); j++)
    //             {
    //                 if (degData[alData[i][j]] > k)
    //                     sum += 1;
    //             };
    //         }
    //     }
        
            
    //     if (n <= 1)
    //     {
    //         rK = 1.0;
    //     }
    //     else
    //     {
    //         double a = 1.0 / (n * (n - 1.0));
    //         rK = a * sum;
    //     }


    //     richCoef.push_back(rK);
    // }
    //Opening the file.
    
    // for(int i=0; i < richCoef.size(); i++)
    // {
        
    // }

    /*
    --------------------------------------------------------------
    Getting how the graph is configured (V and E).
    --------------------------------------------------------------
    
    Note: 
        When we use ">>", the file removes the first element.
        At the end of the next 2 lines, we will have just the
    graph without any of the configuration data (V -conf[0] 
    and E -conf[1]).
    */
    

    //Starting a clock to measure the time the algorithm takes.
    //auto t1 = std::chrono::high_resolution_clock::now();

    /*
    --------------------------------------------------------------
    Creating the adjacency list (aL).
    --------------------------------------------------------------

    Note:
        We decided to create the aL to make it easier to 
    calculate the degree of each node as the degree will be just
    the lenght of the vector at a given position.
        The list is prefered over the matrix because the data is 
    sparce. This saves memory with the "0" entries.
        
    */

    

    /*
    --------------------------------------------------------------
    Creating a vector with the degree of each node and getting the 
    max degree.
    --------------------------------------------------------------
    
    Note:
        Usign function "degree".
    */

    
    

    /*
    --------------------------------------------------------------
    Creating a vector with all nodes of each degree.
    --------------------------------------------------------------
    
    Note:
        Using function "gDegree"
    
    */
    //std::vector<std::vector<int>> graphDegree = gDegree(aL, maxDegree);

    /*
    --------------------------------------------------------------
    Calculating the rich-club coefficient for our graph.
    --------------------------------------------------------------
    
    Note:
        Usign function "faster".
    */

    //omp_set_num_threads(8);
    //std::vector<double> teste = faster(graphDegree, aL, deg);


    

    //Stoping the clock
    /* auto t2 = std::chrono::high_resolution_clock::now();

    auto delta = std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count();
    std::cout << "Time it took: " << delta << " seconds." << std::endl;

    //Writing out the calculated rich-club coeficients.
    std::string nameIn = argv[1];
    std::string nameOut = nameIn.replace(nameIn.end() - 3, nameIn.end(), "rcb");
    std::ofstream out(nameOut);
    for (size_t b = 0; b < teste.size(); b++)
    {

        out << std::fixed;
        out << std::setprecision(5);
        out << teste[b] << std::endl;
    }
    out.close(); */

    MPI_Finalize();
    return 0;
}