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


std::vector<std::vector<double>> calcRC(int rank, std::vector<int> degCount, std::vector<int> degOffset, std::vector<std::vector<int>> aL, std::vector<int> deg)
{
    std::vector<std::vector<double>> richCoef(degCount[rank]);
    for(int k=degOffset[rank]; k<degCount[rank]+degOffset[rank]; k++){
        //std::cout<<k<<" from process " << rank <<std::endl;       

        int n = 0;

        double sum = 0.0;
        double rK=0.0;
            
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

        /* std::vector<double> partial;
        richCoef.push_back(partial); */
        /* partial.push_back(rK);
        partial.push_back(k);
        richCoef.push_back(partial); */

        richCoef[k - degOffset[rank]].push_back(rK);
        richCoef[k - degOffset[rank]].push_back(k);

        // richCoef[0]=rK;
        // richCoef[1]=k;
        
    }
    return richCoef;
}

std::vector<double> calcRCNoK(int rank, std::vector<int> degCount, std::vector<int> degOffset, std::vector<std::vector<int>> aL, std::vector<int> deg)
{
    std::vector<double> richCoef;
    for(int k=degOffset[rank]; k<degCount[rank]+degOffset[rank]; k++){
        //std::cout<<k<<" from process " << rank <<std::endl;       

        int n = 0;

        double sum = 0.0;
        double rK=0.0;
            
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

        /* std::vector<double> partial;
        richCoef.push_back(partial); */
        /* partial.push_back(rK);
        partial.push_back(k);
        richCoef.push_back(partial); */

        richCoef.push_back(rK);
        

        // richCoef[0]=rK;
        // richCoef[1]=k;
        
    }
    return richCoef;
}


/*
Main Function
*/

int main(int argc, char *argv[])
{   
    std::setprecision(5);
    MPI_Init(&argc, &argv);
    int rank, hmProcs;

    MPI_Comm_size(MPI_COMM_WORLD, &hmProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int degSize=-1, vertex=-1, edge=-1, col1Size=-1, col2Size=-1; //conf0 conf1
    
    std::vector<int> col1, col2;

    if(rank == 0)
    {
        std::fstream file(argv[1]);
        
        file >> vertex >> edge;

        col1.reserve(edge);
        col2.reserve(edge);

        for (int i = 0; i < edge; i++)
        {
            //Creating the input graph list.
            int a, b;
            file >> a >> b;
            
            

            //Creating the aL.
            col1.emplace_back(a);
            col2.emplace_back(b);
            
            
            
        }

        col1Size=col1.size();
        col2Size=col2.size();
        // for(int j=0; j<col1Size;j++){
        //     std::cout << "a "  << col1[j] << " " << j <<std::endl;
        // }
        /* deg = degree(aL);       
        degSize=deg.size();*/
        // std::cout << "col1size "  << col1Size <<std::endl;
        // std::cout << "col2size "  << col2Size <<std::endl; 
    }
    

    MPI_Bcast( &col1Size , 1 , MPI_INT , 0, MPI_COMM_WORLD);
    MPI_Bcast( &col2Size , 1 , MPI_INT , 0, MPI_COMM_WORLD);

    MPI_Bcast( &vertex, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast( &edge, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    if(rank!=0){
        col1.reserve(edge);
        col2.reserve(edge);
    }
    //std::cout << "col1Size " << col1Size << " from rank " << rank <<std::endl;
    
    MPI_Bcast(&col1[0], col1Size, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&col2[0], col2Size, MPI_INT, 0, MPI_COMM_WORLD);
    
    std::vector<std::vector<int>> aL(vertex);

    for(int i=0; i<edge; i++){
        aL[col1[i]].push_back(col2[i]);
        aL[col2[i]].push_back(col1[i]);
    }
 
    //std::cout << "Al size "  << aL.size() << " from rank " << rank <<std::endl;



    std::vector<int> deg = degree(aL);

    // std::cout << "deg size "  << deg.size() << " from rank " << rank <<std::endl;

    int maxDegree = *std::max_element(deg.begin(), deg.end());
    // std::cout << "Max deg "  << maxDegree << " from rank " << rank <<std::endl;

    std::vector<int> degCount;
    std::vector<int> degOffset;
    int degOff = 0;


    int degQ = maxDegree/hmProcs;
    int degR = maxDegree%hmProcs; 
    
    
    for (int i=0; i< hmProcs; i++)
    {
        degCount.push_back(i<degR ? degQ+1 : degQ);
        
        degOffset.push_back(degOff);

        degOff += degCount[i];

    }

    //std::vector<std::vector<double>> selfK = calcRC(rank, degCount, degOffset, aL, deg);
    std::vector<double> sKnoK = calcRCNoK(rank, degCount, degOffset, aL, deg);
    /* if(rank==1){
        
        for(int i=0; i<selfK.size(); i++){
                for (int j =0 ; j<selfK[i].size(); j++){
                    std::cout<< selfK[i][j]<< " ";
                }
                std::cout<< std::endl;
            }
    } */
    // if(rank==0){
        
    //     for(int i=0; i<sKnoK.size(); i++){
    //         std::cout<< sKnoK[i] << std::endl;
    //         }
    // }

    std::vector<double> kFinal;


    
    if(rank==0){
        kFinal.reserve(maxDegree); 
        
    }
    
    //MPI_Gatherv( &sKnoK[0] , degCount[rank] , MPI_DOUBLE , &kFinal[0] , degCount[rank] , MPI_DOUBLE , 0 , MPI_COMM_WORLD);
    
    
    

    
    //MPI_Gatherv( const void* sendbuf , int sendcount , MPI_Datatype sendtype , void* recvbuf , const int recvcounts[] , const int displs[] , MPI_Datatype recvtype , int root , MPI_Comm comm);
    std::cout<<"Kfinal "<<kFinal.size() << "from rank " << rank << std::endl;

    // std::cout<<"sKnoK "<<sKnoK.size() << "from rank " << rank << std::endl;
    // std::cout<<"degCount "<<degCount.size() << "from rank " << rank << std::endl;
    // std::cout<<"degoff "<<degOffset.size() << "from rank " << rank << std::endl;
    if(rank==0){
        MPI_Gatherv( &sKnoK[0] , degCount[rank] , MPI_DOUBLE , &kFinal[0] , &degCount[0] , &degOffset[0] , MPI_DOUBLE , 0 ,MPI_COMM_WORLD);
        std::cout<<kFinal.size() << std::endl;
        for(int i=0; i<kFinal.size(); i++){
            std::cout<< kFinal[i] << std::endl;
            }
        
        
        // std::string nameIn = argv[1];
        // std::string nameOut = nameIn.replace(nameIn.end() - 3, nameIn.end(), "rcb");
        // std::ofstream out(nameOut);
        
        // for (size_t b = 0; b < kFinal.size(); b++)
        // {

        //     out << std::fixed;
        //     out << std::setprecision(5);
        //     out << kFinal[b] << std::endl;
        // }
        // out.close();
    }
    
    
    /* std::vector<int> ks;
    ks.reserve(maxDegree);

    if(rank==0){
        
        for(int i=0; i<maxDegree; i++){
                ks[i]=i;
            }
    }
    
    
    
    std::vector<int> selfKs;
    
    selfKs.reserve(degCount[rank]); */
    
    // for(int i=0; i<selfKs.size(); i++){
    //             std::cout<< selfKs[i]<< " from process " << rank << std::endl;
    //         }

    //std::cout<< " from process " << rank << std::endl;
    /* MPI_Scatterv( &ks[0] , &degCount[0] , &degOffset[0] , MPI_INT , &selfKs[0] , degCount[rank] ,MPI_INT , 0 , MPI_COMM_WORLD);
    std::cout<<selfKs.size()<< " from process " << rank << std::endl; */

    // for(int i=0; i<selfKs.size(); i++){
    //             std::cout<< selfKs[i]<< " from process " << rank << std::endl;
    //         }

    //std::cout<< degCount[rank] << " from process " << rank << std::endl;
    // std::cout<< degOffset[rank]<< "from process " << rank << std::endl;



    MPI_Finalize();
    return 0;
}