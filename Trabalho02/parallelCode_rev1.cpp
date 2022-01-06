/*
Bruno Sartorelli Laissener - 10391628
Trabalho 2 - Programação Paralela.

To test:
g++ -O -fopenmp parallelCode_rev1.cpp -o parallelCode_rev1 && ./parallelCode_rev1 large-1-001.net
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
#include <omp.h>

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

        #pragma omp parallel for default(none) shared(k, aL, rK, deg), reduction(+ : sum), reduction(+ : n)

        
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
    //Opening the file.
    std::fstream file(argv[1]);

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
    std::vector<int> conf(2);
    file >> conf[0] >> conf[1];

    //Starting a clock to measure the time the algorithm takes.
    auto t1 = std::chrono::high_resolution_clock::now();

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

    std::vector<std::vector<int>> aL(conf[0]);

    for (int i = 0; i < conf[1]; i++)
    {
        //Creating the input graph list.
        int a, b;
        file >> a >> b;

        //Creating the aL.
        aL[a].push_back(b);
        aL[b].push_back(a);
    }

    /*
    --------------------------------------------------------------
    Creating a vector with the degree of each node and getting the 
    max degree.
    --------------------------------------------------------------
    
    Note:
        Usign function "degree".
    */

    
    std::vector<int> deg = degree(aL);
    int maxDegree = *std::max_element(deg.begin(), deg.end());

    /*
    --------------------------------------------------------------
    Creating a vector with all nodes of each degree.
    --------------------------------------------------------------
    
    Note:
        Using function "gDegree"
    
    */
    std::vector<std::vector<int>> graphDegree = gDegree(aL, maxDegree);

    /*
    --------------------------------------------------------------
    Calculating the rich-club coefficient for our graph.
    --------------------------------------------------------------
    
    Note:
        Usign function "faster".
    */

    omp_set_num_threads(8);
    std::vector<double> teste = faster(graphDegree, aL, deg);


    

    //Stoping the clock
    auto t2 = std::chrono::high_resolution_clock::now();

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
    out.close();

    return 0;
}