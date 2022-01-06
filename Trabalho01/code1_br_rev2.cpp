/*
Bruno Sartorelli Laissener - 10391628
Trabalho 1 - Programação Paralela.
g++ -O code1_br_rev2.cpp -o code1_br_rev2 && ./code1_br_rev2 ./data/in/large-1-001.net
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

std::vector<double> richCCalc(std::vector<std::vector<int>> graphDegree, std::vector<std::vector<int>> aL)
{
    /*
        Function that calculates the rich-club coeficient using the graphDegree list.
        The calculation is based on taking two entries of a certain degree from
        the graphDegree list and searching for a link between them in the aL list.
        This is a second approach to the calculation. The fist one was based only
        in the aL list and it was needed to generate every R(k) set for the calculation.
        
        The current method is a little faster but still not great. For an example, 
        the graph from large-2-009.net took 368 seconds to be calculated. With the old 
        method, the same graph took 374 seconds. As we can see, it is not optimal and
        for larger graphs it is yet impractical.

    */

    int node1;
    int node2;
    std::vector<double> richCoef;

    for (int k = 0; k < graphDegree.size(); k++) //size of graphDegree - 1 is the maximum degree
    {

        double n = graphDegree[k].size();

        double sum = 0.0;

        //setting the condition for cases that our set has less than 2 nodes.
        double rK;
        if (n <= 1.0)
        {
            rK = 1.0;
        }
        else
        {
            double a = 1.0 / (n * (n - 1.0));
            for (int i = 0; i <= n - 1; i++) //Should go from zero to the len of the R(k), which is the same as n.
            {                                //n-1 because j will go from i+1 to the end of the size
                for (int j = i + 1; j < n; j++)
                {
                    node1 = graphDegree[k][i];
                    node2 = graphDegree[k][j];

                    if (std::find(aL[node1].begin(), aL[node1].end(), node2) != aL[node1].end())
                    {
                        //If a link between the two nodes is found, we sum 2 instead of one
                        //because the reverse order of the nodes will also have a link.
                        sum += 2.0;
                    }
                }
            }
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
        Usign function "richCCalc".
    */

    std::vector<double> teste = richCCalc(graphDegree, aL);

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