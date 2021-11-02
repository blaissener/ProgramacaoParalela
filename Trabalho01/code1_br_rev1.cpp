#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <sys/time.h>
#include <string>
#include <vector>
#include <fstream>
#include <algorithm>
#include <iomanip>
/*
    To-Do List:
    [x] - Read file (.net) from terminal
        Done in line 59
    [x] - Get #V and #E from file
        Done in lines 74 and 75
    [x] - Generate graph adijacency list
        Done in lines

*/

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

    for (int i = 0; i < len; i++)
    {
        k[i] = (aL[i].size());
    }

    return k;
}


double rCalc(std::vector<std::vector<int>> aL, std::vector<int> rSet, int n)
{
    /*
        Function that calculates the rich club coefficient of a given aL
        and a set of nodes that have at least the degree "k".

        Parameters:
            aL -> Adjacency List

            rSet -> List of nodes that have at least k degree.

        Return:
            rC -> rK the rich club coefficent for the k degree.

    */

    double sum = 0.0;

    //setting the condition for cases that our set has less than 2 nodes.
    double rK;
    if (n <= 1)
    {
        rK = 1;
    }

    else
    {
        double a = 1.0 / (n * (n - 1.0));
        //Here we run through all nodes that have at least k degree. Note that
        //because we did not especified a direction for a given link, if a
        //node "n" is linked with a node "m", "m" will have the same link.
        //That is important because when a link between the nodes "n" and "m"
        //is observed, we can sum 2 instead of 1 (to the summation of page 1
        //of the pdf), for both directions.
        for (size_t i = 0; i < rSet.size(); i++)
        {   
            for (size_t j = i + 1; j <rSet.size(); j++)
            {
                //Here, we get two numbers, rSet[j] and aL[rSet[i]]. We try to 
                // 
                if (std::find(aL[rSet[i]].begin(), aL[rSet[i]].end(), rSet[j]) != aL[rSet[i]].end())
                {
                    sum += 2.0; 
                }
            }
        }

        rK = a * sum;
    }
    std::cout << "rK " << rK << std::endl;
    return rK;
}

std::vector<double> richClub(std::vector<std::vector<int>> aL, std::vector<int> degree)
{    
    /*
        Function that calculates the rich club coefficient of a given aL
        and a list of degrees.

        Parameters:
            aL -> Adjacency List

            degree -> List with the degrees of all the nodes.

        Return:
            rC -> List with the rich club coefficient of all degrees.

    */
    
    //Finding the major degree
    int maxDegree = *std::max_element(degree.begin(), degree.end());
    
    std::vector<double> rC;
    std::vector<int> rSet;
    std::vector<std::vector<int>> rSetAL; 

    //Creating the set of nodes that are at least k degree. (k goes from 
    // zero to maxDegree - 1)
    
    //run throught the degrees until k is iqual to maxDegree - 1
    for (int k = 0; k < maxDegree; k++)
    { 

        //Starting a fresh set.
        rSet.clear();
        rSetAL.clear();
        
        // search for the degree "k" in the "degree" vector
        for (int l = 0; l < degree.size(); l++)
        { 

            if (degree[l] > k)
            {
                rSet.push_back(l); //---the index of a entry on the "degree" vector is the node
                rSetAL.push_back(aL[l]); //--- the adjacence list of all entries of degree k or more
            }
        }


        
        //getting the cardinality of our set
        int n = rSet.size();

        //Calculating the rich club coefficent with the "rCalc" function.
        double value = rCalc(aL, rSet, n);
        rC.push_back(value);
        
    }

    return rC;
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
    //std::vector<std::vector<int>> inputList(conf[1]);

    for (int i = 0; i < conf[1]; i++)
    {
        //Creating the input graph list.
        int a, b;
        file >> a >> b;
        //inputList[i] = {a, b};

        //Creating the aL.
        aL[a].push_back(b);
        aL[b].push_back(a);
    }

    /*
    --------------------------------------------------------------
    Creating a vector with the degree of each node.
    --------------------------------------------------------------
    
    Note:
        Usign function "degree".
    */

    std::vector<int> deg = degree(aL);

    /*
    --------------------------------------------------------------
    Calculating the rich-club coefficient for our graph.
    --------------------------------------------------------------
    
    Note:
        Usign function "richClub".
    */
    std::vector<double> teste = richClub(aL, deg);
    
    std::string nameIn = argv[1];
    std::string nameOut = nameIn.replace(nameIn.end()-3,nameIn.end(),"rcb");    
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