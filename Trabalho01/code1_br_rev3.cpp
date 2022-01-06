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
    std::vector<std::vector<int>> graphDegree(maxDeg);
     

    int deg = 0;

    for (int i = 0; i < len; i++)
    {
        int index = aL[i].size();
        for(int j =0; j<index; j++)
        {
            graphDegree[j].push_back(i);
        }
        
    }
    
    // for(int i =0; i<graphDegree.size();i++)
    // {
    //     for(int j=0; j<graphDegree[i].size(); j++)
    //     {
    //         std::cout << graphDegree[i][j] << " ";
    //     }
        
    //     std::cout<<std::endl;
    // }

    //max degree is the len of graphDegree
    //std::cout<<"OK"<<std::endl;
    return graphDegree;
}

std::vector<double> richCCalc(std::vector<std::vector<int>> graphDegree,std::vector<std::vector<int>> aL)
{   
    int node1;
    int node2;
    std::vector<double> richCoef;
    std::vector<std::vector<int>> summation(graphDegree.size());
    std::vector<double> ns;
    
    for(int k=0; k<graphDegree.size(); k++) //size of graphDegree - 1 is the maximum degree
    {   
        //std::cout << "Degree: " << k << std::endl <<std::endl;
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
            for(int i=0; i<=n-1; i++)//Should go from zero to the len of the R(k), which is the same as n.
            {                       //n-1 because j will go from i+1 to the end of the size
                for(int j =i+1; j<n; j++)
                {
                    node1 = graphDegree[k][i];
                    node2 = graphDegree[k][j];
                    
                    

                    //std::cout << "há ligação entre " << node1 << " e " << node2 << "??" << std::endl;
                    
                    if(std::find(aL[node1].begin(), aL[node1].end(), node2) != aL[node1].end())
                    {   
                        //std::cout <<  "Encontrado ligação entre " << node1 << " e " << node2 << " em al " << i  <<std::endl;  
                        sum+=2.0;
                        //std::cout << "soma atual: " << sum << std::endl;
                    }
                    else
                    {
                    }
                    
                }
            }
        rK = a*sum;
        

        }
    richCoef.push_back(rK);
    //std::cout << "r" << k << " encontrado: " << rK << std::endl;
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
    auto t1 = std::chrono::high_resolution_clock::now();
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
    
    // for(int i =0; i<aL.size();i++)
    // {
    //     for(int j=0; j<aL[i].size(); j++)
    //     {
    //         std::cout << aL[i][j] << " ";
    //     }
        
    //     std::cout<<std::endl;
    // }
    // std::cout<<std::endl;
    // std::cout<<std::endl;


    std::vector<int> deg = degree(aL);
   
    int maxDegree = *std::max_element(deg.begin(), deg.end());
    //std::cout << maxDegree << std::endl;
    /*
    --------------------------------------------------------------
    Creating a vector with the degree of each node.
    --------------------------------------------------------------
    
    Note:
        Usign function "degree".
    */
    
    std::vector<std::vector<int>> graphDegree = gDegree(aL, maxDegree);

    /*
    --------------------------------------------------------------
    Calculating the rich-club coefficient for our graph.
    --------------------------------------------------------------
    
    Note:
        Usign function "richClub".
    */
    std::vector<double> teste = richCCalc(graphDegree, aL);
    auto t2 = std::chrono::high_resolution_clock::now();

    auto delta = std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count();
    std::cout << "Time it took " << delta << " seconds." << std::endl;
    
    // std::cout << teste.size() << std::endl;

    // for (int i=0; i<12;i++)
    // {
        
    //     std::cout << teste[i] << std::endl;
    // }
    
    // auto t2 = std::chrono::high_resolution_clock::now();

    // auto dif = std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count();
    // std::cout << "Elasped time " << dif << " seconds." << std::endl;
    // std::string nameIn = argv[1];
    // std::string nameOut = nameIn.replace(nameIn.end()-3,nameIn.end(),"rcb");    
    // std::ofstream out(nameOut);
    // for (size_t b = 0; b < teste.size(); b++)
    // {
        
    //     out << std::fixed;
    //     out << std::setprecision(5);
    //     out << teste[b] << std::endl;
        
    // }
    // out.close();
    
    

    return 0;
}