#include <cmath>
#include <algorithm>
#include <iostream>
#include <vector>
#include "prim.hpp"

const int INF = 1000000;
void receiveGraph(int nEdges,int nNodes,int ** Graph){
    for (int i = 0; i < nNodes; i++)
    {
        for (int j = 0; j <= i; j++)
        {
            Graph[i][j] = 0;
            Graph[j][i] = 0;
        }
        
    }
    
    int first,second,weight;
    for (int i = 0; i < nEdges; i++)
    {
        std::cout << "first end of the " << i << "-th edge : ";
        std::cin >> first;
        std::cout << "second end of the " << i << "-th edge : ";
        std::cin >> second;
        std::cout << "weight of the " << i << "-th edge : ";
        std::cin >> weight;
        Graph[first][second] = weight;
        Graph[second][first] = weight;
    }
    
}
int minConnector(int nNodes,int weights[],bool Set[]){

    int min = INF, min_index; 
    for (int i = 0; i < nNodes; i++)
    {
        if (Set[i]==false && weights[i] < min)//weights is the distance of nodes in the Set, Set is the set of all nodes already in the mst
        {
            min = weights[i];
            min_index = i;
        }
        
    }
    return min_index;
}
void prim(int nNodes,int** Graph){
    int *parent = new int[nNodes];
    int weights[nNodes];
    bool Set[nNodes];

    //initialization of these arrays:

    for (int i = 0; i < nNodes; i++)
    {
        Set[i] = false;
        weights[i] = INF;
    }


    parent[0] = -1;//indicating the very first node that we take as the root
    weights[0] = 0;
    /*for (int i = 0; i < nNodes; i++)
    {
        if (Graph[0][i] > 0)//if there is any edge at all
        {
            weights[i] = Graph[0][i];
        }
        
    }*/
    for (int i = 0; i < nNodes; i++)
    {
        parent[i] = -1;
    }
    

    for (int i = 0; i < nNodes; i++)//at each step we add a node to Set and an edge, the tree would be complete in nNodes-1 steps
    {
        int u = minConnector(nNodes,weights,Set);
        Set[u] = true;
        //now probably the weights of u-neighbours have changed, so we update all of them

        for (int j = 0; j < nNodes; j++)
        {

            if (Graph[u][j])//if there is any edge from u to i
            {
                if (Set[j] == false)//if we have not yet the node i in the Set
                {
                    if (weights[j] >= Graph[u][j])//if we can update the distance of the Set from the node i
                    {
                        weights[j] = Graph[u][j];
                        parent[j] = u;
                    }
                    
                }
                
            }            
        }
        
    }

    printGraph(nNodes,parent,Graph);
    return;
}
void printGraph(int nNodes,int *parent,int** Graph){

    
    std::cout << "parent : node : weight" << std::endl;

    for (int i = 0; i < nNodes; i++)
    {
        std::cout << parent[i] << " " << i << std::endl;
    }
    
}


int main(){
    
    int nNodes,nEdges;
    std::cout << "number of nodes? ";
    std::cin >> nNodes;
    std::cout << std::endl << "number of edges? ";
    std::cin >> nEdges;
    std::cout << std::endl;
    int** Graph;
    Graph = new int*[nNodes];
    for (int i = 0; i < nNodes; i++)
    {
        Graph[i] = new int[nNodes];
    }
    
    receiveGraph(nEdges,nNodes,Graph);
    prim(nNodes,Graph);
    for (int i = 0; i < nNodes; i++)
    {
        delete[] Graph[i];
    }
    delete[] Graph;
    return 0;
}