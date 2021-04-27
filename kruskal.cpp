#include <cmath>
#include <algorithm>
#include <iostream>
#include <vector>
#include "kruskal.hpp"


bool edgeCompare(Edge edge1,Edge edge2){
    return !(edge1.weight > edge2.weight);
}
int find(int node,int* parents){
    if (node == parents[node])
    {
        return node;
    }
    parents[node] = find(parents[node],parents);
    return parents[node];
}
void Union(int first,int second,int* parents){
    //the noides in the first and second would be actually the representers of a set
    //so we don't need to perform any find on them
    parents[first] = second;
}
void receiveGraph(Graph* Graph){
    int nNodes,nEdges;
     std::cout << "number of nodes? ";
    std::cin >> nNodes;
    std::cout << std::endl;
    std::cout << "number of edges? ";
    std::cin >> nEdges;
    Graph->edge = new Edge[Graph->nEdges];
    Graph->nNodes = nNodes;
    Graph->nEdges = nEdges;
    int first,second,weight;
    for (int i = 0; i < nEdges; i++)
    {
        std::cout << "first end of the " << i << "-th edge : ";
        std::cin >> first;
        std::cout << "second end of the " << i << "-th edge : ";
        std::cin >> second;
        std::cout << "weight of the " << i << "-th edge : ";
        std::cin >> weight;
        (Graph->edge)[i].start = first;
        (Graph->edge)[i].end = second;
        (Graph->edge)[i].weight = weight;
    }
}

int kruskal(Graph* Graph){
    //first we need to sort all the edges 

    std::sort(Graph->edge,(Graph->edge)+Graph->nEdges,edgeCompare);
    //now we have to go through all edges from the bginning of the sorted array
    //and add as much edge as we can without making loops, until getting a connected graph
    int nNodes = Graph->nNodes;
    int nEdges = Graph->nEdges;
    //int *parents = new int[nNodes];
    int parents[nNodes];
    int nSets = nNodes;
    int weightMST = 0;
    
    for (int i = 0; i < nNodes; i++)
    {
        parents[i] = i;
    }
    int iterate = 0;
    while (nSets > 1)
    {
        int start = (Graph->edge)[iterate].start;
        int end = (Graph->edge)[iterate].end;
        int set1 = find(start,parents);
        int set2 = find(end,parents);

        if (set1 == set2)
        {
            continue;
        }else{
            Union(set1,set2,parents);
            nSets--;
            weightMST += (Graph->edge)[iterate].weight;
        }
        iterate++;
    }
    return weightMST;
}


int main(){
    Graph* Graph = new struct Graph;
    receiveGraph(Graph);
    int minWeight = kruskal(Graph);
    std::cout << std::endl << "weight of the MST tree : " << minWeight << std::endl;
    return 0;
}