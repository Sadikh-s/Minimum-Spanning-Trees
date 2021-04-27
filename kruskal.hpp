#pragma once

struct Edge
{
    int start,end;
    int weight;
};

struct Graph
{
    int nNodes,nEdges;
    Edge* edge;
};


int kruskal(Graph* Graph);
void receiveGraph(Graph* Graph);
void Union(int first,int second,int* parents);
int find(int node,int* parents);
bool edgeCompare(Edge edge1,Edge edge2);