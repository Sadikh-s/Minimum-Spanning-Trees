#pragma once

struct Edge
{
    int start,end,weight;
};


struct Graph
{
    int nNodes,nEdges;
    Edge* edge;
};

struct subset
{
    int parent,rank;
};

int find(struct subset* subsets, int i);
void Union(struct subset* subsets, int x, int y);
void boruvka(Graph* Graph);
void receiveGraph(Graph* Graph);