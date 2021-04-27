#pragma once

void receiveGraph(int nEdges,int nNodes,int** Graph);
int minConnector(int nNodes,int weights[],bool Set[]);
void prim(int nNodes,int** Graph);
void printGraph(int nNodes,int parent[],int** Graph);