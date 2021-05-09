#pragma once
#include <mpi.h>

void start_communication_measure();
void end_communication_measure();
void reduce_global_min(int *, int *, int *, MPI_Datatype *);
int num_process; /* Number of processors */
int rank;                 /* MPI rank */
int root = 0;             /* root process holding the final result */
MPI_Datatype edge;        /* Datatype containing tuple (u,v,d[u]) */
MPI_Op reduce_op;         /* Reduce operation for findimg min edge, when gathering the mins from different processors */
int num_vertice;   /* Number of vertices in the graph, this should be a multiple of number of vertices in each processor */
int num_vertice_process; /* Number of vertices that each process will get*/
int starting_vertex;      /* Initial vertice in MST */
int distance;       /* Total distance of MST (relevant only in root process) */

FILE * f;                 /* Input file containing graph. */

int ** weight;            /* Part of weighted adjacency matrix each process holds */
int * d;                  /* Array holding distances between MST and other vertices */
int * vertice_d;          /* vertice_d[i] holds the index of the vertex that i should be linked to in order to have the distance d[i] */
unsigned char * in_MST;  /* 1 if the vertex is already in the MST; 0 otherwise*/

typedef struct edge_t {   /* struct representing one edge */
    int v;
    int u;
    int weight;
} edge_s;

edge_s * MST;       /* Array of MST edges */
int edges_MST;   /* Number of edges currently added to MST */

double comp_start_time;   /* Start time of computation */
double comp_end_time;     /* End time of computation */

double comm_start_time;   /* Start time of communication */
double comm_end_time;     /* End time of communication */
double total_comm_time;   /* Total time of communication */

double parse_start_time;  /* Start time of parsing */
double parse_end_time;    /* End time of parsing */