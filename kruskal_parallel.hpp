#pragma once

#include <mpi.h>

typedef enum { FALSE, TRUE } boolean;

void parse_input();
void parse_edge_list_input();
int compare_edges(const void *a, const void *b);
boolean should_send();
int get_merge_partner_rank(boolean);
void merge();
void send_recieve_local_msf();
void merge_msf();
void print_received_msf_edges();
void start_sort_measure();
void end_sort_measure();
void start_communication_measure();
void end_communication_measure();

int num_process;  /* Number of processors */
int rank;                  /* MPI rank */
int root = 0;              /* Rank of root process (which holds the final result) */
MPI_Datatype mpi_edge;     /* Datatype containing tuple (u,v,weight) */

FILE * f;                  /* Input file containing graph. */

int num_vertice;    /* Number of vertices in the graph */
int vertice_process;  /* Number of vertices that each process will get */
int num_edge;       /* Number of edges that each processor will get */

int merge_iterate;       /* Current merge iteration */
int max_merge_iterate;   /* Max merge iteration (if there is 2^n processors
                              there will be n merge iterations) */
boolean has_sent;          /* Once processor sends it's data it can terminate */

typedef struct edge_t {    /* struct representing one edge */
    int v;
    int u;
    int weight;
} edge_s;

edge_s * edges;            /* Array of edges belonging to this process */
edge_s * local_msf_edges;  /* Array of edges which form local MSF *///min spanning forest
edge_s * recv_msf_edges;   /* Array of edges received from other process */
edge_s * merged_msf_edges; /* Array of merged local and received edges */
int local_msf_edge_count;  /* Number of edges in array of local MSF edges */
int recv_msf_edge_count;   /* Number of edges in array of received edges */
int merged_msf_edge_count; /* Number of edges in array of merged edges */

typedef struct node_t {    /* Union-Find data structure *///to perform kruskal in each processor
    struct node_t * parent;
    int depth;
} u_node;

u_node * uf_set;           /* Array indicating which set does vertex belong to
                              (used in Union-Find algorithm) */

void uf_make();
u_node * uf_find(u_node * a);
void uf_union(u_node * a, u_node * b);