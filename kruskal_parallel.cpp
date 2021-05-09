/**
 * Parallel implementation of Kruskal's MST algorithm
 *
 * Compile with: mpicc kruskal_mpi.c -o kruskal_mpi -lm
 * Execute with: mpiexec -n 2 kruskal_mpi input_graph
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <mpi.h>

#include "kruskal_parallel.hpp"



void init_MPI(int argc, char** argv) {

	/* Initialize MPI */
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &number_of_processors);//initialization is the same as in the prim_parallel
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	/* Parsing graph size */
	f = fopen(argv[1], "rb");
	fread(&num_vertice, sizeof(num_vertice), 1, f);

	/* Create 'edge' type (u, v, distance) */
	MPI_Type_contiguous(3, MPI_INT, &mpi_edge);//we have already defined mpi_edge as a MPI_Datatype
	MPI_Type_commit(&mpi_edge);

	vertice_process = num_vertice / num_process;


	/*if (rank == root) {
		printf("Number of vertices: %d. Each processor will get %d vertices.\n", number_of_vertices, vertices_per_process);
	}*/

}

void check_conditions() {

	/* Check for error conditions */
	unsigned char isError = 0;
	if (num_vertice < num_process) {
		printf("Number of vertices (%d) is smaller than number of processors (%d)\n", num_vertice, num_process);
		isError = 1;
	}

	if (num_vertice % num_process != 0) {
		puts("Number of vertices % number of processors must be 0\n");
		isError = 1;
	}

	if (!((num_process != 0) && !(num_process & (num_process - 1)))) {/*this is because in each step we have 
                                                                        to merge two different process*/
		puts("Number of processors must be power of 2\n");
		isError = 1;
	}

	if (isError == 1) {
		fclose(f);
		MPI_Type_free(&mpi_edge);
		MPI_Finalize();
		exit(EXIT_FAILURE);
	}

}

void parse_input() {

	/* Start and end index of vertex belonging to this process */
	int rank_range_start = rank * vertice_process;
	int rank_range_end = (rank + 1) * vertice_process;

	num_edge = 0;

	int total_edges, u, v, edge_weight;
	fread(&total_edges, sizeof(int), 1, f);
	
	/* Create arrays for holding edge structures */
	edges = (edge_s *) malloc(vertice_process * num_vertice * sizeof(edge_s));
	local_msf_edges = (edge_s *) malloc((num_vertice - 1) * sizeof(edge_s));
	merged_msf_edges = (edge_s *) malloc(2 * (num_vertice - 1) * sizeof(edge_s));

	int i;
	int f_edge[3];
	for (i=0; i < total_edges; i++) {
		fread(&f_edge, sizeof(int), 3, f);
		
		v = f_edge[0]; u = f_edge[1]; edge_weight = f_edge[2];
		

        //if one end of the edge belongs to this process then we read it in this process
		if (((v >= rank_range_start) && (v < rank_range_end))
		|| ((u >= rank_range_start) && (u < rank_range_end))) {
			edge_s e = { v, u, edge_weight };
			edges[num_edge++] = e;
		}
	}
	
	fclose(f);

	MPI_Barrier(MPI_COMM_WORLD);
}

void finalize() {

	/* Free memory used by dynamic data structures */
	free(uf_set);
	free(edges);
	free(local_msf_edges);
	if (rank % 2 == 0) {//this is because when merging 2k and 2k-1, we keep the result in 2k
		free(recv_msf_edges);
	}
	free(merged_msf_edges);

	/* Shut down MPI */
	MPI_Type_free(&mpi_edge);
	MPI_Finalize();

	//printf("Proc %d: All done. Terminating.\n", rank);
}

int compare_edges(const void *a, const void *b) {
	edge_s *ea = (edge_s *)a;
	edge_s *eb = (edge_s *)b;

	if (ea->weight < eb->weight) {
		return -1;
	} else if (ea->weight > eb->weight) {
		return 1;
	}

	// Never happens if weights are distinct
	return 0;
}

void find_local_msf() {

	local_msf_edge_count = 0;
	merged_msf_edge_count = 0;

	
	qsort(edges, num_edge, sizeof(edge_s), compare_edges);

	uf_make();

	int used_edge_index = 0;
	int i;
	for(i = 1;i<=num_edge;i++) {

		edge_s * min_edge = &edges[used_edge_index++];//that is already sorted

		int v = (*min_edge).v;
		int u = (*min_edge).u;
		u_node * v_node = uf_find(uf_set + v);//finding the set that this node belongs to
		u_node * u_node = uf_find(uf_set + u);

		/* Add edge to MSF if it doesn't form a cycle, that they're not in the same set */
		if(v_node != u_node) {
			local_msf_edges[local_msf_edge_count++] = *min_edge;
			merged_msf_edges[merged_msf_edge_count++] = *min_edge;

			uf_union(v_node, u_node);//putting two nodes in the same set
		}
	}

}

void merge_msf() {

	max_merge_iteration = 0;

	int temp = num_process;
	while(temp > 1) {
		temp = temp / 2;
		max_merge_iteration++;
	}

	/* Skip odd ranks, as they will never receive, everyone else will receive at least once */
	if (rank % 2 == 0) {
		recv_msf_edges = (edge_s *) malloc((number_of_vertices - 1) * sizeof(edge_s));
	}

	has_sent = FALSE;
	merge_iteration = 0;
	while(merge_iteration < max_merge_iteration) {
		send_recieve_local_msf();
		merge_iteration++;
		if (has_sent) { // Process has sent it's data and can terminate now
			break;
		}
	}

}

void send_recieve_local_msf() {

	MPI_Status status;
	boolean sending = should_send();
	int merge_partner_rank = get_merge_partner_rank(sending);//based on receiving or sending
	if (sending) { // should send
		//printf("Proc %d: Sending %d edges to processor #%d\n", rank, local_msf_edge_count, merge_partner_rank);
		MPI_Send(local_msf_edges, local_msf_edge_count, mpi_edge, merge_partner_rank, 0, MPI_COMM_WORLD);
		has_sent = TRUE;
	} else { // should receive
		MPI_Recv(recv_msf_edges, num_vertice - 1, mpi_edge, merge_partner_rank, 0, MPI_COMM_WORLD, &status);
		MPI_Get_count(&status, mpi_edge, &recv_msf_edge_count);
		if (recv_msf_edge_count != MPI_UNDEFINED) {
			//printf("Proc %d: Received %d edges from processor #%d\n", rank, recv_msf_edge_count, merge_partner_rank);
			merge();
		} else {
			MPI_Abort(MPI_COMM_WORLD, 10); // Abort all
		}
	}

}

boolean should_send() {
	if ((rank / (int)pow(2, merge_iteration)) % 2 == 0) {
		return FALSE;
	} else {
		return TRUE;
	}
}

int get_merge_partner_rank(boolean should_send) {
	// Merge partner is rank +/- 2^merge_iteration
	if (should_send) {
		return rank - pow(2, merge_iteration); // whan sending the data it will become closer to the root processor
	} else {
		return rank + pow(2, merge_iteration); // receiving the data from this processor
	}
}

void merge() {

	int i;
	for (i = 0; i < recv_msf_edge_count; i++) {
		edge_s received_edge = recv_msf_edges[i];
		merged_msf_edges[merged_msf_edge_count++] = received_edge;
	}

	/* Sort local and received edges */
	qsort(merged_msf_edges, merged_msf_edge_count, sizeof(edge_s), compare_edges);

	uf_make();

	int used_edge_index = 0;
	local_msf_edge_count = 0;

	for (i = 0; i < merged_msf_edge_count; i++) {

		edge_s * min_edge = &merged_msf_edges[used_edge_index++];

		int v = (*min_edge).v;
		int u = (*min_edge).u;
		u_node * v_node = uf_find(uf_set + v);
		u_node * u_node = uf_find(uf_set + u);

		if(v_node != u_node) {
			local_msf_edges[local_msf_edge_count++] = *min_edge;

			uf_union(v_node, u_node);
		}

	}

	/* Transfer new local MSF edges to merged edges */
	merged_msf_edge_count = 0;
	for (i = 0; i < local_msf_edge_count; i++) {
		edge_s merged_edge = local_msf_edges[i];
		merged_msf_edges[merged_msf_edge_count++] = merged_edge;
	}

}

void uf_make() {
	int size = sizeof(u_node) * num_vertice;// + sizeof(int)*(number_of_vertices - 1);
	uf_set = malloc(size);
	memset(uf_set, 0, size);
}

u_node * uf_find(u_node * a) {
    if (a->parent==NULL) return a;
    else return (a->parent = uf_find(a->parent));
}

void uf_union(u_node * a, u_node * b) {
    if (a->depth > b->depth) {
        b->parent = a;
    } else if (a->depth<b->depth) {
        a->parent = b;
    } else {
        a->parent = b;
        a->depth += 1;
    }
}


int main(int argc, char** argv) {

	if (argc < 2) {
		printf("usage: %s filename\n", argv[0]);
		return 1;
	}

	init_MPI(argc, argv);

	check_conditions();

	parse_input();

	find_local_msf();
	merge_msf();

	if (rank == root) {
		print_merged_msf_edges();
	}

	finalize();

	return 0;
}

void print_local_edges() {
	int i;
	printf("Proc %d local edges:\n", rank);
	for (i = 0; i < number_of_edges; ++i) {
		printf("(%d,%d) = %d\n", edges[i].v, edges[i].u, edges[i].weight);
	}
}

void print_local_msf_edges() {
	int i;
	printf("Proc %d local MSF edges:\n", rank);
	for(i = 0; i < local_msf_edge_count; i++) {
		printf("(%d,%d) = %d\n",local_msf_edges[i].v, local_msf_edges[i].u, local_msf_edges[i].weight);
	}
}

void print_received_msf_edges() {
	int i;
	printf("Proc %d received MSF edges:\n", rank);
	for(i = 0; i < recv_msf_edge_count; i++) {
		printf("(%d,%d) = %d\n",recv_msf_edges[i].v, recv_msf_edges[i].u, recv_msf_edges[i].weight);
	}
}

void print_merged_msf_edges() {
	int i;
	printf("Proc %d merged MSF edges:\n", rank);
	for(i = 0; i < merged_msf_edge_count; i++) {
		printf("(%d,%d) = %d\n",merged_msf_edges[i].v, merged_msf_edges[i].u, merged_msf_edges[i].weight);
	}
}


