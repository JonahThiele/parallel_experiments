#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "../graph.h"
#include "../xxhash.h"

//struct local verticies
typedef struct LocalNode{
    int id;
    uint32_t color;
    int* neighbors;
    int neighbor_count;
} LocalNode;

typedef struct {
    uint32_t hash;
    int original_index;
} HashIndex;

// Determine which id the vertex
int get_owner(int vertex, int verticies, int p) 
{
    int chunk_size = verticies / p;
    int owner = vertex / chunk_size;
    return owner;
}

Graph* load_graph_from_file(const char* filename, int* total_vertices) 
{
    FILE* file = fopen(filename, "r");

    //added some nice error handling while debugging
    if (!file) 
    {
        fprintf(stderr, "Error: Could not open file %s\n", filename);
        return NULL;
    }
    
    //added some more error handling while debugging
    if (fscanf(file, "%d", total_vertices) != 1) 
    {
        fprintf(stderr, "Error: Could not read vertex count\n");
        fclose(file);
        return NULL;
    }
    
    Graph* g = generate_graph(*total_vertices);
    
    int v1, v2;
    while (fscanf(file, "%d,%d", &v1, &v2) == 2) 
    {
        add_edge(g, v1, v2);
    }
    
    fclose(file);
    return g;
}

//used with qsort
int cmp(const void* a, const void* b) {
    //make everything uint64 because otherwise truncation issues arise with the hashses
    uint32_t left = *(uint32_t*)a;
    uint32_t right = *(uint32_t*)b;
    if (left < right)
    {
        return -1;
    }
    if (left > right)
    {
        return 1;
    }

    return 0;
}

int hash_index_cmp(const void* a, const void* b) {
    uint32_t ha = ((HashIndex*)a)->hash;
    uint32_t hb = ((HashIndex*)b)->hash;
    if(ha < hb) return -1;
    if(ha > hb) return 1;
    return 0;
}

//the main program 
uint32_t* color_refinement_mpi(Graph* g, int total_verticies, int id, int p) {
    //determine how many vertices it has
    int chunk_size = 0;
    for (int i = id; i < total_verticies; i += p) 
    {
        chunk_size++;
    }
    
    LocalNode* local_verticies = malloc(chunk_size * sizeof(LocalNode));
    
    int* all_neighbor_counts = NULL; 
    //enumerate all the neighbors for every vertex
    if (id == 0) 
    {
        all_neighbor_counts = (int*)malloc(total_verticies * sizeof(int));
        for (int i = 0; i < total_verticies; i++) 
        {
            all_neighbor_counts[i] = g->list[i]->nodes;
        }
    }
    
    int* local_neighbor_counts = (int*)malloc(chunk_size * sizeof(int));
    
    // count the neighbor counts
    if (id == 0) {
        // fill in its own data
        int local_idx = 0;
        for (int i = 0; i < total_verticies; i += p) {
            local_neighbor_counts[local_idx++] = all_neighbor_counts[i];
        }
        
        // send the counts to other ranks
        for (int rank = 1; rank < p; rank++) 
        {
            int rank_chunk_size = 0;
            for (int i = rank; i < total_verticies; i += p) 
            {
                rank_chunk_size++;
            }
            
            int* temp = (int*)malloc(rank_chunk_size * sizeof(int));
            int temp_idx = 0;
            for (int i = rank; i < total_verticies; i += p) 
            {
                temp[temp_idx++] = all_neighbor_counts[i];
            }
            //sends the counts of the neighbors
            MPI_Send(temp, rank_chunk_size, MPI_INT, rank, 0, MPI_COMM_WORLD);
            free(temp);
        }
        free(all_neighbor_counts);
    } else {
        //get the number of neighbors
        MPI_Recv(local_neighbor_counts, chunk_size, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    
    // setup the local nodes on each process
    int local_idx = 0;
    for (int i = id; i < total_verticies; i += p) 
    {
        //its id is actual number in the global vertices
        local_verticies[local_idx].id = i;
        local_verticies[local_idx].color = 0;
        local_verticies[local_idx].neighbor_count = local_neighbor_counts[local_idx];
        //adjacency list type setup like the global verticies
        local_verticies[local_idx].neighbors = (int*)malloc(local_verticies[local_idx].neighbor_count * sizeof(int));
        local_idx++;
    }
    
    // send the actual global verticies number/id of the neighbors to the processes
    if (id == 0)
    {
        // don't send the neighbors of zero, just add them manually
        for (int i = 0; i < chunk_size; i++) 
        {
            int global_id = local_verticies[i].id;
            //walk along the linked list
            Node* cur = g->list[global_id]->next;
            
            for (int a = 0; a < local_verticies[i].neighbor_count; a++) 
            {
                local_verticies[i].neighbors[a] = cur->v;
                cur = cur->next;
            }
        }
        
        // give the neighbors to the other processes
        for (int rank = 1; rank < p; rank++) 
        {
            int rank_local_idx = 0;
            for (int i = rank; i < total_verticies; i += p)
            {
                //get all the neighbors of the local vertex
                int count = g->list[i]->nodes;

                //go through adajency list and send over global verticie id
                int* neighbors = (int*)malloc(count * sizeof(int));
                Node* cur = g->list[i]->next;
                for (int b = 0; b < count; b++) 
                {
                    neighbors[b] = cur->v;
                    cur = cur->next;
                }
                //give them the neighbors
                MPI_Send(neighbors, count, MPI_INT, rank, rank_local_idx, MPI_COMM_WORLD);
                free(neighbors);
                rank_local_idx++;
            }
        }
    } else {
        // receive neighbors from root
        for (int i = 0; i < chunk_size; i++) 
        {
            if (local_verticies[i].neighbor_count > 0) 
            {
                MPI_Recv(local_verticies[i].neighbors, local_verticies[i].neighbor_count, MPI_INT, 0, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
    }
    
    free(local_neighbor_counts);
    
    //counter variables for debugging
    int iteration = 0;
    int max_iterations = 100; //debug stop if not converging properly

    //is the graph coloring stabilized yet?
    int converged = 0;
    
    //thses persist between rounds
    //use malloc and then fill with zeros
    uint32_t* global_colors = malloc(total_verticies * sizeof(uint32_t));
    uint32_t* old_global_colors = malloc(total_verticies * sizeof(uint32_t));
    //calloc
    for(int i = 0; i < total_verticies; i++)
    {
        global_colors[i] = 0;
        old_global_colors[i] = 0;
    }

    //rounds to check the colors
    while (!converged && iteration < max_iterations) 
    {
        iteration++;
        
        // save old colors
        for(int i = 0; i < total_verticies; i++)
        {
            old_global_colors[i] = global_colors[i];
        }
        
        //start making the new colors
        uint32_t* new_colors = malloc(chunk_size * sizeof(uint32_t));
        
        //for each of the local verticies
        for (int i = 0; i < chunk_size; i++) 
        {
            uint32_t* neighbor_colors = malloc(local_verticies[i].neighbor_count * sizeof(uint32_t));
            
            //get them from the global for all the neighbors
            for (int a = 0; a < local_verticies[i].neighbor_count; a++) 
            {
                neighbor_colors[a] = global_colors[local_verticies[i].neighbors[a]];
            }
            
            //make sure the colors are in order smallest to largest hashs
            qsort(neighbor_colors, local_verticies[i].neighbor_count, sizeof(uint32_t), cmp);
            
            XXH32_state_t* state = XXH32_createState();
            XXH32_reset(state, 5050);
            
            //update hash with the current vertex
            uint32_t current_color = local_verticies[i].color;
            XXH32_update(state, &current_color, sizeof(current_color));
            
            for (int a = 0; a < local_verticies[i].neighbor_count; a++) 
            {
                //update the vertex hash with its neighbor colors
                XXH32_update(state, &neighbor_colors[a], sizeof(neighbor_colors[a]));
            }
            
            new_colors[i] = XXH32_digest(state);
            XXH32_freeState(state);
            free(neighbor_colors);
        }
        
        // Update local vertex colors
        for (int i = 0; i < chunk_size; i++) 
        {
            local_verticies[i].color = new_colors[i];
        }
        free(new_colors);
    
        //we want to send the data right next to each other for the id and the color,
        int* gather_sendbuf = malloc(chunk_size * 2 * sizeof(int));
        for (int i = 0; i < chunk_size; i++) 
        {
            gather_sendbuf[i * 2] = local_verticies[i].id;
            gather_sendbuf[i * 2 + 1] = (int)local_verticies[i].color;
        }
        
        // Use Allgatherv for cyclic distribution
        // https://docs.open-mpi.org/en/main/man-openmpi/man3/MPI_Allgatherv.3.html link for reference guide
        int* recvcounts = (int*)malloc(p * sizeof(int));
        int* recvdispls = (int*)malloc(p * sizeof(int));
        
        int offset = 0;
        for (int rank = 0; rank < p; rank++) 
        {
            int rank_chunk_size = 0;
            for (int i = rank; i < total_verticies; i += p) 
            {
                rank_chunk_size++;
            }
            recvcounts[rank] = rank_chunk_size * 2;
            recvdispls[rank] = offset;
            offset += rank_chunk_size * 2;
        }
        
        int* gather_recvbuf = malloc(offset * sizeof(int));
        //send updated nodes to everyone
        MPI_Allgatherv(gather_sendbuf, chunk_size * 2, MPI_INT, gather_recvbuf, recvcounts, recvdispls, MPI_INT, MPI_COMM_WORLD);
        
        // create global colors
        for (int rank = 0; rank < p; rank++) 
        {
            int disp_idx = recvdispls[rank] / 2;
            int count = recvcounts[rank] / 2;
            for (int i = 0; i < count; i++) 
            {
                int vertex_id = gather_recvbuf[(disp_idx + i) * 2];
                uint32_t color = (uint32_t)gather_recvbuf[(disp_idx + i) * 2 + 1];
                global_colors[vertex_id] = color;
            }
        }
        
        free(gather_sendbuf);
        free(gather_recvbuf);
        free(recvcounts);
        free(recvdispls);
        
        // save hashes before relabeling
        uint32_t* saved_hashes = malloc(total_verticies * sizeof(uint32_t));
        for (int a = 0; a < total_verticies; a++) 
        {
            saved_hashes[a] = global_colors[a];
        }
        
        //create mapping for index
        HashIndex* hash_array = malloc(total_verticies * sizeof(HashIndex));
        for (int a = 0; a < total_verticies; a++) 
        {
            hash_array[a].hash = global_colors[a];
            hash_array[a].original_index = a;
        }
        
        //sort the hashes
        qsort(hash_array, total_verticies, sizeof(HashIndex), hash_index_cmp);
        
        // cannonical relabel
        uint32_t* temp_colors = malloc(total_verticies * sizeof(uint32_t));
        int current_color = 0;
        temp_colors[hash_array[0].original_index] = current_color;
        
        for (int a = 1; a < total_verticies; a++) 
        {
            if (hash_array[a].hash != hash_array[a-1].hash) 
            {
                current_color++;
            }
            temp_colors[hash_array[a].original_index] = current_color;
        }
        
        //copy to global colors
        for (int a = 0; a < total_verticies; a++) 
        {
            global_colors[a] = temp_colors[a];
        }
        
        free(hash_array);
        free(temp_colors);
        
        //update the local vertex with the new colors
        for (int i = 0; i < chunk_size; i++) 
        {
            local_verticies[i].color = global_colors[local_verticies[i].id];
        }
        
        // Check if stable by comparing saved hashes, similar to pthreads version
        int local_changed = 0;
        for (int i = 0; i < chunk_size; i++) 
        {
            int global_id = local_verticies[i].id;
            
            // compute old hash
            XXH32_state_t* state = XXH32_createState();
            XXH32_reset(state, 5050);
            
            uint32_t old_color = old_global_colors[global_id];
            XXH32_update(state, &old_color, sizeof(old_color));
            
            // get neighbor colors for old hash
            uint32_t* old_neighbor_colors = malloc(local_verticies[i].neighbor_count * sizeof(uint32_t));
            for (int a = 0; a < local_verticies[i].neighbor_count; a++) 
            {
                old_neighbor_colors[a] = old_global_colors[local_verticies[i].neighbors[a]];
            }
            qsort(old_neighbor_colors, local_verticies[i].neighbor_count, sizeof(uint32_t), cmp);
            
            for (int a = 0; a < local_verticies[i].neighbor_count; a++) 
            {
                XXH32_update(state, &old_neighbor_colors[a], sizeof(old_neighbor_colors[a]));
            }
            
            uint32_t old_hash = XXH32_digest(state);
            XXH32_freeState(state);
            free(old_neighbor_colors);
            
            // are the old hashes different then the new hashes
            if (old_hash != saved_hashes[global_id]) 
            {
                local_changed = 1;
                break;
            }
        }
        
        free(saved_hashes);
        
        int global_changes = 0;
        //if converged this will return zero
        MPI_Allreduce(&local_changed, &global_changes, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        
        if (id == 0) 
        {
            printf("Round %d: %s\n", iteration, global_changes ? "partition changed" : "CONVERGED");
        }
        
        //we are done
        if (global_changes == 0) 
        {
            converged = 1;
        }
    }
    
    if (id == 0) 
    {
        if (converged) 
        {
            printf("Converged after %d round\n", iteration);
        } else {
            printf("Error didn't converge correctly\n");
        }
        
        // debug final color labels for graph
        int color_count = 0;
        uint32_t* unique_colors = malloc(total_verticies * sizeof(uint32_t));
        int* counts = malloc(total_verticies * sizeof(int));

        for(int i = 0; i < total_verticies; i++)
        {
            counts[i] = 0;
        }
        

        for (int i = 0; i < total_verticies; i++) 
        {
            int found = -1;
            for (int a = 0; a < color_count; a++) 
            {
                if (unique_colors[a] == global_colors[i]) 
                {
                    found = a;
                    break;
                }
            }
            if (found == -1) 
            {
                unique_colors[color_count] = global_colors[i];
                counts[color_count] = 1;
                color_count++;
            } else {
                counts[found]++;
            }
        }
        
        printf("colors for vertices:\n");
        for (int i = 0; i < color_count; i++) 
        {
            printf("  Color %u: %d vertex\n", unique_colors[i], counts[i]);
        }
        
        //clean up for 
        free(unique_colors);
        free(counts);
    }
    

    free(old_global_colors);

    for (int i = 0; i < chunk_size; i++) {
        free(local_verticies[i].neighbors);
    }
    free(local_verticies);
    
    return global_colors;
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    //get the filenames of the two graphs to compare
    if(argc < 3)
    {
        printf("add the two graph files to compare isomorphism\n");
        MPI_Finalize();
        return 1;
    }

    int g1_len = strlen(argv[1]);
    int g2_len = strlen(argv[2]);
    char* filename_g1 = malloc(sizeof(char) * (g1_len + 1));
    char* filename_g2 = malloc(sizeof(char) * (g2_len + 1));

    for(int i = 0; i < g1_len; i++)
    {
        filename_g1[i] = argv[1][i];
    }
    filename_g1[g1_len] = '\0';

    for(int i = 0; i < g2_len; i++)
    {
        filename_g2[i] = argv[2][i];
    }
    filename_g2[g2_len] = '\0';
    
    int id, p, len;
    char hostname[300];
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    //get the hostname like previous assignments
    MPI_Get_processor_name(hostname, &len);
    printf("Hostname: %s, rank: %d, size: %d\n", hostname, id, p);
    fflush(stdout);

    
    int GLOBAL_NODES = 0;
    
    // G1
    Graph* g1 = NULL;
    if (id == 0) 
    {
        g1 = load_graph_from_file(filename_g1, &GLOBAL_NODES);
        if (!g1) {
            fprintf(stderr, "Error loading graph from %s\n", filename_g1);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        printf("Graph 1\n");
        fflush(stdout);
    }
    
    // Broadcast the number of vertices to all processes
    MPI_Bcast(&GLOBAL_NODES, 1, MPI_INT, 0, MPI_COMM_WORLD);

    //timing for the first graph
    double start_graph1 = MPI_Wtime();
    uint32_t* colors1 = color_refinement_mpi(g1, GLOBAL_NODES, id, p);
    double end_graph1 = MPI_Wtime();
    
    // G2 - same graph for testing
    Graph* g2 = NULL;
    if (id == 0) 
    {
        g2 = load_graph_from_file(filename_g2, &GLOBAL_NODES);
        if (!g2) {
            fprintf(stderr, "Error loading graph from %s\n", filename_g2);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        printf("\nGraph 2\n");
        fflush(stdout);
    }
    
    //timing for the second graph
    double start_graph2 = MPI_Wtime();
    uint32_t* colors2 = color_refinement_mpi(g2, GLOBAL_NODES, id, p);
    double end_graph2 = MPI_Wtime();
    
    // Compare the color lists to check for isomorphism
    if (id == 0) 
    {
        int same = 1;

        for (int i = 0; i < GLOBAL_NODES; i++) 
        {
            if (colors1[i] != colors2[i]) 
            {
                same = 0;
                break;
            }
        }
        
        if (same) 
        {
            printf("They are likely isomorphic graphs\n");
            fflush(stdout);
        } else {
            printf("These two graphs can not possibly be the same graphs\n");
            fflush(stdout);
        }

        printf("The elasped time for graph1 is: %f secs\n", end_graph1 - start_graph1);
        printf("The elasped time for graph2 is: %f secs\n", end_graph2 - start_graph2);
        printf("total parallel elasped time for program: %f secs\n", (end_graph1 - start_graph1) + (end_graph2 - start_graph2));
        fflush(stdout);
        
        // Free graph
        for (int i = 0; i < GLOBAL_NODES; i++) 
        {
            Node* cur = g1->list[i];
            while (cur) 
            {
                Node* tmp = cur;
                cur = cur->next;
                free(tmp);
            }
            cur = g2->list[i];
            while (cur) 
            {
                Node* tmp = cur;
                cur = cur->next;
                free(tmp);
            }
        }
        free(g1->list);
        free(g1);
        free(g2->list);
        free(g2);
    }
    
    free(colors1);
    free(colors2);
    free(filename_g1);
    free(filename_g2);
    
    MPI_Finalize();
    return 0;
}