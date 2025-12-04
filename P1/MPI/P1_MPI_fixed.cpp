#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#ifdef __cplusplus
extern "C" {
#endif
#include "../graph.h"
#include "../xxhash.h"
#ifdef __cplusplus
}
#endif

typedef struct LocalNode{
    int id;
    uint64_t color;
    int* neighbors;
    int neighbor_count;
} LocalNode;

// Determine which id the vertex
int get_owner(int vertex, int verticies, int p) {
    int chunk_size = verticies / p;
    int owner = vertex / chunk_size;
    return owner;
}

Graph* create_full_graph(int total_verticies) {
    Graph* g = generate_graph(total_verticies);
    
    add_edge(g, 0, 1);
    add_edge(g, 0, 2);
    add_edge(g, 1, 3);
    add_edge(g, 1, 4);
    add_edge(g, 2, 5);
    add_edge(g, 2, 6);
    add_edge(g, 5, 7);
    add_edge(g, 6, 8);
    add_edge(g, 7, 8);
    add_edge(g, 8, 9);
    
    return g;
}

//used with qsort
int cmp(const void* a, const void* b) {
    //make everything uint64 because otherwise truncation issues arise with the hashses
    uint64_t left = *(uint64_t*)a;
    uint64_t right = *(uint64_t*)b;
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

//the main program 
uint64_t* color_refinement_mpi(Graph* g, int total_verticies, int id, int p) {
    //determine the p of the chunks of verticies
    int chunk_size = total_verticies / p;
    int start = id * chunk_size;
    
    LocalNode* local_verticies = malloc(chunk_size * sizeof(LocalNode));
    
    int* all_neighbor_counts; 
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
    //how do we handle uneven chunk sizes via the scatter
    //this tell each process what the amount of neighbors for its nodes are
    MPI_Scatter(all_neighbor_counts, chunk_size, MPI_INT, local_neighbor_counts, chunk_size, MPI_INT, 0, MPI_COMM_WORLD);
    
    // setup the initial nodes
    for (int i = 0; i < chunk_size; i++) 
    {
        //its id is actual number in the global vertices
        local_verticies[i].id = start + i;
        local_verticies[i].color = 0;
        local_verticies[i].neighbor_count = local_neighbor_counts[i];
        //adjacency list type setup like the global verticies
        local_verticies[i].neighbors = (int*)malloc(local_verticies[i].neighbor_count * sizeof(int));
    }
    
    // send the actual global verticies number/id of the neighbors to the processes
    if (id == 0)
    {
        // don't send the neighbors of zero, just add them manually
        for (int i = 0; i < chunk_size; i++) 
        {
            //walk along the linked list
            Node* cur = g->list[i]->next;
            
            for (int a = 0; a < local_verticies[i].neighbor_count; a++) 
            {
                local_verticies[i].neighbors[a] = cur->v;
                cur = cur->next;
            }
        }
        
        // Send to other processes
        for (int i = 1; i < p; i++) 
        {
            //find the global id start of the vertex for that specific process
            int chunk_start = i * chunk_size;
            for (int a = 0; a < chunk_size; a++)
            {

                int idx = chunk_start + a;
                //get all the neighbors of the local vertex
                int count = g->list[idx]->nodes;

                //this shouldn't happend unless the node is by itself which is not going to happen
                // if (count > 0) 
                // {
                int* neighbors = (int*)malloc(count * sizeof(int));
                Node* cur = g->list[idx]->next;
                for (int b = 0; b < count; b++) 
                {
                    neighbors[b] = cur->v;
                    cur = cur->next;
                }
                MPI_Send(neighbors, count, MPI_INT, i, a, MPI_COMM_WORLD);
                free(neighbors);
                //}
            }
        }
        //this is very strange I will just remove the check and keep the free
        //if (all_neighbor_counts) 
        free(all_neighbor_counts);
    } else {
        // Receive neighbors from root
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
    uint64_t* global_colors = malloc(total_verticies * sizeof(uint64_t));
    uint64_t* old_global_colors = malloc(total_verticies * sizeof(uint64_t));
    for(int i = 0; i < total_verticies; i++)
    {
        global_colors[i] = 0;
        old_global_colors[i] = 0;
    }

    //rounds to check the colors
    while (!converged && iteration < max_iterations) 
    {
        iteration++;
        
        // Save old colors
        //memcpy(old_global_colors, global_colors, total_verticies * sizeof(uint64_t));
        for(int i = 0; i < total_verticies; i++)
        {
            old_global_colors[i] = global_colors[i];
        }
        
        // Compute new colors using current global color state
        uint64_t* new_colors = malloc(chunk_size * sizeof(uint64_t));
        
        //for each of the local verticies
        for (int i = 0; i < chunk_size; i++) 
        {
            uint64_t* neighbor_colors = malloc(local_verticies[i].neighbor_count * sizeof(uint64_t));
            
            //get them from the global for all the neighbors
            for (int a = 0; a < local_verticies[i].neighbor_count; a++) 
            {
                neighbor_colors[a] = global_colors[local_verticies[i].neighbors[a]];
            }
            
            //make sure the colors are in order smallest to largest hashs
            qsort(neighbor_colors, local_verticies[i].neighbor_count, sizeof(uint64_t), cmp);
            
            XXH64_state_t* state = XXH64_createState();
            XXH64_reset(state, 5050);
            
            //update hash with the current vertex
            uint64_t current_color = local_verticies[i].color;
            XXH64_update(state, &current_color, sizeof(current_color));
            
            for (int a = 0; a < local_verticies[i].neighbor_count; a++) 
            {
                //update the vertex hash with its neighbor colors
                XXH64_update(state, &neighbor_colors[a], sizeof(neighbor_colors[a]));
            }
            
            new_colors[i] = XXH64_digest(state);
            XXH64_freeState(state);
            free(neighbor_colors);
        }
        
        // Update local vertex colors
        for (int i = 0; i < chunk_size; i++) 
        {
            local_verticies[i].color = new_colors[i];
        }
        free(new_colors);
        
        // what is actually happening right here?
        //we want to send the data right next to each other for the id and the color,
        //but it seems likes a really strange implementation
        int* gather_sendbuf = malloc(chunk_size * 3 * sizeof(int));
        for (int i = 0; i < chunk_size; i++) 
        {
            gather_sendbuf[i * 3] = local_verticies[i].id;
            //what color shifting strangeness
            gather_sendbuf[i * 3 + 1] = (int)(local_verticies[i].color & 0xFFFFFFFF);
            gather_sendbuf[i * 3 + 2] = (int)(local_verticies[i].color >> 32);
        }
        
        int* gather_recvbuf = malloc(total_verticies * 3 * sizeof(int));
        //send updated nodes to everyone
        MPI_Allgather(gather_sendbuf, chunk_size * 3, MPI_INT, gather_recvbuf, chunk_size * 3, MPI_INT, MPI_COMM_WORLD);
        
        // Build global color array from gathered data
        for (int i = 0; i < total_verticies * 3; i += 3) {
            int id = gather_recvbuf[i];
            //what the sigma is this color shifting 
            uint64_t color = ((uint64_t)gather_recvbuf[i + 2] << 32) | (uint64_t)gather_recvbuf[i + 1];
            global_colors[id] = color;
        }
        
        free(gather_sendbuf);
        free(gather_recvbuf);
        
        // relabel the colors with a smaller number so it is easier to work with
        uint64_t* unique = malloc(total_verticies * sizeof(uint64_t));
        int unique_count = 0;
        
        //
        for (int a = 0; a < total_verticies; a++) 
        {
            int found = -1;
            //check the unique colors and see if the hashed color already has a compressed color associated 
            for (int b = 0; b < unique_count; b++) 
            {
                //give it that associated color
                if (unique[b] == global_colors[a]) 
                {
                    found = b;
                    break;
                }
            }
            //give it the next number as a compressed number
            if (found == -1) 
            {
                unique[unique_count] = global_colors[a];
                found = unique_count;
                unique_count++;
            }
            global_colors[a] = found;
        }
        
        free(unique);
        
        // Update local node colors with new compressed colors
        for (int i = 0; i < chunk_size; i++) 
        {
            local_verticies[i].color = global_colors[local_verticies[i].id];
        }
        
        // Count verticies changed
        int local_changes = 0;
        for (int i = 0; i < chunk_size; i++) 
        {
            if (old_global_colors[local_verticies[i].id] != global_colors[local_verticies[i].id]) 
            {
                local_changes++;
            }
        }
        
        int global_changes = 0;
        //if converged this will return zero
        MPI_Allreduce(&local_changes, &global_changes, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        
        if (id == 0) 
        {
            printf("Round %d: %d nodes changed globally\n", iteration, global_changes);
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
        uint64_t* unique_colors = malloc(total_verticies * sizeof(uint64_t));
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
            printf("  Color %llu: %d vertex\n", unique_colors[i], counts[i]);
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
    
    int id, p;
    MPI_Comm_id(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    
    const int GLOBAL_NODES = 10;
    
    // G1
    Graph* g1 = NULL;
    if (id == 0) {
        g1 = create_full_graph(GLOBAL_NODES);
        printf("Graph 1\n");
    }
    
    uint64_t* colors1 = color_refinement_mpi(g1, GLOBAL_NODES, id, p);
    
    // G2 - same graph for testing
    Graph* g2 = NULL;
    if (id == 0) {
        g2 = create_full_graph(GLOBAL_NODES);
        printf("\nGraph 2\n");
    }
    
    uint64_t* colors2 = color_refinement_mpi(g2, GLOBAL_NODES, id, p);
    
    // Compare the color lists to check for isomorphism
    if (id == 0) {
        int same = 1;

        for (int i = 0; i < GLOBAL_NODES; i++) {
            if (colors1[i] != colors2[i]) {
                same = 0;
                break;
            }
        }
        
        if (same) {
            printf("They are likely isomorphic graphs\n");
        } else {
            printf("These two graphs can not possibly be the same graphs\n");
        }
        
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
    
    MPI_Finalize();
    return 0;
}