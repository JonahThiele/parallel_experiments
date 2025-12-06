#include "../graph.h"
#include "../xxhash.h"
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
//using xxhash64 instead of my own homemade hash because its used in literature
//Here is where I got the hash function from https://github.com/stbrumme/xxhash

//struct for giving args to void pthreads function

typedef struct {
    Graph *g;
    int g_size;
    int t_total;
    int id;
    //share the mutex
    pthread_mutex_t *changed_mutex;
    //share the barrier to be used for the separate parts
    pthread_barrier_t *barrier;
    //share the color list as well
    uint64_t *new_colors;
    //share the changed flag
    int *changed;
    //share the iteration counter
    int *iteration;
    //share the changes counter
    int *global_changes;
} threadArgs;

//so we can use acutal graphs instead of hardcoded ones 
Graph* load_graph_from_file(const char* filename, int* total_vertices) {
    FILE* file = fopen(filename, "r");
    
    Graph* g = generate_graph(*total_vertices);
    
    int v1, v2;
    while (fscanf(file, "%d,%d", &v1, &v2) == 2) 
    {
        add_edge(g, v1, v2);
    }
    
    fclose(file);
    return g;
}

//used to sort the color list because I didn't want to implement my own sort
//https://stackoverflow.com/questions/1787996/c-library-function-to-perform-sort
int comp (const void * el, const void * el2)
{
     uint64_t left = *((uint64_t*)el);
     uint64_t right = *((uint64_t*)el2);
     if(left < right)
     {
          return -1;
     } else if(left > right) {
          return 1;
     } else {
          return 0;
     }
}

//for cannoical graph relabeling
typedef struct {
    uint64_t hash;
    int original_index;
} HashIndex;

int hash_index_comp(const void* a, const void* b) {
    uint64_t ha = ((HashIndex*)a)->hash;
    uint64_t hb = ((HashIndex*)b)->hash;
    if(ha < hb) 
    {
        return -1;
    }
    if(ha > hb)
    {
        return 1;
    }
    return 0;
}

//given to each threads
void *color_refinement_partitioned(void * args)
{
    threadArgs* t_args = (threadArgs*)args;
    Graph *g = t_args->g;
    int local_changed = 1;
    while(local_changed)
    {
          local_changed = 0;
          
          //reset global changed flag at the start of each iteration
          if(t_args->id == 0)
          {
              pthread_mutex_lock(t_args->changed_mutex);
              *(t_args->changed) = 0;
              *(t_args->global_changes) = 0;
              (*(t_args->iteration))++;
              pthread_mutex_unlock(t_args->changed_mutex);
          }
          
          pthread_barrier_wait(t_args->barrier);

          //cyclic allocation
          for(int i = t_args->id; i < t_args->g_size; i += t_args->t_total)
          {  
            //we don't need mutexes because they are only reading these
            int neighbors = g->list[i]->nodes;
            uint64_t *neighbor_colors = malloc(neighbors * sizeof(uint64_t));

            Node *curr = g->list[i]->next;
            for(int a = 0; a < neighbors; a++)
            {
                //start off with the neighbor not the current one
                neighbor_colors[a] = curr->Color;
                curr = curr->next;
            }

            //sort so that different orientations don't give us incorrect or different hashes
            qsort(neighbor_colors, neighbors, sizeof(uint64_t), comp);
            //api that I used for the hash
            //https://xxhash.com/doc/v0.8.3/group___x_x_h64__family.html#gadcfc7fb0e3382dc0b13afb125a0fea5c
               
            XXH64_state_t* state = XXH64_createState();
            XXH64_reset(state, 5050);
            
            //hash the current vertex color first (like MPI version)
            uint64_t current_color = g->list[i]->Color;
            XXH64_update(state, &current_color, sizeof(current_color));
            
            //hash all neighbor colors
            for(int b = 0; b < neighbors; b++)
            {
                XXH64_update(state, &neighbor_colors[b], sizeof(neighbor_colors[b]));
            }  

            //create the hash
            uint64_t hash = XXH64_digest(state);
            XXH64_freeState(state); 	

            t_args->new_colors[i] = hash;
            
            free(neighbor_colors);
          }

          pthread_barrier_wait(t_args->barrier);
          
          //Use a single-threaded relabeling to ensure consistency
          if(t_args->id == 0)
          {
              
              //SAVE the hashes before relabeling so we can check convergence
              uint64_t* saved_hashes = malloc(t_args->g_size * sizeof(uint64_t));
              for(int a = 0; a < t_args->g_size; a++)
              {
                  saved_hashes[a] = t_args->new_colors[a];
              }
              
              //Create a mapping of hash -> color_id
              HashIndex* hash_array = malloc(t_args->g_size * sizeof(HashIndex));
              for(int a = 0; a < t_args->g_size; a++)
              {
                  hash_array[a].hash = t_args->new_colors[a];
                  hash_array[a].original_index = a;
              }
              
              //Sort by hash value
              qsort(hash_array, t_args->g_size, sizeof(HashIndex), hash_index_comp);
              
              //Assign new color IDs: same hash gets same ID
              uint64_t* temp_colors = malloc(t_args->g_size * sizeof(uint64_t));
              int current_color = 0;
              temp_colors[hash_array[0].original_index] = current_color;
              
              for(int a = 1; a < t_args->g_size; a++)
              {
                  if(hash_array[a].hash != hash_array[a-1].hash)
                  {
                      current_color++;
                  }
                  temp_colors[hash_array[a].original_index] = current_color;
              }
              
              //update the colors 
              for(int a = 0; a < t_args->g_size; a++)
              {
                  t_args->new_colors[a] = temp_colors[a];
              }
              
              //check for stabilized colors by comparing the hash values
              *(t_args->changed) = 0;  
              for(int i = 0; i < t_args->g_size; i++)
              {
                  //Compute hash from most recent state of 
                  XXH64_state_t* state = XXH64_createState();
                  XXH64_reset(state, 5050);
                  
                  uint64_t curr_color = g->list[i]->Color;
                  XXH64_update(state, &curr_color, sizeof(curr_color));
                  
                  int neighbors = g->list[i]->nodes;
                  uint64_t *neighbor_colors = malloc(neighbors * sizeof(uint64_t));
                  Node *curr = g->list[i]->next;
                  for(int a = 0; a < neighbors; a++)
                  {
                      neighbor_colors[a] = curr->Color;
                      curr = curr->next;
                  }
                  qsort(neighbor_colors, neighbors, sizeof(uint64_t), comp);
                  
                  for(int b = 0; b < neighbors; b++)
                  {
                      XXH64_update(state, &neighbor_colors[b], sizeof(neighbor_colors[b]));
                  }
                  
                  uint64_t current_hash = XXH64_digest(state);
                  XXH64_freeState(state);
                  free(neighbor_colors);
                  
                  //If this hash differs from the saved hash, it didn't converge 
                  if(current_hash != saved_hashes[i])
                  {
                      *(t_args->changed) = 1;
                      break;  //At least one vertex changed
                  }
              }
              
              free(hash_array);
              free(temp_colors);
              free(saved_hashes);
          }

          pthread_barrier_wait(t_args->barrier);
          
          //Update the colors for the graph structure
          for(int b = t_args->id; b < t_args->g_size; b += t_args->t_total)
          {
               g->list[b]->Color = t_args->new_colors[b];
          }
          
          *(t_args->global_changes) = *(t_args->changed) ? 1 : 0;
          
          //update changed together
          pthread_barrier_wait(t_args->barrier);
          local_changed = *(t_args->changed);
    }
    //should be in shared memory so we don't need to return the color list like the sequential version
    pthread_exit(NULL);
}


//the pthread wrapper setting up everything
uint64_t * color_refinement(Graph *g, int size,  int threads)
{
    //initialize all colors to 0
    for(int i = 0; i < size; i++)
    {
        g->list[i]->Color = 0;
    }
    
    //initialize the threads and the constructs they need 
    pthread_t threadlist[threads];
    pthread_mutex_t changed_mutex = PTHREAD_MUTEX_INITIALIZER;
    pthread_barrier_t barrier;

    pthread_barrier_init(&barrier, NULL, threads);
    uint64_t *new_colors = malloc(sizeof(uint64_t) * size);
    int changed = 1;
    int iteration = 0;
    int global_changes = 0;
    
    threadArgs *args_array = malloc(sizeof(threadArgs) * threads);

    //give each thread a couple vertex in classic cyclic fashion
    for(int i = 0; i < threads; i++)
    {
        args_array[i].g = g;
        args_array[i].id = i;
        args_array[i].g_size = size;
        args_array[i].t_total = threads;
        args_array[i].changed_mutex = &changed_mutex;
        args_array[i].barrier = &barrier;
        args_array[i].new_colors = new_colors;
        args_array[i].changed = &changed;
        args_array[i].iteration = &iteration;
        args_array[i].global_changes = &global_changes;
        pthread_create(&threadlist[i], NULL, color_refinement_partitioned, (void *)&args_array[i]);
    }

    //join them right after they will have internal barriers to handle the different for loops
    for(int i = 0; i < threads; i++)
    {
        pthread_join(threadlist[i], NULL);
    }

    printf("Converged after %d rounds\n", iteration);
    
    // debug final color labels for graph
    int color_count = 0;
    uint64_t* unique_colors = malloc(size * sizeof(uint64_t));
    int* counts = malloc(size * sizeof(int));

    //again should have used calloc
    for(int i = 0; i < size; i++)
    {
        counts[i] = 0;
    }
    
    //generate a frequency list of the the relabeled colors
    for (int i = 0; i < size; i++) 
    {
        int found = -1;
        for (int a = 0; a < color_count; a++) 
        {
            if (unique_colors[a] == g->list[i]->Color) 
            {
                found = a;
                break;
            }
        }
        if (found == -1) 
        {
            unique_colors[color_count] = g->list[i]->Color;
            counts[color_count] = 1;
            color_count++;
        } else {
            counts[found]++;
        }
    }
    
    printf("colors for vertices:\n");
    for (int i = 0; i < color_count; i++) 
    {
        printf("  Color %lu: %d vertices\n", unique_colors[i], counts[i]);
    }
    
    //clean up for 
    free(unique_colors);
    free(counts);

    pthread_barrier_destroy(&barrier);
    pthread_mutex_destroy(&changed_mutex);
    free(args_array);

    return new_colors;
}

int main(int argc, char*argv[])
{
    //load the files and the number of threads as arguments
    if(argc < 4)
    {
        printf("Enter the file of the graphs and amount of threads\n");
        return -1;
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

    int threads = atoi(argv[3]);

    //G1
    int g1_nodes = 0;
    Graph *g1 = load_graph_from_file(filename_g1, &g1_nodes);
    
    if (!g1) {
        fprintf(stderr, "Error loading graph from %s\n", filename_g1);
        return 1;
    }
    printf("Graph 1\n");
    fflush(stdout);

    //G2 because we need to compare to the other one for it to work
    int g2_nodes = 0;
    Graph *g2 = load_graph_from_file(filename_g2, &g2_nodes);
    
    if (!g2) {
        fprintf(stderr, "Error loading graph from %s\n", filename_g2);
        return 1;
    }

    //link to reference I used to write the timing code
    //https://thelinuxcode.com/clock-gettime-c-function/

    //timing for the first graph
    struct timespec start_g1;

    clock_gettime(CLOCK_MONOTONIC, &start_g1);
    uint64_t *colors1 = color_refinement(g1, g1_nodes, threads);

    struct timespec current_g1;  

    clock_gettime(CLOCK_MONOTONIC, &current_g1);
    
    printf("\nGraph 2\n");
    fflush(stdout);
    
    //timing for the second graph
    struct timespec start_g2;

    clock_gettime(CLOCK_MONOTONIC, &start_g2);
    uint64_t *colors2 = color_refinement(g2, g2_nodes, threads);

    struct timespec current_g2;  

    clock_gettime(CLOCK_MONOTONIC, &current_g2);

    int same = 1;
    for(int i = 0; i < g1_nodes; i++)
    {
          //its not possible for them to be the same graph
          if(colors1[i] != colors2[i])
          {
               same = 0;
               break;
          }
    }

    if(same)
    {
      printf("They are likely isomorphic graphs\n");
    }else{
     printf("These two graphs can not possibly be the same graphs\n");
    }

    //printing out elapsed times
    double elapsed_g1 = current_g1.tv_sec - start_g1.tv_sec;
    elapsed_g1 += (current_g1.tv_nsec - start_g1.tv_nsec) / 1000000000.0;
    
    printf("elapsed time for first graph is: %f secs\n", elapsed_g1);

    double elapsed_g2 = current_g2.tv_sec - start_g2.tv_sec;
    elapsed_g2 += (current_g2.tv_nsec - start_g2.tv_nsec) / 1000000000.0;

    printf("elapsed time for the second graph is: %f secs\n", elapsed_g2);

    printf("total elapsed time for both graphs is: %f secs\n", elapsed_g2 + elapsed_g1);

    free(colors1);
    free(colors2);
    free(filename_g1);
    free(filename_g2);
    //should also free graphs here

    return 0;
}