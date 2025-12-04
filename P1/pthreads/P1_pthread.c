#include "../graph.h"
#include "../xxhash.h"
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#define THREADS 4
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
} threadArgs;

//so we can use acutal graphs instead of hardcoded ones 
Graph* load_graph_from_file(const char* filename, int* total_vertices) {
    FILE* file = fopen(filename, "r");
    if (!file) {
        fprintf(stderr, "Error: Could not open file %s\n", filename);
        return NULL;
    }
    
    if (fscanf(file, "%d", total_vertices) != 1) {
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

//used to sort the color list because I didn't want to implement my own sort
//https://stackoverflow.com/questions/1787996/c-library-function-to-perform-sort
int comp (const void * el, const void * el2)
{
     int left = *((int*)el);
     int right = *((int*)el2);
     if(left < right)
     {
          return -1;
     } else if(left > right) {
          return 1;
     } else {
          return 0;
     }

}

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
              pthread_mutex_unlock(t_args->changed_mutex);
          }
          
          pthread_barrier_wait(t_args->barrier);

          //clyclic allocation
          for(int i = t_args->id; i < t_args->g_size; i += t_args->t_total)
          {  
            //we don't need mutexes because they are only reading these
            int neighbors = g->list[i]->nodes;
            int *neighbor_colors = malloc(neighbors * sizeof(int));

            Node *curr = g->list[i]->next;
            for(int a = 0; a < neighbors; a++)
            {
                //start off with the neighbor not the current one
                neighbor_colors[a] = curr->Color;
                curr = curr->next;
                    
            }

            //sort so that different orientations don't give us incorrect or different hashes
            qsort(neighbor_colors, neighbors, sizeof(int), comp);
            //api that I used for the hash
            //https://xxhash.com/doc/v0.8.3/group___x_x_h64__family.html#gadcfc7fb0e3382dc0b13afb125a0fea5c
               
            XXH64_state_t* state = XXH64_createState();
            XXH64_reset(state, 5050);
            //what would the hash function need to be with refinement, I should look into various ones that might be useful
            for(int b = 0; b < neighbors; b++)
            {
                XXH64_update(state, &neighbor_colors[b], sizeof(neighbor_colors[b]));
            }  

            //create the hash
            uint64_t hash = XXH64_digest(state);
            XXH64_freeState(state); 	

            pthread_mutex_lock(t_args->changed_mutex);
            t_args->new_colors[i] = hash;
            pthread_mutex_unlock(t_args->changed_mutex);
            
            free(neighbor_colors);
          }

          uint64_t unique[t_args->g_size];
          int unique_count = 0;

          //simplify the colors to make it easier to handle
          pthread_barrier_wait(t_args->barrier);
          //the cannonical color relabeling has to make the new and unique color list into shared memory
          //More cyclic allocation
          for(int a = t_args->id; a < t_args->g_size; a += t_args->t_total)
          {
               int found = -1;
               for(int j = 0; j < unique_count; j++)
               {
                    if(unique[j] == t_args->new_colors[a])
                    {
                        found = j;
                        break;
                    }
               }
               if(found == -1)
               {
                    unique[unique_count] = t_args->new_colors[a];
                    found = unique_count;
                    unique_count++;
               }
               //guess we only need one mutex because of the barriers, so we will just reuse it
               pthread_mutex_lock(t_args->changed_mutex);
               t_args->new_colors[a] = found;
               pthread_mutex_unlock(t_args->changed_mutex);
          }

          //check if the colors have stabilized/not changed since the last round
          pthread_barrier_wait(t_args->barrier);
          //more cyclic allocation
          for(int b = t_args->id; b < t_args->g_size; b += t_args->t_total)
          {
               if (g->list[b]->Color != t_args->new_colors[b])
               {
                    //update the actual graph colors
                    pthread_mutex_lock(t_args->changed_mutex);
                    //should I wrap these in two separate mutexes
                    g->list[b]->Color = t_args->new_colors[b];
                    *(t_args->changed) = 1;
                    pthread_mutex_unlock(t_args->changed_mutex);
               }
          }
          
          //synchronize to read the changed flag
          pthread_barrier_wait(t_args->barrier);
          local_changed = *(t_args->changed);
    }
    //should be in shared memory so we don't need to return the color list like the sequential version
    pthread_exit(NULL);
}



uint64_t * color_refinement(Graph *g, int size,  int threads)
{
    //initialize the threads
    pthread_t threadlist[threads];
    pthread_mutex_t changed_mutex = PTHREAD_MUTEX_INITIALIZER;
    pthread_barrier_t barrier;

    pthread_barrier_init(&barrier, NULL, threads);
    uint64_t *new_colors = malloc(sizeof(uint64_t) * size);
    int changed = 1;
    
    //allocate thread args in heap so they persist
    threadArgs *args_array = malloc(sizeof(threadArgs) * threads);

    //give each thread a couple nodes
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
        pthread_create(&threadlist[i], NULL, color_refinement_partitioned, (void *)&args_array[i]);
    }

    //join them right after they will have internal barriers to handle the different for loops
    for(int i = 0; i < threads; i++)
    {
        pthread_join(threadlist[i], NULL);
    }

    pthread_barrier_destroy(&barrier);
    pthread_mutex_destroy(&changed_mutex);
    free(args_array);

    return new_colors;
}

int main(int argc, char*argv[])
{
    
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

    //G2 because we need to compare to the other one for it to work
    int g2_nodes = 0;
    Graph *g2 = load_graph_from_file(filename_g2, &g2_nodes);

    //link to reference I used to write the timing code
    //https://thelinuxcode.com/clock-gettime-c-function/

    //timing for the first graph
    struct timespec start_g1;

    clock_gettime(CLOCK_MONOTONIC, &start_g1);
    uint64_t *colors1 = color_refinement(g1, g1_nodes, threads);

    struct timespec current_g1;  

    clock_gettime(CLOCK_MONOTONIC, &current_g1);
    
    
    //timing for the secong graph
    struct timespec start_g2;

    clock_gettime(CLOCK_MONOTONIC, &start_g2);
    uint64_t *colors2 = color_refinement(g2, g2_nodes, threads);

    struct timespec current_g2;  

    clock_gettime(CLOCK_MONOTONIC, &current_g2);

    int same = 1;
    for(int i = 0; i < 10; i++)
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

    //printing out elasped times
    double elapsed_g1 = current_g1.tv_sec - start_g1.tv_sec;
    elapsed_g1 += (current_g1.tv_nsec - start_g1.tv_nsec) / 1000000000.0;
    
    printf("elasped time for first graph is: %f secs\n", elapsed_g1);

    double elapsed_g2 = current_g2.tv_sec - start_g2.tv_sec;
    elapsed_g2 += (current_g2.tv_nsec - start_g2.tv_nsec) / 1000000000.0;

    printf("elasped time for the second graph is: %f secs\n", elapsed_g2);

    printf("total elasped time for both graphs is: %f secs\n", elapsed_g2 + elapsed_g1);

    free(colors1);
    free(colors2);
    //should also free graphs here

    return 0;
}