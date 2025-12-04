#include "../graph.h"
#include "../xxhash.h"
#include <stdio.h>
#include <time.h>
#include <stdlib.h>

//using xxhash64 instead of my own homemade hash because its used in literature
//Here is where I got the hash function from https://github.com/stbrumme/xxhash

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



uint64_t * color_refinement(Graph *g, int size)
{
    uint64_t *new_colors = malloc(sizeof(uint64_t) * size);

    //start algorithm
    int changed = 1;
    while(changed)
    {
          changed = 0;

          //go through every node
          for(int i = 0; i < size; i++)
          {    
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
                    XXH64_update(state, &neighbor_colors[b], sizeof(&neighbor_colors[b]));
               }

               //create the hash
               uint64_t hash = XXH64_digest(state);
               XXH64_freeState(state); 	

               new_colors[i] = hash;
               free(neighbor_colors);
          }

          uint64_t unique[size];
          int unique_count = 0;

          for(int a = 0; a < size; a++)
          {
               int found = -1;
               for(int j = 0; j < unique_count; j++)
               {
                    if(unique[j] == new_colors[a])
                    {
                         found = j;
                         break;
                    }
               }
               if(found == -1)
               {
                    unique[unique_count] = new_colors[a];
                    found = unique_count++;
               }
               new_colors[a] = found;
          }

          for(int b = 0; b < size; b++)
          {
               if (g->list[b]->Color != new_colors[b])
               {
                    g->list[b]->Color = new_colors[b];
                    changed = 1;
               }
          }
    }
    return new_colors;

}

int main(int argc, char **argv)
{
    if(argc < 3)
    {
        printf("Enter the file of the graphs\n");
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

    int threads = atoi(argv[2]);

    //G1
    int g1_nodes = 0;
    Graph *g1 = load_graph_from_file(filename_g1, &g1_nodes);

    //G2 because we need to compare to the other one for it to work
    int g2_nodes = 0;
    Graph *g2 = load_graph_from_file(filename_g2, &g2_nodes);

    
    //timing for graph 1
    struct timespec start_g1;

    clock_gettime(CLOCK_MONOTONIC, &start_g1);
    uint64_t *colors1 = color_refinement(g1, g1_nodes);
    
    struct timespec current_g1;  

    clock_gettime(CLOCK_MONOTONIC, &current_g1);
    
    //timing for graph 2
    struct timespec start_g2;

    clock_gettime(CLOCK_MONOTONIC, &start_g2);
    uint64_t *colors2 = color_refinement(g2, g2_nodes);

    struct timespec current_g2;  

    clock_gettime(CLOCK_MONOTONIC, &current_g2);

    int same = 1;
    for(int i = 0; i < g1_nodes; i++)
    {
          //its not possible for them to be the same graph
          if(colors1[i] != colors2[i])
          {
               break;
          }
    }

    if(same)
    {
      printf("They are likely isomorphic graphs\n");
    }else{
     printf("Thess two graphs can not possibly be the same graphs\n");
    }
    
    //printing out elasped times
    double elapsed_g1 = current_g1.tv_sec - start_g1.tv_sec;
    elapsed_g1 += (current_g1.tv_nsec - start_g1.tv_nsec) / 1000000000.0;
    
    printf("elasped time for first graph is: %f secs\n", elapsed_g1);

    double elapsed_g2 = current_g2.tv_sec - start_g2.tv_sec;
    elapsed_g2 += (current_g2.tv_nsec - start_g2.tv_nsec) / 1000000000.0;

    printf("elasped time for the second graph is: %f secs\n", elapsed_g2);

    printf("total elasped time for both graphs is: %f secs\n", elapsed_g2 + elapsed_g1);

    return 0;
}