#include "../graph.h"
#include "../xxHash64.h"
#include <stdio.h>
#include <stdlib.h>

//using xxhash64 instead of my own homemade hash because its used in literature
//Here is where I got the hash function from https://github.com/stbrumme/xxhash


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

int main()
{
    //Tommorrow I need to figure out a dynamic way to load the graph files into the 
    //hardcode the graph initially but later use some dynamic loading from other sources
    Graph *g = generate_graph(10);
    
    /* Graph looks like this
         0
        / \
       1   2
      /|   |\
     3 4   5 6
          / \
         7---8
              \
               9
    */

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

    int new_colors[10];

    //start algorithm
    int changed = 1;
    while(changed)
    {
          changed = 0;

          //go through every node
          for(int i = 0; i < 10; i++)
          {    
               int neighbors = g->list[i]->nodes;
               int *neighbor_colors = malloc(neighbors * sizeof(int));

               Node *curr = g->list[i];
               for(int a = 0; a < neighbors; a++)
               {
                    //start off with the neighbor not the current one
                    curr = curr->next;
                    neighbor_colors[a] = curr->Color;
               }

               //sort so that different orientations don't give us incorrect or different hashes
               qsort(neighbor_colors, neighbors, sizeof(int), comp);

               //what would the hash function need to be with refinement, I should look into various ones that might be useful
               int hash = g->list[i]->Color * 31 + neighbors;
               for(int b = 0; b < neighbors; b++)
               {
                    hash = hash * 17 + neighbor_colors[b];
                    XXH64_hash_t hash = XXH64(buffer, size, seed);
               }

               new_colors[i] = hash;
               free(neighbor_colors);
          }

          int unique[10];
          int unique_count = 0;

          for(int a = 0; a < 10; a++)
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

          for(int b = 0; b < 10; b++)
          {
               if (g->list[b]->Color != new_colors[b])
               {
                    g->list[b]->Color = new_colors[b];
                    changed = 1;
               }
          }
    }

    return 0;
}