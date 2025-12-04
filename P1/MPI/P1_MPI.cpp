#include "../graph.h"
#include "../xxhash.h"
#include <stdio.h>
//some sorting stuff required to make the many to many work
//without my own bad implementation of linked lists
#include <vector>
#include <algorithm>
#include <map>

#include <stdlib.h>
#include <mpi.h>

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

typedef struct otherNodes {
     int amount;
     std::vector<int> nodes;
     int owner;
} otherNodes;

typedef struct colorRequests {
     int vector;
     int color;
} colorRequests;

int main(int argc, char **argv)
{
    
    //G1
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

    //G2 because we need to compare to the other one for it to work
    Graph *g2 = generate_graph(10);

    add_edge(g2, 0, 1);
    add_edge(g2, 0, 2);
    add_edge(g2, 1, 3);
    add_edge(g2, 1, 4);
    add_edge(g2, 2, 5);
    add_edge(g2, 2, 6);
    add_edge(g2, 5, 7);
    add_edge(g2, 6, 8);
    add_edge(g2, 7, 8);
    add_edge(g2, 8, 9);

    int g_size = 10;

    //Setup for standard MPI things
    MPI_Init(&argc, &argv);
    int id, p;
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    //local nodes and neighbors
    int *amount_neighbors;
    int *local_nodes;
    int *local_colors;
    int **local_neighbors;
    int range, chunk_start, chunk_end;

//     All processes:
//     Load graph partition (each process owns some subset of nodes)
//     For each local node v:
//         neighbors[v] = list of all neighbors (local or remote)
//         Initialize color[v] = 1

     //worker dividing up the nodes
     if(id == 0)
     {    
          //for each processors determine the range of nodes
          //bundle up the processor data and then send it
          for(int i = 1; i < p; i++)
          {
               int start = i * (g_size/p);
               int end = (i + 1) * (g_size/p);
               if (end > g_size)
               {
                    end = g_size;
               }
               int range[2] = {start, end};
               MPI_Send(range, 2, MPI_INT, i, 0, MPI_COMM_WORLD);

               //send adjacency list so that processors have 
               for(int a = start; a < end; a++)
               {
                    int *neighbors = malloc(sizeof(int) * g->list[a]->nodes);
                    Node *curr = g->list[a]->next;
                    int k = 0;
                    while(curr != NULL)
                    {
                         neighbors[k] = curr->v;
                         k++;
                         curr = curr->next;
                    }
                    //convert from a linked list to just an array of ints to send
                    MPI_Send(neighbors, g->list[a]->nodes, MPI_INT, i, a, MPI_COMM_WORLD);
                    free(neighbors);
               }

          }
     //get from the worker
     } else {
          int chunk[2];
          int start, end;
          //get range
          MPI_Recv(chunk, 2, MPI_INT, 0, 0, MPI_COMM_WORLD);
          start = chunk[0];
          end = chunk[1];
          chunk_end = end;
          chunk_start = start;
          range = end - start;

          //fill with range of nodes
          local_nodes = malloc(sizeof(int) * range);
          local_colors = malloc(sizeof(int) * range);
          for(int i = 0; i < range; i++)
          {
               local_nodes[i] = start + i;
               //initial color is no color
               local_colors[i] = 0;
          }

          //fill up the neighbors
          local_neighbors = malloc(sizeof(int*) * range);
          amount_neighbors = malloc(sizeof(int) * range);
          for(int i = 0 ; i < (end - start); i++)
          {
               //we don't know how many neighbors there will be so we will have to use probe
               MPI_Status status;
               MPI_Probe(0, 0, MPI_COMM_WORLD, &status);

               int neighbor_count;
               MPI_Get_count(&status, MPI_INT, &neighbor_count);
               amount_neighbors[i] = neighbor_count;

               int *neighbors = malloc(neighbor_count * sizeof(int));

               MPI_Recv(neighbors, neighbor_count, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
               local_neighbors[i] = neighbors;
          }
     }

     //this is done by each of the MPI nodes until it stabilizes
     int changed_global = 1;
     while(changed_global)
     {
          std::vector<otherNodes*> remote_verticies;
          int non_local_cnt = 0;
          //ask other MPI nodes for the neighbors of 
          for(int i = 0; i < range; i++)
          {
               for(int a = 0; a < amount_neighbors[i] ; a++)
               {
                    //we don't have node locally
                    if(!(local_neighbors[i][a] >= chunk_start && local_neighbors[i][a] <= chunk_end))
                    {
                         //add to send requests, I don't know what I need to make for that
                         //calculate the owner
                         int vertex = local_neighbors[i][a];
                         int block_size = g_size / p;
                         int owner_process = vertex / block_size;
                         if(owner_process >= p)
                         {
                              owner_process = p - 1;
                         }
                         //using linked lists got a bit painful so I replaced them with some nice vectors
                         //https://stackoverflow.com/questions/15517991/search-a-vector-of-objects-by-object-attribute
                         auto it = find_if(remote_verticies.begin(), remote_verticies.end(), [&owner_process](const otherNodes& obj) {return obj->owner == owner_process;});
                         if(it != remote_verticies.end())
                         {
                              //this is not the best way to do because we just need the element at the index
                              auto index = std::distance(remote_verticies.begin(), it);
                              non_local_cnt++;
                              remote_verticies[index]->amount += 1;
                              remote_verticies[index]->nodes.push_back(vertex);
                         }else{
                              //new owner
                              otherNodes *newnode = new otherNodes;
                              std::vector<int> nodes;
                              nodes.push_back(vertex);
                              newnode->amount = 1;
                              newnode->nodes = nodes;
                              newnode->owner = owner_process;
                              remote_verticies.push_back(newnode);
                         }


                         
                    }
               }
          }
//     ------------------------------------------------------------
//     (1) Build arrays of remote neighbors that I must query
//     ------------------------------------------------------------
//     For each local node v:
//         For each neighbor u in neighbors[v]:
//             owner = process_owner(u)
//             If owner != my_rank:
//                 Append u to send_request[owner]

          //order the rank of the remote verticies
          std::sort(remote_verticies.begin(), remote_verticies.end(), [](const otherNodes *a, const otherNodes *b)
          {
               return a->owner < b->owner;
          });
          //all node ids to get
          int * remote_vec_request = malloc(sizeof(int) * non_local_cnt);
          int * remote_count_request = malloc(sizeof(int) * p);
          int n = 0;
          for(int i = 0; i < p; i++)
          {
               int sum = 0;
               if(remote_verticies[i]->owner == i)
               {
                    for( int a = 0; a < remote_verticies[i]->nodes.size(); a++)
                    {
                         remote_vec_request[n++] = remote_verticies[i]->nodes[a];
                         sum += 1;
                    }
               }
               remote_count_request[i] = sum;
          }
//     ------------------------------------------------------------
//     (2) Pack outgoing color requests into send buffers
//     ------------------------------------------------------------
//     For each process q:
//         sendbuf[q] = all node IDs in send_request[q]
//         sendcounts[q] = number of IDs in sendbuf[q]

          int * recv_counts = malloc(sizeof(int) * p);
          MPI_Alltoall(remote_count_request, p, MPI_INT, recv_counts, p, MPI_INT, MPI_COMM_WORLD);


          //create the total vertices that are requested by other processes
          int other_requested_recvs = 0;
          for(int i = 0; i < p; i++)
          {
               other_requested_recvs += recv_counts[i];
          }

          //size recv properly
          int * other_requested_verts = malloc(sizeof(int) * other_requested_recvs);

          //setup the offsets
          int * send_offsets = malloc(sizeof(int) * p);
          int * recv_offsets = malloc(sizeof(int) * p);

          send_offsets[0] = 0;
          for(int i = 0; i < p; i++)
          {
               send_offsets[i] = send_offsets[i - 1] + remove_count_request[i - 1];
          }

          recv_offsets[0] = 0;
          for(int i = 0; i < p; i++)
          {
               recv_offsets[i] = recv_offsets[i - 1] + recv_counts[i - 1];
          }

          MPI_Alltoallv(remote_vec_request, remote_count_request, send_offsets, MPI_INT, other_requested_verts, recv_counts, recv_offsets, MPI_INT, MPI_COMM_WORLD);


//     ------------------------------------------------------------
//     (3) Exchange variable-sized neighbor-requests (MPI_Alltoallv)
//         Everyone learns what colors it must send back
//     ------------------------------------------------------------
//     MPI_Alltoall(sendcounts => recvcounts)
//         # recvcounts[q] = how many nodes q wants colors for

//     Compute rdispls[] and sdispls[] from counts

//     MPI_Alltoallv(sendbuf  => recvbuf) 
//         # recvbuf contains the node IDs remote processes want
       
          std::vector<colorRequests> replies;
          int *replycounts = malloc(sizeof(int) * p);
          //fill with zeros
          for(int i = 0; i < p; i++)
          {    
               replycounts[i] = 0;
          }

          for(int i = 0; i < other_requested_recvs; i++)
          {
               colorRequests cr;
               cr.vector = other_requested_verts[i];
               cr.color = local_colors[i - chunk_start];
               replies.push_back(cr);

               int block_size = g_size / p;
               int owner_process = cr.vector / block_size;
               if(owner_process >= p)
               {
                    owner_process = p - 1;
               }
               replycounts[owner_process] += 1;
          }



//     ------------------------------------------------------------
//     (4) Prepare buffer of color replies
//     ------------------------------------------------------------
//     replies = empty list
//     For each node u in recvbuf:        # someone asked for these
//         Append (u, color[u]) to replies

//     replycounts[q] = how many pairs for each destination q
          int *incoming_replycounts = malloc(sizeof(int) * p);

          MPI_Alltoall(replycounts, 1, MPI_INT, incoming_replycounts, 1, MPI_INT, MPI_COMM_WORLD);

          //total reply size
          int reply_size = 0;
          for(int i = 0; i < p; i++)
          {

          }

          int *send_offsets2 = malloc(sizeof(int) * p);
          int *recv_offsets2 = malloc(sizeof(int) * p);

          send_offsets2[0] = 0;
          recv_offsets2[0] = 0;

          for(int i = 1; i < p; i++)
          {
               send_offsets2[i] = send_offsets2[i-1] + replycounts[i-1];
               recv_offsets2[i] = recv_offsets2[i-1] + incoming_replycounts[i-1];
          }

          int total_recv = recv_offsets2[p-1] + incoming_replycounts[p-1];
          MPI_Alltoallv(replies, replycounts, send_offsets2, MPI_2INT, incoming_replies, incoming_replycounts, recv_offsets2, MPI_2INT, MPI_COMM_WORLD);



//     ------------------------------------------------------------
//     (5) Exchange color information (MPI_Alltoallv)
//     ------------------------------------------------------------
//     MPI_Alltoall(replycounts => incoming_replycounts)
//     Compute displacement arrays

//     MPI_Alltoallv(replies => incoming_replies)
//         # incoming_replies contains:
//         #    list of (node_id, color) for all remote neighbors

          //hashmap for local node and colors of neighbors
          //don't want to implement a hashmap so im using the STL version
          std::map<int, std::vector<int>> neighbor_colors;

          //go through all the verticies
          for(int i = chunk_start; i < chunk_end; i++)
          {
               std::vector<int> v;
               neighbor_colors[i] = v;

               for(int a = 0; a < amount_neighbors[i]; a++)
               {
                    int neighbor = local_neighbors[i][a];
                    int c;
                    if(neighbor <= chunk_end && neighbor >= chunk_start)
                    {
                         //local so just grab color
                         c = local_colors[neighbor - chunk_start];
                    }else{
                         //get from the all to all MPI communication
                         //how would I get this from the simple arrays transfered between processes
                    }
                    v.push_back(c);
               }
          }

//     ------------------------------------------------------------
//     (6) Build full neighbor-color map for this iteration
//     ------------------------------------------------------------
//     For each local node v:
//         neighbor_colors[v] = empty multiset

//         For each neighbor u in neighbors[v]:
//             if u is local:
//                 Insert color[u] into neighbor_colors[v]
//             else:
//                 Look up color[u] from incoming_replies
//                 Insert that into neighbor_colors[v]

          //hashing time
          XXH64_state_t *state = XXH64_createState();
          XXH64_reset(state, 5050);

          int changed_local = 0;

          //new color list to update for next iteration
          int * new_color = malloc(sizeof(int) * range);
          for(int i = chunk_start; i < chunk_end; i++)
          {
               XXH64_update(state, &local_colors[i], sizeof(int));
               for(int a = 0; a < neighbor_colors[i]; a++)
               {
                   XXH64_update(state, &neighbor_colors[i][a], sizeof(int)); 
               }

               new_colors[i] = XXH64_digest(state);

               if (new_color[i] != local_colors[i])
               {
                    changed_local = 1;
               }
          }


//     ------------------------------------------------------------
//     (7) Compute new colors using hashing
//     ------------------------------------------------------------
//     changed_local = false

//     For each local node v:
//         new_color[v] =
//             hash( color[v], sorted(neighbor_colors[v]) )

//         If new_color[v] != color[v]:
//             changed_local = true

          //update the local verticies
          for(int i = 0; i < range; i++)
          {
               local_colors[i] = new_colors[i];
          }

//     ------------------------------------------------------------
//     (8) Commit new colors
//     ------------------------------------------------------------
//     For each local node v:
//         color[v] = new_color[v]

          //OR all the changed, so a single false will keep it going, all processes 
          MPI_Allreduce(&changed_local, &changed_global, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);
//     ------------------------------------------------------------
//     (9) Check global convergence
//     ------------------------------------------------------------
//     changed_global = MPI_Allreduce( changed_local, OR )

     }


    //Distribute vertices of graph G1 among ranks
    //local_vertices = subset of vertices owned by this rank

    uint64_t *colors1 = color_refinement(g);

    MPI_Barrier(MPI_COMM_WORLD);
    //Distribute vertices of graph G2 among ranks
    //local_vertices = subset of vertices owned by this rank


    uint64_t *colors2 = color_refinement(g2);

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
     printf("Thess two graphs can not possibly be the same graphs\n");
    }

    MPI_Finalize();
    return 0;
}