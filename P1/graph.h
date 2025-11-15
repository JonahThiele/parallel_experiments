#include <stdio.h>
#include <stdlib.h>

typedef struct Node{
    //so we don't have to traverse the list to count the nodes
    int nodes;
    int Color;
    int v;
    struct Node* next;
} Node;

typedef struct Graph{
    int size;
    struct Node** list;
} Graph;

Graph* generate_graph(int v);

//also sets the color because we load the full graph before running the algorithm
void add_edge(Graph* g, int a, int b);

//debugging functions
void printGraph(Graph* g);