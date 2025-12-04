#include <stdio.h>
#include <stdlib.h>
#include "./graph.h"

Graph* generate_graph(int v)
{
    struct Graph* g = malloc(sizeof(Graph));
    g->list = malloc(sizeof(Node*) * v);
    g->size = v;
    for(int i = 0; i < v; i++)
    {
        g->list[i] = NULL;
    }
    return g;
}

//also sets the color because we load the full graph before running the algorithm
void add_edge(Graph* g, int a, int b)
{
    //if first one create head otherwise get to end of the list and add a new node
    if(g->list[a] == NULL)
    {
        Node* head = malloc(sizeof(Node));
        head->Color = 0;
        head->next = NULL;
        head->v = b;
        head->nodes = 0;
        g->list[a] = head;
    }else{
        Node *cur = g->list[a];
        cur->nodes += 1;
        while(cur->next != NULL)
        {
            cur = cur->next;
        }

        //create new node 
        Node* node = malloc(sizeof(Node));
        node->Color = 0;
        node->next = NULL;
        node->v = b;
        cur->next = node;
    }

    //if first one create head otherwise get to the end of the list and add a new node
    if(g->list[b] == NULL)
    {
        Node* head = malloc(sizeof(Node));
        head->Color = 0;
        head->next = NULL;
        head->v = a;
        head->nodes = 0;
        g->list[b] = head;

    }else{
        Node *cur = g->list[b];
        cur->nodes += 1;
        while(cur->next != NULL)
        {
            cur = cur->next;
        }

        //create new node 
        Node* node = malloc(sizeof(Node));
        node->Color = 0;
        node->next = NULL;
        node->v = a;
        cur->next = node;
    }
}

//debugging functions
void printGraph(Graph* g)
{
    printf("V: Adjacency List\n");
    for(int v = 0; v < g->size;v++)
    {
        printf("V: %d -> ", v);
        Node *cur = g->list[v];
        while(cur)
        {
            printf("%d ", cur->v);
            cur = cur->next;
        }
        printf("\n");
    }
}