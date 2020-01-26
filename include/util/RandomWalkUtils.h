//
// Created by Suman Kalyan Bera on 2020-01-23.
//

#ifndef SUBGRAPHCOUNT_RANDOMWALKUTILS_H
#define SUBGRAPHCOUNT_RANDOMWALKUTILS_H

#include <random>

#include "../Graph.h"

VertexIdx AlternateSeed (CGraph* cg, std::mt19937 mt)
{
    int attempt = 0;
    VertexIdx deg_of_parent = 0, parent = 0;
    std::uniform_int_distribution<VertexIdx> dist_seed_vertex(0, cg->nVertices - 1); // Initialize a uniform distribution for generating seed vertices
    while (deg_of_parent == 0 && attempt < 1000) {
        printf("BAd luck! Trying again %d\n", attempt);
        parent = dist_seed_vertex(mt);
        deg_of_parent = cg->degree(parent);
        attempt++;
    }
    return parent;
}

#endif //SUBGRAPHCOUNT_RANDOMWALKUTILS_H
