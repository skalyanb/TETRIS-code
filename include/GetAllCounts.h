//
// Created by Suman Kalyan Bera on 2019-06-17.
//

#ifndef SUBGRAPHCOUNT_GETALLCOUNTS_H
#define SUBGRAPHCOUNT_GETALLCOUNTS_H

#include <algorithm>
#include "GraphIO.h"
#include "Triadic.h"
#include "Graph.h"
#include "Digraph.h"

using namespace Escape;

// This function generates all non-induced counts for 3-vertex patterns.
// It is a wrapper function that calls the main algorithmic parts, and finally
// calls conversion functions to get induced counts.
//
// Input: pointer to CGraph, corresponding DAG, empty array nonInd with 4 entries
// No output: nonInd will have non-induced counts

void getAllThree(CGraph *cg, CDAG *dag, double (&nonInd)[4]) {
    double n, m, w;
    TriangleInfo info;

    n = cg->nVertices;
    m = 0;
    w = 0;

    for (VertexIdx i = 0; i < n; i++) {
        VertexIdx deg = cg->offsets[i + 1] - cg->offsets[i]; // degree of i
        m = m + deg;
        w = w + (deg * (deg - 1)) / 2; // updating total wedge count
    }
    m = m / 2;

    nonInd[0] = (n * (n - 1) * (n - 2)) / 6;  // number of independent sets
    nonInd[1] = m * (n - 2);    // number of plain edges
    nonInd[2] = w; // number of plain wedges

    info = betterWedgeEnumerator(&(dag->outlist));
    nonInd[3] = info.total;
}


#endif //SUBGRAPHCOUNT_GETALLCOUNTS_H
