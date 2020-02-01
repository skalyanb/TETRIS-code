//
// Created by Suman Kalyan Bera on 2019-06-17.
//

#ifndef SUBGRAPHCOUNT_TRIANGLEESTIMATORS_H
#define SUBGRAPHCOUNT_TRIANGLEESTIMATORS_H

#include <random>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <iostream>
#include <unordered_set>
#include <utility>
#include <string>
#include <ctime>
#include <iterator>


#include <algorithm>
#include "GraphIO.h"
#include "Triadic.h"
#include "Graph.h"
#include "Digraph.h"
#include "EstimatorUtilStruct.h"
//#include "EstimateEdgeCount.h"


using namespace Escape;

/**
 * A structure to keep all the input parameters to various
 triangle estimators.
 * filename: the path to the datafile in escape format
 * seed_count: the number of seed vertices to start random walks from (only for random walk based estimators)
 * walk_length: the length of each random walk (only for random walk based estimators)
 * subsample_size:
 * no_of_repeat: The number of time the estimators are going to be run.
 * sparsification_prob: The sampling probability of the edges in graph sparsification based estimators
 */


/**
 * Exactly counting the number of triangles in a graph.
 * This subroutine can be used to count the number of triangles in sparsified graphs
 * @param cg
 * @return
 */

Estimates CountExactTriangles (CGraph *cg)
{
//    printf ("Exactly counting the number of triangles...\n");

    Estimates output;
    cg->sortById();
    CGraph cg_relabel = cg->renameByDegreeOrder();
    cg_relabel.sortById();
    CDAG dag = degreeOrdered(&cg_relabel);
    (dag.outlist).sortById();

    //printf("Conversion comple. Going to count triangles now. \n");
    TriangleInfo info;
    info = betterWedgeEnumerator(&(dag.outlist));
    output.estimate = info.total;

    // Free up the memory
    delCGraph(dag.inlist);
    delCGraph(dag.outlist);
    delCGraph(cg_relabel);
    //printf ("  Triangle=%lld ",info.total);
    return output;
}

/**
 * Algorithm:
 * 1. Sample a vertex u.a.r. Call this vertex u.
 * 2. Repeat p many times: // Start a random walk from u.
 * 3.   Sample a neighbor of u u.a.r. Call this vertex v.
 * 4.   Let e= (u,v) and Assume u < v according to degree ordering.
 * 5.   Sample a neighbor w of u u.a.r.
 * 6.   Check for triangle (e,w) where u < v < w according to the degree ordering.
 * 7.   If triangle found, then set Z = d_u, else set Z = 0.
 * 8.   Y += scaling * Z
 * 9.   set u=v (move on to the next step)
 *10. X = Y/p (average over all the steps)
 */




#endif //SUBGRAPHCOUNT_TRIANGLEESTIMATORS_H
