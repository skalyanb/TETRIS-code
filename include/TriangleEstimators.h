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
#include <set>
#include <utility>
#include <string>
#include <ctime>
#include <iterator>


#include <algorithm>
#include "GraphIO.h"
#include "Triadic.h"
#include "Graph.h"
#include "Digraph.h"

using namespace Escape;

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

/**
 * There needs to be three functions:
 * 1.   First one does a random walk and collect a long path (collection of edges)
 *      This should be stored in a structure with two elements: edge index and edge degree.
 * 2. Now there would be few variations for comparison.
 *      First variation: all edges collected so far are used for triangle estimation.
 *      Second variation: sample some edges proportional to their degrees and then
 *                        use them for triangle estimation.
 */

// A structure to store an edge with all the relavant information
struct EdgeStruct {
    VertexIdx u;  // the lower degree end point
    VertexIdx v; // the higher degree end point
    EdgeIdx index;   // nbors[index] represents (src, dest)
    VertexIdx degree; // deg(u)
};

// A collection of EdgeStruct
struct EdgeStructCollection {
    EdgeIdx no_of_edges; // No of edges in the collection
    std::vector<EdgeStruct> edge_list; // The list of edges
    std::set<EdgeIdx> visited_edge_set; // The set of unique edges in the collection
    std::set<VertexIdx> visited_vertex_set; // The set of unique vertices in the collection
};

// A structure to keep all the input parameters to various version of our algorithms
struct Parameters {
    VertexIdx seed_count;
    EdgeIdx walk_length;
    EdgeIdx subsample_size;
    int no_of_repeat;
    std::string filename;
};

// A structure to hold the estimated triangle count and the fraction of the graph seen in the process.
struct Estimates {
    double triangle_estimate;
    double fraction_of_vertices_seen;
    double fraction_of_edges_seen; // This fraction is with respect to twice the number of edges
};

struct ExactCount {
    double triangle_count;
};

bool ComparatorByEdgesSeen(Estimates a, Estimates b) {
    return (a.fraction_of_edges_seen < b.fraction_of_edges_seen);
}

bool ComparatorByVerticesSeen(Estimates a, Estimates b) {
    return (a.fraction_of_vertices_seen < b.fraction_of_vertices_seen);
}

// A structure to store the details of vertices and edges observed during the algorithm
struct ObservedGraphStats {
    VertexIdx VerticesSeenInRandomWalk;
    VertexIdx VerticesSeenAsNbors;
    EdgeIdx EdgesSeenInRandomWalk;
    EdgeIdx EdgesSeenAsNbors;
};


ExactCount CountExactTriangles (CGraph *cg)
{
    printf ("Exactly counting triangles...\n");

    ExactCount output;
    cg->sortById();
    CGraph cg_relabel = cg->renameByDegreeOrder();
    cg_relabel.sortById();
    CDAG dag = degreeOrdered(&cg_relabel);
    (dag.outlist).sortById();
//  (dag.inlist).sortById();

    TriangleInfo info;
    info = betterWedgeEnumerator(&(dag.outlist));
    output.triangle_count = info.total;

    printf ("Triangle=%lld \n",info.total);
    return output;
}

EdgeStructCollection GetEdgesByRandomWalk(CGraph *cg, Parameters params, std::mt19937 mt) {
    // Get the number of vertices, number of edges, Initialize scaling to number of edges.
    VertexIdx n = cg->nVertices;
    EdgeIdx m = cg->nEdges; // TODO note that m is double of the number of edges in the graph
    VertexIdx seed_count = params.seed_count;
    EdgeIdx walk_length = params.walk_length;

    // The list of edges found by the random walk process
    std::vector<EdgeStruct> edge_list;

    // Keep track of the vertices and edges seen so far by the random walk
    std::set<EdgeIdx> visited_edge_set;
    std::set<VertexIdx> visited_vertex_set;

    // Set up random number generator
    //std::random_device rd;
    // std::mt19937 mt(rd());

    // Initialize a uniform distribution for generating seed vertices
    std::uniform_int_distribution<VertexIdx> dist_seed_vertex(0, n - 1);

    // Perform a random walk for seed_count many times
    for (VertexIdx sC = 0; sC < seed_count; sC++) {
        // Pick a random seed vertex
        VertexIdx seed = dist_seed_vertex(mt); // TODO: verify randomness
        VertexIdx parent, child;
        parent = seed;
        // Perform a random walk of length walk_length from the seed vertex
        for (EdgeIdx wL = 0; wL < walk_length; wL++) {
            VertexIdx deg_of_parent = cg->offsets[parent + 1] - cg->offsets[parent]; // degree of parent vertex
            // TODO assumption: deg_of_parent is non-zero, add boundary check
            // If the degree of parent vertex is 0, then we attempt seed another vertex
            // We continue this process for certain no of times before giving up and continuing with the run
            int attempt = 0;
            while (deg_of_parent == 0 && attempt < 1000) {
                printf("BAd luck! Trying again %d\n", attempt);
                parent = dist_seed_vertex(mt);
                deg_of_parent = cg->offsets[parent + 1] - cg->offsets[parent];
                attempt++;
            }

            // Take a step of the random walk
            std::uniform_int_distribution<VertexIdx> dist_parent_nbor(0, deg_of_parent - 1);
            EdgeIdx random_nbor_edge =
                    cg->offsets[parent] + dist_parent_nbor(mt); // TODO: check randomness. same seeding ok?
            child = cg->nbors[random_nbor_edge]; // child is the next vertex on the random walk
            VertexIdx deg_of_child = cg->offsets[child + 1] - cg->offsets[child]; // degree of the child vertex

            // Create the edge info structure
            VertexIdx u, v, low_deg;
            if (deg_of_child < deg_of_parent || (deg_of_child == deg_of_parent && child < parent)) {
                u = child;
                v = parent;
                low_deg = deg_of_child;
            } else {
                u = parent;
                v = child;
                low_deg = deg_of_parent;
            }
            edge_list.push_back(EdgeStruct{u, v, random_nbor_edge, low_deg});

            // Collect the distinct edges and distinct vertices visited so far in a set
            visited_edge_set.insert(random_nbor_edge);
            visited_vertex_set.insert(parent);

            // The random walk proceeds with the vertex child
            parent = child;
            // If {u,v} is an isolated edge, then we will stuck here: so find a different seed and continue from there
            if (deg_of_child == 1 && deg_of_parent == 1)
                parent = dist_seed_vertex(mt);
        }
    }

    // Create return structure for this function.
    // The number of edges produced by the random walk is walk_length;
    EdgeStructCollection returnEdgeCollection = {walk_length, edge_list, visited_edge_set, visited_vertex_set};
//    printf("Random walk: walk length = %lld ",walk_length);
//    printf("Edges Seen=%lu, Vertices Seen=%lu\n",visited_edge_set.size(),visited_vertex_set.size());

    return returnEdgeCollection;
}

Estimates SampleByEdgeDegree(CGraph *cg, EdgeStructCollection &edge_collection, Parameters params, std::mt19937 mt) {

    // Set up random number generator
    //std::random_device rd;
    //std::mt19937 mt(rd());

    // Collect the degrees of all the edges in a vector so that we can iterate over them
    // The iterator is useful in setting up the discrete_distribution to pick an edge
    // proportionate to its degree
    std::vector<VertexIdx> edge_degree_list(edge_collection.no_of_edges);
    double sum_of_degree = 0.0;
    for (EdgeIdx i = 0; i < edge_collection.no_of_edges; i++) {
        edge_degree_list[i] = edge_collection.edge_list[i].degree;
        sum_of_degree += edge_collection.edge_list[i].degree;
    }

    // std::discrete_distribution is ideal for sampling an edge proportionate to its degree
    // Refer to third answer here: https://stackoverflow.com/questions/1761626/weighted-random-numbers
    std::discrete_distribution<EdgeIdx> wghted_dist_on_edges(edge_degree_list.begin(), edge_degree_list.end());
    //for (double x:wghted_dist_on_edges.probabilities()) std::cout << x << " ";

    // Variables used for estimating the triangle count
    EdgeIdx m = cg->nEdges; // TODO note that m is double of the number of edges in the graph
    double X = 0, Y = 0, Z = 0, scaling = edge_collection.no_of_edges;
    Count raw_count = 0; // TODO remove, does not have any purpose

    // Stats about the algorithm: how many new distinct vertices and edges are
    // seen during neighbor sampling process
    ObservedGraphStats obs_graph_stats = {};
    obs_graph_stats.VerticesSeenInRandomWalk = edge_collection.visited_vertex_set.size();
    obs_graph_stats.EdgesSeenInRandomWalk = edge_collection.visited_edge_set.size();

    for (EdgeIdx i = 0; i < params.subsample_size; i++) {

        // Sample an edge proportionate to its degree
        EdgeIdx index = wghted_dist_on_edges(mt);
        EdgeStruct edge = edge_collection.edge_list[index];

        // Sample u.a.r a neighbor w of the edge (of lower degree end point)
        std::uniform_int_distribution<VertexIdx> distNbor(0, edge.degree - 1);
        EdgeIdx random_nbor_edge = cg->offsets[edge.u] + distNbor(mt);
        VertexIdx w = cg->nbors[random_nbor_edge];

        // Update the edges and vertices seen so far data structure
        // TODO to decide whether to pass a pointer or pass a reference for updating the sets
        edge_collection.visited_edge_set.insert(random_nbor_edge);
        edge_collection.visited_vertex_set.insert(w);

        // Now check for a triangle: a triangle is only found if u < v < w and {u,v,w} forms a triangle.
        VertexIdx deg_of_w = cg->offsets[w + 1] - cg->offsets[w];
        VertexIdx deg_of_v = cg->offsets[edge.v + 1] - cg->offsets[edge.v];
        if (cg->isEdgeBinary(w, edge.v) && (deg_of_w > deg_of_v || (deg_of_w == deg_of_v && w > edge.v))) {
            EdgeIdx third_edge = cg->getEdgeBinary(edge.v, w);
            edge_collection.visited_edge_set.insert(third_edge);
            //Z = edge.degree;
            Z = 1;
            raw_count++; // Found a triangle
        } else
            Z = 0;
        //Y += scaling * Z;
        Y += Z;
    }
    X = Y / params.subsample_size;
    double sum_degree_estimate = (cg->nEdges / 2.0) * sum_of_degree / edge_collection.no_of_edges;
    X = sum_degree_estimate * X;
    // Create return object and store relevant stats
    Estimates return_estimate = {};
    return_estimate.triangle_estimate = X;
    return_estimate.fraction_of_vertices_seen = edge_collection.visited_vertex_set.size() * 100.0 / cg->nVertices;
    return_estimate.fraction_of_edges_seen = edge_collection.visited_edge_set.size() * 100.0 / m;

    // Some other stats about the algorithm TODO decide what to do with them
    obs_graph_stats.VerticesSeenAsNbors =
            edge_collection.visited_vertex_set.size() - obs_graph_stats.VerticesSeenInRandomWalk;
    obs_graph_stats.EdgesSeenAsNbors = edge_collection.visited_edge_set.size() - obs_graph_stats.EdgesSeenInRandomWalk;

//    printf("**** Sample weighted *****\n");
//    printf("Sum of degree estimate = %lf,actual edges = %lld, err =%lf\n",
//           sum_degree_estimate,cg->nEdges,std::abs(sum_degree_estimate-cg->nEdges)*100/cg->nEdges);
//    printf("Edges Seen=%lu, Vertices Seen=%lu\n",edge_collection.visited_edge_set.size(),edge_collection.visited_vertex_set.size());
//    printf ("%lf,%lf,%lf\n",X,vertex_fraction,edge_fraction);
    return return_estimate;
}

Estimates SampleAllEdges(CGraph *cg, EdgeStructCollection &edge_collection, std::mt19937 mt) {

    // Retrieve relevant information from the edge_collection structure
    std::vector<EdgeStruct> edge_list;
    edge_list = edge_collection.edge_list;
    EdgeIdx no_of_edges = edge_collection.no_of_edges;

    // Set up random number generator
    //std::random_device rd;
    //std::mt19937 mt(rd());

    // Variables used for estimating the triangle count
    EdgeIdx m = cg->nEdges; // TODO note that m is double of the number of edges in the graph
    double X = 0, Y = 0, Z = 0, scaling = m / 2;
    Count raw_count = 0; // TODO remove, does not have any purpose

    // Stats about the algorithm: how many new distinct vertices and edges are
    // seen during neighbor sampling process
    ObservedGraphStats obs_graph_stats = {};
    obs_graph_stats.VerticesSeenInRandomWalk = edge_collection.visited_vertex_set.size();
    obs_graph_stats.EdgesSeenInRandomWalk = edge_collection.visited_edge_set.size();

    for (EdgeIdx i = 0; i < no_of_edges; i++) {
        // e = (u,v) is the current with u being the lower degree end-point
        // degree of the edge is degree of u.
        // Sample a neighbor w of u u.a.r.
        EdgeStruct edge = edge_list[i];

        std::uniform_int_distribution<VertexIdx> dist_nbor(0, edge.degree - 1);
        EdgeIdx random_nbor_edge = cg->offsets[edge.u] + dist_nbor(mt);
        VertexIdx w = cg->nbors[random_nbor_edge];

        // Update the edges and vertices seen so far data structure
        // TODO to decide whether to pass a pointer or pass a reference for updating the sets
        edge_collection.visited_edge_set.insert(random_nbor_edge);
        edge_collection.visited_vertex_set.insert(w);

        // Now check for a triangle: a triangle is only found if u < v < w and {u,v,w} forms a triangle.
        VertexIdx deg_of_w = cg->offsets[w + 1] - cg->offsets[w];
        VertexIdx deg_of_v = cg->offsets[edge.v + 1] - cg->offsets[edge.v];
        if (cg->isEdgeBinary(w, edge.v) && (deg_of_w > deg_of_v || (deg_of_w == deg_of_v && w > edge.v))) {
            EdgeIdx third_edge = cg->getEdgeBinary(edge.v, w);
            edge_collection.visited_edge_set.insert(third_edge);
            Z = edge.degree;
            raw_count++; // Found a triangle
        } else
            Z = 0;
        Y += scaling * Z;
    }
    X = Y * 1.0 / no_of_edges;
    // Create return object and store relevant stats
    Estimates return_estimate = {};
    return_estimate.triangle_estimate = X;
    return_estimate.fraction_of_vertices_seen = edge_collection.visited_vertex_set.size() * 100.0 / cg->nVertices;
    return_estimate.fraction_of_edges_seen = edge_collection.visited_edge_set.size() * 100.0 / m;

    // Some other stats about the algorithm TODO decide what to do with them
    obs_graph_stats.VerticesSeenAsNbors =
            edge_collection.visited_vertex_set.size() - obs_graph_stats.VerticesSeenInRandomWalk;
    obs_graph_stats.EdgesSeenAsNbors = edge_collection.visited_edge_set.size() - obs_graph_stats.EdgesSeenInRandomWalk;
//    printf("**** Sample all *****\n");
//    printf("Edges Seen=%lu, Vertices Seen=%lu\n",edge_collection.visited_edge_set.size(),edge_collection.visited_vertex_set.size());
//    printf ("%lf,%lf,%lf\n",X,vertex_fraction,edge_fraction);
    return return_estimate;
}

//EsTRaW
Estimates EstTriByRWandWghtedSampling(CGraph *cg, Parameters params) {
    Estimates output;
    EdgeStructCollection randomEdgeCollection;

    // Set up random number generator
    std::random_device rd;
    // Using this random number generator initializize a PRNG: this PRNG is passed along to
    // draw an element from various distribution
    std::mt19937 mt(rd());

    randomEdgeCollection = GetEdgesByRandomWalk(cg, params, mt);
    output = SampleByEdgeDegree(cg, randomEdgeCollection, params, mt);
    return output;
}

Estimates EstTriByRWandSimpleSampling(CGraph *cg, Parameters params) {

    Estimates output;
    EdgeStructCollection randomEdgeCollection;

    // Set up random number generator
    std::random_device rd;
    // Using this random number generator initializize a PRNG: this PRNG is passed along to
    // draw an element from various distributions
    std::mt19937 mt(rd());

    randomEdgeCollection = GetEdgesByRandomWalk(cg, params, mt);
    output = SampleAllEdges(cg, randomEdgeCollection, mt);
    return output;
}




Estimates EstTriByRW(CGraph *cg, Parameters params) {}

/**
 * Algorithm:
 * 1. Sample every edge with probability p
 * 2. Let the subsampled graph be G_p.
 * 3. Count the number of triangles in G_p and scale the count up by 1/p^3
 * @param cg
 * @param params
 * @return
 */
Estimates EstTriByEdgeSparsification(CGraph *cg, Parameters params)
{
    // This algorithm is implemented by a three step process.
    // First, we construct a new graph object from the subsampled edges
    // Then, we convert it into a CGraph object
    // Finally, we call exact triangle count routine on CGraph object

    double p = 0.1;
    VertexIdx n = cg->nVertices;
    EdgeIdx m = cg->nEdges; // Note that m is twice the number of edges
    // Set up random number generator
    std::random_device rd;
    // Using this random number generator initializize a PRNG: this PRNG is passed along to
    // draw an element from various distribution0-
    std::mt19937 mt(rd());
    std::bernoulli_distribution sparsifier(p);

    // Iterate over the graph to sparsify it
    VertexIdx src = 0, dst = 0;
    EdgeIdx nEdges = 0;
    std::vector<VertexIdx > srcs;
    std::vector<VertexIdx > dsts;
    for (EdgeIdx i = 0; i < m; i++) {
        if ( i == cg->offsets[src+1])
            ++src;
        dst = cg->nbors[i]; // (src,dst) is the i-th edge in the nbor list
        if (src < dst) {    // This check ensures that each edge {src,dst} is considered only once
            if (sparsifier(mt)) { // The edge {src,dst} will be part of the sparsified graph
                srcs.emplace_back(src);
                dsts.emplace_back(dst);
                ++nEdges;
                srcs.emplace_back(dst);  // The graph representation (adjacency list) requires both the edges to be present
                dsts.emplace_back(src);
                ++nEdges;
            }
        }
    }
    Graph g_p;
    g_p.nVertices = n;
    g_p.nEdges = nEdges;
    g_p.srcs = new VertexIdx[g_p.nEdges];
    g_p.dsts = new VertexIdx[g_p.nEdges];

    std::copy(std::begin(srcs), std::end(srcs), g_p.srcs);
    std::copy(std::begin(dsts), std::end(dsts), g_p.dsts);

    CGraph cg_p = makeCSR(g_p);

    ExactCount triangles_in_Gp = CountExactTriangles(&cg_p);
    double triangleEstimate = (triangles_in_Gp.triangle_count * 1.0)/ (p*p*p);
    std::cout << nEdges <<"  " << triangleEstimate << std::endl;
    printf ("%lld,  %lf",nEdges,triangleEstimate);

    Estimates return_estimate = {};
    return_estimate.triangle_estimate = triangleEstimate;
    return_estimate.fraction_of_vertices_seen = 1; // TODO fix this
    return_estimate.fraction_of_edges_seen = nEdges * 100.0 / m;

    return return_estimate;
}


#endif //SUBGRAPHCOUNT_TRIANGLEESTIMATORS_H
