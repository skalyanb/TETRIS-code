//
// Created by Suman Kalyan Bera on 2019-09-19.
//

#ifndef SUBGRAPHCOUNT_BASELINEESTIMATORS_H
#define SUBGRAPHCOUNT_BASELINEESTIMATORS_H

#include "TriangleEstimators.h"


/**
 * Algorithm:
 * 1. Sample every edge with probability p
 * 2. Let the subsampled graph be G_p.
 * 3. Count the number of triangles in G_p and scale the count up by 1/p^3
 * @param cg
 * @param params
 * @return
 */
Estimates EstTriBySparsification(CGraph *cg, Parameters params)
{
    // This algorithm is implemented by a three step process.
    // First, we construct a new graph object from the subsampled edges; call it G_p
    // Then, we convert it into a CGraph object
    // Finally, we call exact triangle count routine on CGraph object

    double p = params.sparsification_prob;
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
    std::set <VertexIdx > visited_vertex_set; // This will be used to find the number of vertices in G_p
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
                visited_vertex_set.insert(src);
                visited_vertex_set.insert(dst);
            }
        }
    }
    // Construct a graph from the sprsified edges
    Graph G_p;
    G_p.nVertices = n;
    // Note that we do not use visited_vertex_set.size() here. The graph representation requires the
    // indices to be in the range 0, number_of_vertices-1.
    G_p.nEdges = nEdges;
    G_p.srcs = new VertexIdx[G_p.nEdges];
    G_p.dsts = new VertexIdx[G_p.nEdges];

    std::copy(std::begin(srcs), std::end(srcs), G_p.srcs);
    std::copy(std::begin(dsts), std::end(dsts), G_p.dsts);

    CGraph CG_p = makeCSR(G_p);

    // Count the exact number of triangles in the sparsified graph G_p and scale it up.
    Estimates triangles_in_Gp = CountExactTriangles(&CG_p);
    double triangleEstimate = (triangles_in_Gp.triangle_estimate * 1.0)/ (p*p*p);
//    std::cout << nEdges <<"  " << triangleEstimate << std::endl;
//    printf ("%lld,  %lf",nEdges,triangleEstimate);

    Estimates return_estimate = {};
    return_estimate.triangle_estimate = triangleEstimate;
    return_estimate.fraction_of_vertices_seen = visited_vertex_set.size() * 100.0 / n;
    return_estimate.fraction_of_edges_seen = G_p.nEdges * 100.0 / m;

    return return_estimate;
}


/**
 * Algorithm:
 * 0. Let p be the sparsification parameter. The goal is to sample mp many edges of the graph.
 * 1. Sample mp many edges from the graph, independently and uniformly at random.
 * 2. Let the subsampled graph be G_p.
 * 3. Count the number of triangles in G_p and scale the count up by 1/p^3
 * @param cg
 * @param params
 * @return
 */
Estimates EstTriByUniformSampling(CGraph *cg, Parameters params)
{
    // This algorithm is implemented by a three step process.
    // First, we construct a new graph object from the subsampled edges; call it G_p
    // Then, we convert it into a CGraph object
    // Finally, we call exact triangle count routine on CGraph object

    double p = params.sparsification_prob;
    VertexIdx n = cg->nVertices;
    EdgeIdx m = cg->nEdges; // Note that m is twice the number of edges
    EdgeIdx edges_in_G_p = floor(m*p)/2; // We sample edges_in_G_p many edges independently anf uniformly at random
    // Why we divide by 2? Because the number of edges in cg ig m/2. So, number of edges in G_p = mp/2.
    // This will ensure that we only see p fraction of the total edges. Note that while constructing the
    // graph G_p, we will have mp many entries in the adjacency list representation.

    // Set up random number generator
    std::random_device rd;
    // Using this random number generator initializize a PRNG: this PRNG is passed along to
    // draw an element from various distribution0-
    std::mt19937 mt(rd());
    std::uniform_int_distribution<EdgeIdx> unif_rand_edge(0, m - 1);

    // Sample edges_in_G_p many edges independently and u.a.r from G
    VertexIdx src = 0, dst = 0;
    EdgeIdx nEdges = 0;
    std::vector<VertexIdx > srcs;
    std::vector<VertexIdx > dsts;

    std::set <VertexIdx > visited_vertex_set; // This will be used to find the number of vertices in G_p
    std::set <EdgeIdx > visited_edge_set; // This will be used to find the number of vertices in G_p

    for (EdgeIdx i = 0; i < edges_in_G_p ; i++) {
        EdgeIdx e = unif_rand_edge(mt);
        dst = cg->nbors[e];
        // binary search the array cg->offset to find the index v, such that edges[e] lies between cg->offset[v]
        // and cg->offset[v+1]. We achieve this by calling std::upper_bound
        src = std::upper_bound (cg->offsets, cg->offsets+n, e) - cg->offsets -1;

        // sanity check: the edge {src,dst} must have index edges[e].
        if (cg->getEdgeBinary(src,dst) != e)
            printf("Bug in the code!! Run!!!! \n\n");

        srcs.emplace_back(src);
        dsts.emplace_back(dst);
        ++nEdges;
        srcs.emplace_back(dst);  // The graph representation (adjacency list) requires both the edges to be present
        dsts.emplace_back(src);
        ++nEdges;

        visited_vertex_set.insert(src);
        visited_vertex_set.insert(dst);
        visited_edge_set.insert(e);
    }

    // Construct the induced graph G_p
    // Note that G_p is a multi-graph. When counting triangles exactly in G_p, we consider this multiplicity.
    Graph G_p;
    G_p.nVertices = n;
    G_p.nEdges = nEdges;
    G_p.srcs = new VertexIdx[G_p.nEdges];
    G_p.dsts = new VertexIdx[G_p.nEdges];

    std::copy(std::begin(srcs), std::end(srcs), G_p.srcs);
    std::copy(std::begin(dsts), std::end(dsts), G_p.dsts);

    // Convert the graph into CSR representation
    CGraph CG_p = makeCSR(G_p);

    // Count the exact number of triangles in the sparsified graph G_p and scale it up.
    Estimates triangles_in_Gp = CountExactTriangles(&CG_p);
    double scale = 1.0 / pow(p,3);
    double triangleEstimate = triangles_in_Gp.triangle_estimate * scale ;
//    std::cout << nEdges <<"  " << triangleEstimate << std::endl;
    printf ("%lld,  %lf",nEdges,triangleEstimate);

    Estimates return_estimate = {};
    return_estimate.triangle_estimate = triangleEstimate;
    return_estimate.fraction_of_vertices_seen = visited_vertex_set.size() * 100.0 / n;
    return_estimate.fraction_of_edges_seen = 2 *visited_edge_set.size() * 100.0 / m;
    // Why the multiplication by 2? The fraction of edges seen is effectively
    // compared against all the entries in the adjacency list, which is 2m.

    return return_estimate;
}


Estimates EstTriByRW(CGraph *cg, Parameters params) {

    VertexIdx n = cg->nVertices;
    EdgeIdx m = cg->nEdges; // Note that m is double of the number of edges in the graph
    VertexIdx seed_count = params.seed_count;
    //EdgeIdx walk_length = params.walk_length;

    double p = params.sparsification_prob;
    EdgeIdx walk_length = floor(m*p)/2;

    // Set up random number generator
    std::random_device rd;
    // Using this random number generator initializize a PRNG: this PRNG is passed along to
    // draw an element from various distribution0-
    std::mt19937 mt(rd());

    // Keep track of the vertices and edges seen so far by the random walk
    std::set<EdgeIdx> visited_edge_set;
    std::set<VertexIdx> visited_vertex_set;
    EdgeIdx nEdges = 0;
    std::vector<VertexIdx > srcs;
    std::vector<VertexIdx > dsts;
    std::set<std::pair <VertexIdx,VertexIdx>> edgeList; // This list will be used to track distinct edges in the random walk.

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
            VertexIdx deg_of_parent = cg->degree(parent); // degree of parent vertex
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
            VertexIdx deg_of_child = cg->degree(child);// degree of the child vertex

            srcs.emplace_back(parent);
            dsts.emplace_back(child);
            ++nEdges;
            srcs.emplace_back(child);  // The graph representation (adjacency list) requires both the edges to be present
            dsts.emplace_back(parent);
            ++nEdges;

            edgeList.insert(std::make_pair(parent,child));
            edgeList.insert(std::make_pair(child,parent));

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


    // Construct a graph from the edges of the random walk
    // Construct the induced graph G_p
    // Note that G_p is a multi-graph. When counting triangles exactly in G_p, we consider this multiplicity.
    Graph G_p;
    G_p.nVertices = n;
    G_p.nEdges = nEdges;
    G_p.srcs = new VertexIdx[G_p.nEdges];
    G_p.dsts = new VertexIdx[G_p.nEdges];

    std::copy(std::begin(srcs), std::end(srcs), G_p.srcs);
    std::copy(std::begin(dsts), std::end(dsts), G_p.dsts);

    // Convert the graph into CSR representation
    CGraph CG_p = makeCSR(G_p);

    // Construct a simple graph from the edges of the random walk
    Graph G_p_s;
    G_p_s.nVertices = n;
    G_p_s.nEdges = edgeList.size();
    G_p_s.srcs = new VertexIdx[G_p.nEdges];
    G_p_s.dsts = new VertexIdx[G_p.nEdges];

    std::copy(std::begin(srcs), std::end(srcs), G_p.srcs);
    std::copy(std::begin(dsts), std::end(dsts), G_p.dsts);

    // Convert the graph into CSR representation
    CGraph CG_p_s = makeCSR(G_p_s);


    // Count the exact number of triangles in the sparsified graph G_p and scale it up.
    Estimates triangles_in_Gp = CountExactTriangles(&CG_p);
    //double p = walk_length * 1.0/m;
    double scale = 1.0 / pow(p,3);
    printf("p=%lf,  scale = %lf \n",p,scale);
    double triangleEstimate = triangles_in_Gp.triangle_estimate * scale ;
    printf ("%lld,  %lf",nEdges,triangleEstimate);

    Estimates return_estimate = {};
    return_estimate.triangle_estimate = triangleEstimate;
    return_estimate.fraction_of_vertices_seen = visited_vertex_set.size() * 100.0 / n;
    return_estimate.fraction_of_edges_seen = 2*visited_edge_set.size() * 100.0 / m;
    // Why the multiplication by 2? The fraction of edges seen is effectively
    // compared against all the entries in the adjacency list, which is 2m.

    return return_estimate;
}


#endif //SUBGRAPHCOUNT_BASELINEESTIMATORS_H
