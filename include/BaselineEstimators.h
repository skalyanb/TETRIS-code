//
// Created by Suman Kalyan Bera on 2019-09-19.
//

#ifndef SUBGRAPHCOUNT_BASELINEESTIMATORS_H
#define SUBGRAPHCOUNT_BASELINEESTIMATORS_H

#include <set>
#include <unordered_map>

#include "EstimatorUtilStruct.h"
#include "TriangleEstimators.h"
#include "TETRIS.h"


bool IsTrue (bool b) {return b;}

CGraph MakeMultiGraph (VertexIdx n, std::vector<VertexIdx > srcs, std::vector<VertexIdx > dsts) {

    // First rename the vertices in srcs and dsts to figure out how many vertices are actually
    // and remove those vertice.

    // Create a map that will map every vertex id the in the input file to an integer in the output file
    std::unordered_map <VertexIdx, VertexIdx> dict;
    dict.reserve(100000000); // 100M
    VertexIdx index = 0;

    for (VertexIdx idx=0; idx < srcs.size(); idx++) {
        VertexIdx node1 = srcs[idx];
        VertexIdx node2 = dsts[idx];
        // Find if the vertices are already present in the dictionary
        auto it1 = dict.find (node1);
        auto it2 = dict.find (node2);

        // If the first vertex is present then replace the node1 with the corresponding node id;
        // Otherwise create a new node for this vertex.
        if (it1!= dict.end())
            node1 = it1->second;
        else {
            dict[node1] = index;
            node1 = index;
            index++;
        }
        // Repeat the same for the second node
        if (it2!= dict.end())
            node2 = it2->second;
        else {
            dict[node2] = index;
            node2 = index;
            index++;
        }
        //Add the edge (node1,node2)
        srcs[idx] = node1;
        dsts[idx] = node2;
    }

    //Sanity check
    if (dict.size()!= index)
        printf("Dict size=%lu, index = %lld, Something went wrong.",dict.size(),index);

    Graph G_p;
    G_p.nVertices = dict.size();
    G_p.nEdges = srcs.size();
    G_p.srcs = new VertexIdx[G_p.nEdges];
    G_p.dsts = new VertexIdx[G_p.nEdges];

    std::copy(std::begin(srcs), std::end(srcs), G_p.srcs);
    std::copy(std::begin(dsts), std::end(dsts), G_p.dsts);

    // Convert the graph into CSR representation. The second parameters requests for in-pace operation.
    CGraph CG_p = makeCSR(G_p,true);
    return CG_p;
}

CGraph MakeSimpleGraph (VertexIdx n, std::set<std::pair<VertexIdx ,VertexIdx >> edge_list) {

    // Construct a simple graph from the edges of the random walk

    Graph G_p;
    G_p.nVertices = n;
    G_p.nEdges = edge_list.size();
    G_p.srcs = new VertexIdx[G_p.nEdges];
    G_p.dsts = new VertexIdx[G_p.nEdges];

    std::vector<VertexIdx > srcs_s, dsts_s;

    // The next step is to extract the elements of edgeList and copy them into srcs and dsts array
    std::transform(std::begin(edge_list),std::end(edge_list),
                   std::back_inserter(srcs_s),
                   [](auto & pair) {return pair.first;});

    std::transform(std::begin(edge_list),std::end(edge_list),
                   std::back_inserter(dsts_s),
                   [](auto & pair) {return pair.second;});


    std::copy(std::begin(srcs_s), std::end(srcs_s), G_p.srcs);
    std::copy(std::begin(dsts_s), std::end(dsts_s), G_p.dsts);

    CGraph CG_p = makeCSR(G_p);
    return CG_p;
}

/** Count the number of triangles incident on the edge e.
 *  The input is given as the Edgeinfo format.
 * @param cg
 * @param e_struct
 * @return
 */
EdgeIdx TriangleByEdge (CGraph* cg, EdgeInfo e_struct) {
    VertexIdx u = e_struct.src;
    VertexIdx v = e_struct.dest;
    // Populate the neighbor list for both u and v
    std::vector<VertexIdx > N_u, N_v;
    N_u.assign(cg->nbors+cg->offsets[u],cg->nbors+cg->offsets[u]+cg->degree(u));
    N_v.assign(cg->nbors+cg->offsets[v],cg->nbors+cg->offsets[v]+cg->degree(v));
    std::sort (N_u.begin(),N_u.end());
    std::sort (N_v.begin(),N_v.end());
    std::vector <VertexIdx > common_neighbors;
    std::set_intersection(N_u.begin(),N_u.end(), N_v.begin(),N_v.end(),std::back_inserter(common_neighbors));
    return common_neighbors.size();
}

/** Count the number of triangles incident on the edge e.
 *  The input is given as the Edgeinfo format.
 * @param cg
 * @param e_struct
 * @return
 */
EdgeIdx TriangleByVertex (CGraph* cg, VertexIdx u) {
    Count v_tri= 0;
    for ( EdgeIdx i=cg->offsets[u]; i < cg->offsets[u+1];i++)
        for ( EdgeIdx j=i+1; j < cg->offsets[u+1];j++) {
            VertexIdx v = cg->nbors[i];
            VertexIdx w = cg->nbors[j];
            if (cg->isEdgeBinary(v,w))
                v_tri ++;
        }
    return v_tri;
}


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
    // Note that the graph will be a simple graph since each edge is considered exactly once.

    CGraph CG_p = MakeMultiGraph(n,srcs,dsts);

    // Count the exact number of triangles in the sparsified graph G_p and scale it up.
    Estimates triangles_in_Gp = CountExactTriangles(&CG_p);
    double triangleEstimate = (triangles_in_Gp.estimate * 1.0)/ (p*p*p);
//    std::cout << nEdges <<"  " << triangleEstimate << std::endl;
//    printf ("%lld,  %lf",nEdges,triangleEstimate);

    Estimates return_estimate = {};
    return_estimate.estimate = triangleEstimate;
    return_estimate.fraction_of_vertices_seen = visited_vertex_set.size() * 100.0 / n;
    return_estimate.fraction_of_edges_seen = CG_p.nEdges * 100.0 / m;

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

    std::unordered_set <VertexIdx > visited_vertex_set; // This will be used to find the number of vertices in G_p
    std::unordered_set <EdgeIdx > visited_edge_set; // This will be used to find the number of vertices in G_p

//    std::set<std::pair <VertexIdx,VertexIdx>> edge_list; // This list will be used to track distinct edges in the random walk.

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

//        edge_list.insert(std::make_pair(src,dst));
//        edge_list.insert(std::make_pair(dst,src));

        visited_vertex_set.insert(src);
        visited_vertex_set.insert(dst);
        visited_edge_set.insert(e);
    }

    // Construct the induced graph G_p
    // Note that G_p is a multi-graph. When counting triangles exactly in G_p, we consider this multiplicity.
    CGraph CG_p = MakeMultiGraph(n,srcs,dsts);

    // Now construct a simple version of this graph by ignoring the multiplicity
//    CGraph CG_p_s = MakeSimpleGraph(n,edge_list);

    // Count the exact number of triangles in the sparsified graph CG_p and scale it up.
    Estimates triangles_in_Gp = CountExactTriangles(&CG_p);
    delCGraph(CG_p);
    double scale = 1.0 / pow(p,3);
    double triangleEstimate = triangles_in_Gp.estimate * scale ;
//    printf ("Multi: %lld,  %lf\n",CG_p.nEdges,triangleEstimate);

    // Count the exact number of triangles in the sparsified simple graph CG_p_s and scale it up.
//    Estimates triangles_in_Gp_s = CountExactTriangles(&CG_p_s);
//    double scale_s = 1.0 / pow(p,3);
//    double triangleEstimate_s = triangles_in_Gp_s.triangle_estimate * scale_s ;
////    printf ("Simple : %lld,  %lf\n",CG_p_s.nEdges,triangleEstimate_s);

    Estimates return_estimate = {};
    return_estimate.estimate = triangleEstimate;
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
    std::unordered_set<EdgeIdx> visited_edge_set;
    std::unordered_set<VertexIdx> visited_vertex_set;
    EdgeIdx nEdges = 0;
    std::vector<VertexIdx > srcs;
    std::vector<VertexIdx > dsts;
//    std::set<std::pair <VertexIdx,VertexIdx>> edgeList; // This list will be used to track distinct edges in the random walk.

    // Initialize a uniform distribution for generating seed vertices
    std::uniform_int_distribution<VertexIdx> dist_seed_vertex(0, n - 1);

    // Perform a random walk for seed_count many times
    for (VertexIdx sC = 0; sC < seed_count; sC++) {
        // Pick a random seed vertex
        //VertexIdx seed = dist_seed_vertex(mt); // TODO: verify randomness
        VertexIdx seed = params.seed_vertices[sC];
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

            //if (wL > 5000) {// ignore the edges till tmix
            srcs.emplace_back(parent);
            dsts.emplace_back(child);
            ++nEdges;
            srcs.emplace_back(child);  // The graph representation (adjacency list) requires both the edges to be present
            dsts.emplace_back(parent);
            ++nEdges;

//            edgeList.insert(std::make_pair(parent, child));
//            edgeList.insert(std::make_pair(child, parent));

            // Collect the distinct edges and distinct vertices visited so far in a set
            visited_edge_set.insert(random_nbor_edge);
            visited_vertex_set.insert(parent);
            //}
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
    // Create a multi-graph using the edges from the random walk
    CGraph CG_p = MakeMultiGraph(n, srcs,dsts);

    // Create a simple graph from the edges in the random walk
//    CGraph CG_p_s = MakeSimpleGraph(n, edgeList);


    //p = 2.0*visited_edge_set.size() / m;
    // Count the exact number of triangles in the sparsified graph CG_p and scale it up.
    Estimates triangles_in_Gp = CountExactTriangles(&CG_p);
    delCGraph(CG_p);
    double scale = 1.0 / pow(p,3);
    double triangleEstimate = triangles_in_Gp.estimate * scale ;
    printf("p=%lf,  scale = %lf, desired scale = %lf \n",p,scale,548658705.0/triangles_in_Gp.estimate *1.0);
//    printf ("Multi: %lld,  %lf\n",CG_p.nEdges,triangleEstimate);

    // Count the exact number of triangles in the sparsified graph CG_p_s and scale it up.
//    Estimates triangles_in_Gp_s = CountExactTriangles(&CG_p_s);
//    double scale_s = 1.0 / pow(p,3);
//    double triangleEstimate_s = triangles_in_Gp_s.triangle_estimate * scale_s ;
//    double another_scale = 1.0* CG_p_s.nEdges / m;
//    printf("p=%lf,  scale = %lf , desired scale = %lf, another scale=%lf\n",
//            p,scale_s,58771288.0/triangles_in_Gp_s.triangle_estimate *1.0, 1.0/pow(another_scale,3));
//    printf ("Simple: %lld,  %lf\n",CG_p_s.nEdges,triangleEstimate_s);


    Estimates return_estimate = {};
    return_estimate.estimate = triangleEstimate;
    return_estimate.fraction_of_vertices_seen = visited_vertex_set.size() * 100.0 / n;
    return_estimate.fraction_of_edges_seen = 2*visited_edge_set.size() * 100.0 / m;
    // Why the multiplication by 2? The fraction of edges seen is effectively
    // compared against all the entries in the adjacency list, which is 2m.

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
Estimates EstTriByEdgeSampleAndCount(CGraph *cg, Parameters params)
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
    EdgeIdx count_match = 0;
    std::vector<VertexIdx > srcs;
    std::vector<VertexIdx > dsts;
    double running_count = 0;

    std::unordered_set <VertexIdx > visited_vertex_set; // This will be used to find the number of vertices in G_p
    std::unordered_set <EdgeIdx > visited_edge_set; // This will be used to find the number of vertices in G_p
    visited_vertex_set.reserve(edges_in_G_p*50);
    visited_edge_set.reserve(edges_in_G_p*50);


    for (EdgeIdx i = 0; i < edges_in_G_p ; i++) {
        EdgeIdx e = unif_rand_edge(mt);
        dst = cg->nbors[e];
        // binary search the array cg->offset to find the index v, such that edges[e] lies between cg->offset[v]
        // and cg->offset[v+1]. We achieve this by calling std::upper_bound
        src = std::upper_bound (cg->offsets, cg->offsets+n, e) - cg->offsets -1;

        // sanity check: the edge {src,dst} must have index edges[e].
        if (cg->getEdgeBinary(src,dst) != e)
            printf("Bug in the code!! Run!!!! \n\n");

        EdgeInfo e_struct;
        e_struct.dest = dst;
        e_struct.src = src;
        e_struct.index = e;

        running_count += TriangleByEdge(cg,e_struct);

        visited_vertex_set.insert(src);
        visited_vertex_set.insert(dst);
        visited_edge_set.insert(e);
        // Now for the vertex with lesser degree, we visit all of its neighbors
        // and all the incident edges on it
        // we copy the lower degree vertex in src.

        if (cg->degree(src) > cg-> degree(dst)) {
            VertexIdx temp = src;
            src = dst;
            dst = temp;
        }
        for (EdgeIdx idx = cg->offsets[src]; idx <cg->offsets[src+1];idx++ ) {
            visited_vertex_set.insert(cg->nbors[idx]);
            visited_edge_set.insert(idx);
            EdgeIdx e1 = cg->getEdgeBinary(cg->nbors[idx],dst);
            EdgeIdx e2 = cg->getEdgeBinary(dst,cg->nbors[idx]);
            visited_edge_set.insert(e1);
            visited_edge_set.insert(e2);
            if (e2!= -1) {
                count_match++;
                visited_edge_set.insert(e2);
            }
        }
    }
//    printf("nEdges = %lld\n",count_match);
    double scale = 1.0*m / (6.0*edges_in_G_p);
    double triangleEstimate = running_count * scale;
//    printf ("Multi: %lld,  %lf\n",edges_in_G_p,triangleEstimate);

    Estimates return_estimate = {};
    return_estimate.estimate = triangleEstimate;
    VertexIdx edges_seen = visited_edge_set.size();
    VertexIdx vertices_seen = visited_vertex_set.size();
    return_estimate.fraction_of_vertices_seen = vertices_seen * 100.0 / n;
    return_estimate.fraction_of_edges_seen = edges_seen * 100.0 / m;
    // Why the multiplication by 2? The fraction of edges seen is effectively
    // compared against all the entries in the adjacency list, which is 2m.

    return return_estimate;
}

Estimates EstTriByRWAndCountPerEdge(CGraph *cg, Parameters params) {

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

    std::unordered_set <VertexIdx > visited_vertex_set; // This will be used to find the number of vertices in G_p
    std::unordered_set <EdgeIdx > visited_edge_set; // This will be used to find the number of vertices in G_p
    visited_vertex_set.reserve(walk_length*50);
    visited_edge_set.reserve(walk_length*50);

    EdgeIdx nEdges = 0;
    std::vector<VertexIdx > srcs;
    std::vector<VertexIdx > dsts;
    double running_count = 0;
    EdgeIdx repeats = 0;
    std::set<std::pair <VertexIdx,VertexIdx>> edgeList; // This list will be used to track distinct edges in the random walk.

    // Initialize a uniform distribution for generating seed vertices
    std::uniform_int_distribution<VertexIdx> dist_seed_vertex(0, n - 1);

    // Perform a random walk for seed_count many times
    for (VertexIdx sC = 0; sC < seed_count; sC++) {
        // Pick a random seed vertex
        //VertexIdx seed = dist_seed_vertex(mt); // TODO: verify randomness
        VertexIdx seed = params.seed_vertices[sC]; // TODO: verify randomness
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

            if (wL >= 0) {// ignore the edges till tmix
                repeats++;
                EdgeInfo e_struct;
                e_struct.index = random_nbor_edge;
                e_struct.src = parent;
                e_struct.dest = child;
                running_count += TriangleByEdge(cg,e_struct);

                visited_vertex_set.insert(parent);
                visited_vertex_set.insert(child);
                visited_edge_set.insert(random_nbor_edge);
                // Now for the vertex with lesser degree, we visit all of its neighbors
                // and all the incident edges on it
                // we copy the lower degree vertex in src.
                VertexIdx src,dst;
                if (cg->degree(parent) > cg-> degree(child)) {
                    src = child;
                    dst = parent;
                }
                else {
                    src = parent;
                    dst = child;
                }
                for (EdgeIdx idx = cg->offsets[src]; idx <cg->offsets[src+1];idx++ ) {
                    visited_vertex_set.insert(cg->nbors[idx]);
                    visited_edge_set.insert(idx);
                    EdgeIdx e1 = cg->getEdgeBinary(dst,cg->nbors[idx]);
                    EdgeIdx e2 = cg->getEdgeBinary(cg->nbors[idx],dst);
                    if (e2!= -1) {
                        visited_edge_set.insert(e2);
                        //visited_edge_set.insert(e1);
                    }
                }

            }
            // The random walk proceeds with the vertex child
            parent = child;
            // If {u,v} is an isolated edge, then we will stuck here: so find a different seed and continue from there
            if (deg_of_child == 1 && deg_of_parent == 1)
                parent = dist_seed_vertex(mt);
        }
    }


    // Count the exact number of triangles in the sparsified graph CG_p and scale it up.
    double scale = 1.0*m / (6.0*repeats);
    double triangleEstimate = running_count * scale ;
    //printf("p=%lf,  scale = %lf \n",p,scale);
//    printf ("Multi: %lld,  %lf\n",visited_edge_set.size(),triangleEstimate);

    Estimates return_estimate = {};
    return_estimate.estimate = triangleEstimate;
    VertexIdx edges_seen = visited_edge_set.size();
    VertexIdx vertices_seen = visited_vertex_set.size();
    //printf("V=%lld,E=%lld",vertices_seen,edges_seen);
    return_estimate.fraction_of_vertices_seen = vertices_seen * 100.0 / n;
    return_estimate.fraction_of_edges_seen = edges_seen * 100.0 / m;
    //    return_estimate.fraction_of_vertices_seen = visited_with_nbor_vertex_set.size() * 100.0 / n;
//    return_estimate.fraction_of_edges_seen = visited_with_nbor_edge_set.size() * 100.0 / m;
    // Why the multiplication by 2? The fraction of edges seen is effectively
    // compared against all the entries in the adjacency list, which is 2m.

    return return_estimate;
}

Estimates NborSampling(CGraph *cg, OrderedEdgeCollection &edge_collection, std::mt19937 mt) {

    // Retrieve relevant information from the edge_collection structure
    std::vector<OrderedEdge> edge_list;
    edge_list = edge_collection.edge_list;
    EdgeIdx no_of_edges = edge_collection.no_of_edges;

    // Set up random number generator
    //std::random_device rd;
    //std::mt19937 mt(rd());

    // Variables used for estimating the triangle count
    EdgeIdx m = cg->nEdges; // TODO note that m is double of the number of edges in the graph
    double X = 0, Y = 0, Z = 0, scaling = m / 6;
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
        OrderedEdge edge = edge_list[i];

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
        //if (cg->isEdgeBinary(w, edge.v) && (deg_of_w > deg_of_v || (deg_of_w == deg_of_v && w > edge.v))) {
        if (cg->isEdgeBinary(w, edge.v)) {
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
    return_estimate.estimate = X;
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

Estimates EstTriByRWandNborSampling(CGraph *cg, Parameters params) {

    Estimates output;
    OrderedEdgeCollection randomEdgeCollection;

    // Set up random number generator
    std::random_device rd;
    // Using this random number generator initializize a PRNG: this PRNG is passed along to
    // draw an element from various distributions
    std::mt19937 mt(rd());

    randomEdgeCollection = GetEdgesByRandomWalk(cg, params, mt);
    output = NborSampling(cg, randomEdgeCollection, mt);
    return output;
}

#endif //SUBGRAPHCOUNT_BASELINEESTIMATORS_H
