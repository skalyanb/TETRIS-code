//
// Created by Suman Kalyan Bera on 2020-01-28.
//

#ifndef SUBGRAPHCOUNT_BASELINEUTIL_H
#define SUBGRAPHCOUNT_BASELINEUTIL_H

#include <set>
#include <unordered_map>

#include "../EstimatorUtilStruct.h"
#include "../TriangleEstimators.h"
#include "../EstimatorUtilStruct.h"
#include "RandomWalkUtils.h"

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

CGraph MakeMultiGraph (VertexIdx n, vector<VertexIdx > srcs, vector<VertexIdx > dsts) {

    // First rename the vertices in srcs and dsts to figure out how many vertices are actually
    // and remove those vertice.

    // Create a map that will map every vertex id the in the input file to an integer in the output file
    unordered_map <VertexIdx, VertexIdx> dict;
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

    copy(begin(srcs), end(srcs), G_p.srcs);
    copy(begin(dsts), end(dsts), G_p.dsts);

    // Convert the graph into CSR representation. The second parameters requests for in-pace operation.
    CGraph CG_p = makeCSR(G_p,true);
    return CG_p;
}

#endif //SUBGRAPHCOUNT_BASELINEUTIL_H
