//
// Created by Suman Kalyan Bera on 2020-01-28.
//

#ifndef SUBGRAPHCOUNT_TRIANGLLECOUNTUTIL_H
#define SUBGRAPHCOUNT_TRIANGLLECOUNTUTIL_H

#include <set>
#include <unordered_map>

#include "../EstimatorUtilStruct.h"
#include "../TriangleEstimators.h"

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

#endif //SUBGRAPHCOUNT_TRIANGLLECOUNTUTIL_H
