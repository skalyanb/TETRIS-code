//
// Created by Suman Kalyan Bera on 2019-06-17.
//

#ifndef SUBGRAPHCOUNT_DIGRAPH_H
#define SUBGRAPHCOUNT_DIGRAPH_H

#include "ErrorCode.h"
#include "Graph.h"
#include <algorithm>

using namespace Escape;

// DAG structure has two pointers, one to the
// adjacency list of outedge, one to the adjacency list of inedges
struct CDAG
{
    CGraph outlist;
    CGraph inlist;
};

// Structure for comparing nodes according to their degree.
// So u < v if degree of u less than that of v in graph g.

struct DegreeComp
{
    CGraph *g;
    DegreeComp(CGraph *g) { this->g = g;}

    bool operator () (VertexIdx u, VertexIdx v)
    {
        if (u<0 || v<0) // something wrong, so print error message, and exit
        {
            printf("Something wrong in DegreeComp, negative vertices. u is %lld, v is %lld\n",u,v);
            exit(EXIT_FAILURE);
        }
        VertexIdx degu = g->offsets[u+1] - g->offsets[u];  // Degree of u
        VertexIdx degv = g->offsets[v+1] - g->offsets[v];  // Degree of v

        if (degu < degv || (degu == degv && u < v))    // Comparing degrees and breaking ties by id
            return true;
        else
            return false;
    }
};

// degDist: A convenient function to get the degree distribution of graph, as an array
// Input:
//     g: A pointer to CGraph for which degree distribution needs to be computed
//     degdistarray: A pointer to array where degree distribution should be output. This is where output will be placed.
//     It must have length at least g->nVertices+1.
// Output:
//     Count: the maximum degree, so only needs to read degdistarray up to that point. But the array degdistarray
//            has the complete degree distribution, padded appropriately with zeros.
//

Count degDist(CGraph* g, VertexIdx* degdistarray)
{
    VertexIdx deg, maxdeg = 0; // variable for maximum degree
    VertexIdx i; // index for looping

    for(i=0; i<g->nVertices; ++i) //looping to initialize degdistarray to all zeroes
        degdistarray[i] = 0;

    for(i=0; i<g->nVertices; ++i)
    {
        deg = g->offsets[i+1] - g->offsets[i]; // computing degree of vertex i
        if (deg > maxdeg)
            maxdeg = deg; //update maximum degree

        degdistarray[deg]++; //increment number of degree deg vertices
    }

    return maxdeg; //return maximum degree
}

//Construct DAG based on degree ordering
//
// Input: Pointer for CGraph g
// Output: CDAG for degree ordering in g
//         This is the CDAG for the DAG where each edge points from lower degree endpoint to higher degree endpoint.
//
//
//         The outlist in CDAG is guaranteed to be sorted by degrees. This means that the neighbors
//         of every vertex in the outlist are sorted by their degrees in g. This is quite useful in
//         further processing.
CDAG degreeOrdered(CGraph *g)
{
    CDAG ret;     // CDAG to be returned
    CGraph outdag = {g->nVertices, 0, new EdgeIdx[g->nVertices+1], new VertexIdx[g->nEdges+1]};  // Initialize DAG of out-edges
    CGraph indag = {g->nVertices, 0, new EdgeIdx[g->nVertices+1], new VertexIdx[g->nEdges+1]};   // Initialize DAG of in-edges
    EdgeIdx outcur = 0;
    EdgeIdx incur = 0;
    VertexIdx dest;
    VertexIdx degi;
    VertexIdx degdest;

    outdag.offsets[0] = 0;
    indag.offsets[0] = 0;
    for (VertexIdx i=0; i < g->nVertices; ++i)   // Looping over all vertices in g
    {
        for (EdgeIdx j = g->offsets[i]; j < g->offsets[i+1]; ++j)   // Looping over neighbors of i in g
        {
            dest = g->nbors[j];     // We are now looking at edge (i,dest)
            degi = g->offsets[i+1] - g->offsets[i];   // Degree of i
            degdest = g->offsets[dest+1]- g->offsets[dest];   // Degree of dest
            //printf("i=%lld dest=%lld degi=%lld degdest=%lld\n",i,dest,degi,degdest);

            //We now orient the edge depending of degi vs degdest.
            // We break ties according to vertex id.
            // In the output, the g-edge (i,dest) is either pointing to dest (in if condition) or pointing to i (in else condition).

            if (degi < degdest || (degi == degdest && i < dest))
            {
                outdag.nbors[outcur] = dest;   // We want point edge from i to dest. So this directed edge is added to outdag.
                ++outcur;                      // Increment pointer in outdag.nbors and the number of edges in outdag.
                ++outdag.nEdges;
            }
            else
            {
                indag.nbors[incur] = dest;     // We point edge from dest to i. So this edge goes into indag.
                ++incur;                       // Pointer and number of edges incremented
                ++indag.nEdges;
            }
        }
        outdag.offsets[i+1] = outcur;         // We have finished all edges incident to i, so we can update offsets in DAGs.
        indag.offsets[i+1] = incur;
    }

    for (VertexIdx i=0; i < g->nVertices;++i)  // Loops over vertices
        std::sort(outdag.nbors+outdag.offsets[i], outdag.nbors+outdag.offsets[i+1], DegreeComp(g)); // In outdag, sort all neighbors of i according to their degree. Note that DegreeComp gives the desired comparator.

    ret.outlist = outdag;
    ret.inlist = indag;

    return ret;
}


// degenOrdered: This function produces the degeneracy ordered DAG.
// Input:
//        g: A pointer to a CGraph, for which we we desire the corresponding DAG
// Output:
//       CDAG: This stores the degeneracy ordered DAG for g.
//
// The algorithm used is the standard Matula-Beck algorithm. We iteratively remove the lowest degree
// vertex, to get the removal times of each vertex. Then, we orient all the edges according to removal times.
//

CDAG degenOrdered(CGraph* g)
{
    CDAG ret;     // CDAG to be returned
    CGraph outdag = {g->nVertices, 0, new EdgeIdx[g->nVertices+1], new VertexIdx[g->nEdges+1]};  // Initialize DAG of out-edges
    CGraph indag = {g->nVertices, 0, new EdgeIdx[g->nVertices+1], new VertexIdx[g->nEdges+1]};   // Initialize DAG of in-edges
    EdgeIdx outcur = 0;
    EdgeIdx incur = 0;
    VertexIdx *removal_time = new VertexIdx[g->nVertices]; // array of removal times
    VertexIdx *verts_by_deg = new VertexIdx[g->nVertices]; // array of vertices, sorted by current degree
    VertexIdx *degrees = new VertexIdx[g->nVertices]; // array of degrees
    VertexIdx *vert_ptrs = new VertexIdx[g->nVertices]; // array of vertex pointers, that point into verts_by_deg
    VertexIdx *deg_ptrs = new VertexIdx[g->nVertices+1]; // array of degree pointers, pointing into verts_by_deg. deg_ptrs[deg] will point to the first vertex in vert_by_deg with degree deg
    bool *is_removed = new bool[g->nVertices]; // to keep track of whether vertex was removed

    //first, set up a sorted list of vertices by degree
    for(VertexIdx i=0; i < g->nVertices; i++) // looping over vertices
    {
        verts_by_deg[i] = i; //initializing verts_by_deg to store all vertex labels
        is_removed[i] = false; //intiailizing is_removed to all false
        degrees[i] = g->offsets[i+1] - g->offsets[i]; // storing the degree of vertex i
    }

    std::sort(verts_by_deg, verts_by_deg + g->nVertices, DegreeComp(g)); // sort all the vertices by degree
    // set up vertex pointers
    for (VertexIdx i=0; i < g->nVertices; i++) //looping over verts_by_deg
        vert_ptrs[verts_by_deg[i]] = i; // the vertex verts_by_deg[i] is a position i

    //now set up degree pointers
    VertexIdx current_vert = verts_by_deg[0]; //initialize current vertex to smallest degree vertex
    VertexIdx min_deg = g->offsets[current_vert+1] - g->offsets[current_vert]; //initialize running degree to lowest degree
    for (VertexIdx deg=0; deg < g->nVertices; deg++) //looping over degrees
        deg_ptrs[deg] = -1; // initialize these pointers to -1
    deg_ptrs[min_deg] = 0; // initialize the pointer for minimum degree to beginning on array

    VertexIdx current_deg, running_deg; // we will need this variables while looping
    running_deg = min_deg; //initialize the running degree to the minimum degree
    for(VertexIdx i=0; i < g->nVertices; i++) //looping over vertices
    {
        current_vert = verts_by_deg[i]; // the ith vertex in sorted list by degree
        current_deg = g->offsets[current_vert+1] - g->offsets[current_vert]; // get degree of current vertex
        if (current_deg == running_deg) // the current degree is same as running degree, so we move to next vertex
            continue;
        // current_deg is larger, so we need to set pointer in deg_ptrs
        deg_ptrs[current_deg] = i; // the first vertex with current_deg is at location i
        running_deg = current_deg; //update the running degree
    }

    // at this point, we have set up verts_by_deg and deg_ptrs. Now we begin the actual removal of vertices
    VertexIdx to_remove, nbr, swap, next; // convenient variables used in next loop
    for(VertexIdx time=0; time < g->nVertices; time++) // a loop to remove all vertices, with loop variable being the "time"
    {
        to_remove = verts_by_deg[deg_ptrs[min_deg]]; // remove the vertex of minimum degree
        removal_time[to_remove] = time; // set removal time of vertex
        is_removed[to_remove] = true; // remove this vertex

//         printf("Removed %lld, min deg = %lld, position = %lld\n",to_remove,min_deg,deg_ptrs[min_deg]); //(for debug)
//         if (time == g->nVertices-1) // last vertex removed
//             break; // just exit. This is to prevent possible pointer/off-by-1 issues in the remaining code

        next = verts_by_deg[deg_ptrs[min_deg]+1]; // the next vertex in the current ordering of degree
        verts_by_deg[deg_ptrs[min_deg]] = -1; // remove vertex from array
        if (degrees[next] != degrees[to_remove]) // next has higher degree
        {
            deg_ptrs[min_deg] = -1; // there is no vertex of degree min_deg
            min_deg = degrees[next]; //update minimum degree
        }
        else //next vertex has same degree, so only need to increment pointer
            deg_ptrs[min_deg]++; // increment degree pointer for minimum degree, since first vertex is removed

        for (VertexIdx j = g->offsets[to_remove]; j < g->offsets[to_remove+1]; j++) // loop over neighbors of removed vertex
        {
            nbr = g->nbors[j]; // looking at the neighbor
            if (is_removed[nbr]) // if neighbor is removed, move to next neighbor
                continue;
//             for (VertexIdx ind=0; ind < g->nVertices; ind++) FOR DEBUG
//                 printf("%lld ",verts_by_deg[ind]);
//             printf("Moving %lld from %lld to %lld\n",nbr,vert_ptrs[nbr],deg_ptrs[degrees[nbr]]);

            swap = verts_by_deg[deg_ptrs[degrees[nbr]]]; //vertex to be swapped with nbr
            if (swap == -1) // some problem, swapping with deleted vertex
            {
                printf("Error in degeneracy removal. Swapping with removed vertex\n nbr = %lld, deg_ptrs[degrees[nbr]] = %lld",nbr,deg_ptrs[degrees[nbr]]);
                exit(EXIT_FAILURE);
            }
            verts_by_deg[deg_ptrs[degrees[nbr]]] = nbr; //place nbr at beginning of list of vertices of its degree
            verts_by_deg[vert_ptrs[nbr]] = swap; //place swapped vertex at location of nbr

            vert_ptrs[swap] = vert_ptrs[nbr]; // update position of swapped vertex
            vert_ptrs[nbr] = deg_ptrs[degrees[nbr]]; // update position of nbr

            if (degrees[verts_by_deg[deg_ptrs[degrees[nbr]]+1]] == degrees[nbr]) // if next vertex in verts_by_deg after nbr has same degree
                deg_ptrs[degrees[nbr]]++; //increment degree pointer for vertices of degree of nbr, since nbr's degree will be decremented
            else // next vertex has higher degree
                deg_ptrs[degrees[nbr]] = -1; // no more vertices of that degree
            degrees[nbr]--; //decrement degree of neighbor
            if (deg_ptrs[degrees[nbr]] == -1) // if there was no vertex with the same degree as nbr
                deg_ptrs[degrees[nbr]] = vert_ptrs[nbr]; // set the degree pointer to current position of nbr

            if (min_deg > degrees[nbr]) // update minimum degree, if needed
                min_deg = degrees[nbr];
        }
    }

    //sanity check; make sure all vertices are removed
    for (VertexIdx i=0; i < g->nVertices; i++)
    {
        if (!is_removed[i]) // i was not removed
        {
            printf("Error in degeneracy removal. Vertex %lld not removed\n",i);
            exit(EXIT_FAILURE);
        }
    }

    //now we build the DAG
    outdag.offsets[0] = 0;
    indag.offsets[0] = 0;
    VertexIdx dest;
    outcur = 0;
    incur = 0;

    for (VertexIdx i=0; i < g->nVertices; ++i)   // Looping over all vertices in g
    {
        for (EdgeIdx j = g->offsets[i]; j < g->offsets[i+1]; ++j)   // Looping over neighbors of i in g
        {
            dest = g->nbors[j];     // We are now looking at edge (i,dest)

            //We now orient the edge depending of removal times
            // In the output, the g-edge (i,dest) is either pointing to dest (in if condition) or pointing to i (in else condition).

            if (removal_time[i] < removal_time[dest]) // i was removed before dest
            {
                outdag.nbors[outcur] = dest;   // We want point edge from i to dest. So this directed edge is added to outdag.
                ++outcur;                      // Increment pointer in outdag.nbors and the number of edges in outdag.
                ++outdag.nEdges;
            }
            else
            {
                indag.nbors[incur] = dest;     // We point edge from dest to i. So this edge goes into indag.
                ++incur;                       // Pointer and number of edges incremented
                ++indag.nEdges;
            }
        }
        outdag.offsets[i+1] = outcur;         // We have finished all edges incident to i, so we can update offsets in DAGs.
        indag.offsets[i+1] = incur;
    }

    ret.outlist = outdag;
    ret.inlist = indag;

    // free memory
    delete[] removal_time;
    delete[] verts_by_deg;
    delete[] degrees;
    delete[] vert_ptrs;
    delete[] deg_ptrs;
    delete[] is_removed;

    return ret;
}

#endif //SUBGRAPHCOUNT_DIGRAPH_H
