#include "include/GraphIO.h"


using namespace Escape;


static bool isBlankLine(const char *line) {
    while (*line) {
        if (!isspace(*line))
            return false;
        ++line;
    }
    return true;
}


static ErrorCode loadGraph_Escape(const char *path, Graph &graph, int undirected) {
    FILE *f = fopen(path, "r");
    if (!f) {
        fprintf(stderr, "could not open file %s\n", path);
        return ecInvalidInput;
    }

    graph.nVertices = 0;
    EdgeIdx iEdge = 0;
    char line[1024];
    while (fgets(line, sizeof(line), f)) {
        //Ignore comment lines.
        if (line[0] == '#')
            continue;

        if (isBlankLine(line))
            continue;

        int64_t i1, i2;
        sscanf(line, "%lld%lld", &i1, &i2);
        if (graph.nVertices == 0) {
            graph.nVertices = i1;
            if (undirected)
                graph.nEdges = 2 * i2;
            else
                graph.nEdges = i2;
            graph.srcs = new VertexIdx[graph.nEdges];
            graph.dsts = new VertexIdx[graph.nEdges];
        } else {
            graph.srcs[iEdge] = i1;
            graph.dsts[iEdge] = i2;
            ++iEdge;
            if (undirected) {
                graph.srcs[iEdge] = i2;
                graph.dsts[iEdge] = i1;
                ++iEdge;
            }
        }
    }
    fclose(f);

    if (iEdge < graph.nEdges) {
        fprintf(stderr, "expected %lld edges, only got %lld\n", graph.nEdges, iEdge);
        return ecIOError;
    }

    return ecNone;
}

ErrorCode Escape::loadGraphCSR(const char *path, CGraph &cg, int undirected)
{
    VertexIdx nVertices;
    EdgeIdx nEdges;
    auto in_file = std::ifstream(path, std::ios::in | std::ios::binary);
    if(!in_file) {
        std::cout << "Cannot open file " << path << std::endl;
        return ecIOError;
    }
    in_file.read(reinterpret_cast<char*> (&nVertices), sizeof(nVertices));
    in_file.read(reinterpret_cast<char*>(&nEdges), sizeof(nEdges));

    cg.nVertices = nVertices;
    cg.nEdges = nEdges;

    cg.offsets = new EdgeIdx[cg.nVertices+1];
    cg.nbors = new VertexIdx[cg.nEdges];

    in_file.read(reinterpret_cast<char*>(cg.offsets), (nVertices+1)* sizeof(VertexIdx));
    in_file.read(reinterpret_cast<char*>(cg.nbors), nEdges*sizeof(EdgeIdx));
    in_file.close();
    if(!in_file.good()) {
        std::cout << "Error occurred at writing time!" << std::endl;
    }
    return ecNone;
}

ErrorCode Escape::loadGraph(const char *path, Graph &graph, int undirected, IOFormat fmt) {
    switch (fmt) {
        case IOFormat::escape:
            return loadGraph_Escape(path, graph, undirected);
        default:
            return ecUnsupportedFormat;
    }
}
