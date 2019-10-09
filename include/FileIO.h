//
// Created by Suman Kalyan Bera on 2019-10-08.
//

#ifndef SUBGRAPHCOUNT_FILEIO_H
#define SUBGRAPHCOUNT_FILEIO_H

#include <algorithm>
#include <iostream>
#include <cstring>

// for mmap:
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

#include "GraphIO.h"
#include "Graph.h"

using namespace Escape;

const char* map_file(const char* fname, size_t& length);

static ErrorCode loadGraph_Fast(const char *path, Graph &graph) {
    size_t length;
    auto f = map_file(path, length);
    auto l = f + length;

    uintmax_t m_numLines = 0;

}
int main()
{
    size_t length;
    auto f = map_file("test.cpp", length);
    auto l = f + length;

    uintmax_t m_numLines = 0;
    while (f && f!=l)
        if ((f = static_cast<const char*>(memchr(f, '\n', l-f))))
            m_numLines++, f++;

    std::cout << "m_numLines = " << m_numLines << "\n";
}

void handle_error(const char* msg) {
    perror(msg);
    exit(255);
}

const char* map_file(const char* fname, size_t& length)
{
    int fd = open(fname, O_RDONLY);
    if (fd == -1)
        handle_error("open");

    // obtain file size
    struct stat sb;
    if (fstat(fd, &sb) == -1)
        handle_error("fstat");

    length = sb.st_size;

    const char* addr = static_cast<const char*>(mmap(NULL, length, PROT_READ, MAP_PRIVATE, fd, 0u));
    if (addr == MAP_FAILED)
        handle_error("mmap");

    close(fd);
    // TODO close fd at some point in time, call munmap(...)
    return addr;
}

#endif //SUBGRAPHCOUNT_FILEIO_H
