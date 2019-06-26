//
// Created by Suman Kalyan Bera on 2019-05-30.
//

#ifndef SUBGRAPHCOUNT_ERRORCODE_H
#define SUBGRAPHCOUNT_ERRORCODE_H

namespace Escape
{
    enum ErrorCode : int
    {
        ecNone
        , ecInvalidInput
        , ecSystemError
        , ecUnsupportedFormat
        , ecIOError
    };
}
#endif //SUBGRAPHCOUNT_ERRORCODE_H
