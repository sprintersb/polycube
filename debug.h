// -*- c++ -*-
#ifndef DEBUG_H
#define DEBUG_H

#include "diagnostic.h"

#ifdef DEBUG
#define Assert(cond, fmt, x...) \
    do { if (! (cond)) fatal (fmt, ##x); } while (0)
#else
#define Assert(...) \
    do {} while (0)
#endif // DEBUG

#define Unreachable(fmt, x...) fatal ("unreachable code: " fmt, ##x)

#endif // DEBUG_H
