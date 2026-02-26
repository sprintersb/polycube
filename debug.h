// -*- c++ -*-
#ifndef DEBUG_H
#define DEBUG_H

#include "diagnostic.h"
#include "config-debug.h" // auto-generated

#ifdef DEBUG
#define Assert(cond, fmt, x...)                                         \
    do {                                                                \
        if (! (cond))                                                   \
            fatal_at (__FILE__, __LINE__,                               \
                      "Assert \"" #cond "\" failed", fmt, ##x);         \
    } while (0)
#else
#define Assert(...) \
    do {} while (0)
#endif // DEBUG

#endif // DEBUG_H
