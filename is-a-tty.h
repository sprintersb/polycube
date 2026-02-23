#ifndef IS_A_TTY_H
#define IS_A_TTY_H

#include <cstdio>

#ifdef __has_include
#if __has_include (<unistd.h>)
#include <unistd.h>
static inline bool is_a_tty (FILE *stream)
{
    return isatty (fileno (stream));
}
#else
static inline bool is_a_tty (FILE*)
{
    return false;
}
#endif // has <unistd.h> ?
#else
static inline bool is_a_tty (FILE*)
{
    return false;
}
#endif // has __has_include ?

#endif // IS_A_TTY_H
