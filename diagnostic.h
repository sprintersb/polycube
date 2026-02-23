#ifndef DIAGNOSTIC_H
#define DIAGNOSTIC_H

#include <stdio.h>
#include <stdarg.h>

/* Set to a function like gmp_vfprintf from <gmp.h> to support "%Zd"
   to print mpz_t.  Or set to mpfr_vfprintf from <mpfr.h> to support "%RNf"
   to print mpfr_t.  Notice that gmp.h and mpfr.h declare these functions
   if <stdarg.h> or <cstdarg> was included prior to <gmp.h> resp. <mpfr.h>.
   Note that mpfr_printf also supports the print modifiers from gmp_printf.
   If not set or nullptr, use vfprintf. */
extern int (*_diagnostic_vfprintf) (FILE*, const char*, va_list);

#ifdef DIAGNOSTIC_NO_FORMAT_CHECK
[[noreturn]]
extern void fatal_at (const char*, int line, const char*, const char*, ...);
[[noreturn]]
extern void error_at (const char *file, int line, const char *fmt, ...);
extern void warning_at (const char *file, int line, const char *fmt, ...);
extern void info_at (const char *file, int line, const char *func,
                     const char *fmt, ...);
extern void out_at (const char *file, int line, const char *func,
                    const char *fmt, ...);
#else
[[noreturn]] [[gnu::format(__printf__,4,5)]]
extern void fatal_at (const char*, int line, const char*, const char*, ...);
[[noreturn]] [[gnu::format(__printf__,3,4)]]
extern void error_at (const char *file, int line, const char *fmt, ...);
[[gnu::format(__printf__,3,4)]]
extern void warning_at (const char *file, int line, const char *fmt, ...);
[[gnu::format(__printf__,4,5)]]
extern void info_at (const char *file, int line, const char *func,
                     const char *fmt, ...);
[[gnu::format(__printf__,4,5)]]
extern void out_at (const char *file, int line, const char *func,
                    const char *fmt, ...);
#endif /* DIAGNOSTIC_NO_FORMAT */

extern bool info_function;
extern bool out_context;

#define error(x...) \
    error_at (__FILE__, __LINE__, ##x)

#define fatal(x...) \
    fatal_at (__FILE__, __LINE__, "fatal", ##x)

#define unreachable(x...) \
    fatal_at (__FILE__, __LINE__, "unreachable code", ##x)

#define warning(x...) \
    warning_at (__FILE__, __LINE__, ##x)

#define info(active, x...)                                             \
    do {                                                               \
        if ((active))                                                  \
            info_at (__FILE__, __LINE__, __PRETTY_FUNCTION__, ##x);    \
    } while (0)

#define out(x...)                                               \
    do {                                                        \
        out_at (__FILE__, __LINE__, __PRETTY_FUNCTION__, ##x);  \
    } while (0)

#endif /* DIAGNOSTIC_H */
