#include "diagnostic.h"

#include <cstdarg>
#include <cstdlib>
#include <cstdio>
#include <cstring>

int (*_diagnostic_vfprintf) (FILE*, const char*, va_list);

bool info_function;
bool out_context;

enum { x_out, x_info, x_warning, x_error, x_fatal, x_None };

static constexpr int exit_code[] =
{
    [x_out] = EXIT_SUCCESS,
    [x_info]  = EXIT_SUCCESS,
    [x_warning] = EXIT_SUCCESS,
    [x_error] = -2,
    [x_fatal] = -42,
    [x_None] = 0,
};

static FILE* const stream[] =
{
    [x_out] = stdout,
    [x_info]  = stderr,
    [x_warning] = stderr,
    [x_error] = stderr,
    [x_fatal] = stderr,
    [x_None] = nullptr
};

static int _cont = x_None;

static void
diagnostic (int what, const char *name, const char *file, int line,
            const char *fmt, va_list args)
{
    FILE *s = stream[what];
    fflush (stdout);
    fflush (stderr);

    auto vfp = _diagnostic_vfprintf ? _diagnostic_vfprintf : vfprintf;

    if (_cont != what)
        if (out_context || what != x_out)
            fprintf (s, "%s%s:%d: %s: ", _cont != x_None ? "\n" : "",
                     file, line, name);

    size_t len = std::strlen (fmt);
    _cont = what != x_out && len && fmt[len - 1] == ' ' ? what : x_None;

    vfp (s, fmt, args);
    if (what != x_out && _cont == x_None)
        fprintf (s, "\n");
    fflush (s);

#ifdef __SANITIZE_ADDRESS__
    if (what == x_fatal)
    {
        // Raise a signal, which is nice with, say -fsanitize=address to
        // show a stack backtrace.
        fprintf (s, "%s:%d: raising SIGNAL for address sanitizer...\n",
                 __FILE__, __LINE__);
        fflush(s);
        volatile int * volatile p = 0;
        *p = 0;
    }
#endif // -fsanitize=address

    if (exit_code[what] != EXIT_SUCCESS)
        std::exit (exit_code[what]);
}

void out_at (const char *file, int line, const char *func,
             const char *fmt, ...)
{
    va_list args;
    va_start (args, fmt);
    static char name[300];

    *name = '\0';
    if (out_context)
        sprintf (name, "<%s>", func);

    diagnostic (x_out, name, file, line, fmt, args);

    va_end (args);
}

void info_at (const char *file, int line, const char *func,
              const char *fmt, ...)
{
    va_list args;
    va_start (args, fmt);
    static char name[300];

    if (_cont != x_info)
    {
        if (info_function)
            sprintf (name, "info <%s>", func);
        else
            sprintf (name, "info");
    }

    diagnostic (x_info, name, file, line, fmt, args);

    va_end (args);
}

void warning_at (const char *file, int line, const char *fmt, ...)
{
    va_list args;
    va_start (args, fmt);

    diagnostic (x_warning, "warning", file, line, fmt, args);

    va_end (args);
}

void error_at (const char *file, int line, const char *fmt, ...)
{
    va_list args;
    va_start (args, fmt);

    diagnostic (x_error, "error", file, line, fmt, args);

    va_end (args);

    std::exit (-2); // Get rid of "noreturn function does return".
}

void fatal_at (const char *file, int line, const char *tag, const char *fmt,...)
{
    va_list args;
    va_start (args, fmt);

    diagnostic (x_fatal, tag, file, line, fmt, args);

    va_end (args);

    std::exit (-42); // Get rid of "noreturn function does return".
}
