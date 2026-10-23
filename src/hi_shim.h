/* hi_shim.h -- CRAN I/O compliance shim for the bundled HiGHS sources.
 *
 * Force-included (via -include, see inst/build_highs.sh) into EVERY
 * translation unit of the cmake-built HiGHS static library, C and C++
 * alike.  It remaps, at preprocessing time, every token whose compiled
 * symbol CRAN's check flags (tools:::so_symbol_names_table: printf, puts,
 * putchar, vprintf, abort, exit, _exit/_Exit, stdout, stderr, ...) to
 * package-local hi_shim_* functions implemented in src/hi_shim.cpp over
 * R's sanctioned output/error API (Rvprintf/REvprintf; Writing R
 * Extensions, "Printing").  Because the remap happens before compilation,
 * the banned symbols are never generated -- which is exactly what the
 * symbol-table scan of the installed .so looks for.
 *
 * Function-like macros expand only at `name(` call sites, so identifiers
 * that merely contain a banned name (gk_errexit, atexit, `bool abort`)
 * are untouched.
 *
 * src/hi_shim.cpp defines HI_SHIM_NO_REMAP before including this header:
 * it gets the declarations but keeps the real libc names, since it must
 * call through to them for genuine (non-console) FILE* streams.
 *
 * This header is NEVER included by any vendored HiGHS file.
 */
#ifndef HI_SHIM_H
#define HI_SHIM_H

/* Real headers FIRST, so their own declarations are processed before any
 * remapping and later (guarded) re-inclusion inside HiGHS TUs is a no-op.
 * The full iostream family is pre-included because the C++ standard library
 * headers themselves call std::fwrite/std::fflush/... internally (e.g.
 * libc++ <fstream>); if a vendored TU included them after this shim, those
 * qualified calls would hit the macros and fail to compile. */
#ifdef __cplusplus
  #include <cstdio>
  #include <cstdarg>
  #include <cstdlib>
  #include <iostream>
  #include <fstream>
  #include <sstream>
#else
  #include <stdio.h>
  #include <stdarg.h>
  #include <stdlib.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* Sentinel console streams: unique non-NULL FILE* values, never
 * dereferenced.  Pointer identity is preserved, so vendored logic like
 * `if (file != stdout) fclose(file);` keeps its exact meaning. */
FILE* hi_shim_stdout(void);
FILE* hi_shim_stderr(void);
FILE* hi_shim_stdin(void);

int   hi_shim_printf(const char* fmt, ...);
int   hi_shim_vprintf(const char* fmt, va_list ap);
int   hi_shim_fprintf(FILE* f, const char* fmt, ...);
int   hi_shim_vfprintf(FILE* f, const char* fmt, va_list ap);
int   hi_shim_puts(const char* s);
int   hi_shim_putchar(int c);
int   hi_shim_fputs(const char* s, FILE* f);
int   hi_shim_fputc(int c, FILE* f);
size_t hi_shim_fwrite(const void* p, size_t sz, size_t n, FILE* f);
int   hi_shim_fflush(FILE* f);
int   hi_shim_fclose(FILE* f);
char* hi_shim_fgets(char* s, int n, FILE* f);
/* noreturn (they throw): without this, replacing libc's noreturn exit()/
 * abort() makes gcc emit -Wreturn-type at vendored call sites where the
 * original terminator was the last statement (observed on win-builder:
 * HighsCliqueTable.cpp separateCliques, whose body starts with exit(0)). */
#if defined(__cplusplus)
#  define HI_SHIM_NORETURN [[noreturn]]
#elif defined(__GNUC__)
#  define HI_SHIM_NORETURN __attribute__((noreturn))
#else
#  define HI_SHIM_NORETURN
#endif
HI_SHIM_NORETURN void hi_shim_exit(int status);
HI_SHIM_NORETURN void hi_shim_abort(void);
/* Flush queued worker-thread messages; called from the package interface
 * when control returns to R (R API is main-thread-only, R-exts 2591-2596). */
void  hi_shim_drain(void);
/* Thread-local console mute for the CALLING thread only: set on = 1 before
 * operations whose construction unavoidably emits (the Highs constructor
 * banner), restore afterwards.  Returns the previous state (for nesting /
 * RAII).  Worker-thread messages are unaffected -- they queue as usual and
 * appear at the next unmuted main-thread emission or hi_shim_drain(). */
int   hi_shim_set_quiet(int on);

#ifdef __cplusplus
} /* extern "C" */

/* Replacement target for `std::cout` sites, rewritten by inst/patch_min.R
 * (macros cannot reach qualified names).  Backed by a streambuf over
 * hi_shim_printf -- no std::cout/std::cerr symbol is ever referenced. */
std::ostream& hi_shim_cout(void);
#endif

#ifndef HI_SHIM_NO_REMAP

/* stdout/stderr/stdin are macros on some platforms (e.g. Windows UCRT
 * expands stdout to (__acrt_iob_func(1))) -- #undef first. */
#undef  stdout
#undef  stderr
#undef  stdin
#define stdout (hi_shim_stdout())
#define stderr (hi_shim_stderr())
#define stdin  (hi_shim_stdin())

#undef  printf
#undef  vprintf
#undef  fprintf
#undef  vfprintf
#undef  puts
#undef  putchar
#undef  fputs
#undef  fputc
#undef  fwrite
#undef  fflush
#undef  fclose
#undef  fgets
#define printf(...)          hi_shim_printf(__VA_ARGS__)
#define vprintf(fmt, ap)     hi_shim_vprintf(fmt, ap)
#define fprintf(...)         hi_shim_fprintf(__VA_ARGS__)
#define vfprintf(f, fmt, ap) hi_shim_vfprintf(f, fmt, ap)
#define puts(s)              hi_shim_puts(s)
#define putchar(c)           hi_shim_putchar(c)
#define fputs(s, f)          hi_shim_fputs(s, f)
#define fputc(c, f)          hi_shim_fputc(c, f)
#define fwrite(p, sz, n, f)  hi_shim_fwrite(p, sz, n, f)
#define fflush(f)            hi_shim_fflush(f)
#define fclose(f)            hi_shim_fclose(f)
#define fgets(s, n, f)       hi_shim_fgets(s, n, f)
#define exit(status)         hi_shim_exit(status)
#define abort()              hi_shim_abort()

#endif /* HI_SHIM_NO_REMAP */
#endif /* HI_SHIM_H */
