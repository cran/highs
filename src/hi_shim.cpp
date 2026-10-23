/* hi_shim.cpp -- implementation of the CRAN I/O compliance shim declared in
 * hi_shim.h.  Compiled by R CMD INSTALL as part of src/ (with R's headers),
 * NOT by the cmake build; defines HI_SHIM_NO_REMAP so the real libc names
 * remain callable for genuine FILE* streams.
 *
 * Design points:
 *  - The only output sinks are Rvprintf/REvprintf (R_ext/Print.h), the
 *    va_list forms the manual provides for exactly this purpose.
 *  - The R API may be used from the main thread only (R-exts 2591-2596).
 *    The main thread id is captured when this object file is loaded
 *    (package load happens on the R main thread); messages arriving on
 *    other threads are formatted R-free with vsnprintf and queued, then
 *    drained on the next main-thread emission or hi_shim_drain() call.
 *  - exit()/abort() become C++ exceptions: destructors unwind, and Rcpp's
 *    .Call boundary (BEGIN_RCPP/END_RCPP in RcppExports.cpp) converts them
 *    to an R error.  This avoids Rf_error's longjmp over C++ destructors.
 */
#define HI_SHIM_NO_REMAP
#include "hi_shim.h"

#define R_USE_C99_IN_CXX
#include <R_ext/Print.h>

#include <cstring>
#include <mutex>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>
#include <vector>

namespace {

char tag_out, tag_err, tag_in; /* unique addresses; never dereferenced */

/* Captured at package load (dlopen runs on the R main thread). */
const std::thread::id main_thread_id = std::this_thread::get_id();

std::mutex queue_mutex;
std::vector<std::pair<bool, std::string>> queued; /* (is_stderr, text) */

/* Per-thread console mute (see hi_shim_set_quiet in hi_shim.h).  thread_local:
 * no synchronization needed, and the scope is exactly the muting thread --
 * concurrent worker-thread messages keep queueing untouched. */
thread_local bool quiet = false;

bool on_main_thread() { return std::this_thread::get_id() == main_thread_id; }

void drain_locked_free() {
  std::vector<std::pair<bool, std::string>> q;
  {
    std::lock_guard<std::mutex> lock(queue_mutex);
    q.swap(queued);
  }
  for (const auto& m : q) {
    if (m.first)
      REprintf("%s", m.second.c_str());
    else
      Rprintf("%s", m.second.c_str());
  }
}

int emit(bool is_err, const char* fmt, va_list ap) {
  if (on_main_thread()) {
    /* Muted: skip the print AND the drain -- draining would surface queued
     * worker messages in the middle of a deliberately quiet call.  Deferred,
     * not lost: they appear at the next unmuted emission or drain. */
    if (quiet) return 0;
    drain_locked_free();
    if (is_err)
      REvprintf(fmt, ap);
    else
      Rvprintf(fmt, ap);
  } else {
    char buf[4096];
    vsnprintf(buf, sizeof buf, fmt, ap);
    std::lock_guard<std::mutex> lock(queue_mutex);
    queued.emplace_back(is_err, std::string(buf));
  }
  return 0;
}

bool is_sentinel(const FILE* f) {
  return f == (const FILE*)&tag_out || f == (const FILE*)&tag_err ||
         f == (const FILE*)&tag_in;
}

} // namespace

extern "C" {

FILE* hi_shim_stdout(void) { return (FILE*)&tag_out; }
FILE* hi_shim_stderr(void) { return (FILE*)&tag_err; }
FILE* hi_shim_stdin(void)  { return (FILE*)&tag_in; }

void hi_shim_drain(void) {
  if (on_main_thread() && !quiet) drain_locked_free();
}

int hi_shim_set_quiet(int on) {
  int prev = quiet ? 1 : 0;
  quiet = (on != 0);
  return prev;
}

int hi_shim_vfprintf(FILE* f, const char* fmt, va_list ap) {
  if (f == (FILE*)&tag_out) return emit(false, fmt, ap);
  if (f == (FILE*)&tag_err) return emit(true, fmt, ap);
  return vfprintf(f, fmt, ap); /* genuine stream: pass through */
}

int hi_shim_printf(const char* fmt, ...) {
  va_list ap;
  va_start(ap, fmt);
  int r = emit(false, fmt, ap);
  va_end(ap);
  return r;
}

int hi_shim_vprintf(const char* fmt, va_list ap) { return emit(false, fmt, ap); }

int hi_shim_fprintf(FILE* f, const char* fmt, ...) {
  va_list ap;
  va_start(ap, fmt);
  int r = hi_shim_vfprintf(f, fmt, ap);
  va_end(ap);
  return r;
}

int hi_shim_puts(const char* s) { return hi_shim_printf("%s\n", s); }

int hi_shim_putchar(int c) { return hi_shim_printf("%c", c); }

int hi_shim_fputs(const char* s, FILE* f) { return hi_shim_fprintf(f, "%s", s); }

int hi_shim_fputc(int c, FILE* f) { return hi_shim_fprintf(f, "%c", c); }

size_t hi_shim_fwrite(const void* p, size_t sz, size_t n, FILE* f) {
  if (is_sentinel(f)) {
    hi_shim_fprintf(f, "%.*s", (int)(sz * n), (const char*)p);
    return n;
  }
  return fwrite(p, sz, n, f);
}

int hi_shim_fflush(FILE* f) {
  if (f == NULL || is_sentinel(f)) return 0;
  return fflush(f);
}

int hi_shim_fclose(FILE* f) {
  if (is_sentinel(f)) return 0;
  return fclose(f);
}

char* hi_shim_fgets(char* s, int n, FILE* f) {
  if (is_sentinel(f)) return NULL; /* console input: report EOF */
  return fgets(s, n, f);
}

void hi_shim_exit(int status) {
  hi_shim_drain();
  throw std::runtime_error("HiGHS attempted exit(" + std::to_string(status) +
                           ")");
}

void hi_shim_abort(void) {
  hi_shim_drain();
  throw std::runtime_error(
      "HiGHS internal error (abort/highsAssert reached)");
}

} // extern "C"

namespace {
struct hi_shim_streambuf : std::streambuf {
  std::streamsize xsputn(const char* s, std::streamsize n) override {
    hi_shim_printf("%.*s", (int)n, s);
    return n;
  }
  int overflow(int c) override {
    if (c != EOF) hi_shim_printf("%c", (char)c);
    return c;
  }
};
} // namespace

std::ostream& hi_shim_cout(void) {
  static hi_shim_streambuf buf;
  static std::ostream os(&buf);
  return os;
}
