/* A minimal consumer of inst/include/highs_ccallable.h, compiled as C by
 * inst/tinytest/test_ccallable.R (R_HIGHS_LINK_TEST=1) with R CMD SHLIB and
 * nothing but -I<highs>/include: no HiGHS header, no HiGHS library.
 *
 * Solves  min x + y  s.t.  x + y >= 1, 0 <= x, y <= 1  with a logging
 * callback attached (output_flag on, log_to_console off, so the log goes to
 * the callback only) and returns
 * list(version, sizeof_int, status_pass, status_run, model_status,
 *      objective, col_value, callback_calls, callback_item_ok).
 */

#include <highs_ccallable.h>

typedef struct {
  int calls;
  int item_ok;
} callback_state;

/* Counts logging callbacks; item_ok records whether the C API's accessor for
 * the log_type item returns the address of the struct member, which is what
 * the library does when its HighsCallbackDataOut layout matches ours. */
static void logging_callback(int callback_type, const char* message,
                             const HighsCallbackDataOut* data_out,
                             HighsCallbackDataIn* data_in, void* user_data) {
  callback_state* state = (callback_state*) user_data;
  (void) message;
  (void) data_in;
  if (callback_type != kHighsCallbackLogging) return;
  state->calls++;
  if (Highs_getCallbackDataOutItem(data_out, kHighsCallbackDataOutLogTypeName) ==
      (const void*) &data_out->log_type)
    state->item_ok = 1;
}

SEXP ccallable_consumer(void) {
  void* highs = Highs_create();
  const double inf = Highs_getInfinity(highs);

  const HighsInt num_col = 2, num_row = 1, num_nz = 2;
  const double col_cost[2] = {1.0, 1.0};
  const double col_lower[2] = {0.0, 0.0};
  const double col_upper[2] = {1.0, 1.0};
  const double row_lower[1] = {1.0};
  const double row_upper[1] = {inf};
  const HighsInt a_start[3] = {0, 1, 2};
  const HighsInt a_index[2] = {0, 0};
  const double a_value[2] = {1.0, 1.0};

  callback_state state = {0, 0};
  Highs_setCallback(highs, logging_callback, &state);
  Highs_startCallback(highs, kHighsCallbackLogging);
  Highs_setBoolOptionValue(highs, "output_flag", 1);
  Highs_setBoolOptionValue(highs, "log_to_console", 0);

  const HighsInt status_pass = Highs_passModel(
      highs, num_col, num_row, num_nz, 0, kHighsMatrixFormatColwise,
      kHighsHessianFormatTriangular, kHighsObjSenseMinimize, 0.0, col_cost,
      col_lower, col_upper, row_lower, row_upper, a_start, a_index, a_value,
      NULL, NULL, NULL, NULL);
  const HighsInt status_run = Highs_run(highs);
  Highs_stopCallback(highs, kHighsCallbackLogging);
  const HighsInt model_status = Highs_getModelStatus(highs);
  const double objective = Highs_getObjectiveValue(highs);

  SEXP col_value = PROTECT(Rf_allocVector(REALSXP, num_col));
  Highs_getSolution(highs, REAL(col_value), NULL, NULL, NULL);

  SEXP out = PROTECT(Rf_allocVector(VECSXP, 9));
  SET_VECTOR_ELT(out, 0, Rf_mkString(Highs_version()));
  SET_VECTOR_ELT(out, 1, Rf_ScalarInteger(Highs_getSizeofHighsInt(highs)));
  SET_VECTOR_ELT(out, 2, Rf_ScalarInteger(status_pass));
  SET_VECTOR_ELT(out, 3, Rf_ScalarInteger(status_run));
  SET_VECTOR_ELT(out, 4, Rf_ScalarInteger(model_status));
  SET_VECTOR_ELT(out, 5, Rf_ScalarReal(objective));
  SET_VECTOR_ELT(out, 6, col_value);
  SET_VECTOR_ELT(out, 7, Rf_ScalarInteger(state.calls));
  SET_VECTOR_ELT(out, 8, Rf_ScalarInteger(state.item_ok));

  SEXP names = PROTECT(Rf_allocVector(STRSXP, 9));
  SET_STRING_ELT(names, 0, Rf_mkChar("version"));
  SET_STRING_ELT(names, 1, Rf_mkChar("sizeof_int"));
  SET_STRING_ELT(names, 2, Rf_mkChar("status_pass"));
  SET_STRING_ELT(names, 3, Rf_mkChar("status_run"));
  SET_STRING_ELT(names, 4, Rf_mkChar("model_status"));
  SET_STRING_ELT(names, 5, Rf_mkChar("objective"));
  SET_STRING_ELT(names, 6, Rf_mkChar("col_value"));
  SET_STRING_ELT(names, 7, Rf_mkChar("callback_calls"));
  SET_STRING_ELT(names, 8, Rf_mkChar("callback_item_ok"));
  Rf_setAttrib(out, R_NamesSymbol, names);

  Highs_destroy(highs);
  UNPROTECT(3);
  return out;
}
