
#include "highs_interface.h"


using namespace Rcpp;


static void R_message_handler(HighsLogType type, const char* message, void* log_callback_data) {
  Rcpp::Rcout << message << std::endl;
}


// [[Rcpp::export]]
SEXP new_model() {
    Rcpp::XPtr<HighsModel> highs_model(new HighsModel(), true);
    return highs_model;
}

// [[Rcpp::export]]
SEXP model_set_ncol(SEXP mpt, int32_t ncol) {
    Rcpp::XPtr<HighsModel>model(mpt);
    model->lp_.num_col_ = ncol;
    return R_NilValue;
}

// [[Rcpp::export]]
SEXP model_set_nrow(SEXP mpt, int32_t nrow) {
    Rcpp::XPtr<HighsModel>model(mpt);
    model->lp_.num_row_ = nrow;
    return R_NilValue;
}

// [[Rcpp::export]]
SEXP model_set_sense(SEXP mpt, bool maximum) {
    Rcpp::XPtr<HighsModel>model(mpt);
    if (maximum) {
        model->lp_.sense_ = ObjSense::kMaximize;
    } else {
        model->lp_.sense_ = ObjSense::kMinimize;
    }
    return R_NilValue;
}

// [[Rcpp::export]]
SEXP model_set_offset(SEXP mpt, double_t offset) {
    Rcpp::XPtr<HighsModel>model(mpt);
    model->lp_.offset_ = offset;
    return R_NilValue;
}

// [[Rcpp::export]]
SEXP model_set_objective(SEXP mpt, std::vector<double> objective) {
    Rcpp::XPtr<HighsModel>model(mpt);
    model->lp_.col_cost_ = objective;
    return R_NilValue;
}

// [[Rcpp::export]]
SEXP model_set_lower(SEXP mpt, std::vector<double> lower) {
    Rcpp::XPtr<HighsModel>model(mpt);
    model->lp_.col_lower_ = lower;
    return R_NilValue;
}

// [[Rcpp::export]]
SEXP model_set_upper(SEXP mpt, std::vector<double> upper) {
    Rcpp::XPtr<HighsModel>model(mpt);
    model->lp_.col_upper_ = upper;
    return R_NilValue;
}

// [[Rcpp::export]]
SEXP model_set_constraint_matrix(SEXP mpt, std::string format,
                                 std::vector<int32_t> start,
                                 std::vector<int32_t> index,
                                 std::vector<double> value) {
    Rcpp::XPtr<HighsModel>model(mpt);
    if (format == "colwise") {
        model->lp_.a_matrix_.format_ = MatrixFormat::kColwise;
    } else if (format == "rowwise") {
        model->lp_.a_matrix_.format_ = MatrixFormat::kRowwise;
    } else if (format == "rowwise_partitioned") {
        model->lp_.a_matrix_.format_ = MatrixFormat::kRowwisePartitioned;
    } else {
        Rcpp::stop("unkown format!");
    }
    model->lp_.a_matrix_.start_ = start;
    model->lp_.a_matrix_.index_ = index;
    model->lp_.a_matrix_.value_ = value;
    return R_NilValue;
}

// [[Rcpp::export]]
SEXP model_set_lhs(SEXP mpt, std::vector<double> lower) {
    Rcpp::XPtr<HighsModel>model(mpt);
    model->lp_.row_lower_ = lower;
    return R_NilValue;
}

// [[Rcpp::export]]
SEXP model_set_rhs(SEXP mpt, std::vector<double> upper) {
    Rcpp::XPtr<HighsModel>model(mpt);
    model->lp_.row_upper_ = upper;
    return R_NilValue;
}

// [[Rcpp::export]]
SEXP model_set_hessian(SEXP mpt, std::string format, int32_t dim,
                       std::vector<int32_t> start,
                       std::vector<int32_t> index,
                       std::vector<double> value) {
    Rcpp::XPtr<HighsModel>model(mpt);
    model->hessian_.dim_ = dim;
    if (format == "triangular") {
        model->hessian_.format_ = HessianFormat::kTriangular;
    } else if (format == "square") {
        model->hessian_.format_ = HessianFormat::kSquare;
    } else {
        Rcpp::stop("unkown format!");
    }
    model->hessian_.start_ = start;
    model->hessian_.index_ = index;
    model->hessian_.value_ = value;
    return R_NilValue;
}

// [[Rcpp::export]]
SEXP model_set_vartype(SEXP mpt, std::vector<int32_t> type) {
    Rcpp::XPtr<HighsModel>model(mpt);
    if (model->lp_.integrality_.size() < type.size()) {
        model->lp_.integrality_.resize(type.size());
    }
    std::vector<HighsVarType> variable_types = {
        HighsVarType::kContinuous, HighsVarType::kInteger, HighsVarType::kSemiContinuous,
        HighsVarType::kSemiInteger, HighsVarType::kImplicitInteger};
    for (std::size_t i = 0; i < type.size(); ++i) {
        model->lp_.integrality_[i] = variable_types[type[i]];
    }
    return R_NilValue;
}

// [[Rcpp::export]]
int32_t model_get_nvars(SEXP mpt) {
    Rcpp::XPtr<HighsModel>model(mpt);
    return static_cast<int32_t>(model->lp_.num_col_);
}

// [[Rcpp::export]]
int32_t model_get_ncons(SEXP mpt) {
    Rcpp::XPtr<HighsModel>model(mpt);
    return static_cast<int32_t>(model->lp_.num_row_);
}

// [[Rcpp::export]]
Rcpp::IntegerVector model_get_vartype(SEXP mpt) {
    Rcpp::XPtr<HighsModel>model(mpt);
    IntegerVector type(model->lp_.integrality_.size());
    for (R_xlen_t i = 0; i < type.size(); ++i) {
        type[i] = static_cast<int32_t>(model->lp_.integrality_[i]);
    }
    return type;
}

// [[Rcpp::export]]
SEXP new_solver(SEXP mpt) {
    Rcpp::XPtr<HighsModel>model(mpt);
    Rcpp::XPtr<Highs> highs(new Highs(), true);
    highs->setLogCallback(R_message_handler);
    HighsStatus return_status = highs->passModel(*model.get());
    if (return_status != HighsStatus::kOk) {
        return R_NilValue;
    }
    return highs;
}

// [[Rcpp::export]]
int32_t solver_set_option(SEXP hi, std::string key, SEXP value) {
    Rcpp::XPtr<Highs>highs(hi);
    HighsStatus status;
    if (Rf_isLogical(value)) {
        bool logval = Rcpp::as<bool>(value);
        status = highs->setOptionValue(key, logval);
    } else if (Rf_isInteger(value)) {
        HighsInt intval = Rcpp::as<int32_t>(value);
        status = highs->setOptionValue(key, intval);
    } else if (Rf_isNumeric(value)) {
        double numval = Rcpp::as<double>(value);
        status = highs->setOptionValue(key, numval);
    } else if (Rf_isString(value)) {
        std::string strval = Rcpp::as<std::string>(value);
        status = highs->setOptionValue(key, strval);
    } else {
        Rcpp::stop("unkown type of value.");
    }
    return static_cast<int32_t>(status);
}


// [[Rcpp::export]]
int32_t solver_clear(SEXP hi) {
    Rcpp::XPtr<Highs>highs(hi);
    HighsStatus return_status = highs->clear();
    return static_cast<int32_t>(return_status);
}


// [[Rcpp::export]]
int32_t solver_clear_model(SEXP hi) {
    Rcpp::XPtr<Highs>highs(hi);
    HighsStatus return_status = highs->clearModel();
    return static_cast<int32_t>(return_status);
}


// [[Rcpp::export]]
int32_t solver_clear_solver(SEXP hi) {
    Rcpp::XPtr<Highs>highs(hi);
    HighsStatus return_status = highs->clearSolver();
    return static_cast<int32_t>(return_status);
}



// [[Rcpp::export]]
int32_t solver_run(SEXP hi) {
    Rcpp::XPtr<Highs>highs(hi);
    HighsStatus return_status = highs->run();
    return static_cast<int32_t>(return_status);
}


//
// int32_t solver_presolve(SEXP hi) {
//     Rcpp::XPtr<Highs>highs(hi);
//     HighsStatus return_status = highs->presolve();
//     return static_cast<int32_t>(return_status);
// }


// TODO
// HighsStatus postsolve(const HighsSolution& solution, const HighsBasis& basis);


// TODO
// const HighsLp& getPresolvedLp()


// TODO
// const HighsModel& getPresolvedModel()


// [[Rcpp::export]]
int32_t solver_write_model(SEXP hi, const std::string filename) {
    Rcpp::XPtr<Highs>highs(hi);
    HighsStatus return_status = highs->writeModel(filename);
    return static_cast<int32_t>(return_status);
}


// [[Rcpp::export]]
int32_t solver_write_basis(SEXP hi, const std::string filename) {
    Rcpp::XPtr<Highs>highs(hi);
    HighsStatus return_status = highs->writeBasis(filename);
    return static_cast<int32_t>(return_status);
}


// [[Rcpp::export]]
std::string solver_status_message(SEXP hi) {
    Rcpp::XPtr<Highs>highs(hi);
    const HighsModelStatus& model_status = highs->getModelStatus();
    return highs->modelStatusToString(model_status);
}


// [[Rcpp::export]]
int32_t solver_status(SEXP hi) {
    Rcpp::XPtr<Highs>highs(hi);
    const HighsModelStatus& model_status = highs->getModelStatus();
    return static_cast<int32_t>(model_status);
}


// [[Rcpp::export]]
double_t  solver_infinity(SEXP hi) {
    Rcpp::XPtr<Highs>highs(hi);
    return highs->getInfinity();
}


// [[Rcpp::export]]
Rcpp::List solver_info(SEXP hi) {
    Rcpp::XPtr<Highs>highs(hi);
    const HighsInfo& info = highs->getInfo();
    List z = List::create(
        Named("valid") = info.valid,
        Named("mip_node_count") = info.mip_node_count,
        Named("simplex_iteration_count") = info.simplex_iteration_count,
        Named("ipm_iteration_count") = info.ipm_iteration_count,
        Named("qp_iteration_count") = info.qp_iteration_count,
        Named("crossover_iteration_count") = info.crossover_iteration_count,
        Named("primal_solution_status") = highs->solutionStatusToString(info.primal_solution_status),
        Named("dual_solution_status") = highs->solutionStatusToString(info.dual_solution_status),
        Named("basis_validity") = info.basis_validity,
        Named("objective_function_value") = info.objective_function_value,
        Named("mip_dual_bound") = info.mip_dual_bound,
        Named("mip_gap") = info.mip_gap,
        Named("num_primal_infeasibilities") = info.num_primal_infeasibilities,
        Named("max_primal_infeasibility") = info.max_primal_infeasibility,
        Named("sum_primal_infeasibilities") = info.sum_primal_infeasibilities,
        Named("num_dual_infeasibilities") = info.num_dual_infeasibilities,
        Named("max_dual_infeasibility") = info.max_dual_infeasibility,
        Named("sum_dual_infeasibilities") = info.sum_dual_infeasibilities
    );
    return z;
}

// [[Rcpp::export]]
Rcpp::List solver_solution(SEXP hi) {
    Rcpp::XPtr<Highs>highs(hi);
    const HighsSolution& solution = highs->getSolution();
    List z = List::create(
        Named("value_valid") = solution.value_valid,
        Named("dual_valid") = solution.dual_valid,
        Named("col_value") = solution.col_value,
        Named("col_dual") = solution.col_dual,
        Named("row_value") = solution.row_value,
        Named("row_dual") = solution.row_dual
    );
    return z;
}

// [[Rcpp::export]]
bool solver_get_bool_option(SEXP hi, std::string key) {
    Rcpp::XPtr<Highs>highs(hi);
    bool value;
    highs->getOptionValue(key, value);
    // HighsInt    
    return value;
}

// [[Rcpp::export]]
int32_t solver_get_int_option(SEXP hi, std::string key) {
    Rcpp::XPtr<Highs>highs(hi);
    HighsInt value;
    highs->getOptionValue(key, value);
    return value;
}

// [[Rcpp::export]]
double_t solver_get_dbl_option(SEXP hi, std::string key) {
    Rcpp::XPtr<Highs>highs(hi);
    double_t value;
    highs->getOptionValue(key, value);
    return value;
}

// [[Rcpp::export]]
std::string solver_get_str_option(SEXP hi, std::string key) {
    Rcpp::XPtr<Highs>highs(hi);
    std::string value;
    highs->getOptionValue(key, value);   
    return value;
}


// [[Rcpp::export]]
int32_t solver_change_bounds(SEXP hi, IntegerVector idx, NumericVector lower, NumericVector upper) {
    Rcpp::XPtr<Highs>highs(hi);
    HighsStatus return_status = highs->changeColsBounds(idx.size(), &(idx[0]), &(lower[0]), &(upper[0]));
    return static_cast<int32_t>(return_status);
}


// [[Rcpp::export]]
int32_t solver_change_lrhs(SEXP hi, IntegerVector idx, NumericVector lhs, NumericVector rhs) {
    Rcpp::XPtr<Highs>highs(hi);
    HighsStatus return_status = highs->changeRowsBounds(idx.size(), &(idx[0]), &(lhs[0]), &(rhs[0]));
    return static_cast<int32_t>(return_status);
}


// [[Rcpp::export]]
int32_t solver_add_rows(SEXP hi, NumericVector lhs, NumericVector rhs,
    IntegerVector start, IntegerVector index, NumericVector value) {
    Rcpp::XPtr<Highs>highs(hi);
    HighsStatus return_status = highs->addRows(
        lhs.size(), //!< num_new_row = Number of new rows
        &(lhs[0]),  //!< lower = Array of size num_new_row with lower bounds
        &(rhs[0]),  //!< upper = Array of size num_new_row with upper bounds
        value.size(), //!< num_new_nz = Number of new nonzeros
        &(start[0]),   //!< starts = Array of size num_new_row with start indices of the rows
        &(index[0]),   //!< indices = Array of size num_new_nz with column indices for all rows
        &(value[0])    //!< values = Array of size num_new_nz with column values for all rows
    );
    return static_cast<int32_t>(return_status);
}


// [[Rcpp::export]]
int32_t solver_add_cols(SEXP hi, NumericVector costs,
    NumericVector lower, NumericVector upper,
    IntegerVector start, IntegerVector index, NumericVector value) {
    Rcpp::XPtr<Highs>highs(hi);
    HighsStatus return_status = highs->addCols(
        lower.size(), //!< num_new_col = Number of new columns
        &(costs[0]),   //!< costs = Array of size num_new_col with costs
        &(lower[0]),   //!< lower = Array of size num_new_col with lower bounds
        &(upper[0]),   //!< upper = Array of size num_new_col with upper bounds
        value.size(), //!< num_new_nz = Number of new nonzeros
        &(start[0]),   //!< starts = Array of size num_new_row with start indices of the columns
        &(index[0]),   //!< indices = Array of size num_new_nz with row indices for all columns
        &(value[0])    //!< values = Array of size num_new_nz with row values for all columns
    );
    return static_cast<int32_t>(return_status);
}

